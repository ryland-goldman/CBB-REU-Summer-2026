"""
CESR gun in WarpX (RZ): accelerate the cathode-emitted electrons through the
scaled CESR_gun.gdf electrode field with self-consistent (electromagnetostatic)
space charge, then reconstruct the timed exit beam for the injector handoff.

Drives lume-warpx from config/gun.yaml (which holds every constant); this module reads those
back, overrides only the runtime-computed values (dt, step count, the seed arrays, diagnostic
periods), and builds up the time-release beam with a beforestep callback. The 1 nC bunch is
released over the 2 ns grid pulse (the physical low-line-density representation), so no single
volumetric dump holds the ~0.4 m drifting beam; build_exit_handoff() reconstructs the full exit
beam by particle id across the dumps. See docs/gun.md for physics, the beam representation, the
exit handoff, and gotchas.

main() runs ONLY the simulation; sim/plot/gun.py produces the figures.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import shutil

import numpy as np

from sim.helpers.tools import C_LIGHT as c, M_E as m_e, E_CHARGE as q_e, prepare_env
from sim.helpers.loadparticles import (
    make_particle_group, downsample, open_particle_series, write_openpmd_particles)
from sim.helpers.buildfields import build_gun_field

CONFIG = "config/gun.yaml"
CATHODE_DIAG = "logs/diags/cathode/particles"
DIAG_DIR = "logs/diags/gun"
HANDOFF_DIR = "logs/diags/gun/handoff"   # reconstruct the exit beam here for the injector


def load_cathode_bunch(rmax, zmax, bunch_charge, rng_seed, max_part, pulse_width):
    """Import the last cathode snapshot and remap the (x, z) slab into RZ.

    Returns dict of x, y, z, ux, uy, uz, w, t arrays for the seed + time-release callback
    (ux/uy/uz are proper velocity u = γβc [m/s]); t is the per-macroparticle emission time,
    released over pulse_width and t-sorted so the per-step injection callback walks them in
    one pass.
    """
    ts = open_particle_series(CATHODE_DIAG, "cathode")
    it = ts.iterations[-1]
    x, z, ux, uy, uz, w = ts.get_particle(
        ["x", "z", "ux", "uy", "uz", "w"], species="electrons", iteration=it)

    rng = np.random.default_rng(rng_seed)
    r = np.abs(x)
    keep = r < rmax
    if not keep.any():
        raise RuntimeError(f"no cathode particles with r < rmax={rmax} m")
    # Keep the masked signed x (`xk`) so the radial-momentum sign survives the downsample.
    xk = x[keep]
    r, z, ux, uy, uz, w = (a[keep] for a in (r, z, ux, uy, uz, w))

    (xk, r, z, ux, uy, uz), w = downsample((xk, r, z, ux, uy, uz), w, max_part, rng)

    # slab(x) → RZ disc: importance-resample (with replacement) by r·w to supply the 2πr
    # revolution Jacobian the naive r=|x| map omits (else n(r) ∝ 1/r on-axis cusp).
    if r.max() > 0.0:
        rw = r * w
        sel = rng.choice(r.size, r.size, replace=True, p=rw / rw.sum())
        xk, r, z, ux, uy, uz, w = (a[sel] for a in (xk, r, z, ux, uy, uz, w))

    theta = rng.uniform(0.0, 2.0 * np.pi, size=r.size)
    ct, st = np.cos(theta), np.sin(theta)
    # slab x -> radius; transverse momentum: radial = ux·sign(x), azimuthal = uy
    ur = ux * np.sign(np.where(xk == 0.0, 1.0, xk))
    xpos, ypos = r * ct, r * st
    uxn = ur * ct - uy * st
    uyn = ur * st + uy * ct
    zpos = np.clip(z, 0.0, zmax)
    w = w * (bunch_charge / (w.sum() * q_e))            # renormalize to BUNCH_CHARGE

    # Emission time per macroparticle: released over PULSE_WIDTH, t-sorted so the per-step
    # injection callback walks them in one pass.
    t = rng.uniform(0.0, pulse_width, size=r.size)
    order = np.argsort(t)
    t = t[order]
    xpos, ypos, zpos, uxn, uyn, uz, w = (
        a[order] for a in (xpos, ypos, zpos, uxn, uyn, uz, w))

    print(f"Imported {r.size} macroparticles from cathode (iter {it}); "
          f"renormalized to {bunch_charge*1e9:.3f} nC, r ≤ {r.max()*1e3:.2f} mm; "
          f"timed release over {pulse_width*1e9:.1f} ns", flush=True)
    # openPMD ux/uy/uz are dimensionless γβ → proper velocity u = γβc [m/s], so ×c.
    return dict(x=xpos, y=ypos, z=zpos, ux=uxn * c, uy=uyn * c, uz=uz * c, w=w, t=t)


def build_exit_handoff(zmax_field):
    """Reconstruct the full time-released exit beam and write it for the injector.

    Reconstruct from the volumetric dumps by particle id: take each id's first appearance in
    the FIELD-FREE pad (z ≥ zmax_field) as its exit-plane phase space, then drift all to a
    common reference time. Sampling in the pad (not the last in-field dump) preserves εn,x.
    Writes an openPMD dump to HANDOFF_DIR. See docs/gun.md → Exit-beam handoff.
    """
    ts = open_particle_series(os.path.join(DIAG_DIR, "particles"), "gun")

    cols = {k: [] for k in ("id", "t", "x", "y", "z", "ux", "uy", "uz", "w")}
    for it, t_dump in zip(ts.iterations, ts.t):
        idp, x, y, z, ux, uy, uz, w = ts.get_particle(
            ["id", "x", "y", "z", "ux", "uy", "uz", "w"], species="electrons", iteration=it)
        if len(idp) == 0:
            continue
        cols["id"].append(idp.astype(np.int64)); cols["t"].append(np.full(len(idp), t_dump))
        cols["x"].append(x); cols["y"].append(y); cols["z"].append(z)
        cols["ux"].append(ux); cols["uy"].append(uy); cols["uz"].append(uz); cols["w"].append(w)
    if not cols["id"]:
        print("  handoff: no particles tracked — skipped", flush=True)
        return
    cat = {k: np.concatenate(v) for k, v in cols.items()}
    n_ids = np.unique(cat["id"]).size

    # First appearance in the field-free pad (z ≥ zmax_field) per id.
    pad = np.where(cat["z"] >= zmax_field)[0]
    if pad.size == 0:
        print(f"  handoff: no particle reached the field-free pad (z ≥ {zmax_field*1e3:.1f} mm) "
              "— run longer? skipped", flush=True)
        return
    order = np.lexsort((cat["t"][pad], cat["id"][pad]))
    sid = cat["id"][pad][order]
    first = np.empty(sid.size, dtype=bool)
    first[0] = True
    first[1:] = sid[1:] != sid[:-1]                     # True at each id-group start (earliest t)
    pick = pad[order[first]]
    t_s = cat["t"][pick]
    x, y, z = cat["x"][pick], cat["y"][pick], cat["z"][pick]
    ux, uy, uz, w = cat["ux"][pick], cat["uy"][pick], cat["uz"][pick], cat["w"][pick]
    n_exit = pick.size

    # Classify the ids that never reached the pad: still in the final dump → un-flushed tail;
    # else absorbed before the pad → radial (r=RMAX) loss.
    exited_ids = set(cat["id"][pick].tolist())
    final_ids = set(cat["id"][cat["t"] == ts.t[-1]].tolist()) if len(ts.t) else set()
    not_exited = [i for i in np.unique(cat["id"]).tolist() if i not in exited_ids]
    n_unflushed = sum(1 for i in not_exited if i in final_ids)
    n_radial = len(not_exited) - n_unflushed

    # Ballistic drift of the field-free samples to a common reference time t_ref = max(t_s).
    gamma = np.sqrt(1.0 + ux**2 + uy**2 + uz**2)       # γ from u = γβ
    vx, vy, vz = (c * ux / gamma, c * uy / gamma, c * uz / gamma)
    dtau = t_s.max() - t_s
    xh, yh, zh = x + vx * dtau, y + vy * dtau, z + vz * dtau

    pg = make_particle_group(xh, yh, zh, ux, uy, uz, w)
    if os.path.isdir(HANDOFF_DIR):
        shutil.rmtree(HANDOFF_DIR)
    write_openpmd_particles(pg, HANDOFF_DIR, iteration=0, time=0.0)
    print(f"  handoff: {n_exit}/{n_ids} macroparticles cleared the field "
          f"({pg.charge*1e9:.3f} nC, transmission {100*n_exit/n_ids:.0f}%; "
          f"{n_radial} radial loss, {n_unflushed} un-flushed), "
          f"z-extent {(zh.max()-zh.min())*1e3:.0f} mm → {HANDOFF_DIR}", flush=True)


def main():
    prepare_env()
    from warpx import WarpX

    w = WarpX(input_file=CONFIG, path=DIAG_DIR)
    nr, nz = w.get("grid/number_of_cells")
    rmax, zmax = w.get("grid/upper_bound")
    p = w.get("params")
    gun_voltage = p["GUN_VOLTAGE"]

    # Build the applied gun field from the GDF (single-source: field-map scale + the kinematics
    # below agree on GUN_VOLTAGE). Idempotent — a fresh checkout rebuilds fieldmaps/h5/gun_E.h5.
    build_gun_field(gun_voltage)

    bunch = load_cathode_bunch(rmax, zmax, p["BUNCH_CHARGE"], p["RNG_SEED"], p["MAX_PART"],
                               p["PULSE_WIDTH"])

    # ── Time step / duration ──────────────────────────────────────────────────
    gamma = 1.0 + q_e * gun_voltage / (m_e * c**2)
    v_exit = c * np.sqrt(1.0 - 1.0 / gamma**2)
    dt = p["CFL"] * (zmax / nz) / v_exit
    # Run length is sized on the FIELD transit (ZMAX_FIELD), not the padded domain, so the run
    # stops while the beam is in the pad (the padded domain drains and aborts MLMG if over-run).
    transit_field = p["ZMAX_FIELD"] / (p["AVG_SPEED_FRAC"] * v_exit)
    run_time = p["PULSE_WIDTH"] + p["TRANSIT_MARGIN"] * transit_field
    max_steps = p["MAX_STEPS"] or int(run_time / dt)
    period = max(1, max_steps // p["N_DIAGS"])

    print(f"Gun: {gun_voltage/1e3:.0f} kV  ->  γ={gamma:.3f}, β={v_exit/c:.3f}, "
          f"v_exit={v_exit:.2e} m/s", flush=True)
    print(f"dt = {dt:.3e} s, max_steps = {max_steps} "
          f"(release over {p['PULSE_WIDTH']*1e9:.1f} ns + transit)", flush=True)

    # Fresh diags: the h5 backend appends one file per dump, so stale files would corrupt plots.
    if os.path.isdir(os.path.join(DIAG_DIR, "fields")):
        shutil.rmtree(os.path.join(DIAG_DIR, "fields"))
    if os.path.isdir(os.path.join(DIAG_DIR, "particles")):
        shutil.rmtree(os.path.join(DIAG_DIR, "particles"))

    # Seed only the earliest macroparticle (PICMI needs a non-empty initial distribution);
    # inject the rest over the pulse via the callback below.
    w.update({
        "simulation/max_steps": max_steps,
        "simulation/time_step_size": dt,
        "species/0/x": bunch["x"][0:1], "species/0/y": bunch["y"][0:1], "species/0/z": bunch["z"][0:1],
        "species/0/ux": bunch["ux"][0:1], "species/0/uy": bunch["uy"][0:1], "species/0/uz": bunch["uz"][0:1],
        "species/0/weight": bunch["w"][0:1],
        "diagnostics/0/period": period,
        "diagnostics/1/period": period,
    })

    # ── Time-release injection ────────────────────────────────────────────────
    # Build up the bunch over PULSE_WIDTH: each step inject the macroparticles whose emission
    # time falls in the step window. add_particles uses the WRITE path (the broken RZ accessor
    # is the READ path). ux/uy/uz are proper velocity γβc [m/s].
    from pywarpx import particle_containers
    pc = [None]
    bt = bunch["t"]
    state = {"next": 1, "step": 0}

    def _inject():
        if pc[0] is None:
            pc[0] = particle_containers.ParticleContainerWrapper("electrons")
        t_hi = (state["step"] + 1) * dt
        j = k = state["next"]
        while k < bt.size and bt[k] < t_hi:
            k += 1
        if k > j:
            pc[0].add_particles(
                x=bunch["x"][j:k], y=bunch["y"][j:k], z=bunch["z"][j:k],
                ux=bunch["ux"][j:k], uy=bunch["uy"][j:k], uz=bunch["uz"][j:k], w=bunch["w"][j:k])
            state["next"] = k
        state["step"] += 1

    w.install_callback("beforestep", _inject)

    print(f"\nRunning {max_steps} steps (diag every {period}) …")
    w.run(progress="gun")
    print(f"Injected {state['next']}/{bunch['t'].size} macroparticles over the pulse", flush=True)
    build_exit_handoff(p["ZMAX_FIELD"])
    print(f"\nDone. openPMD output → {DIAG_DIR}/{{fields,particles}}/")


if __name__ == "__main__":
    main()
