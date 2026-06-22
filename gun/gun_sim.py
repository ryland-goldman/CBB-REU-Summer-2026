"""
CESR gun in WarpX (RZ): accelerate the cathode-emitted electrons through the
scaled CESR_gun.gdf electrode field with self-consistent (electromagnetostatic)
space charge, then reconstruct the timed exit beam for the injector handoff.

Drives lume-warpx from gun/gun.yaml (which holds every constant); this module reads those
back, overrides only the runtime-computed values (dt, step count, the seed arrays, diagnostic
periods), and builds up the time-release beam with a beforestep callback. See gun/README.md
for physics, the beam representation, the exit handoff, and gotchas.
"""

import os
import shutil
import sys

# Self-insert the repo root so `python gun/gun_sim.py` resolves pipeline.* standalone.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from openpmd_viewer import OpenPMDTimeSeries

from pipeline.constants import C_LIGHT as c, M_E as m_e, E_CHARGE as q_e, MC2_EV
from gun.build_gun_field import GUN_VOLTAGE   # single-source: field-map scale + kinematics agree

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG = os.path.join(HERE, "gun.yaml")
CATHODE_DIAG = "cathode/diags/particles"
DIAG_DIR = "gun/diags"
HANDOFF_DIR = "gun/diags/handoff"      # timed mode reconstructs the exit beam here for the injector


def load_cathode_bunch(rmax, zmax, bunch_charge, rng_seed, max_part, beam_release, pulse_width):
    """Import the last cathode snapshot and remap the (x, z) slab into RZ.

    Returns dict of x, y, z, ux, uy, uz, w, t arrays for the seed + time-release callback
    (ux/uy/uz are proper velocity u = γβc [m/s]).
    """
    ts = OpenPMDTimeSeries(CATHODE_DIAG)
    if len(ts.iterations) == 0:
        raise RuntimeError(f"{CATHODE_DIAG} has no iterations — did the cathode stage run?")
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

    if max_part and r.size > max_part:
        sel = rng.choice(r.size, max_part, replace=False)
        scale_w = r.size / max_part
        xk, r, z, ux, uy, uz, w = (a[sel] for a in (xk, r, z, ux, uy, uz, w))
        w = w * scale_w

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

    # Emission time per macroparticle. timed: released over PULSE_WIDTH, t-sorted so the
    # per-step injection callback walks them in one pass. snapshot: all emitted at t=0.
    if beam_release == "timed":
        t = rng.uniform(0.0, pulse_width, size=r.size)
        order = np.argsort(t)
        t = t[order]
        xpos, ypos, zpos, uxn, uyn, uz, w = (
            a[order] for a in (xpos, ypos, zpos, uxn, uyn, uz, w))
    else:
        t = np.zeros(r.size)

    print(f"Imported {r.size} macroparticles from cathode (iter {it}); "
          f"renormalized to {bunch_charge*1e9:.3f} nC, r ≤ {r.max()*1e3:.2f} mm; "
          f"release={beam_release}"
          + (f" over {pulse_width*1e9:.1f} ns" if beam_release == "timed" else ""), flush=True)
    # openPMD ux/uy/uz are dimensionless γβ → proper velocity u = γβc [m/s], so ×c.
    return dict(x=xpos, y=ypos, z=zpos, ux=uxn * c, uy=uyn * c, uz=uz * c, w=w, t=t)


def build_exit_handoff(zmax_field):
    """Reconstruct the full time-released exit beam and write it for the injector.

    Reconstruct from the volumetric dumps by particle id: take each id's first appearance in
    the FIELD-FREE pad (z ≥ zmax_field) as its exit-plane phase space, then drift all to a
    common reference time. Sampling in the pad (not the last in-field dump) preserves εn,x.
    Writes an openPMD dump to HANDOFF_DIR via pipeline.impact_io. See README → Exit-beam handoff.
    """
    from pmd_beamphysics import ParticleGroup
    from pipeline.impact_io import write_openpmd_particles

    ts = OpenPMDTimeSeries(os.path.join(DIAG_DIR, "particles"))
    if len(ts.iterations) == 0:
        print("  handoff: no volumetric dumps to reconstruct from — skipped", flush=True)
        return

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

    pg = ParticleGroup(data=dict(
        x=xh, y=yh, z=zh,
        px=ux * MC2_EV, py=uy * MC2_EV, pz=uz * MC2_EV,   # γβ·(mc² in eV) [eV/c]
        t=np.zeros(xh.size), weight=w * q_e,              # macro count → charge [C]
        status=np.ones(xh.size, dtype=np.int64), species="electron"))
    if os.path.isdir(HANDOFF_DIR):
        shutil.rmtree(HANDOFF_DIR)
    write_openpmd_particles(pg, HANDOFF_DIR, iteration=0, time=0.0)
    print(f"  handoff: {n_exit}/{n_ids} macroparticles cleared the field "
          f"({pg.charge*1e9:.3f} nC, transmission {100*n_exit/n_ids:.0f}%; "
          f"{n_radial} radial loss, {n_unflushed} un-flushed), "
          f"z-extent {(zh.max()-zh.min())*1e3:.0f} mm → {HANDOFF_DIR}", flush=True)


def main():
    from warpx import WarpX

    w = WarpX(input_file=CONFIG, path="gun")
    nr, nz = w.get("grid/number_of_cells")
    rmax, zmax = w.get("grid/upper_bound")
    p = w.get("params")
    timed = (p["BEAM_RELEASE"] == "timed")

    bunch = load_cathode_bunch(rmax, zmax, p["BUNCH_CHARGE"], p["RNG_SEED"], p["MAX_PART"],
                               p["BEAM_RELEASE"], p["PULSE_WIDTH"])

    # ── Time step / duration ──────────────────────────────────────────────────
    gamma = 1.0 + q_e * GUN_VOLTAGE / (m_e * c**2)
    v_exit = c * np.sqrt(1.0 - 1.0 / gamma**2)
    dt = p["CFL"] * (zmax / nz) / v_exit
    # Run length is sized on the FIELD transit (ZMAX_FIELD), not the padded domain, so the run
    # stops while the beam is in the pad (the padded domain drains and aborts MLMG if over-run).
    transit_field = p["ZMAX_FIELD"] / (p["AVG_SPEED_FRAC"] * v_exit)
    run_time = (p["PULSE_WIDTH"] if timed else 0.0) + p["TRANSIT_MARGIN"] * transit_field
    max_steps = p["MAX_STEPS"] or int(run_time / dt)
    period = max(1, max_steps // p["N_DIAGS"])

    print(f"Gun: {GUN_VOLTAGE/1e3:.0f} kV  ->  γ={gamma:.3f}, β={v_exit/c:.3f}, "
          f"v_exit={v_exit:.2e} m/s", flush=True)
    print(f"dt = {dt:.3e} s, max_steps = {max_steps}"
          + (f" (release over {p['PULSE_WIDTH']*1e9:.1f} ns + transit)" if timed else ""), flush=True)

    if os.path.isdir(DIAG_DIR):
        shutil.rmtree(DIAG_DIR)

    # snapshot: seed the whole bunch. timed: seed only the earliest macroparticle, inject the
    # rest over the pulse via the callback below.
    sl = slice(0, 1) if timed else slice(None)
    w.update({
        "simulation/max_steps": max_steps,
        "simulation/time_step_size": dt,
        "species/0/x": bunch["x"][sl], "species/0/y": bunch["y"][sl], "species/0/z": bunch["z"][sl],
        "species/0/ux": bunch["ux"][sl], "species/0/uy": bunch["uy"][sl], "species/0/uz": bunch["uz"][sl],
        "species/0/weight": bunch["w"][sl],
        "diagnostics/0/period": period,
        "diagnostics/1/period": period,
    })

    # ── Time-release injection ────────────────────────────────────────────────
    # Build up the bunch over PULSE_WIDTH: each step inject the macroparticles whose emission
    # time falls in the step window. add_particles uses the WRITE path (the broken RZ accessor
    # is the READ path). ux/uy/uz are proper velocity γβc [m/s].
    state = {"next": 1, "step": 0}
    if timed:
        from pywarpx import particle_containers
        pc = [None]
        bt = bunch["t"]

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
    if timed:
        print(f"Injected {state['next']}/{bunch['t'].size} macroparticles over the pulse", flush=True)
        build_exit_handoff(p["ZMAX_FIELD"])
    print(f"\nDone. openPMD output → {DIAG_DIR}/{{fields,particles}}/")


if __name__ == "__main__":
    main()
