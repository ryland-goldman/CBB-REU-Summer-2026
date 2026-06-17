"""
CESR gun in WarpX (RZ): accelerate the cathode-emitted electrons through the
scaled CESR_gun.gdf electrode field with self-consistent (electromagnetostatic)
space charge, then reconstruct the timed exit beam for the injector handoff.

Second stage of the WarpX Linac chain. See gun/README.md for physics, the
beam-source slab→RZ remap, the timed/snapshot representation, the exit handoff,
and parameters. Use the gun facade (`import gun; gun.run()`); direct
`python gun/gun_sim.py` fails (needs the repo root on sys.path for pipeline.*).
"""

import os
import shutil

import numpy as np
import pywarpx
from pywarpx import picmi, callbacks, particle_containers
from openpmd_viewer import OpenPMDTimeSeries

from pipeline._runner import run_step

c = picmi.constants.c
m_e = picmi.constants.m_e
q_e = picmi.constants.q_e
ep0 = picmi.constants.ep0

# ── Gun / field-map parameters (must match build_gun_field.py) ────────────────
GUN_FIELD = "gun/gun_field/gun_E.h5"
GUN_VOLTAGE = 150.0e3        # [V]
RMAX = 0.015                 # field-map R extent [m]
ZMAX_FIELD = 0.051765        # field-map Z extent [m]; Ez→0 at the map edge (the exit plane)
# Field-free drift pad past the map (WarpX zero-fills applied field beyond it). Required so
# the timed handoff samples each particle in field-free space (drift-to-common-time preserves
# εn,x; sampling still-in-field manufactures a spurious x–u correlation, ~8× inflation).
ZPAD = 0.020                 # field-free drift pad past the field map [m]
ZMAX = ZMAX_FIELD + ZPAD     # full RZ domain z-extent [m]

CATHODE_DIAG = "cathode/diags/particles"
BUNCH_CHARGE = 1.0e-9        # renormalized gun bunch charge [C] (raw cathode snapshot is ~82 nC)
RNG_SEED = 0

# ── Beam representation (snapshot vs time-release); see README → "Beam source" ──
BEAM_RELEASE = "timed"       # "timed" → release over PULSE_WIDTH (default); "snapshot" → all at t=0
PULSE_WIDTH = 2.0e-9         # grid-pulse emission window [s] (uniform/flat-top)
HANDOFF_DIR = "gun/diags/handoff"   # timed mode reconstructs the exit beam here for the injector

# ── Grid (RZ, single azimuthal mode — the gun field is m = 0) ─────────────────
nr, nz = 128, 712            # divisible by the blocking factor (8)

# ── Diagnostics output directory ──────────────────────────────────────────────
DIAG_DIR = "gun/diags"

# ── Performance knobs (tunable via gun.config(...); see README → "Performance knobs") ─
REQUIRED_PRECISION = 1e-5            # MLMG Poisson solve relative tolerance
SPACE_CHARGE = True                  # beam self-field on/off (False → warpx_do_not_deposit)
MAX_ITERS = None                     # MLMG iteration cap (None → PICMI default)
CFL = 0.4                            # dt = CFL · dz / v_exit
TRANSIT_MARGIN = 1.15                # run length = TRANSIT_MARGIN × gun-transit time
AVG_SPEED_FRAC = 0.6                 # bunch average speed as a fraction of v_exit
MAX_STEPS = 0                        # 0 → auto-derive from CFL/margins; >0 → fixed
N_DIAGS = 40                         # openPMD dumps over the run (≥20 keeps space_charge.png snapshot)
MAX_PART = 0                         # 0/None → keep all cathode particles; >0 → cap


def load_cathode_bunch():
    """Import the last cathode snapshot and remap the (x, z) slab into RZ.

    Returns dict of x, y, z, ux, uy, uz, w arrays for ParticleListDistribution.
    """
    ts = OpenPMDTimeSeries(CATHODE_DIAG)
    if len(ts.iterations) == 0:
        raise RuntimeError(
            f"{CATHODE_DIAG} has no iterations — did the cathode stage run and "
            f"produce particles?")
    it = ts.iterations[-1]
    x, z, ux, uy, uz, w = ts.get_particle(
        ["x", "z", "ux", "uy", "uz", "w"], species="electrons", iteration=it,
    )

    rng = np.random.default_rng(RNG_SEED)
    r = np.abs(x)
    keep = r < RMAX
    if not keep.any():
        raise RuntimeError(
            f"no cathode particles with r < RMAX={RMAX} m; check RMAX or the "
            f"upstream cathode output")
    # Keep the masked signed x (`xk`) so the radial-momentum sign below survives the downsample.
    xk = x[keep]
    r, z, ux, uy, uz, w = (a[keep] for a in (r, z, ux, uy, uz, w))

    # Optionally downsample (reweighted to preserve total charge) to cap the cost.
    if MAX_PART and r.size > MAX_PART:
        sel = rng.choice(r.size, MAX_PART, replace=False)
        scale_w = r.size / MAX_PART
        xk, r, z, ux, uy, uz, w = (a[sel] for a in (xk, r, z, ux, uy, uz, w))
        w = w * scale_w

    # slab(x) → RZ disc: importance-resample (with replacement) by r·w to supply the 2πr
    # revolution Jacobian the naive r=|x| map omits (else n(r) ∝ 1/r on-axis cusp).
    # See README → "Beam source".
    if r.max() > 0.0:
        rw = r * w
        sel = rng.choice(r.size, r.size, replace=True, p=rw / rw.sum())
        xk, r, z, ux, uy, uz, w = (a[sel] for a in (xk, r, z, ux, uy, uz, w))

    theta = rng.uniform(0.0, 2.0 * np.pi, size=r.size)
    ct, st = np.cos(theta), np.sin(theta)

    # slab x -> radius; transverse momentum: radial = ux·sign(x), azimuthal = uy
    ur = ux * np.sign(np.where(xk == 0.0, 1.0, xk))
    ut = uy
    xpos = r * ct
    ypos = r * st
    uxn = ur * ct - ut * st
    uyn = ur * st + ut * ct
    zpos = np.clip(z, 0.0, ZMAX)

    # Renormalize weights so the imported distribution carries BUNCH_CHARGE.
    w = w * (BUNCH_CHARGE / (w.sum() * q_e))

    # Emission time per macroparticle. timed: released over PULSE_WIDTH, t-sorted so the
    # per-step injection callback walks them in one pass. snapshot: all emitted at t=0.
    if BEAM_RELEASE == "timed":
        t = rng.uniform(0.0, PULSE_WIDTH, size=r.size)
        order = np.argsort(t)
        t = t[order]
        xpos, ypos, zpos, uxn, uyn, uz, w = (
            a[order] for a in (xpos, ypos, zpos, uxn, uyn, uz, w))
    else:
        t = np.zeros(r.size)

    print(f"Imported {r.size} macroparticles from cathode (iter {it}); "
          f"renormalized to {BUNCH_CHARGE*1e9:.3f} nC, r ≤ {r.max()*1e3:.2f} mm; "
          f"release={BEAM_RELEASE}"
          + (f" over {PULSE_WIDTH*1e9:.1f} ns" if BEAM_RELEASE == "timed" else ""),
          flush=True)
    # openPMD ux/uy/uz are dimensionless γβ; PICMI ParticleListDistribution wants proper
    # velocity u = γβc [m/s], so ×c (else the beam is injected essentially at rest).
    return dict(x=xpos, y=ypos, z=zpos, ux=uxn * c, uy=uyn * c, uz=uz * c, w=w, t=t)


def build_exit_handoff():
    """Reconstruct the full time-released exit beam and write it for the injector.

    The ~2 ns stream's ballistic z-extent exceeds the gun domain, so reconstruct it from the
    volumetric dumps by particle id (the `pipeline/collimator.py` idiom): take each id's first
    appearance in the FIELD-FREE pad (z ≥ ZMAX_FIELD) as its exit-plane phase space, then drift
    all to a common reference time. Sampling in the pad — NOT the last in-field dump — preserves
    εn,x (a still-in-field sample manufactures a spurious x–u correlation). Ids never reaching
    the pad are classified r=RMAX radial loss vs un-flushed tail. See README → "Exit-beam handoff".

    Writes an openPMD dump to HANDOFF_DIR via `pipeline.impact_io`; the injector reads it when present.
    """
    from pmd_beamphysics import ParticleGroup
    from pipeline.impact_io import write_openpmd_particles

    MC2_EV = 510998.95069
    pdir = os.path.join(DIAG_DIR, "particles")
    ts = OpenPMDTimeSeries(pdir)
    if len(ts.iterations) == 0:
        print("  handoff: no volumetric dumps to reconstruct from — skipped", flush=True)
        return

    # Stack every (id, dump) row.
    cols = {k: [] for k in ("id", "t", "x", "y", "z", "ux", "uy", "uz", "w")}
    final_it = ts.iterations[-1]
    for it, t_dump in zip(ts.iterations, ts.t):
        idp, x, y, z, ux, uy, uz, w = ts.get_particle(
            ["id", "x", "y", "z", "ux", "uy", "uz", "w"],
            species="electrons", iteration=it)
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

    # First appearance in the field-free pad (z ≥ ZMAX_FIELD) per id: sort pad rows by (id, t),
    # take the earliest-t row of each id-group.
    pad = np.where(cat["z"] >= ZMAX_FIELD)[0]
    if pad.size == 0:
        print("  handoff: no particle reached the field-free pad (z ≥ "
              f"{ZMAX_FIELD*1e3:.1f} mm) — run longer? skipped", flush=True)
        return
    order = np.lexsort((cat["t"][pad], cat["id"][pad]))
    sid = cat["id"][pad][order]
    first = np.empty(sid.size, dtype=bool)
    first[0] = True
    first[1:] = sid[1:] != sid[:-1]                # True at each id-group start (earliest t)
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
    gamma = np.sqrt(1.0 + ux**2 + uy**2 + uz**2)  # γ from u = γβ
    vx, vy, vz = (c * ux / gamma, c * uy / gamma, c * uz / gamma)
    dtau = t_s.max() - t_s
    xh, yh, zh = x + vx * dtau, y + vy * dtau, z + vz * dtau

    pg = ParticleGroup(data=dict(
        x=xh, y=yh, z=zh,
        px=ux * MC2_EV, py=uy * MC2_EV, pz=uz * MC2_EV,   # γβ·m_ec² /c = γβ·(mc² in eV) [eV/c]
        t=np.zeros(xh.size),
        weight=w * q_e,                                   # macro count → charge [C]
        status=np.ones(xh.size, dtype=np.int64),
        species="electron",
    ))
    if os.path.isdir(HANDOFF_DIR):
        shutil.rmtree(HANDOFF_DIR)
    write_openpmd_particles(pg, HANDOFF_DIR, iteration=0, time=0.0)
    print(f"  handoff: {n_exit}/{n_ids} macroparticles cleared the field "
          f"({pg.charge*1e9:.3f} nC, transmission {100*n_exit/n_ids:.0f}%; "
          f"{n_radial} radial loss, {n_unflushed} un-flushed), "
          f"z-extent {(zh.max()-zh.min())*1e3:.0f} mm → {HANDOFF_DIR}", flush=True)


def main():
    grid = picmi.CylindricalGrid(
        number_of_cells=[nr, nz],
        n_azimuthal_modes=1,
        lower_bound=[0.0, 0.0],
        upper_bound=[RMAX, ZMAX],
        # r=0 is "none" (axis). The outer radial wall MUST be dirichlet (not neumann): the
        # magnetostatic A_z solve (driven by beam j_z) is otherwise an all-Neumann singular
        # operator and MLMG diverges. See README → "Space-charge model".
        lower_boundary_conditions=["none", "dirichlet"],
        upper_boundary_conditions=["dirichlet", "dirichlet"],
        lower_boundary_conditions_particles=["none", "absorbing"],
        upper_boundary_conditions_particles=["absorbing", "absorbing"],
        warpx_blocking_factor=8,
    )

    # Electromagnetostatic solver: warpx_magnetostatic=True adds the self-B (qβ×B pinch ⇒ net
    # qE_r/γ²) on top of the ES Poisson solve. See README → "Space-charge model".
    solver_kw = dict(grid=grid, method="Multigrid",
                     required_precision=REQUIRED_PRECISION,
                     warpx_magnetostatic=True,
                     warpx_magnetostatic_required_precision=REQUIRED_PRECISION,
                     warpx_self_fields_verbosity=0,
                     warpx_magnetostatic_verbosity=0)
    if MAX_ITERS:                                     # omit when None → PICMI default
        solver_kw["maximum_iterations"] = MAX_ITERS
        solver_kw["warpx_magnetostatic_max_iters"] = MAX_ITERS
    solver = picmi.ElectrostaticSolver(**solver_kw)

    # Applied gun electrode field (scaled CESR_gun.gdf), read from file via raw WarpX inputs
    # (PICMI has no tabulated particle-applied-field class). See README → "The gun field map".
    pywarpx.particles.E_ext_particle_init_style = "read_from_file"
    pywarpx.particles.read_fields_from_path = GUN_FIELD
    pywarpx.particles.B_ext_particle_init_style = "none"

    bunch = load_cathode_bunch()
    timed = (BEAM_RELEASE == "timed")
    # snapshot: seed the whole bunch. timed: seed only the earliest macroparticle (PICMI needs
    # a non-empty initial distribution), inject the rest over the pulse via the callback below.
    sl = slice(0, 1) if timed else slice(None)
    electrons = picmi.Species(
        particle_type="electron",
        name="electrons",
        initial_distribution=picmi.ParticleListDistribution(
            x=bunch["x"][sl], y=bunch["y"][sl], z=bunch["z"][sl],
            ux=bunch["ux"][sl], uy=bunch["uy"][sl], uz=bunch["uz"][sl],
            weight=bunch["w"][sl],
        ),
        warpx_do_not_deposit=not SPACE_CHARGE,   # SPACE_CHARGE=False → no beam self-field
    )

    # ── Time step / duration ──────────────────────────────────────────────────
    gamma = 1.0 + q_e * GUN_VOLTAGE / (m_e * c**2)
    v_exit = c * np.sqrt(1.0 - 1.0 / gamma**2)
    dz = ZMAX / nz
    dt = CFL * dz / v_exit
    # Run length is sized on the time to clear the FIELD region (ZMAX_FIELD), NOT the padded
    # domain: the run must stop while the beam is still in the pad, before it drains out ZMAX —
    # an empty domain aborts the MLMG solve (`MLMG failed`). MAX_STEPS (0 = auto) overrides.
    transit_field = ZMAX_FIELD / (AVG_SPEED_FRAC * v_exit)
    run_time = (PULSE_WIDTH if timed else 0.0) + TRANSIT_MARGIN * transit_field
    max_steps = MAX_STEPS or int(run_time / dt)

    print(f"Gun: {GUN_VOLTAGE/1e3:.0f} kV  ->  γ={gamma:.3f}, β={v_exit/c:.3f}, "
          f"v_exit={v_exit:.2e} m/s", flush=True)
    print(f"dt = {dt:.3e} s, max_steps = {max_steps}"
          + (f" (release over {PULSE_WIDTH*1e9:.1f} ns + transit)" if timed else ""),
          flush=True)

    # ── Diagnostics (openPMD, HDF5) ───────────────────────────────────────────
    # Fresh diags: WarpX appends one openPMD file per dump, so a stale series from a prior
    # grid/step count would interleave with the new one. See README → "Fresh diags on rerun".
    if os.path.isdir(DIAG_DIR):
        shutil.rmtree(DIAG_DIR)

    period = max(1, max_steps // N_DIAGS)
    field_diag = picmi.FieldDiagnostic(
        name="fields",
        grid=grid,
        period=period,
        data_list=["phi", "rho", "E"],
        write_dir=DIAG_DIR,
        warpx_format="openpmd",
        warpx_openpmd_backend="h5",
    )
    part_diag = picmi.ParticleDiagnostic(
        name="particles",
        period=period,
        species=[electrons],
        data_list=["position", "momentum", "weighting"],
        write_dir=DIAG_DIR,
        warpx_format="openpmd",
        warpx_openpmd_backend="h5",
    )

    sim = picmi.Simulation(
        solver=solver,
        max_steps=max_steps,
        time_step_size=dt,
        verbose=0,                     # the tqdm bar is the progress display
        particle_shape="linear",
    )
    # ParticleListDistribution supplies the macroparticles explicitly, so this layout is inert.
    sim.add_species(
        electrons,
        layout=picmi.PseudoRandomLayout(n_macroparticles_per_cell=1, grid=grid),
    )
    sim.add_diagnostic(field_diag)
    sim.add_diagnostic(part_diag)

    # ── Time-release injection ────────────────────────────────────────────────
    # Build up the bunch over PULSE_WIDTH: each step inject the macroparticles whose emission
    # time falls in the step window (particle 0 already seeded, so start at index 1).
    # add_particles writes into the live container — the broken RZ accessor is the READ path
    # (`Component x does not exist`), not this WRITE path. ux/uy/uz are proper velocity γβc [m/s].
    if timed:
        state = {"next": 1, "step": 0}
        pc = [None]
        bt = bunch["t"]

        def _inject():
            if pc[0] is None:
                pc[0] = particle_containers.ParticleContainerWrapper("electrons")
            t_hi = (state["step"] + 1) * dt
            j = state["next"]
            k = j
            while k < bt.size and bt[k] < t_hi:
                k += 1
            if k > j:
                pc[0].add_particles(
                    x=bunch["x"][j:k], y=bunch["y"][j:k], z=bunch["z"][j:k],
                    ux=bunch["ux"][j:k], uy=bunch["uy"][j:k], uz=bunch["uz"][j:k],
                    w=bunch["w"][j:k])
                state["next"] = k
            state["step"] += 1

        callbacks.installbeforestep(_inject)

    print(f"\nRunning {max_steps} steps (diag every {period}) …")
    run_step(sim, max_steps, desc="gun")
    if timed:
        print(f"Injected {state['next']}/{bunch['t'].size} macroparticles over the pulse",
              flush=True)
        build_exit_handoff()
    print(f"\nDone. openPMD output → {DIAG_DIR}/{{fields,particles}}/")


if __name__ == "__main__":
    main()
