"""
CESR gun in WarpX (RZ): accelerate the cathode-emitted electrons through the
Poisson–Superfish gun field, with self-consistent space charge.

This is the second stage of the Cornell Linac chain modelled in WarpX. Stage 1
(`cathode/`) is the thermionic cathode operating at the Child–Langmuir
limit; here we take its emitted electrons and track them through the gun's
electrostatic accelerating field — the `CESR_gun.gdf` map scaled to 150 kV by
`build_gun_field.py` and applied as an external (electrode) field on the
particles, while WarpX's electrostatic solver supplies the beam self-field.

Geometry is RZ (cylindrical), matching the gun field map's native symmetry.

Run with (from the repo root, with `conda activate CBB`):
    python -c "import gun; gun.run()"               # build field map + sim + plots
    python -c "import gun; gun.run(plots=False)"    # build + sim only
    python -c "import gun; gun.plot()"              # plots only
Direct script invocation (`python gun/gun_sim.py`) does NOT work — this module
imports `pipeline._runner`, which is only on sys.path when launched from the
repo root (either via the facade above or `python -m gun.gun_sim`).

Beam source — see README: the cathode run is a continuous (DC) emitter, so the
weights in its last particle snapshot encode the steady-state population in
transit through the diode (~82 nC), not a bunch charge. We import the emitted
**phase-space distribution** (positions + momenta), remap the 2D (x, z) slab
into RZ by treating |x| as the radius r and smearing the particles uniformly in
azimuth — importance-resampling by r so the revolution supplies its 2πr Jacobian
(a uniform-in-x slab → a uniform-density disc, not a spurious 1/r on-axis cusp) —
and renormalize the total weight to a physical gun bunch charge
`BUNCH_CHARGE` (the CESR gun is pulse-grid gated). The full 82 nC injected as
one instantaneous bunch is unphysical — its radial space-charge field (~50 MV/m)
dwarfs the gun field and blows the beam apart before it accelerates.
"""

import os
import shutil

import numpy as np
import pywarpx
from pywarpx import picmi
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
ZMAX = 0.051765              # field-map Z extent [m]

CATHODE_DIAG = "cathode/diags/particles"
BUNCH_CHARGE = 1.0e-9        # renormalized gun bunch charge [C] = 1 nC, matching the
                             # original LinacSim gpt_master.in total_charge = -1e-9;
                             # raw cathode snapshot is ~82 nC
RNG_SEED = 0

# ── Grid (RZ, single azimuthal mode — the gun field is m = 0) ─────────────────
nr, nz = 96, 384             # divisible by the blocking factor (8)

# ── Diagnostics output directory ──────────────────────────────────────────────
DIAG_DIR = "gun/diags"

# ── Performance knobs (tunable via gun.config(...); see pipeline/run_pipeline.py) ─
# Defaults reproduce the original run exactly; lower them to trade accuracy for speed.
# Runtime ≈ nz² (per-step cost ∝ cells, and dz=ZMAX/nz ⇒ fewer derived steps as nz drops).
REQUIRED_PRECISION = 1e-5            # MLMG Poisson solve relative tolerance
SPACE_CHARGE = True                  # beam self-field (space charge) on/off. False →
                                     # warpx_do_not_deposit: the beam deposits no charge,
                                     # so only the applied gun field acts (no self-repulsion).
MAX_ITERS = None                     # MLMG iteration cap (None → PICMI default)
CFL = 0.4                            # dt = CFL · dz / v_exit
TRANSIT_MARGIN = 1.15                # run length = TRANSIT_MARGIN × gun-transit time
AVG_SPEED_FRAC = 0.6                 # bunch average speed as a fraction of v_exit
MAX_STEPS = 0                        # 0 → auto-derive from CFL/margins; >0 → fixed
N_DIAGS = 40                         # number of openPMD dumps over the run (≥20 keeps
                                     # space_charge.png's near-launch field snapshot)
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
    # Keep `xk` (the masked signed x) alongside the kept arrays so the radial-momentum
    # sign below survives the optional downsample (x[keep] would re-mask the full set).
    xk = x[keep]
    r, z, ux, uy, uz, w = (a[keep] for a in (r, z, ux, uy, uz, w))

    # Optionally downsample (reweighted to preserve total charge) to cap the cost.
    if MAX_PART and r.size > MAX_PART:
        sel = rng.choice(r.size, MAX_PART, replace=False)
        scale_w = r.size / MAX_PART
        xk, r, z, ux, uy, uz, w = (a[sel] for a in (xk, r, z, ux, uy, uz, w))
        w = w * scale_w

    # slab(x) → RZ disc: supply the 2πr revolution Jacobian that the naive r=|x|
    # map omits. A 2D Cartesian slab uniform in x has a flat dN/dr; revolving it
    # with r=|x| and unchanged weight yields areal density n(r) ∝ 1/r — a spurious
    # on-axis charge cusp that gives a radially-flat (nonlinear) self-field and
    # corrupts the σ_r, φ-well, and emittance this stage is meant to deliver.
    # Importance-resample (with replacement) with probability ∝ r·w (charge-correct;
    # ≡ ∝ r for the cathode's uniform weights), so dN/dr → r·dN/dr and the areal
    # density matches the cathode's true radial profile (a flat-top emitting strip
    # → a uniform-density disc). Drawing from the actual particles keeps weights
    # uniform (no weight-variance, so downstream RMS/emittance stay
    # unweighted-valid) and preserves the cathode-edge position–momentum correlations.
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
    print(f"Imported {r.size} macroparticles from cathode (iter {it}); "
          f"renormalized to {BUNCH_CHARGE*1e9:.3f} nC, r ≤ {r.max()*1e3:.2f} mm",
          flush=True)
    # openPMD ux/uy/uz are the dimensionless normalized momenta γβ; PICMI's
    # ParticleListDistribution wants proper velocity u = γβc in m/s, so ×c.
    # (Without this the beam is injected essentially at rest and the cathode's
    # thermal transverse momentum — hence its emittance — is lost; the energy
    # gain is insensitive because the cathode KE ≪ 150 keV gun voltage.)
    return dict(x=xpos, y=ypos, z=zpos, ux=uxn * c, uy=uyn * c, uz=uz * c, w=w)


def main():
    grid = picmi.CylindricalGrid(
        number_of_cells=[nr, nz],
        n_azimuthal_modes=1,
        lower_bound=[0.0, 0.0],
        upper_bound=[RMAX, ZMAX],
        # r=0 must be "none" (axis); the electrode field is applied externally, so
        # the self-field Poisson solve just needs grounded z plates (dirichlet) with
        # a zero-normal-field (neumann) outer radial wall.
        lower_boundary_conditions=["none", "dirichlet"],
        upper_boundary_conditions=["neumann", "dirichlet"],
        lower_boundary_conditions_particles=["none", "absorbing"],
        upper_boundary_conditions_particles=["absorbing", "absorbing"],
        warpx_blocking_factor=8,
    )

    # Electrostatic solver for the beam self-field only.
    solver_kw = dict(grid=grid, method="Multigrid",
                     required_precision=REQUIRED_PRECISION,
                     warpx_self_fields_verbosity=0)   # silence MLMG per-iteration chatter
    if MAX_ITERS:                                     # omit when None → PICMI default
        solver_kw["maximum_iterations"] = MAX_ITERS
    solver = picmi.ElectrostaticSolver(**solver_kw)

    # ── Applied gun field: the scaled CESR_gun.gdf map, read from file ────────
    # Applied directly to particles every step (the electrode field), on top of the
    # self-consistent space-charge field from the Poisson solve. PICMI has no class
    # for a tabulated particle-applied field, so set the raw WarpX inputs.
    pywarpx.particles.E_ext_particle_init_style = "read_from_file"
    pywarpx.particles.read_fields_from_path = GUN_FIELD
    pywarpx.particles.B_ext_particle_init_style = "none"

    bunch = load_cathode_bunch()
    electrons = picmi.Species(
        particle_type="electron",
        name="electrons",
        initial_distribution=picmi.ParticleListDistribution(
            x=bunch["x"], y=bunch["y"], z=bunch["z"],
            ux=bunch["ux"], uy=bunch["uy"], uz=bunch["uz"],
            weight=bunch["w"],
        ),
        warpx_do_not_deposit=not SPACE_CHARGE,   # SPACE_CHARGE=False → no beam self-field
    )

    # ── Time step / duration ──────────────────────────────────────────────────
    # Exit kinetic energy ≈ 150 keV -> γ ≈ 1.29, β ≈ 0.63, v_exit ≈ 1.9e8 m/s.
    gamma = 1.0 + q_e * GUN_VOLTAGE / (m_e * c**2)
    v_exit = c * np.sqrt(1.0 - 1.0 / gamma**2)
    dz = ZMAX / nz
    dt = CFL * dz / v_exit
    # Steps for the bunch to just cross the full gun (average speed ~AVG_SPEED_FRAC·v_exit).
    # We stop as the beam reaches the exit: running longer empties the domain, and
    # the Multigrid self-field solve aborts when there is essentially no charge left.
    # MAX_STEPS (module constant, 0 = auto) overrides the derived value when set.
    max_steps = MAX_STEPS or int(
        TRANSIT_MARGIN * ZMAX / (AVG_SPEED_FRAC * v_exit) / dt)

    print(f"Gun: {GUN_VOLTAGE/1e3:.0f} kV  ->  γ={gamma:.3f}, β={v_exit/c:.3f}, "
          f"v_exit={v_exit:.2e} m/s", flush=True)
    print(f"dt = {dt:.3e} s, max_steps = {max_steps}", flush=True)

    # ── Diagnostics (openPMD, HDF5) ───────────────────────────────────────────
    # Fresh diags: WarpX appends one openPMD file per dump, so re-running with a
    # different grid/step count would otherwise mix old and new iterations (whose
    # diag steps interleave) into one series — the plots then show a fan of
    # overlapping curves. diags are git-ignored and regenerated, so clearing is
    # safe. (Mirrors injector_sim.py / linac_sec1_sim.py.)
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
        verbose=0,                     # silence per-step "STEP N starts" — the tqdm bar is the progress display
        particle_shape="linear",
    )
    # ParticleListDistribution supplies the macroparticles explicitly, so this layout
    # (and its n_macroparticles_per_cell) is inert — the count is the imported list size.
    sim.add_species(
        electrons,
        layout=picmi.PseudoRandomLayout(n_macroparticles_per_cell=1, grid=grid),
    )
    sim.add_diagnostic(field_diag)
    sim.add_diagnostic(part_diag)

    print(f"\nRunning {max_steps} steps (diag every {period}) …")
    run_step(sim, max_steps, desc="gun")
    print(f"\nDone. openPMD output → {DIAG_DIR}/{{fields,particles}}/")


if __name__ == "__main__":
    main()
