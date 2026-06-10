"""
Finite thermionic cathode at the space-charge (Child–Langmuir) limit — WarpX 2D.

This demo models the *electron source* of the Cornell Linac: Adam Bartnik's
"Region 1" is a thermionic cathode held a short distance from a positively biased
grid/anode, operating in the **space-charge-limited (SCL)** regime where the
emitted current saturates at the Child–Langmuir value regardless of how many
electrons the hot cathode could supply (see
`reference/Linac Simulation Documentation/details.md`).

Unlike the canonical 1D Pierce-diode example
(`reference/WarpX Documentation/usage/examples/pierce_diode/`), here the cathode
has a **finite transverse extent**, simulated in 2D (x–z). This lets us see two
things at once:

  1. On axis (center of the cathode, far from the edges) the potential and field
     reproduce the 1D Child–Langmuir laws  φ(z) = V (z/d)^{4/3},
     Ez(z) = -(4V/3d)(z/d)^{1/3}.
  2. Near the cathode *edges* the equipotentials crowd together, enhancing the
     local field — a genuinely 2D effect absent from the planar theory.

Physics setup (electrostatic lab frame):
  - cathode plane at z = 0 held at 0 V, anode plane at z = d held at +V
  - electrons emitted from the cathode patch |x| < R via continuous flux injection
    (PICMI `UniformFluxDistribution`).  We deliberately *over-inject* at 2 × the Child–Langmuir
    current; WarpX's self-consistent self-fields build a virtual cathode that
    reflects the excess, so the transmitted current self-limits to J_CL.

Run with (from the repo root, with `conda activate CBB`):
    python -c "import cathode; cathode.run()"            # sim + plots
    python -c "import cathode; cathode.run(plots=False)" # sim only
Direct script invocation (`python cathode/cathode_diode.py`) does NOT work —
this module imports `pipeline._runner`, which is only on sys.path when launched
from the repo root (either via the facade above or `python -m cathode.cathode_diode`).
"""

import os
import shutil

import numpy as np
from pywarpx import picmi

from pipeline._runner import run_step

c    = picmi.constants.c
m_e  = picmi.constants.m_e
q_e  = picmi.constants.q_e          # = +1.602e-19 C (elementary charge)
ep0  = picmi.constants.ep0
kb   = 1.380649e-23                  # Boltzmann constant [J/K]

# ── Diode parameters (Adam's Region-1 thermionic cathode) — matched to the
#    original LinacSim inputs (reference/.../input_files/cathode_master.in +
#    gpt_master.in): cathode-grid distance l=0.2 mm, pulse voltage Vpulse=60 V,
#    cathode diameter egun_cath_diam=16 mm, cathode temperature egun_cath_T=1425 K.
V_anode   = 60.0         # anode (grid) bias [V] — cathode at 0 V (Vpulse)
gap_d     = 200.0e-6     # cathode→anode gap [m] (l = 0.2 mm)
R_cathode = 8.0e-3       # cathode half-width [m] (16 mm diameter)
T_cathode = 1425.0       # cathode temperature [K] → small thermal emittance

over_inject = 2.0        # inject this multiple of the Child–Langmuir current

# ── Grid (2D x–z) ──────────────────────────────────────────────────────────────
W = 16.0e-3              # transverse half-width of the domain [m] (2× cathode half-width)
nx, nz = 128, 64         # both divisible by the blocking factor (8)

# ── Diagnostics output directory ───────────────────────────────────────────────
DIAG_DIR = "cathode/diags"

# ── Run length ─────────────────────────────────────────────────────────────────
MAX_STEPS = 2000                     # ~4× the gap-fill time → reaches steady state

# ── Performance knobs (tunable via cathode.config(...); see pipeline/run_pipeline.py) ─
# Defaults reproduce the original run exactly; lower them to trade accuracy for speed.
REQUIRED_PRECISION = 1e-5            # MLMG Poisson solve relative tolerance
MAX_ITERS = None                     # MLMG iteration cap (None → PICMI default)
SPACE_CHARGE = True                  # beam self-field (space charge) on/off. KEEP TRUE — unlike the
                                     # downstream stages (where SC is a correction on an applied-field
                                     # beam), space charge is the SOLE current-limiting mechanism of
                                     # THIS stage. False → warpx_do_not_deposit zeros the beam ρ, so
                                     # there is no virtual cathode, no field suppression at the cathode
                                     # surface, and NO Child–Langmuir limiting: the diode then passes
                                     # the full 2×J_CL over-injection unreflected (≈double the physical
                                     # current) and the validation figures (child_langmuir.png,
                                     # current_saturation.png, rho_z_time.png — the latter two read the
                                     # now-zeroed ρ/j) are invalid. So SPACE_CHARGE=False is NOT a
                                     # meaningful cathode operating point — it removes the very effect
                                     # the stage exists to demonstrate (provided only for parity with
                                     # the other stages' toggle / a forces-off sanity check).
PPC = 10                             # macroparticles per cell (PseudoRandomLayout)
CFL = 0.4                            # dt = CFL · dz / v_final
DIAG_PERIOD = None                   # None → dense-early union slice (keeps figs 3,4);
                                     # an int → uniform period for both diagnostics


def main():
    # Child–Langmuir current density for electrons across a planar gap:
    #   J_CL = (4/9) eps0 sqrt(2 e / m_e) V^{3/2} / d^2
    J_CL = (4.0 / 9.0) * ep0 * np.sqrt(2.0 * q_e / m_e) * V_anode**1.5 / gap_d**2
    flux = over_inject * J_CL / q_e      # particle flux [# / m^2 / s]

    # Thermal velocity spread of emitted electrons (sets the intrinsic emittance):
    v_th = np.sqrt(kb * T_cathode / m_e)

    # Characteristic final velocity (cold electron falling through the full bias):
    v_final = np.sqrt(2.0 * q_e * V_anode / m_e)

    print(f"Diode : V = {V_anode:.0f} V, gap d = {gap_d*1e3:.1f} mm, "
          f"cathode 2R = {2*R_cathode*1e3:.1f} mm")
    print(f"Child–Langmuir J_CL = {J_CL:.1f} A/m^2  "
          f"(injecting {over_inject:.0f}× = {over_inject*J_CL:.1f} A/m^2)")
    print(f"v_th = {v_th:.2e} m/s, v_final = {v_final:.2e} m/s")

    grid = picmi.Cartesian2DGrid(
        number_of_cells=[nx, nz],
        lower_bound=[-W, 0.0],
        upper_bound=[ W, gap_d],
        # x walls: insulating (zero normal field); z plates: fixed potentials
        lower_boundary_conditions=["neumann", "dirichlet"],
        upper_boundary_conditions=["neumann", "dirichlet"],
        lower_boundary_conditions_particles=["absorbing", "absorbing"],
        upper_boundary_conditions_particles=["absorbing", "absorbing"],
        warpx_potential_lo_z=0.0,        # cathode
        warpx_potential_hi_z=V_anode,    # anode / grid
        warpx_blocking_factor=8,
    )

    solver_kw = dict(grid=grid, method="Multigrid",
                     required_precision=REQUIRED_PRECISION,
                     warpx_self_fields_verbosity=0)   # silence MLMG per-iteration chatter
    if MAX_ITERS:                                     # omit when None → PICMI default
        solver_kw["maximum_iterations"] = MAX_ITERS
    solver = picmi.ElectrostaticSolver(**solver_kw)

    # ── Cathode emission (continuous space-charge-limited flux) ─────────────────
    dz = gap_d / nz
    emission = picmi.UniformFluxDistribution(
        flux=flux,
        flux_normal_axis="z",
        surface_flux_position=0.0,                 # emit from the z = 0 plane
        flux_direction=+1,                          # into the gap (+z)
        lower_bound=[-R_cathode, None, None],       # finite cathode patch in x …
        upper_bound=[ R_cathode, None, None],       # … emission is zero outside |x|<R
        rms_velocity=[v_th, v_th, v_th],            # thermal spread (the y component is
                                                    # inert in this 2D x–z run, but the gun
                                                    # reuses cathode uy as the RZ azimuthal
                                                    # thermal momentum — keep it)
        directed_velocity=[0.0, 0.0, 0.0],          # emitted ~at rest, field-accelerated
        gaussian_flux_momentum_distribution=True,   # half-Maxwellian normal to surface
    )
    if not SPACE_CHARGE:
        print("WARNING: cathode SPACE_CHARGE=False — beam self-field deposition is OFF, so the "
              "space-charge-limited (Child–Langmuir) mechanism is disabled. The diode will pass the "
              "full 2×J_CL over-injection unlimited (~2× the physical current) and the validation "
              "figures (child_langmuir / current_saturation / rho_z_time) are NOT valid. This is a "
              "forces-off diagnostic only, not a meaningful cathode operating point.", flush=True)
    electrons = picmi.Species(
        particle_type="electron",
        name="electrons",
        initial_distribution=emission,
        warpx_do_not_deposit=not SPACE_CHARGE,   # SPACE_CHARGE=False → no beam self-field
    )

    # ── Time step / duration ────────────────────────────────────────────────────
    dt = CFL * dz / v_final

    # ── Diagnostics (openPMD) ───────────────────────────────────────────────────
    # Fresh diags: the h5 backend appends one file per dump, so re-running with a
    # different step count/period would leave stale files that interleave with the
    # new ones and corrupt the plots (a fan of overlapping curves). diags are
    # git-ignored and regenerated, so clearing is safe. (Mirrors the other stages.)
    if os.path.isdir(DIAG_DIR):
        shutil.rmtree(DIAG_DIR)

    # Sample densely through the gap-fill transient (≤ 0.13 ns ≈ step 470, since
    # dt ≈ 2.7e-13 s) and sparsely once the diode reaches steady state.  WarpX's
    # interval syntax unions the slices: every 5 steps to 470, then every 80.
    # DIAG_PERIOD=None keeps the dense-early union slice (figs 3/4 need it); an int
    # override applies one uniform period (the field-diag period must be a string).
    # max(MAX_STEPS, 471) guards the second slice against inverting if a caller sets
    # MAX_STEPS ≤ 470 (the default is 2000, so this is defensive only).
    field_period = (str(DIAG_PERIOD) if DIAG_PERIOD
                    else f"0:470:5, 470:{max(MAX_STEPS, 471)}:80")
    field_diag = picmi.FieldDiagnostic(
        name="fields",
        grid=grid,
        period=field_period,
        data_list=["phi", "rho", "E", "J"],
        write_dir=DIAG_DIR,
        warpx_format="openpmd",
        # HDF5 writes one clean file per iteration; the ADIOS2 BP5 default clobbers
        # files under the rapid successive flushes our dense early sampling produces.
        warpx_openpmd_backend="h5",
    )
    part_diag = picmi.ParticleDiagnostic(
        name="particles",
        period=(DIAG_PERIOD or 200),
        species=[electrons],
        data_list=["position", "momentum", "weighting"],
        write_dir=DIAG_DIR,
        warpx_format="openpmd",
        warpx_openpmd_backend="h5",     # pin h5 like the field diag + downstream stages
    )

    sim = picmi.Simulation(
        solver=solver,
        max_steps=MAX_STEPS,
        time_step_size=dt,
        verbose=0,                     # silence per-step "STEP N starts" — the tqdm bar is the progress display
        particle_shape="linear",
    )
    sim.add_species(
        electrons,
        layout=picmi.PseudoRandomLayout(n_macroparticles_per_cell=PPC, grid=grid),
    )
    sim.add_diagnostic(field_diag)
    sim.add_diagnostic(part_diag)

    print(f"\nRunning {MAX_STEPS} steps  dt = {dt:.3e} s  "
          f"(gap-fill ≈ {int(3*gap_d/v_final/dt)} steps)")
    run_step(sim, MAX_STEPS, desc="cathode")

    print(f"\nDone. openPMD output → {DIAG_DIR}/{{fields,particles}}/")


if __name__ == "__main__":
    main()
