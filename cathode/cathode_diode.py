"""
Finite thermionic cathode at the space-charge (Child–Langmuir) limit — WarpX 2D.

Over-injects at 2×J_CL and lets WarpX's self-consistent fields build a virtual
cathode that self-limits the transmitted current to J_CL. See cathode/README.md
for physics, parameters, and gotchas.

Run from the repo root (the module imports pipeline._runner, only on sys.path
from there): python -c "import cathode; cathode.run()".
"""

import os
import shutil

import numpy as np

from pipeline._runner import run_step
from pipeline.constants import E_CHARGE as q_e, M_E as m_e   # = +1.602e-19 C, electron mass
from pipeline.emission import child_langmuir_current_density, thermal_velocity_sigma

V_anode   = 30.0         # peak anode (grid) bias [V] — cathode at 0 V (= Voff+Vpulse, see README)
gap_d     = 200.0e-6     # cathode→anode gap [m]
R_cathode = 8.0e-3       # cathode half-width [m]
T_cathode = 1425.0       # cathode temperature [K]

over_inject = 2.0        # inject this multiple of the Child–Langmuir current

W = 16.0e-3              # transverse half-width of the domain [m]
nx, nz = 128, 64         # both divisible by the blocking factor (8)

DIAG_DIR = "cathode/diags"

MAX_STEPS = 2000

REQUIRED_PRECISION = 1e-5            # MLMG Poisson solve relative tolerance
MAX_ITERS = None                     # MLMG iteration cap (None → PICMI default)
SPACE_CHARGE = True                  # KEEP TRUE — space charge is the SOLE current-limiting
                                     # mechanism here; False disables Child–Langmuir (see README).
PPC = 10                             # macroparticles per cell (PseudoRandomLayout)
CFL = 0.4                            # dt = CFL · dz / v_final
DIAG_PERIOD = None                   # None → dense-early union slice (keeps figs 3,4);
                                     # an int → uniform period for both diagnostics


def main():
    from pywarpx import picmi          # lazy: keeps the module pywarpx-free to import (plot reuse)

    # Child–Langmuir current density (electrons, planar gap); shared with plot_cathode
    J_CL = float(child_langmuir_current_density(V_anode, gap_d))
    flux = over_inject * J_CL / q_e      # particle flux [# / m^2 / s]

    v_th = thermal_velocity_sigma(T_cathode)        # thermal velocity spread
    v_final = np.sqrt(2.0 * q_e * V_anode / m_e)    # cold final velocity through full bias

    print(f"Diode : V = {V_anode:.0f} V, gap d = {gap_d*1e3:.1f} mm, "
          f"cathode 2R = {2*R_cathode*1e3:.1f} mm")
    print(f"Child–Langmuir J_CL = {J_CL:.1f} A/m^2  "
          f"(injecting {over_inject:.0f}× = {over_inject*J_CL:.1f} A/m^2)")
    print(f"v_th = {v_th:.2e} m/s, v_final = {v_final:.2e} m/s")

    grid = picmi.Cartesian2DGrid(
        number_of_cells=[nx, nz],
        lower_bound=[-W, 0.0],
        upper_bound=[ W, gap_d],
        # x walls neumann (insulating); z plates dirichlet (fixed potential)
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
        lower_bound=[-R_cathode, None, None],       # finite cathode patch: emission zero outside |x|<R
        upper_bound=[ R_cathode, None, None],
        rms_velocity=[v_th, v_th, v_th],            # y component inert in 2D, but gun reuses uy as
                                                    # the RZ azimuthal thermal momentum — keep it
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

    dt = CFL * dz / v_final

    # Fresh diags: the h5 backend appends one file per dump, so stale files from a
    # prior run with a different step count/period would interleave and corrupt the plots.
    if os.path.isdir(DIAG_DIR):
        shutil.rmtree(DIAG_DIR)

    # Union slice: dense through the gap-fill transient (every 5 steps to 470), sparse after.
    # The field-diag period must be a string. max(MAX_STEPS, 471) guards the second slice
    # against inverting if a caller sets MAX_STEPS ≤ 470 (defensive; default is 2000).
    field_period = (str(DIAG_PERIOD) if DIAG_PERIOD
                    else f"0:470:5, 470:{max(MAX_STEPS, 471)}:80")
    field_diag = picmi.FieldDiagnostic(
        name="fields",
        grid=grid,
        period=field_period,
        data_list=["phi", "rho", "E", "J"],
        write_dir=DIAG_DIR,
        warpx_format="openpmd",
        # h5: one clean file per iteration; ADIOS2 BP5 default clobbers files under the
        # rapid successive flushes the dense early sampling produces.
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
