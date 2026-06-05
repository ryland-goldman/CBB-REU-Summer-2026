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

Run with:
    conda run -n CBB python cathode/cathode_diode.py
"""

import numpy as np
from pywarpx import picmi

from pipeline._runner import run_step

c    = picmi.constants.c
m_e  = picmi.constants.m_e
q_e  = picmi.constants.q_e          # = +1.602e-19 C (elementary charge)
ep0  = picmi.constants.ep0
kb   = 1.380649e-23                  # Boltzmann constant [J/K]

# ── Diode parameters (Adam's Region-1 thermionic cathode, scaled for a 2D demo) ─
V_anode   = 50.0         # anode (grid) bias [V] — cathode at 0 V
gap_d     = 100.0e-6     # cathode→anode gap [m]
R_cathode = 6.0e-3       # cathode half-width [m]
T_cathode = 1200.0       # cathode temperature [K] → small thermal emittance

over_inject = 2.0        # inject this multiple of the Child–Langmuir current

# ── Grid (2D x–z) ──────────────────────────────────────────────────────────────
W = 12.0e-3              # transverse half-width of the domain [m]
nx, nz = 128, 64         # both divisible by the blocking factor (8)

# ── Diagnostics output directory ───────────────────────────────────────────────
DIAG_DIR = "cathode/diags"

# ── Run length ─────────────────────────────────────────────────────────────────
MAX_STEPS = 2000                     # ~4× the gap-fill time → reaches steady state


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

    solver = picmi.ElectrostaticSolver(
        grid=grid, method="Multigrid", required_precision=1e-5,
        warpx_self_fields_verbosity=0,                # silence MLMG per-iteration chatter
    )

    # ── Cathode emission (continuous space-charge-limited flux) ─────────────────
    dz = gap_d / nz
    emission = picmi.UniformFluxDistribution(
        flux=flux,
        flux_normal_axis="z",
        surface_flux_position=0.0,                 # emit from the z = 0 plane
        flux_direction=+1,                          # into the gap (+z)
        lower_bound=[-R_cathode, None, None],       # finite cathode patch in x …
        upper_bound=[ R_cathode, None, None],       # … emission is zero outside |x|<R
        rms_velocity=[v_th, v_th, v_th],            # thermal spread
        directed_velocity=[0.0, 0.0, 0.0],          # emitted ~at rest, field-accelerated
        gaussian_flux_momentum_distribution=True,   # half-Maxwellian normal to surface
    )
    electrons = picmi.Species(
        particle_type="electron",
        name="electrons",
        initial_distribution=emission,
    )

    # ── Time step / duration ────────────────────────────────────────────────────
    dt = 0.4 * dz / v_final

    # ── Diagnostics (openPMD) ───────────────────────────────────────────────────
    # Sample densely through the gap-fill transient (≤ 0.07 ns ≈ step 470, since
    # dt ≈ 1.49e-13 s) and sparsely once the diode reaches steady state.  WarpX's
    # interval syntax unions the slices: every 5 steps to 470, then every 80.
    field_diag = picmi.FieldDiagnostic(
        name="fields",
        grid=grid,
        period="0:470:5, 470:2000:80",
        data_list=["phi", "rho", "E", "J"],
        write_dir=DIAG_DIR,
        warpx_format="openpmd",
        # HDF5 writes one clean file per iteration; the ADIOS2 BP5 default clobbers
        # files under the rapid successive flushes our dense early sampling produces.
        warpx_openpmd_backend="h5",
    )
    part_diag = picmi.ParticleDiagnostic(
        name="particles",
        period=200,
        species=[electrons],
        data_list=["position", "momentum", "weighting"],
        write_dir=DIAG_DIR,
        warpx_format="openpmd",
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
        layout=picmi.PseudoRandomLayout(n_macroparticles_per_cell=10, grid=grid),
    )
    sim.add_diagnostic(field_diag)
    sim.add_diagnostic(part_diag)

    print(f"\nRunning {MAX_STEPS} steps  dt = {dt:.3e} s  "
          f"(gap-fill ≈ {int(3*gap_d/v_final/dt)} steps)")
    run_step(sim, MAX_STEPS, desc="cathode")

    print(f"\nDone. openPMD output → {DIAG_DIR}/{{fields,particles}}/")


if __name__ == "__main__":
    main()
