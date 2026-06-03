"""
SPEBS gun, self-consistent: finite thermionic cathode in WarpX at the SPEBS
Pierce-gun operating point, in RZ (cylindrical) geometry to match their gun.

The SPEBS deck (Maxson & Andorf, "Modelling SPEBS + Spectrometer Design") models
the gun field as an *analytic vacuum* Pierce diode and injects ~1 A. Here we
instead solve the cathode region **self-consistently** — the beam's own space
charge depresses the field near the cathode — so we can quantify how far the real
on-axis field departs from the vacuum value at their operating point.

Parameters (from the SPEBS "Gun" / "Gun Fieldmap" slides):
    V = 50 kV, cathode→anode gap d = 23.1 mm, 8 mm-dia (R = 4 mm) flat cathode,
    emission ~1 A  →  J ≈ 2.0e4 A/m²  ≈ 0.41 × the Child–Langmuir limit.

Usage (run twice, once per emission level):
    python warpx_cathode/spebs_gun.py 0.41 spebs   # SPEBS 1 A operating point
    python warpx_cathode/spebs_gun.py 1.50 cl      # over-injected → CL-limited
"""

import sys

import numpy as np
from pywarpx import picmi

c    = picmi.constants.c
m_e  = picmi.constants.m_e
q_e  = picmi.constants.q_e
ep0  = picmi.constants.ep0
kb   = 1.380649e-23

# ── CLI: emission level (× J_CL) and an output tag ──────────────────────────────
inject_factor = float(sys.argv[1]) if len(sys.argv) > 1 else 0.41
tag           = sys.argv[2] if len(sys.argv) > 2 else "spebs"

# ── SPEBS gun parameters ────────────────────────────────────────────────────────
V_anode   = 50.0e3       # 50 kV bias
gap_d     = 23.1e-3      # cathode→anode gap [m]
R_cathode = 4.0e-3       # 8 mm-diameter cathode → 4 mm radius
T_cathode = 1500.0

J_CL = (4.0 / 9.0) * ep0 * np.sqrt(2.0 * q_e / m_e) * V_anode**1.5 / gap_d**2
flux = inject_factor * J_CL / q_e
v_th = np.sqrt(kb * T_cathode / m_e)
# relativistic final velocity at 50 keV (gamma ~ 1.10)
gamma_f = 1.0 + q_e * V_anode / (m_e * c**2)
v_final = c * np.sqrt(1.0 - 1.0 / gamma_f**2)

print(f"[{tag}] V = {V_anode/1e3:.0f} kV, gap = {gap_d*1e3:.1f} mm, "
      f"R = {R_cathode*1e3:.0f} mm")
print(f"[{tag}] J_CL = {J_CL:.3e} A/m²,  injecting {inject_factor:.2f}× "
      f"= {inject_factor*J_CL:.3e} A/m²  ({inject_factor*J_CL*np.pi*R_cathode**2:.3f} A)")

# ── RZ grid (cylindrical, axisymmetric) ─────────────────────────────────────────
rmax   = 8.0e-3
nr, nz = 64, 160         # divisible by the blocking factor (8)

grid = picmi.CylindricalGrid(
    number_of_cells=[nr, nz],
    lower_bound=[0.0, 0.0],
    upper_bound=[rmax, gap_d],
    n_azimuthal_modes=1,
    lower_boundary_conditions=["none", "dirichlet"],      # r=0 axis, z=0 cathode
    upper_boundary_conditions=["neumann", "dirichlet"],   # open wall, z=d anode
    lower_boundary_conditions_particles=["none", "absorbing"],
    upper_boundary_conditions_particles=["absorbing", "absorbing"],
    warpx_potential_lo_z=0.0,        # cathode
    warpx_potential_hi_z=V_anode,    # anode
    warpx_blocking_factor=8,
)

solver = picmi.ElectrostaticSolver(grid=grid, method="Multigrid",
                                   required_precision=1e-5)

# ── Cathode emission (finite disk r < R at z = 0) ───────────────────────────────
dz = gap_d / nz
emission = picmi.UniformFluxDistribution(
    flux=flux,
    flux_normal_axis="z",
    surface_flux_position=0.0,
    flux_direction=+1,
    lower_bound=[0.0,        None, None],
    upper_bound=[R_cathode,  None, None],     # finite cathode disk
    rms_velocity=[v_th, v_th, v_th],
    directed_velocity=[0.0, 0.0, 0.0],
    gaussian_flux_momentum_distribution=True,
)
electrons = picmi.Species(particle_type="electron", name="electrons",
                          initial_distribution=emission)

dt        = 0.4 * dz / v_final
max_steps = 3000                     # ~3× the gap-fill time

field_diag = picmi.FieldDiagnostic(
    name="fields", grid=grid, period=500,
    data_list=["phi", "rho", "E", "J"],
    write_dir=f"warpx_cathode/diags_{tag}", warpx_format="openpmd",
)

sim = picmi.Simulation(solver=solver, max_steps=max_steps, time_step_size=dt,
                       verbose=1, particle_shape="linear")
sim.add_species(electrons,
                layout=picmi.PseudoRandomLayout(n_macroparticles_per_cell=8, grid=grid))
sim.add_diagnostic(field_diag)

print(f"[{tag}] running {max_steps} steps  dt = {dt:.3e} s")
sim.step(max_steps)
print(f"[{tag}] done → warpx_cathode/diags_{tag}/")
