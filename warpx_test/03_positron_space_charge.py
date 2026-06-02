"""
Gaussian positron bunch at 5 MeV — strong space charge (1 nC).

Identical setup to 02_positron_bunch.py except the bunch charge is 1 nC
(1000× larger).  At 5 MeV (gamma ≈ 11) the space-charge suppression factor
1/gamma² ≈ 1/120 is small enough that the Coulomb self-repulsion drives visible
transverse beam expansion within 200 time steps, in stark contrast to the
stable 1 pC reference beam.

Physics parameters match those of a Cornell ERL-style photoinjector:
  - MeV-class energy where space charge dominates emittance growth
  - mm-scale transverse size, nC-level charge
  - Relativistic Poisson solver accounts for the Lorentz-contracted fields

Run with:
    conda run -n CBB python warpx_test/03_positron_space_charge.py
"""

import numpy as np
from pywarpx import picmi

c   = picmi.constants.c
m_e = picmi.constants.m_e
q_e = picmi.constants.q_e

# ── Beam parameters (identical to 02 except charge) ──────────────────────────
KE_MeV     = 5.0
rest_MeV   = 0.5109989
gamma      = 1.0 + KE_MeV / rest_MeV
beta       = np.sqrt(1.0 - 1.0 / gamma**2)
uz0        = gamma * beta * c

total_charge_C = 1e-9              # 1 nC — strong space charge
n_real         = total_charge_C / q_e

sigma_r = 0.5e-3   # 0.5 mm transverse (4× smaller → ~16× stronger space charge)
sigma_z = 1.0e-3

n_macroparticles = 20000

print(f"Beam : KE = {KE_MeV} MeV, gamma = {gamma:.2f}")
print(f"Bunch: Q = {total_charge_C*1e9:.1f} nC  "
      f"sigma_r = {sigma_r*1e3:.1f} mm  sigma_z = {sigma_z*1e3:.1f} mm")

# ── Simulation parameters ─────────────────────────────────────────────────────
dz = 0.25e-3
dt = 0.99 * dz / c
max_steps = 500

# ── Grid ──────────────────────────────────────────────────────────────────────
grid = picmi.Cartesian3DGrid(
    number_of_cells=[32, 32, 32],
    lower_bound=[-4e-3, -4e-3, -4e-3],
    upper_bound=[ 4e-3,  4e-3,  4e-3],
    lower_boundary_conditions=["open", "open", "open"],
    upper_boundary_conditions=["open", "open", "open"],
    lower_boundary_conditions_particles=["absorbing", "absorbing", "absorbing"],
    upper_boundary_conditions_particles=["absorbing", "absorbing", "absorbing"],
    moving_window_velocity=[0.0, 0.0, beta * c],
    warpx_blocking_factor=8,
)

solver = picmi.ElectrostaticSolver(grid=grid, method="FFT", warpx_relativistic=True)

# ── Species ───────────────────────────────────────────────────────────────────
dist = picmi.GaussianBunchDistribution(
    n_physical_particles=n_real,
    rms_bunch_size=[sigma_r, sigma_r, sigma_z],
    centroid_position=[0.0, 0.0, 0.0],
    centroid_velocity=[0.0, 0.0, uz0],
    rms_velocity=[0.0, 0.0, 0.0],  # zero divergence: expansion = space charge only
)
positrons = picmi.Species(
    particle_type="positron",
    name="positrons",
    initial_distribution=dist,
)

# ── Diagnostics ───────────────────────────────────────────────────────────────
beam_diag = picmi.ReducedDiagnostic(
    diag_type="BeamRelevant",
    name="beam_stats_sc",
    period=1,
    path="warpx_test/diags_space_charge/",
    species=positrons,
)
part_diag = picmi.ParticleDiagnostic(
    name="part",
    period=50,
    species=[positrons],
    data_list=["position", "momentum", "weighting"],
    write_dir="warpx_test/diags_space_charge",
    warpx_format="plotfile",
)

# ── Simulation ────────────────────────────────────────────────────────────────
sim = picmi.Simulation(
    solver=solver,
    max_steps=max_steps,
    time_step_size=dt,
    verbose=1,
    particle_shape="linear",
    warpx_particle_pusher_algo="vay",
)
sim.add_species(positrons, layout=picmi.PseudoRandomLayout(n_macroparticles=n_macroparticles))
sim.add_diagnostic(beam_diag)
sim.add_diagnostic(part_diag)

print(f"\nRunning {max_steps} steps with 1 nC space charge  dt = {dt:.3e} s")
sim.step(max_steps)

print("\nDone. Diagnostics → warpx_test/diags_space_charge/beam_stats_sc.txt")
