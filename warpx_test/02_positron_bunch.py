"""
Gaussian positron bunch at 5 MeV — negligible space charge (1 pC).

A 1 pC, 5 MeV Gaussian positron bunch propagates through a moving-window
grid using the relativistic electrostatic Poisson solver.  With only 1 pC the
self-field is negligible, so the beam stays nearly frozen (the initial emittance
is also set to zero to isolate the effect).  This is the baseline for comparing
against 03_positron_space_charge.py (1 nC).

Why 5 MeV?
  At 5 MeV (gamma ≈ 11) the relativistic space-charge suppression factor is
  1/gamma² ≈ 1/120, small enough that 1 nC of charge causes visible beam blow-up
  within ~200 time steps — unlike the 1 GeV case where 1/gamma² ≈ 1/3.8M.

Why electrostatic solver?
  The Yee FDTD solver suffers from the Numerical Cherenkov Instability (NCI)
  for relativistic beams (gamma >> 1), causing unphysical energy gain.  The
  FFT-based (Integrated Green Function) Poisson solver is NCI-free.

Run with:
    conda run -n CBB python warpx_test/02_positron_bunch.py
"""

import numpy as np
from pywarpx import picmi

c   = picmi.constants.c
m_e = picmi.constants.m_e
q_e = picmi.constants.q_e

# ── Beam parameters ───────────────────────────────────────────────────────────
KE_MeV     = 5.0
rest_MeV   = 0.5109989
gamma      = 1.0 + KE_MeV / rest_MeV  # ≈ 10.78
beta       = np.sqrt(1.0 - 1.0 / gamma**2)
uz0        = gamma * beta * c           # normalized z-momentum

total_charge_C = 1e-12                 # 1 pC → negligible space charge
n_real         = total_charge_C / q_e

sigma_r = 0.5e-3   # 0.5 mm transverse RMS
sigma_z = 1.0e-3   # 1 mm longitudinal RMS

n_macroparticles = 20000

print(f"Beam : KE = {KE_MeV} MeV, gamma = {gamma:.2f}, beta = {beta:.5f}")
print(f"Bunch: Q = {total_charge_C*1e12:.1f} pC  "
      f"sigma_r = {sigma_r*1e3:.1f} mm  sigma_z = {sigma_z*1e3:.1f} mm")

# ── Simulation parameters ─────────────────────────────────────────────────────
dx = dy = dz = 0.25e-3             # 0.25 mm cells → sigma_r / dx = 2
dt = 0.99 * dz / c                 # ≈ 0.825 ps
max_steps = 500

# ── Grid with moving window at beam velocity ──────────────────────────────────
# ±8 mm transverse (8 sigma), ±10 mm longitudinal (5 sigma).
# Window follows beam at beta*c (not c) since beam is mildly relativistic.
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

# FFT-based relativistic Poisson — stable for any gamma, handles open BCs
solver = picmi.ElectrostaticSolver(grid=grid, method="FFT", warpx_relativistic=True)

# ── Species ───────────────────────────────────────────────────────────────────
dist = picmi.GaussianBunchDistribution(
    n_physical_particles=n_real,
    rms_bunch_size=[sigma_r, sigma_r, sigma_z],
    centroid_position=[0.0, 0.0, 0.0],
    centroid_velocity=[0.0, 0.0, uz0],
    rms_velocity=[0.0, 0.0, 0.0],   # zero divergence: expansion = space charge only
)
positrons = picmi.Species(
    particle_type="positron",
    name="positrons",
    initial_distribution=dist,
)

# ── Diagnostics ───────────────────────────────────────────────────────────────
beam_diag = picmi.ReducedDiagnostic(
    diag_type="BeamRelevant",
    name="beam_stats",
    period=1,
    path="warpx_test/diags_bunch/",
    species=positrons,
)
part_diag = picmi.ParticleDiagnostic(
    name="part",
    period=50,
    species=[positrons],
    data_list=["position", "momentum", "weighting"],
    write_dir="warpx_test/diags_bunch",
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

print(f"\nRunning {max_steps} steps  dt = {dt:.3e} s  "
      f"beam travel ≈ {beta*c*max_steps*dt*100:.1f} cm")
sim.step(max_steps)

print("\nDone. Diagnostics → warpx_test/diags_bunch/beam_stats.txt")
