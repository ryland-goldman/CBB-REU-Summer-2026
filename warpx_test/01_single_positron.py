"""
Single-positron trajectory in free space using WarpX/PICMI.

A single 100 MeV positron is injected at the origin and propagates
in the +z direction. This is the simplest sanity check for WarpX:
confirm the particle travels at the correct relativistic velocity and
that diagnostic output is written correctly.

Run with:
    conda run -n CBB python warpx_test/01_single_positron.py
"""

import numpy as np
from pywarpx import picmi

# ── Physical constants ────────────────────────────────────────────────────────
c = picmi.constants.c      # 2.998e8 m/s
m_e = picmi.constants.m_e  # 9.109e-31 kg
q_e = picmi.constants.q_e  # 1.602e-19 C

# ── Beam parameters ───────────────────────────────────────────────────────────
KE_MeV = 100.0                           # kinetic energy [MeV]
rest_energy_MeV = 0.5109989             # m_e c^2 in MeV
gamma = 1.0 + KE_MeV / rest_energy_MeV  # Lorentz factor ≈ 196.8
beta = np.sqrt(1.0 - 1.0 / gamma**2)    # ≈ 1
uz0 = gamma * beta * c                   # normalized momentum in z [m/s]

print(f"Positron: KE = {KE_MeV} MeV, gamma = {gamma:.2f}, uz = {uz0:.3e} m/s")

# ── Simulation parameters ─────────────────────────────────────────────────────
max_steps = 50
# CFL time step: dt = 0.99 * dz / c (dominated by longitudinal cell size)
dz = 1.92 / 96    # 2 cm cell size in z
dt = 0.99 * dz / c  # ≈ 6.6e-11 s
# At nearly c, the positron travels ~dz per step → ~1 m in 50 steps

# ── Grid ──────────────────────────────────────────────────────────────────────
# Thin (10 × 10 cells) in x,y; long (100 cells) in z.
# Positron starts at z = 0.1 m, travels to z ≈ 1.1 m before exiting.
grid = picmi.Cartesian3DGrid(
    number_of_cells=[8, 8, 96],
    lower_bound=[-0.1, -0.1, 0.0],
    upper_bound=[ 0.1,  0.1, 1.92],
    lower_boundary_conditions=["open", "open", "open"],
    upper_boundary_conditions=["open", "open", "open"],
    lower_boundary_conditions_particles=["absorbing", "absorbing", "absorbing"],
    upper_boundary_conditions_particles=["absorbing", "absorbing", "absorbing"],
    warpx_blocking_factor=8,
)

# ── Field solver ──────────────────────────────────────────────────────────────
# Standard Yee FDTD; single particle so fields are trivial.
solver = picmi.ElectromagneticSolver(grid=grid, cfl=0.99)

# ── Particle species: single positron ────────────────────────────────────────
initial_dist = picmi.ParticleListDistribution(
    x=0.0, y=0.0, z=0.1,   # start 10 cm from lower z boundary
    ux=0.0, uy=0.0, uz=uz0,
    weight=1.0,             # one real positron per macroparticle
)
positrons = picmi.Species(
    particle_type="positron",
    name="positrons",
    initial_distribution=initial_dist,
)

# ── Diagnostics ───────────────────────────────────────────────────────────────
# Save full particle state (position + momentum) at every time step.
part_diag = picmi.ParticleDiagnostic(
    name="part",
    period=1,
    species=[positrons],
    data_list=["position", "momentum", "weighting"],
    write_dir="warpx_test/diags_single",
    warpx_format="plotfile",
)

# ── Simulation ────────────────────────────────────────────────────────────────
sim = picmi.Simulation(
    solver=solver,
    max_steps=max_steps,
    time_step_size=dt,
    verbose=1,
    particle_shape="linear",
)

sim.add_species(positrons, layout=picmi.PseudoRandomLayout(n_macroparticles=1))
sim.add_diagnostic(part_diag)

print(f"Running {max_steps} steps (dt = {dt:.3e} s)")
print(f"Expected final z ≈ {0.1 + beta * c * max_steps * dt:.3f} m")
sim.step(max_steps)

print("Done. Particle diagnostics written to warpx_test/diags_single/")
print("Use yt or openpmd_viewer to inspect the plotfile output.")
