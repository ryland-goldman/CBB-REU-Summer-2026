"""
CESR gun in WarpX (RZ): accelerate the cathode-emitted electrons through the
Poisson–Superfish gun field, with self-consistent space charge.

This is the second stage of the Cornell Linac chain modelled in WarpX. Stage 1
(`warpx_cathode/`) is the thermionic cathode operating at the Child–Langmuir
limit; here we take its emitted electrons and track them through the gun's
electrostatic accelerating field — the `CESR_gun.gdf` map scaled to 150 kV by
`build_gun_field.py` and applied as an external (electrode) field on the
particles, while WarpX's electrostatic solver supplies the beam self-field.

Geometry is RZ (cylindrical), matching the gun field map's native symmetry.

Pipeline:
    conda run -n CBB python warpx_gun/build_gun_field.py   # writes gun_field/gun_E.h5
    conda run -n CBB python warpx_gun/gun_sim.py           # this script
    conda run -n CBB python warpx_gun/plot_gun.py

Beam source — see README: the cathode run is a continuous (DC) emitter, so the
weights in its last particle snapshot encode the steady-state population in
transit through the diode (~102 nC), not a bunch charge. We import the emitted
**phase-space distribution** (positions + momenta), remap the 2D (x, z) slab
into RZ by treating |x| as the radius r and smearing the particles uniformly in
azimuth, and renormalize the total weight to a physical gun bunch charge
`BUNCH_CHARGE` (the CESR gun is pulse-grid gated). The full 102 nC injected as
one instantaneous bunch is unphysical — its radial space-charge field (~50 MV/m)
dwarfs the gun field and blows the beam apart before it accelerates.
"""

import numpy as np
import pywarpx
from pywarpx import picmi
from openpmd_viewer import OpenPMDTimeSeries

c = picmi.constants.c
m_e = picmi.constants.m_e
q_e = picmi.constants.q_e
ep0 = picmi.constants.ep0

# ── Gun / field-map parameters (must match build_gun_field.py) ────────────────
GUN_FIELD = "warpx_gun/gun_field/gun_E.h5"
GUN_VOLTAGE = 150.0e3        # [V]
RMAX = 0.015                 # field-map R extent [m]
ZMAX = 0.051765              # field-map Z extent [m]

CATHODE_DIAG = "warpx_cathode/diags/particles"
BUNCH_CHARGE = 0.1e-9        # renormalized gun bunch charge [C] (matches the 0.1 nC
                             # IMPACT-T gun model); raw cathode snapshot is ~102 nC
RNG_SEED = 0

# ── Grid (RZ, single azimuthal mode — the gun field is m = 0) ─────────────────
nr, nz = 96, 384             # divisible by the blocking factor (8)

grid = picmi.CylindricalGrid(
    number_of_cells=[nr, nz],
    n_azimuthal_modes=1,
    lower_bound=[0.0, 0.0],
    upper_bound=[RMAX, ZMAX],
    # r=0 must be "none" (axis); the electrode field is applied externally, so
    # the self-field Poisson solve just needs grounded outer walls.
    lower_boundary_conditions=["none", "dirichlet"],
    upper_boundary_conditions=["neumann", "dirichlet"],
    lower_boundary_conditions_particles=["none", "absorbing"],
    upper_boundary_conditions_particles=["absorbing", "absorbing"],
    warpx_blocking_factor=8,
)

# Electrostatic solver for the beam self-field only.
solver = picmi.ElectrostaticSolver(
    grid=grid, method="Multigrid", required_precision=1e-5,
)

# ── Applied gun field: the scaled CESR_gun.gdf map, read from file ────────────
# Applied directly to particles every step (the electrode field), on top of the
# self-consistent space-charge field from the Poisson solve. PICMI has no class
# for a tabulated particle-applied field, so set the raw WarpX inputs.
pywarpx.particles.E_ext_particle_init_style = "read_from_file"
pywarpx.particles.read_fields_from_path = GUN_FIELD
pywarpx.particles.B_ext_particle_init_style = "none"


def load_cathode_bunch():
    """Import the last cathode snapshot and remap the (x, z) slab into RZ.

    Returns dict of x, y, z, ux, uy, uz, w arrays for ParticleListDistribution.
    """
    ts = OpenPMDTimeSeries(CATHODE_DIAG)
    it = ts.iterations[-1]
    x, z, ux, uy, uz, w = ts.get_particle(
        ["x", "z", "ux", "uy", "uz", "w"], species="electrons", iteration=it,
    )

    rng = np.random.default_rng(RNG_SEED)
    r = np.abs(x)
    keep = r < RMAX
    r, z, ux, uy, uz, w = (a[keep] for a in (r, z, ux, uy, uz, w))

    theta = rng.uniform(0.0, 2.0 * np.pi, size=r.size)
    ct, st = np.cos(theta), np.sin(theta)

    # slab x -> radius; transverse momentum: radial = ux·sign(x), azimuthal = uy
    ur = ux * np.sign(np.where(x[keep] == 0.0, 1.0, x[keep]))
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


bunch = load_cathode_bunch()
electrons = picmi.Species(
    particle_type="electron",
    name="electrons",
    initial_distribution=picmi.ParticleListDistribution(
        x=bunch["x"], y=bunch["y"], z=bunch["z"],
        ux=bunch["ux"], uy=bunch["uy"], uz=bunch["uz"],
        weight=bunch["w"],
    ),
)

# ── Time step / duration ──────────────────────────────────────────────────────
# Exit kinetic energy ≈ 150 keV -> γ ≈ 1.29, β ≈ 0.63, v_exit ≈ 1.9e8 m/s.
gamma = 1.0 + q_e * GUN_VOLTAGE / (m_e * c**2)
v_exit = c * np.sqrt(1.0 - 1.0 / gamma**2)
dz = ZMAX / nz
dt = 0.4 * dz / v_exit
# Steps for the bunch to just cross the full gun (average speed ~0.6·v_exit).
# We stop as the beam reaches the exit: running longer empties the domain, and
# the Multigrid self-field solve aborts when there is essentially no charge left.
max_steps = int(1.15 * ZMAX / (0.6 * v_exit) / dt)

print(f"Gun: {GUN_VOLTAGE/1e3:.0f} kV  ->  γ={gamma:.3f}, β={v_exit/c:.3f}, "
      f"v_exit={v_exit:.2e} m/s", flush=True)
print(f"dt = {dt:.3e} s, max_steps = {max_steps}", flush=True)

# ── Diagnostics (openPMD, HDF5) ───────────────────────────────────────────────
# Sample finely enough to resolve the near-cathode launch and the acceleration.
period = max(1, max_steps // 40)
field_diag = picmi.FieldDiagnostic(
    name="fields",
    grid=grid,
    period=period,
    data_list=["phi", "rho", "E"],
    write_dir="warpx_gun/diags",
    warpx_format="openpmd",
    warpx_openpmd_backend="h5",
)
part_diag = picmi.ParticleDiagnostic(
    name="particles",
    period=period,
    species=[electrons],
    data_list=["position", "momentum", "weighting"],
    write_dir="warpx_gun/diags",
    warpx_format="openpmd",
    warpx_openpmd_backend="h5",
)

# ── Simulation ────────────────────────────────────────────────────────────────
sim = picmi.Simulation(
    solver=solver,
    max_steps=max_steps,
    time_step_size=dt,
    verbose=1,
    particle_shape="linear",
)
sim.add_species(
    electrons,
    layout=picmi.PseudoRandomLayout(n_macroparticles_per_cell=1, grid=grid),
)
sim.add_diagnostic(field_diag)
sim.add_diagnostic(part_diag)

print(f"\nRunning {max_steps} steps (diag every {period}) …")
sim.step(max_steps)
print("\nDone. openPMD output → warpx_gun/diags/{fields,particles}/")
