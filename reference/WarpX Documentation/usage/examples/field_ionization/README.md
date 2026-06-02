<a id="examples-tests-field-ionization"></a>

# Field Ionization

## Run Test

For [MPI-parallel](https://www.mpi-forum.org) runs, prefix these lines with `mpiexec -n 4 ...` or `srun -n 4 ...`, depending on the system.

### lab frame

This example can be run **either** as:

* **Python** script: `python3 inputs_test_2d_ionization_picmi.py` or
* WarpX **executable** using an input file: `warpx.2d inputs_test_2d_ionization_lab max_step=1600`

### Python: Script

```python3
#!/usr/bin/env python3

from pywarpx import picmi

# Physical constants
c = picmi.constants.c

# Number of time steps
max_steps = 1600

# Number of cells
nx = 16
nz = 800

# Physical domain
xmin = -5e-06
xmax = 5e-06
zmin = 0e-06
zmax = 20e-06

# Domain decomposition
max_grid_size = 64
blocking_factor = 16

# Create grid
grid = picmi.Cartesian2DGrid(
    number_of_cells=[nx, nz],
    lower_bound=[xmin, zmin],
    upper_bound=[xmax, zmax],
    lower_boundary_conditions=["periodic", "open"],
    upper_boundary_conditions=["periodic", "open"],
    warpx_max_grid_size=max_grid_size,
    warpx_blocking_factor=blocking_factor,
)

# Particles: electrons and ions
ions_density = 1
ions_xmin = None
ions_ymin = None
ions_zmin = 5e-06
ions_xmax = None
ions_ymax = None
ions_zmax = 15e-06
uniform_distribution = picmi.UniformDistribution(
    density=ions_density,
    lower_bound=[ions_xmin, ions_ymin, ions_zmin],
    upper_bound=[ions_xmax, ions_ymax, ions_zmax],
    fill_in=True,
)
electrons = picmi.Species(
    particle_type="electron",
    name="electrons",
    warpx_add_real_attributes={"orig_z": "z"},
)
ions = picmi.Species(
    particle_type="N",
    name="ions",
    charge_state=2,
    initial_distribution=uniform_distribution,
    warpx_add_real_attributes={"orig_z": "z"},
)

# Field ionization
nitrogen_ionization = picmi.FieldIonization(
    model="ADK",  # Ammosov-Delone-Krainov model
    ionized_species=ions,
    product_species=electrons,
)

# Laser
position_z = 3e-06
profile_t_peak = 60.0e-15
laser = picmi.GaussianLaser(
    wavelength=0.8e-06,
    waist=1e10,
    duration=26.685e-15,
    focal_position=[0, 0, position_z],
    centroid_position=[0, 0, position_z - c * profile_t_peak],
    propagation_direction=[0, 0, 1],
    polarization_direction=[1, 0, 0],
    a0=1.8,
    fill_in=False,
)
laser_antenna = picmi.LaserAntenna(
    position=[0.0, 0.0, position_z], normal_vector=[0, 0, 1]
)

# Electromagnetic solver
solver = picmi.ElectromagneticSolver(grid=grid, method="CKC", cfl=0.999)

# Diagnostics
particle_diag = picmi.ParticleDiagnostic(
    name="diag1",
    period=10000,
    species=[electrons, ions],
    data_list=["ux", "uy", "uz", "x", "z", "weighting", "orig_z"],
)
field_diag = picmi.FieldDiagnostic(
    name="diag1",
    grid=grid,
    period=10000,
    data_list=["Bx", "By", "Bz", "Ex", "Ey", "Ez", "Jx", "Jy", "Jz"],
)

# Set up simulation
sim = picmi.Simulation(
    solver=solver, max_steps=max_steps, particle_shape="linear", warpx_use_filter=0
)

# Add electrons and ions
sim.add_species(
    electrons, layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=[0, 0])
)
sim.add_species(
    ions, layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=[2, 1])
)

# Add field ionization
sim.add_interaction(nitrogen_ionization)

# Add laser
sim.add_laser(laser, injection_method=laser_antenna)

# Add diagnostics
sim.add_diagnostic(particle_diag)
sim.add_diagnostic(field_diag)

# Write input file that can be used to run with the compiled version
sim.write_input_file(file_name="inputs_2d_picmi")

# Initialize inputs and WarpX instance
sim.initialize_inputs()
sim.initialize_warpx()

# Advance simulation until last time step
sim.step(max_steps)
```

### Executable: Input File

```none
max_step = 1600
amr.n_cell =  16 800
amr.max_grid_size = 64
amr.blocking_factor = 16
geometry.dims = 2
geometry.prob_lo     = -5.e-6   0.e-6
geometry.prob_hi     =  5.e-6  20.e-6
amr.max_level = 0

boundary.field_lo = periodic pml
boundary.field_hi = periodic pml

algo.maxwell_solver = ckc
warpx.cfl = .999
warpx.use_filter = 0

# Order of particle shape factors
algo.particle_shape = 1

particles.species_names = electrons ions

ions.mass = 2.3428415e-26
ions.charge = q_e
ions.injection_style = nuniformpercell
ions.num_particles_per_cell_each_dim = 2 1
ions.zmin =  5.e-6
ions.zmax = 15.e-6
ions.profile = constant
ions.density = 1.
ions.momentum_distribution_type = at_rest
ions.do_field_ionization = 1
ions.ionization_initial_level = 2
ions.ionization_product_species = electrons
ions.physical_element = N

electrons.mass = m_e
electrons.charge = -q_e
electrons.injection_style = none
electrons.addRealAttributes = orig_z
electrons.attribute.orig_z(x,y,z,ux,uy,uz,t) = z

lasers.names        = laser1
laser1.profile      = Gaussian
laser1.position     = 0. 0. 3.e-6
laser1.direction    = 0. 0. 1.
laser1.polarization = 1. 0. 0.
laser1.a0           = 1.8
laser1.profile_waist = 1.e10
laser1.profile_duration = 26.685e-15
laser1.profile_t_peak = 60.e-15
laser1.profile_focal_distance = 0
laser1.wavelength = 0.8e-6

# Diagnostics
diagnostics.diags_names = diag1 chk

diag1.intervals = 10000
diag1.diag_type = Full

chk.intervals = 1000
chk.diag_type = Full
chk.format = checkpoint
```

### boosted frame

This example can be run as:

* WarpX **executable** using an input file: `warpx.2d inputs_test_2d_ionization_boost max_step=420`

```none
max_step = 420
amr.n_cell =  16 800
amr.max_grid_size = 64
amr.blocking_factor = 16
geometry.dims = 2
geometry.prob_lo     = -5.e-6 -40.e-6
geometry.prob_hi     =  5.e-6   0.e-6
amr.max_level = 0

boundary.field_lo = periodic pml
boundary.field_hi = periodic pml

algo.maxwell_solver = ckc
warpx.cfl = .999
warpx.do_moving_window = 1
warpx.moving_window_dir = z
warpx.moving_window_v = 1.0
warpx.gamma_boost = 2.
warpx.boost_direction = z
warpx.use_filter = 0

# Order of particle shape factors
algo.particle_shape = 1

particles.species_names = electrons ions

ions.mass = 2.3428415e-26
ions.charge = q_e
ions.injection_style = nuniformpercell
ions.num_particles_per_cell_each_dim = 2 2
ions.zmin =  0.
ions.zmax =  50.e-6
ions.profile = constant
ions.density = 1.
ions.momentum_distribution_type = at_rest
ions.do_field_ionization = 1
ions.ionization_initial_level = 2
ions.ionization_product_species = electrons
ions.physical_element = N
ions.do_continuous_injection=1

electrons.mass = m_e
electrons.charge = -q_e
electrons.injection_style = nuniformpercell
electrons.num_particles_per_cell_each_dim = 2 2
electrons.zmin = 0.
electrons.zmax =  50.e-6
electrons.profile = constant
electrons.density = 2.
electrons.momentum_distribution_type = at_rest
electrons.do_continuous_injection = 1

lasers.names        = laser1
laser1.profile      = Gaussian
laser1.position     = 0. 0. -1.e-6
laser1.direction    = 0. 0. 1.
laser1.polarization = 1. 0. 0.
laser1.a0           = 1.8
laser1.profile_waist = 1.e10
laser1.profile_duration = 26.685e-15
laser1.profile_t_peak = 60.e-15
laser1.profile_focal_distance = 0
laser1.wavelength = 0.8e-6

# Diagnostics
diagnostics.diags_names = diag1
diag1.intervals = 10000
diag1.diag_type = Full
```

## Analyze

### Script `analysis.py`

```python3
#!/usr/bin/env python3

# Copyright 2019-2020 Luca Fedeli, Maxence Thevenet
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


"""
This script tests the result of the ionization module in WarpX.

Input files inputs.rt and inputs.bf.rt are used to reproduce the test from
Chen, JCP, 2013, figure 2 (in the lab frame and in a boosted frame,
respectively): a plane-wave laser pulse propagates through a
uniform N2+ neutral plasma and further ionizes the Nitrogen atoms. This test
checks that, after the laser went through the plasma, ~32% of Nitrogen
ions are N5+, in agreement with theory from Chen's article.
"""

import sys

import numpy as np
import yt

yt.funcs.mylog.setLevel(0)

# Open plotfile specified in command line, and get ion's ionization level.
filename = sys.argv[1]
ds = yt.load(filename)
ad = ds.all_data()
ilev = ad["ions", "particle_ionizationLevel"].v

# Fraction of Nitrogen ions that are N5+.
N5_fraction = ilev[ilev == 5].size / float(ilev.size)

print("Number of ions: " + str(ilev.size))
print("Number of N5+ : " + str(ilev[ilev == 5].size))
print("N5_fraction   : " + str(N5_fraction))

do_plot = False
if do_plot:
    import matplotlib.pyplot as plt

    all_data_level_0 = ds.covering_grid(
        level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions
    )
    F = all_data_level_0["boxlib", "Ex"].v.squeeze()
    extent = [
        ds.domain_left_edge[1],
        ds.domain_right_edge[1],
        ds.domain_left_edge[0],
        ds.domain_right_edge[0],
    ]
    ad = ds.all_data()

    # Plot ions with ionization levels
    species = "ions"
    xi = ad[species, "particle_position_x"].v
    zi = ad[species, "particle_position_y"].v
    ii = ad[species, "particle_ionizationLevel"].v
    plt.figure(figsize=(10, 10))
    plt.subplot(211)
    plt.imshow(np.abs(F), extent=extent, aspect="auto", cmap="magma", origin="default")
    plt.colorbar()
    for lev in range(int(np.max(ii) + 1)):
        select = ii == lev
        plt.scatter(
            zi[select], xi[select], s=0.2, label="ionization level: " + str(lev)
        )
    plt.legend()
    plt.title("abs(Ex) (V/m) and ions")
    plt.xlabel("z (m)")
    plt.ylabel("x (m)")
    plt.subplot(212)
    plt.imshow(np.abs(F), extent=extent, aspect="auto", cmap="magma", origin="default")
    plt.colorbar()

    # Plot electrons
    species = "electrons"
    if species in [x[0] for x in ds.field_list]:
        xe = ad[species, "particle_position_x"].v
        ze = ad[species, "particle_position_y"].v
        plt.scatter(ze, xe, s=0.1, c="r", label="electrons")
    plt.title("abs(Ex) (V/m) and electrons")
    plt.xlabel("z (m)")
    plt.ylabel("x (m)")
    plt.savefig("image_ionization.pdf", bbox_inches="tight")

error_rel = abs(N5_fraction - 0.32) / 0.32
tolerance_rel = 0.07

print("error_rel    : " + str(error_rel))
print("tolerance_rel: " + str(tolerance_rel))

assert error_rel < tolerance_rel

# Check that the user runtime component (if it exists) worked as expected
try:
    orig_z = ad["electrons", "particle_orig_z"].v
    print(f"orig_z: min = {np.min(orig_z)}, max = {np.max(orig_z)}")
    assert np.all((orig_z > 0.0) & (orig_z < 1.5e-5))
    print("particle_orig_z has reasonable values")
except yt.utilities.exceptions.YTFieldNotFound:
    pass  # The backtransformed diagnostic version of the test does not have orig_z
```

## Visualize

![Electric field of the laser pulse with (top) ions with ionization levels and (bottom) ionized electrons.](https://gist.githubusercontent.com/johvandewetering/48d092c003915f1d1689b507caa2865b/raw/29f5d12ed77831047ca12f456a07dbf3b99770d5/image_ionization.png)
