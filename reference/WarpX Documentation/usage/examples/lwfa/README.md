<a id="examples-lwfa"></a>

# Laser-Wakefield Acceleration of Electrons

This example shows how to model a laser-wakefield accelerator (LWFA) [[2](../../examples.md#id3), [3](../../examples.md#id4)].

Laser-wakefield acceleration is best performed in 3D or quasi-cylindrical (RZ) geometry, in order to correctly capture some of the key physics (laser diffraction, beamloading, shape of the accelerating bubble in the blowout regime, etc.).
For physical situations that have close-to-cylindrical symmetry, simulations in RZ geometry capture the relevant physics at a fraction of the computational cost of a 3D simulation.
On the other hand, for physical situation with strong asymmetries (e.g., non-round laser driver, strong hosing of the accelerated beam, etc.), only 3D simulations are suitable.

For LWFA scenarios with long propagation lengths, use the [boosted frame method](../../../theory/boosted_frame.md#theory-boostedframe).
An example can be seen in the [PWFA example](../pwfa/README.md#examples-pwfa).

## Run

For [MPI-parallel](https://www.mpi-forum.org) runs, prefix these lines with `mpiexec -n 4 ...` or `srun -n 4 ...`, depending on the system.

### 3D

This example can be run **either** as:

* **Python** script: `python3 inputs_test_3d_laser_acceleration_picmi.py` or
* WarpX **executable** using an input file: `warpx.3d inputs_test_3d_laser_acceleration max_step=400`

### Python: Script

```python3
#!/usr/bin/env python3

from pywarpx import picmi

# Physical constants
c = picmi.constants.c
q_e = picmi.constants.q_e

# Number of time steps
max_steps = 100

# Number of cells
nx = 32
ny = 32
nz = 256

# Physical domain
xmin = -30e-06
xmax = 30e-06
ymin = -30e-06
ymax = 30e-06
zmin = -56e-06
zmax = 12e-06

# Domain decomposition
max_grid_size = 64
blocking_factor = 32

# Create grid
grid = picmi.Cartesian3DGrid(
    number_of_cells=[nx, ny, nz],
    lower_bound=[xmin, ymin, zmin],
    upper_bound=[xmax, ymax, zmax],
    lower_boundary_conditions=["periodic", "periodic", "dirichlet"],
    upper_boundary_conditions=["periodic", "periodic", "dirichlet"],
    lower_boundary_conditions_particles=["periodic", "periodic", "absorbing"],
    upper_boundary_conditions_particles=["periodic", "periodic", "absorbing"],
    moving_window_velocity=[0.0, 0.0, c],
    warpx_max_grid_size=max_grid_size,
    warpx_blocking_factor=blocking_factor,
)

# Particles: plasma electrons
plasma_density = 2e23
plasma_xmin = -20e-06
plasma_ymin = -20e-06
plasma_zmin = 0
plasma_xmax = 20e-06
plasma_ymax = 20e-06
plasma_zmax = None
uniform_distribution = picmi.UniformDistribution(
    density=plasma_density,
    lower_bound=[plasma_xmin, plasma_ymin, plasma_zmin],
    upper_bound=[plasma_xmax, plasma_ymax, plasma_zmax],
    fill_in=True,
)
electrons = picmi.Species(
    particle_type="electron",
    name="electrons",
    initial_distribution=uniform_distribution,
    warpx_add_int_attributes={"regionofinterest": "(z>12.0e-6) * (z<13.0e-6)"},
    warpx_add_real_attributes={"initialenergy": "ux*ux + uy*uy + uz*uz"},
)

# Particles: beam electrons
q_tot = 1e-12
x_m = 0.0
y_m = 0.0
z_m = -28e-06
x_rms = 0.5e-06
y_rms = 0.5e-06
z_rms = 0.5e-06
ux_m = 0.0
uy_m = 0.0
uz_m = 500.0
ux_th = 2.0
uy_th = 2.0
uz_th = 50.0
gaussian_bunch_distribution = picmi.GaussianBunchDistribution(
    n_physical_particles=q_tot / q_e,
    rms_bunch_size=[x_rms, y_rms, z_rms],
    rms_velocity=[c * ux_th, c * uy_th, c * uz_th],
    centroid_position=[x_m, y_m, z_m],
    centroid_velocity=[c * ux_m, c * uy_m, c * uz_m],
)
beam = picmi.Species(
    particle_type="electron",
    name="beam",
    initial_distribution=gaussian_bunch_distribution,
)

# Laser
e_max = 16e12
position_z = 9e-06
profile_t_peak = 30.0e-15
profile_focal_distance = 100e-06
laser = picmi.GaussianLaser(
    wavelength=0.8e-06,
    waist=5e-06,
    duration=15e-15,
    focal_position=[0, 0, profile_focal_distance + position_z],
    centroid_position=[0, 0, position_z - c * profile_t_peak],
    propagation_direction=[0, 0, 1],
    polarization_direction=[0, 1, 0],
    E0=e_max,
    fill_in=False,
)
laser_antenna = picmi.LaserAntenna(
    position=[0.0, 0.0, position_z], normal_vector=[0, 0, 1]
)

# Electromagnetic solver
solver = picmi.ElectromagneticSolver(grid=grid, method="Yee", cfl=1.0, divE_cleaning=0)

# Diagnostics
diag_field_list = ["B", "E", "J", "rho"]
particle_diag = picmi.ParticleDiagnostic(
    name="diag1",
    period=100,
)
field_diag = picmi.FieldDiagnostic(
    name="diag1",
    grid=grid,
    period=100,
    data_list=diag_field_list,
)

# Set up simulation
sim = picmi.Simulation(
    solver=solver,
    max_steps=max_steps,
    verbose=1,
    particle_shape="cubic",
    warpx_use_filter=1,
    warpx_serialize_initial_conditions=1,
    warpx_do_dynamic_scheduling=0,
)

# Add plasma electrons
sim.add_species(
    electrons, layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=[1, 1, 1])
)

# Add beam electrons
sim.add_species(beam, layout=picmi.PseudoRandomLayout(grid=grid, n_macroparticles=100))

# Add laser
sim.add_laser(laser, injection_method=laser_antenna)

# Add diagnostics
sim.add_diagnostic(particle_diag)
sim.add_diagnostic(field_diag)

# Write input file that can be used to run with the compiled version
sim.write_input_file(file_name="inputs_3d_picmi")

# Initialize inputs and WarpX instance
sim.initialize_inputs()
sim.initialize_warpx()

# Advance simulation until last time step
sim.step(max_steps)
```

### Executable: Input File

```none
#################################
####### GENERAL PARAMETERS ######
#################################
max_step = 100           # for production, run for longer time, e.g. max_step = 1000
amr.n_cell = 32 32 256   # for production, run with finer mesh, e.g. amr.n_cell = 64 64 512
amr.max_grid_size = 64   # maximum size of each AMReX box, used to decompose the domain
amr.blocking_factor = 32 # minimum size of each AMReX box, used to decompose the domain
geometry.dims = 3
geometry.prob_lo     = -30.e-6   -30.e-6   -56.e-6    # physical domain
geometry.prob_hi     =  30.e-6    30.e-6    12.e-6
amr.max_level = 0 # Maximum level in hierarchy (1 might be unstable, >1 is not supported)
# warpx.fine_tag_lo = -5.e-6   -5.e-6   -50.e-6
# warpx.fine_tag_hi =  5.e-6    5.e-6   -30.e-6

#################################
####### Boundary condition ######
#################################
boundary.field_lo = periodic periodic pec
boundary.field_hi = periodic periodic pec

#################################
############ NUMERICS ###########
#################################
warpx.verbose = 1
warpx.do_dive_cleaning = 0
warpx.use_filter = 1
warpx.cfl = 1. # if 1., the time step is set to its CFL limit
warpx.do_moving_window = 1
warpx.moving_window_dir = z
warpx.moving_window_v = 1.0 # units of speed of light
warpx.do_dynamic_scheduling = 0 # for production, set this to 1 (default)
warpx.serialize_initial_conditions = 1         # for production, set this to 0 (default)

# Order of particle shape factors
algo.particle_shape = 3

#################################
############ PLASMA #############
#################################
particles.species_names = electrons

electrons.charge = -q_e
electrons.mass = m_e
electrons.injection_style = "NUniformPerCell"
electrons.num_particles_per_cell_each_dim = 1 1 1
electrons.xmin = -20.e-6
electrons.xmax =  20.e-6
electrons.ymin = -20.e-6
electrons.ymax =  20.e-6
electrons.zmin =  0
electrons.profile = constant
electrons.density = 2.e23  # number of electrons per m^3
electrons.momentum_distribution_type = "at_rest"
electrons.do_continuous_injection = 1
electrons.addIntegerAttributes = regionofinterest
electrons.attribute.regionofinterest(x,y,z,ux,uy,uz,t) = "(z>12.0e-6) * (z<13.0e-6)"
electrons.addRealAttributes = initialenergy
electrons.attribute.initialenergy(x,y,z,ux,uy,uz,t) = "ux*ux + uy*uy + uz*uz"

#################################
############ LASER  #############
#################################
lasers.names        = laser1
laser1.profile      = Gaussian
laser1.position     = 0. 0. 9.e-6        # This point is on the laser plane
laser1.direction    = 0. 0. 1.           # The plane normal direction
laser1.polarization = 0. 1. 0.           # The main polarization vector
laser1.e_max        = 16.e12             # Maximum amplitude of the laser field (in V/m)
laser1.profile_waist = 5.e-6             # The waist of the laser (in m)
laser1.profile_duration = 15.e-15        # The duration of the laser (in s)
laser1.profile_t_peak = 30.e-15          # Time at which the laser reaches its peak (in s)
laser1.profile_focal_distance = 100.e-6  # Focal distance from the antenna (in m)
laser1.wavelength = 0.8e-6               # The wavelength of the laser (in m)

# Diagnostics
diagnostics.diags_names = diag1
diag1.intervals = 100
diag1.diag_type = Full
diag1.fields_to_plot = Ex Ey Ez Bx By Bz jx jy jz rho
diag1.format = openpmd

# Reduced Diagnostics
warpx.reduced_diags_names               = FP

FP.type = FieldProbe
FP.intervals = 10
FP.integrate = 0
FP.probe_geometry = Line
FP.x_probe = 0
FP.y_probe = 0
FP.z_probe = -56e-6
FP.x1_probe = 0
FP.y1_probe = 0
FP.z1_probe = 12e-6
FP.resolution = 300
FP.do_moving_window_FP = 1
```

### RZ

This example can be run **either** as:

* **Python** script: `python3 inputs_test_rz_laser_acceleration_picmi.py` or
* WarpX **executable** using an input file: `warpx.rz inputs_test_rz_laser_acceleration max_step=400`

### Python: Script

```python3
#!/usr/bin/env python3

from pywarpx import picmi

# Physical constants
c = picmi.constants.c
q_e = picmi.constants.q_e

# Number of time steps
max_steps = 10

# Number of cells
nr = 64
nz = 512

# Physical domain
rmin = 0
rmax = 30e-06
zmin = -56e-06
zmax = 12e-06

# Domain decomposition
max_grid_size = 64
blocking_factor = 32

# Create grid
grid = picmi.CylindricalGrid(
    number_of_cells=[nr, nz],
    n_azimuthal_modes=2,
    lower_bound=[rmin, zmin],
    upper_bound=[rmax, zmax],
    lower_boundary_conditions=["none", "dirichlet"],
    upper_boundary_conditions=["dirichlet", "dirichlet"],
    lower_boundary_conditions_particles=["none", "absorbing"],
    upper_boundary_conditions_particles=["absorbing", "absorbing"],
    moving_window_velocity=[0.0, c],
    warpx_max_grid_size=max_grid_size,
    warpx_blocking_factor=blocking_factor,
)

# Particles: plasma electrons
plasma_density = 2e23
plasma_xmin = -20e-06
plasma_ymin = None
plasma_zmin = 10e-06
plasma_xmax = 20e-06
plasma_ymax = None
plasma_zmax = None
uniform_distribution = picmi.UniformDistribution(
    density=plasma_density,
    lower_bound=[plasma_xmin, plasma_ymin, plasma_zmin],
    upper_bound=[plasma_xmax, plasma_ymax, plasma_zmax],
    fill_in=True,
)
electrons = picmi.Species(
    particle_type="electron",
    name="electrons",
    initial_distribution=uniform_distribution,
)

# Particles: beam electrons
q_tot = 1e-12
x_m = 0.0
y_m = 0.0
z_m = -28e-06
x_rms = 0.5e-06
y_rms = 0.5e-06
z_rms = 0.5e-06
ux_m = 0.0
uy_m = 0.0
uz_m = 500.0
ux_th = 2.0
uy_th = 2.0
uz_th = 50.0
gaussian_bunch_distribution = picmi.GaussianBunchDistribution(
    n_physical_particles=q_tot / q_e,
    rms_bunch_size=[x_rms, y_rms, z_rms],
    rms_velocity=[c * ux_th, c * uy_th, c * uz_th],
    centroid_position=[x_m, y_m, z_m],
    centroid_velocity=[c * ux_m, c * uy_m, c * uz_m],
)
beam = picmi.Species(
    particle_type="electron",
    name="beam",
    initial_distribution=gaussian_bunch_distribution,
)

# Laser
e_max = 16e12
position_z = 9e-06
profile_t_peak = 30.0e-15
profile_focal_distance = 100e-06
laser = picmi.GaussianLaser(
    wavelength=0.8e-06,
    waist=5e-06,
    duration=15e-15,
    focal_position=[0, 0, profile_focal_distance + position_z],
    centroid_position=[0, 0, position_z - c * profile_t_peak],
    propagation_direction=[0, 0, 1],
    polarization_direction=[0, 1, 0],
    E0=e_max,
    fill_in=False,
)
laser_antenna = picmi.LaserAntenna(
    position=[0.0, 0.0, position_z], normal_vector=[0, 0, 1]
)

# Electromagnetic solver
solver = picmi.ElectromagneticSolver(grid=grid, method="Yee", cfl=1.0, divE_cleaning=0)

# Diagnostics
diag_field_list = ["B", "E", "J", "rho"]
field_diag = picmi.FieldDiagnostic(
    name="diag1",
    grid=grid,
    period=10,
    data_list=diag_field_list,
    warpx_dump_rz_modes=1,
)
diag_particle_list = ["weighting", "momentum"]
particle_diag = picmi.ParticleDiagnostic(
    name="diag1",
    period=10,
    species=[electrons, beam],
    data_list=diag_particle_list,
)

# Set up simulation
sim = picmi.Simulation(
    solver=solver,
    max_steps=max_steps,
    verbose=1,
    particle_shape="cubic",
    warpx_use_filter=0,
)

# Add plasma electrons
sim.add_species(
    electrons, layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=[1, 4, 1])
)

# Add beam electrons
sim.add_species(beam, layout=picmi.PseudoRandomLayout(grid=grid, n_macroparticles=100))

# Add laser
sim.add_laser(laser, injection_method=laser_antenna)

# Add diagnostics
sim.add_diagnostic(field_diag)
sim.add_diagnostic(particle_diag)

# Write input file that can be used to run with the compiled version
sim.write_input_file(file_name="inputs_rz_picmi")

# Initialize inputs and WarpX instance
sim.initialize_inputs()
sim.initialize_warpx()

# Advance simulation until last time step
sim.step(max_steps)
```

### Executable: Input File

```none
#################################
####### GENERAL PARAMETERS ######
#################################
max_step = 10
amr.n_cell =  64  512
amr.max_grid_size = 64   # maximum size of each AMReX box, used to decompose the domain
amr.blocking_factor = 32 # minimum size of each AMReX box, used to decompose the domain
geometry.dims = RZ
geometry.prob_lo     =   0.   -56.e-6    # physical domain
geometry.prob_hi     =  30.e-6    12.e-6
amr.max_level = 0 # Maximum level in hierarchy (1 might be unstable, >1 is not supported)

warpx.n_rz_azimuthal_modes = 2

boundary.field_lo = none pec
boundary.field_hi = pec pec

#################################
############ NUMERICS ###########
#################################
warpx.verbose = 1
warpx.do_dive_cleaning = 0
warpx.use_filter = 1
warpx.filter_npass_each_dir = 0 1
warpx.cfl = 1. # if 1., the time step is set to its CFL limit
warpx.do_moving_window = 1
warpx.moving_window_dir = z
warpx.moving_window_v = 1.0 # units of speed of light

# Order of particle shape factors
algo.particle_shape = 3

#################################
############ PLASMA #############
#################################
particles.species_names = electrons beam

electrons.charge = -q_e
electrons.mass = m_e
electrons.injection_style = "NUniformPerCell"
electrons.num_particles_per_cell_each_dim = 1 4 1
electrons.xmin = -20.e-6
electrons.xmax =  20.e-6
electrons.zmin =  10.e-6
electrons.profile = constant
electrons.density = 2.e23  # number of electrons per m^3
electrons.momentum_distribution_type = "at_rest"
electrons.do_continuous_injection = 1
electrons.addRealAttributes = orig_x orig_z
electrons.attribute.orig_x(x,y,z,ux,uy,uz,t) = "x"
electrons.attribute.orig_z(x,y,z,ux,uy,uz,t) = "z"

beam.charge = -q_e
beam.mass = m_e
beam.injection_style = "gaussian_beam"
beam.x_rms = .5e-6
beam.y_rms = .5e-6
beam.z_rms = .5e-6
beam.x_m = 0.
beam.y_m = 0.
beam.z_m = -28.e-6
beam.npart = 100
beam.q_tot = -1.e-12
beam.momentum_distribution_type = "gaussian"
beam.ux_m = 0.0
beam.uy_m = 0.0
beam.uz_m = 500.
beam.ux_th = 2.
beam.uy_th = 2.
beam.uz_th = 50.

#################################
############ LASER ##############
#################################
lasers.names        = laser1
laser1.profile      = Gaussian
laser1.position     = 0. 0. 9.e-6        # This point is on the laser plane
laser1.direction    = 0. 0. 1.           # The plane normal direction
laser1.polarization = 0. 1. 0.           # The main polarization vector
laser1.e_max        = 16.e12             # Maximum amplitude of the laser field (in V/m)
laser1.profile_waist = 5.e-6             # The waist of the laser (in m)
laser1.profile_duration = 15.e-15        # The duration of the laser (in s)
laser1.profile_t_peak = 30.e-15          # Time at which the laser reaches its peak (in s)
laser1.profile_focal_distance = 100.e-6  # Focal distance from the antenna (in m)
laser1.wavelength = 0.8e-6               # The wavelength of the laser (in m)

# Diagnostics
diagnostics.diags_names = diag1
diag1.intervals = 10
diag1.diag_type = Full
diag1.fields_to_plot = Er Et Ez Br Bt Bz jr jt jz rho
diag1.electrons.variables = x y z w ux uy uz orig_x orig_z
diag1.beam.variables = x y z w ux uy uz
```

## Analyze

#### NOTE
This section is TODO.

## Visualize

You can run the following script to visualize the beam evolution over time:

### Script `plot_3d.py`

```python3
#!/usr/bin/env python3

# Copyright 2023 The WarpX Community
#
# This file is part of WarpX.
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# This is a script plots the wakefield of an LWFA simulation.

import sys

import matplotlib.pyplot as plt
import yt

yt.funcs.mylog.setLevel(50)


def plot_lwfa():
    # this will be the name of the plot file
    fn = sys.argv[1]

    # Read the file
    ds = yt.load(fn)

    # plot the laser field and absolute density
    fields = ["Ey", "rho"]
    normal = "y"
    sl = yt.SlicePlot(ds, normal=normal, fields=fields)
    for field in fields:
        sl.set_log(field, False)

    sl.set_figure_size((4, 8))
    fig = sl.export_to_mpl_figure(nrows_ncols=(2, 1))
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    plot_lwfa()
```

![(top) Electric field of the laser pulse and (bottom) absolute density.](https://user-images.githubusercontent.com/1353258/287800852-f994a020-4ecc-4987-bffc-2cb7df6144a9.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTEiLCJleHAiOjE3MDE3MTYyNjksIm5iZiI6MTcwMTcxNTk2OSwicGF0aCI6Ii8xMzUzMjU4LzI4NzgwMDg1Mi1mOTk0YTAyMC00ZWNjLTQ5ODctYmZmYy0yY2I3ZGY2MTQ0YTkucG5nP1gtQW16LUFsZ29yaXRobT1BV1M0LUhNQUMtU0hBMjU2JlgtQW16LUNyZWRlbnRpYWw9QUtJQUlXTkpZQVg0Q1NWRUg1M0ElMkYyMDIzMTIwNCUyRnVzLWVhc3QtMSUyRnMzJTJGYXdzNF9yZXF1ZXN0JlgtQW16LURhdGU9MjAyMzEyMDRUMTg1MjQ5WiZYLUFtei1FeHBpcmVzPTMwMCZYLUFtei1TaWduYXR1cmU9NDkyNWJkMTg2NWM3ZjcwZjVkMjlmNDE1NmRjNWEyZWM5MzgxMWJhZTVjMGMxNjdkZDg1Zjk0NmQ1NGEwMjNiMiZYLUFtei1TaWduZWRIZWFkZXJzPWhvc3QmYWN0b3JfaWQ9MCZrZXlfaWQ9MCZyZXBvX2lkPTAifQ.C_NQceQcqiCDzBoSzIjm3c8QdTLNDdtjJmkQjkhW4c8)
