<a id="examples-pwfa"></a>

# Beam-Driven Wakefield Acceleration of Electrons

This example shows how to model a beam-driven plasma-wakefield accelerator (PWFA) [[2](../../examples.md#id3), [3](../../examples.md#id4)].

PWFA is best performed in 3D or quasi-cylindrical (RZ) geometry, in order to correctly capture some of the key physics (structure of the space-charge fields, beamloading, shape of the accelerating bubble in the blowout regime, etc.).
For physical situations that have close-to-cylindrical symmetry, simulations in RZ geometry capture the relevant physics at a fraction of the computational cost of a 3D simulation.
On the other hand, for physical situation with strong asymmetries (e.g., non-round driver, strong hosing of the accelerated beam, etc.), only 3D simulations are suitable.

Additionally, to speed up computation, this example uses the [boosted frame method](../../../theory/boosted_frame.md#theory-boostedframe) to effectively model long acceleration lengths.

Alternatively, an other common approximation for PWFAs is quasi-static modeling, e.g., if effects such as self-injection can be ignored.
In the Beam, Plasma & Accelerator Simulation Toolkit (BLAST), [HiPACE++](https://hipace.readthedocs.io) provides such methods.

#### NOTE
TODO: The Python (PICMI) input file should use the boosted frame method, like the `inputs_test_3d_plasma_acceleration_boosted` file.

## Run

This example can be run **either** as:

* **Python** script: `python3 inputs_test_3d_plasma_acceleration_picmi.py` or
* WarpX **executable** using an input file: `warpx.3d inputs_test_3d_plasma_acceleration_boosted`

For [MPI-parallel](https://www.mpi-forum.org) runs, prefix these lines with `mpiexec -n 4 ...` or `srun -n 4 ...`, depending on the system.

### Python: Script

#### NOTE
TODO: This input file should use the boosted frame method, like the `inputs_test_3d_plasma_acceleration_boosted` file.

```python3
#!/usr/bin/env python3

from pywarpx import picmi

# from warp import picmi

constants = picmi.constants

nx = 64
ny = 64
nz = 64

xmin = -200.0e-6
xmax = +200.0e-6
ymin = -200.0e-6
ymax = +200.0e-6
zmin = -200.0e-6
zmax = +200.0e-6

moving_window_velocity = [0.0, 0.0, constants.c]

number_per_cell_each_dim = [2, 2, 1]

max_steps = 10

grid = picmi.Cartesian3DGrid(
    number_of_cells=[nx, ny, nz],
    lower_bound=[xmin, ymin, zmin],
    upper_bound=[xmax, ymax, zmax],
    lower_boundary_conditions=["periodic", "periodic", "open"],
    upper_boundary_conditions=["periodic", "periodic", "open"],
    lower_boundary_conditions_particles=["periodic", "periodic", "absorbing"],
    upper_boundary_conditions_particles=["periodic", "periodic", "absorbing"],
    moving_window_velocity=moving_window_velocity,
    warpx_max_grid_size=32,
)

solver = picmi.ElectromagneticSolver(grid=grid, cfl=1)

beam_distribution = picmi.UniformDistribution(
    density=1.0e23,
    lower_bound=[-20.0e-6, -20.0e-6, -150.0e-6],
    upper_bound=[+20.0e-6, +20.0e-6, -100.0e-6],
    directed_velocity=[0.0, 0.0, 1.0e9],
)

plasma_distribution = picmi.UniformDistribution(
    density=1.0e22,
    lower_bound=[-200.0e-6, -200.0e-6, 0.0],
    upper_bound=[+200.0e-6, +200.0e-6, None],
    fill_in=True,
)

beam = picmi.Species(
    particle_type="electron", name="beam", initial_distribution=beam_distribution
)
plasma = picmi.Species(
    particle_type="electron", name="plasma", initial_distribution=plasma_distribution
)

sim = picmi.Simulation(
    solver=solver,
    max_steps=max_steps,
    verbose=1,
    warpx_current_deposition_algo="esirkepov",
    warpx_use_filter=0,
)

sim.add_species(
    beam,
    layout=picmi.GriddedLayout(
        grid=grid, n_macroparticle_per_cell=number_per_cell_each_dim
    ),
)
sim.add_species(
    plasma,
    layout=picmi.GriddedLayout(
        grid=grid, n_macroparticle_per_cell=number_per_cell_each_dim
    ),
)

field_diag = picmi.FieldDiagnostic(
    name="diag1",
    grid=grid,
    period=max_steps,
    data_list=["Ex", "Ey", "Ez", "Jx", "Jy", "Jz", "part_per_cell"],
)

part_diag = picmi.ParticleDiagnostic(
    name="diag1",
    period=max_steps,
    species=[beam, plasma],
    data_list=["ux", "uy", "uz", "weighting"],
)

sim.add_diagnostic(field_diag)
sim.add_diagnostic(part_diag)

# write_inputs will create an inputs file that can be used to run
# with the compiled version.
# sim.write_input_file(file_name = 'inputs_from_PICMI')

# Alternatively, sim.step will run WarpX, controlling it from Python
sim.step()
```

### Executable: Input File

```none
#################################
####### GENERAL PARAMETERS ######
#################################
stop_time = 3.93151387287e-11
amr.n_cell = 64 64 128 #32 32 256
amr.max_grid_size = 64
amr.blocking_factor = 32
amr.max_level = 0
geometry.dims = 3
geometry.prob_lo = -0.00015 -0.00015 -0.00012
geometry.prob_hi = 0.00015 0.00015 1.e-06

#################################
####### Boundary condition ######
#################################
boundary.field_lo = periodic periodic pml
boundary.field_hi = periodic periodic pml

#################################
############ NUMERICS ###########
#################################
algo.maxwell_solver = ckc
warpx.verbose = 1
warpx.do_dive_cleaning = 0
warpx.use_filter = 1
warpx.cfl = .99
warpx.do_moving_window = 1
warpx.moving_window_dir = z
warpx.moving_window_v = 1. # in units of the speed of light
my_constants.lramp = 8.e-3
my_constants.dens  = 1e+23

# Order of particle shape factors
algo.particle_shape = 3

#################################
######### BOOSTED FRAME #########
#################################
warpx.gamma_boost = 10.0
warpx.boost_direction = z

#################################
############ PLASMA #############
#################################
particles.species_names = driver plasma_e plasma_p beam driverback
particles.use_fdtd_nci_corr = 1
particles.rigid_injected_species = driver beam

driver.charge = -q_e
driver.mass = 1.e10
driver.injection_style = "gaussian_beam"
driver.x_rms = 2.e-6
driver.y_rms = 2.e-6
driver.z_rms = 4.e-6
driver.x_m = 0.
driver.y_m = 0.
driver.z_m = -20.e-6
driver.npart = 1000
driver.q_tot = -1.e-9
driver.momentum_distribution_type = "gaussian"
driver.ux_m = 0.0
driver.uy_m = 0.0
driver.uz_m = 200000.
driver.ux_th = 2.
driver.uy_th = 2.
driver.uz_th = 20000.
driver.zinject_plane = 0.
driver.rigid_advance = true

driverback.charge = q_e
driverback.mass = 1.e10
driverback.injection_style = "gaussian_beam"
driverback.x_rms = 2.e-6
driverback.y_rms = 2.e-6
driverback.z_rms = 4.e-6
driverback.x_m = 0.
driverback.y_m = 0.
driverback.z_m = -20.e-6
driverback.npart = 1000
driverback.q_tot = 1.e-9
driverback.momentum_distribution_type = "gaussian"
driverback.ux_m = 0.0
driverback.uy_m = 0.0
driverback.uz_m = 200000.
driverback.ux_th = 2.
driverback.uy_th = 2.
driverback.uz_th = 20000.
driverback.do_backward_propagation = true

plasma_e.charge = -q_e
plasma_e.mass = m_e
plasma_e.injection_style = "NUniformPerCell"
plasma_e.zmin = -100.e-6 # 0.e-6
plasma_e.zmax = 0.2
plasma_e.xmin = -70.e-6
plasma_e.xmax =  70.e-6
plasma_e.ymin = -70.e-6
plasma_e.ymax =  70.e-6
# plasma_e.profile = constant
# plasma_e.density = 1.e23
plasma_e.profile = parse_density_function
plasma_e.density_function(x,y,z) = "(z<lramp)*0.5*(1-cos(pi*z/lramp))*dens+(z>lramp)*dens"
plasma_e.num_particles_per_cell_each_dim = 1 1 1
plasma_e.momentum_distribution_type = "at_rest"
plasma_e.do_continuous_injection = 1

plasma_p.charge = q_e
plasma_p.mass = m_p
plasma_p.injection_style = "NUniformPerCell"
plasma_p.zmin = -100.e-6 # 0.e-6
plasma_p.zmax = 0.2
# plasma_p.profile = "constant"
# plasma_p.density = 1.e23
plasma_p.profile = parse_density_function
plasma_p.density_function(x,y,z) = "(z<lramp)*0.5*(1-cos(pi*z/lramp))*dens+(z>lramp)*dens"
plasma_p.xmin = -70.e-6
plasma_p.xmax =  70.e-6
plasma_p.ymin = -70.e-6
plasma_p.ymax =  70.e-6
plasma_p.num_particles_per_cell_each_dim = 1 1 1
plasma_p.momentum_distribution_type = "at_rest"
plasma_p.do_continuous_injection = 1

beam.charge = -q_e
beam.mass = m_e
beam.injection_style = "gaussian_beam"
beam.x_rms = .5e-6
beam.y_rms = .5e-6
beam.z_rms = 1.e-6
beam.x_m = 0.
beam.y_m = 0.
beam.z_m = -100.e-6
beam.npart = 1000
beam.q_tot = -5.e-10
beam.momentum_distribution_type = "gaussian"
beam.ux_m = 0.0
beam.uy_m = 0.0
beam.uz_m = 2000.
beam.ux_th = 2.
beam.uy_th = 2.
beam.uz_th = 200.
beam.zinject_plane = .8e-3
beam.rigid_advance = true

# Diagnostics
diagnostics.diags_names = diag1
diag1.intervals = 10000
diag1.diag_type = Full
```

## Analyze

#### NOTE
This section is TODO.

## Visualize

#### NOTE
This section is TODO.
