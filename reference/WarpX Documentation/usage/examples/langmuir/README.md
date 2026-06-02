<a id="examples-langmuir"></a>

# Langmuir Waves

These are examples of Plasma oscillations ([Langmuir waves](https://en.wikipedia.org/wiki/Plasma_oscillation)) in a uniform plasma in 1D, 2D, 3D, and RZ.

In each case, a uniform plasma is setup with a sinusoidal perturbation in the electron momentum along each axis.
The plasma is followed for a short period of time, long enough so that E fields develop.
The resulting fields can be compared to the analytic solutions.

## Run

For [MPI-parallel](https://www.mpi-forum.org) runs, prefix these lines with `mpiexec -n 4 ...` or `srun -n 4 ...`, depending on the system.

### 3D

### Python: Script

This example can be run as a **Python** script: `python3 inputs_test_3d_langmuir_multi_picmi.py`.

```python3
#!/usr/bin/env python3
#
# --- Simple example of Langmuir oscillations in a uniform plasma

from pywarpx import picmi

constants = picmi.constants

##########################
# physics parameters
##########################

plasma_density = 1.0e25
plasma_xmin = 0.0
plasma_x_velocity = 0.1 * constants.c

##########################
# numerics parameters
##########################

# --- Number of time steps
max_steps = 40
diagnostic_interval = 10

# --- Grid
nx = 64
ny = 64
nz = 64

xmin = -20.0e-6
ymin = -20.0e-6
zmin = -20.0e-6
xmax = +20.0e-6
ymax = +20.0e-6
zmax = +20.0e-6

number_per_cell_each_dim = [2, 2, 2]

##########################
# physics components
##########################

uniform_plasma = picmi.UniformDistribution(
    density=1.0e25,
    upper_bound=[0.0, None, None],
    directed_velocity=[0.1 * constants.c, 0.0, 0.0],
)

electrons = picmi.Species(
    particle_type="electron", name="electrons", initial_distribution=uniform_plasma
)

##########################
# numerics components
##########################

grid = picmi.Cartesian3DGrid(
    number_of_cells=[nx, ny, nz],
    lower_bound=[xmin, ymin, zmin],
    upper_bound=[xmax, ymax, zmax],
    lower_boundary_conditions=["periodic", "periodic", "periodic"],
    upper_boundary_conditions=["periodic", "periodic", "periodic"],
    moving_window_velocity=[0.0, 0.0, 0.0],
    warpx_max_grid_size=32,
)

solver = picmi.ElectromagneticSolver(grid=grid, cfl=1.0)

##########################
# diagnostics
##########################

field_diag1 = picmi.FieldDiagnostic(
    name="diag1",
    grid=grid,
    period=diagnostic_interval,
    data_list=["Ex", "Jx"],
)

part_diag1 = picmi.ParticleDiagnostic(
    name="diag1",
    period=diagnostic_interval,
    species=[electrons],
    data_list=["weighting", "ux"],
)

##########################
# simulation setup
##########################

sim = picmi.Simulation(
    solver=solver,
    max_steps=max_steps,
    verbose=1,
    warpx_current_deposition_algo="direct",
)

sim.add_species(
    electrons,
    layout=picmi.GriddedLayout(
        n_macroparticle_per_cell=number_per_cell_each_dim, grid=grid
    ),
)

sim.add_diagnostic(field_diag1)
sim.add_diagnostic(part_diag1)

##########################
# simulation run
##########################

# write_inputs will create an inputs file that can be used to run
# with the compiled version.
# sim.write_input_file(file_name = 'inputs_from_PICMI')

# Alternatively, sim.step will run WarpX, controlling it from Python
sim.step()
```

### Executable: Input File

This example can be run as WarpX **executable** using an input file: `warpx.3d inputs_base_3d`.
Check out `Examples/Tests/langmuir/inputs_test_3d_langmuir_multi` for additional input parameters.

```none
# Parameters for the plasma wave
my_constants.max_step = 40
my_constants.lx = 40.e-6 # length of sides
my_constants.dx = 6.25e-07 # grid cell size
my_constants.nx = lx/dx # number of cells in each dimension
my_constants.epsilon = 0.01
my_constants.n0 = 2.e24  # electron and positron densities, #/m^3
my_constants.wp = sqrt(2.*n0*q_e**2/(epsilon0*m_e))  # plasma frequency
my_constants.kp = wp/clight  # plasma wavenumber
my_constants.k = 2.*2.*pi/lx  # perturbation wavenumber
# Note: kp is calculated in SI for a density of 4e24 (i.e. 2e24 electrons + 2e24 positrons)
# k is calculated so as to have 2 periods within the 40e-6 wide box.

# Maximum number of time steps
max_step = max_step

# number of grid points
amr.n_cell =  nx nx nx

# Maximum allowable size of each subdomain in the problem domain;
#    this is used to decompose the domain for parallel calculations.
amr.max_grid_size = nx nx nx

# Maximum level in hierarchy (for now must be 0, i.e., one level in total)
amr.max_level = 0

# Geometry
geometry.dims = 3
geometry.prob_lo     = -lx/2.   -lx/2.   -lx/2.    # physical domain
geometry.prob_hi     =  lx/2.    lx/2.    lx/2.

# Boundary condition
boundary.field_lo = periodic periodic periodic
boundary.field_hi = periodic periodic periodic

warpx.serialize_initial_conditions = 1

# Verbosity
warpx.verbose = 1

# Algorithms
algo.current_deposition = esirkepov
algo.field_gathering = energy-conserving
warpx.use_filter = 0

# Order of particle shape factors
algo.particle_shape = 1

# CFL
warpx.cfl = 1.0

# Particles
particles.species_names = electrons positrons

electrons.charge = -q_e
electrons.mass = m_e
electrons.injection_style = "NUniformPerCell"
electrons.num_particles_per_cell_each_dim = 1 1 1
electrons.xmin = -20.e-6
electrons.xmax =  20.e-6
electrons.ymin = -20.e-6
electrons.ymax = 20.e-6
electrons.zmin = -20.e-6
electrons.zmax = 20.e-6

electrons.profile = constant
electrons.density = n0   # number of electrons per m^3
electrons.momentum_distribution_type = parse_momentum_function
electrons.momentum_function_ux(x,y,z) = "epsilon * k/kp * sin(k*x) * cos(k*y) * cos(k*z)"
electrons.momentum_function_uy(x,y,z) = "epsilon * k/kp * cos(k*x) * sin(k*y) * cos(k*z)"
electrons.momentum_function_uz(x,y,z) = "epsilon * k/kp * cos(k*x) * cos(k*y) * sin(k*z)"

positrons.charge = q_e
positrons.mass = m_e
positrons.injection_style = "NUniformPerCell"
positrons.num_particles_per_cell_each_dim = 1 1 1
positrons.xmin = -20.e-6
positrons.xmax =  20.e-6
positrons.ymin = -20.e-6
positrons.ymax = 20.e-6
positrons.zmin = -20.e-6
positrons.zmax = 20.e-6

positrons.profile = constant
positrons.density = n0   # number of positrons per m^3
positrons.momentum_distribution_type = parse_momentum_function
positrons.momentum_function_ux(x,y,z) = "-epsilon * k/kp * sin(k*x) * cos(k*y) * cos(k*z)"
positrons.momentum_function_uy(x,y,z) = "-epsilon * k/kp * cos(k*x) * sin(k*y) * cos(k*z)"
positrons.momentum_function_uz(x,y,z) = "-epsilon * k/kp * cos(k*x) * cos(k*y) * sin(k*z)"

# Diagnostics
warpx.synchronize_velocity_for_diagnostics = 1
diagnostics.diags_names = diag1 openpmd
diag1.intervals = max_step
diag1.diag_type = Full
diag1.fields_to_plot = Ex Ey Ez Bx By Bz jx jy jz part_per_cell rho divE
diag1.electrons.variables = x y z w ux
diag1.positrons.variables = x y z uz

openpmd.intervals = 40
openpmd.diag_type = Full
openpmd.format = openpmd
openpmd.electrons.additional_variables = Ex Ey Ez
```

### 2D

### Python: Script

This example can be run as a **Python** script: `python3 inputs_test_2d_langmuir_multi_picmi.py`.

```python3
#!/usr/bin/env python3
#
# --- Simple example of Langmuir oscillations in a uniform plasma
# --- in two dimensions

from pywarpx import picmi

constants = picmi.constants

##########################
# physics parameters
##########################

plasma_density = 1.0e25
plasma_xmin = 0.0
plasma_x_velocity = 0.1 * constants.c

##########################
# numerics parameters
##########################

# --- Number of time steps
max_steps = 40
diagnostic_intervals = "::10"

# --- Grid
nx = 64
nz = 64

xmin = -20.0e-6
zmin = -20.0e-6
xmax = +20.0e-6
zmax = +20.0e-6

number_per_cell_each_dim = [2, 2]

##########################
# physics components
##########################

uniform_plasma = picmi.UniformDistribution(
    density=1.0e25,
    upper_bound=[0.0, None, None],
    directed_velocity=[0.1 * constants.c, 0.0, 0.0],
)

electrons = picmi.Species(
    particle_type="electron", name="electrons", initial_distribution=uniform_plasma
)

##########################
# numerics components
##########################

grid = picmi.Cartesian2DGrid(
    number_of_cells=[nx, nz],
    lower_bound=[xmin, zmin],
    upper_bound=[xmax, zmax],
    lower_boundary_conditions=["periodic", "periodic"],
    upper_boundary_conditions=["periodic", "periodic"],
    moving_window_velocity=[0.0, 0.0, 0.0],
    warpx_max_grid_size=32,
)

solver = picmi.ElectromagneticSolver(grid=grid, cfl=1.0)

##########################
# diagnostics
##########################

field_diag1 = picmi.FieldDiagnostic(
    name="diag1",
    grid=grid,
    period=diagnostic_intervals,
    data_list=["Ex", "Jx"],
)

part_diag1 = picmi.ParticleDiagnostic(
    name="diag1",
    period=diagnostic_intervals,
    species=[electrons],
    data_list=["weighting", "ux"],
)

##########################
# simulation setup
##########################

sim = picmi.Simulation(
    solver=solver,
    max_steps=max_steps,
    verbose=1,
    warpx_current_deposition_algo="direct",
    warpx_use_filter=0,
)

sim.add_species(
    electrons,
    layout=picmi.GriddedLayout(
        n_macroparticle_per_cell=number_per_cell_each_dim, grid=grid
    ),
)

sim.add_diagnostic(field_diag1)
sim.add_diagnostic(part_diag1)

##########################
# simulation run
##########################

# write_inputs will create an inputs file that can be used to run
# with the compiled version.
sim.write_input_file(file_name="inputs2d_from_PICMI")

# Alternatively, sim.step will run WarpX, controlling it from Python
sim.step()
```

### Executable: Input File

This example can be run as WarpX **executable** using an input file: `warpx.2d inputs_base_2d`
Check out `Examples/Tests/langmuir/inputs_test_2d_langmuir_multi` for additional input parameters.

```none
# Maximum number of time steps
max_step = 80

# number of grid points
amr.n_cell =   128  128

# Maximum allowable size of each subdomain in the problem domain;
#    this is used to decompose the domain for parallel calculations.
amr.max_grid_size = 64

# Maximum level in hierarchy (for now must be 0, i.e., one level in total)
amr.max_level = 0

# Geometry
geometry.dims = 2
geometry.prob_lo     = -20.e-6   -20.e-6    # physical domain
geometry.prob_hi     =  20.e-6    20.e-6

# Boundary condition
boundary.field_lo = periodic periodic
boundary.field_hi = periodic periodic

warpx.serialize_initial_conditions = 1

# Verbosity
warpx.verbose = 1

# Algorithms
algo.field_gathering = energy-conserving
warpx.use_filter = 0

# Order of particle shape factors
algo.particle_shape = 1

# CFL
warpx.cfl = 1.0

# Parameters for the plasma wave
my_constants.epsilon = 0.01
my_constants.n0 = 2.e24  # electron and positron densities, #/m^3
my_constants.wp = sqrt(2.*n0*q_e**2/(epsilon0*m_e))  # plasma frequency
my_constants.kp = wp/clight  # plasma wavenumber
my_constants.k = 2.*pi/20.e-6  # perturbation wavenumber
# Note: kp is calculated in SI for a density of 4e24 (i.e. 2e24 electrons + 2e24 positrons)
# k is calculated so as to have 2 periods within the 40e-6 wide box.

# Particles
particles.species_names = electrons positrons

electrons.charge = -q_e
electrons.mass = m_e
electrons.injection_style = "NUniformPerCell"
electrons.num_particles_per_cell_each_dim = 2 2
electrons.xmin = -20.e-6
electrons.xmax =  20.e-6
electrons.ymin = -20.e-6
electrons.ymax = 20.e-6
electrons.zmin = -20.e-6
electrons.zmax = 20.e-6

electrons.profile = constant
electrons.density = n0   # number of electrons per m^3
electrons.momentum_distribution_type = parse_momentum_function
electrons.momentum_function_ux(x,y,z) = "epsilon * k/kp * sin(k*x) * cos(k*y) * cos(k*z)"
electrons.momentum_function_uy(x,y,z) = "epsilon * k/kp * cos(k*x) * sin(k*y) * cos(k*z)"
electrons.momentum_function_uz(x,y,z) = "epsilon * k/kp * cos(k*x) * cos(k*y) * sin(k*z)"

positrons.charge = q_e
positrons.mass = m_e
positrons.injection_style = "NUniformPerCell"
positrons.num_particles_per_cell_each_dim = 2 2
positrons.xmin = -20.e-6
positrons.xmax =  20.e-6
positrons.ymin = -20.e-6
positrons.ymax = 20.e-6
positrons.zmin = -20.e-6
positrons.zmax = 20.e-6

positrons.profile = constant
positrons.density = n0   # number of positrons per m^3
positrons.momentum_distribution_type = parse_momentum_function
positrons.momentum_function_ux(x,y,z) = "-epsilon * k/kp * sin(k*x) * cos(k*y) * cos(k*z)"
positrons.momentum_function_uy(x,y,z) = "-epsilon * k/kp * cos(k*x) * sin(k*y) * cos(k*z)"
positrons.momentum_function_uz(x,y,z) = "-epsilon * k/kp * cos(k*x) * cos(k*y) * sin(k*z)"

# Diagnostics
warpx.synchronize_velocity_for_diagnostics = 1
diagnostics.diags_names = diag1
diag1.intervals = 40
diag1.diag_type = Full
```

### RZ

### Python: Script

This example can be run as a **Python** script: `python3 inputs_test_rz_langmuir_multi_picmi.py`.

```python3
#!/usr/bin/env python3
#
# This is a script that analyses the multimode simulation results.
# This simulates a RZ multimode periodic plasma wave.
# The electric field from the simulation is compared to the analytic value

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pywarpx import picmi

constants = picmi.constants

##########################
# physics parameters
##########################

density = 2.0e24
epsilon0 = 0.001 * constants.c
epsilon1 = 0.001 * constants.c
epsilon2 = 0.001 * constants.c
w0 = 5.0e-6
n_osc_z = 3

# Plasma frequency
wp = np.sqrt((density * constants.q_e**2) / (constants.m_e * constants.ep0))
kp = wp / constants.c

##########################
# numerics parameters
##########################

nr = 64
nz = 200

rmin = 0.0e0
zmin = 0.0e0
rmax = +20.0e-6
zmax = +40.0e-6

# Wave vector of the wave
k0 = 2.0 * np.pi * n_osc_z / (zmax - zmin)

diagnostic_intervals = 40

##########################
# physics components
##########################

uniform_plasma = picmi.UniformDistribution(
    density=density,
    upper_bound=[+18e-6, None, +40e-6],
    directed_velocity=[0.0, 0.0, 0.0],
)

momentum_expressions = [
    """+ epsilon0/kp*2*x/w0**2*exp(-(x**2+y**2)/w0**2)*sin(k0*z)
                           - epsilon1/kp*2/w0*exp(-(x**2+y**2)/w0**2)*sin(k0*z)
                           + epsilon1/kp*4*x**2/w0**3*exp(-(x**2+y**2)/w0**2)*sin(k0*z)
                           - epsilon2/kp*8*x/w0**2*exp(-(x**2+y**2)/w0**2)*sin(k0*z)
                           + epsilon2/kp*8*x*(x**2-y**2)/w0**4*exp(-(x**2+y**2)/w0**2)*sin(k0*z)""",
    """+ epsilon0/kp*2*y/w0**2*exp(-(x**2+y**2)/w0**2)*sin(k0*z)
                           + epsilon1/kp*4*x*y/w0**3*exp(-(x**2+y**2)/w0**2)*sin(k0*z)
                           + epsilon2/kp*8*y/w0**2*exp(-(x**2+y**2)/w0**2)*sin(k0*z)
                           + epsilon2/kp*8*y*(x**2-y**2)/w0**4*exp(-(x**2+y**2)/w0**2)*sin(k0*z)""",
    """- epsilon0/kp*k0*exp(-(x**2+y**2)/w0**2)*cos(k0*z)
                           - epsilon1/kp*k0*2*x/w0*exp(-(x**2+y**2)/w0**2)*cos(k0*z)
                           - epsilon2/kp*k0*4*(x**2-y**2)/w0**2*exp(-(x**2+y**2)/w0**2)*cos(k0*z)""",
]

analytic_plasma = picmi.AnalyticDistribution(
    density_expression=density,
    upper_bound=[+18e-6, None, +40e-6],
    epsilon0=epsilon0,
    epsilon1=epsilon1,
    epsilon2=epsilon2,
    kp=kp,
    k0=k0,
    w0=w0,
    momentum_expressions=momentum_expressions,
)

electrons = picmi.Species(
    particle_type="electron", name="electrons", initial_distribution=analytic_plasma
)
protons = picmi.Species(
    particle_type="proton", name="protons", initial_distribution=uniform_plasma
)

##########################
# numerics components
##########################

grid = picmi.CylindricalGrid(
    number_of_cells=[nr, nz],
    n_azimuthal_modes=3,
    lower_bound=[rmin, zmin],
    upper_bound=[rmax, zmax],
    lower_boundary_conditions=["none", "periodic"],
    upper_boundary_conditions=["none", "periodic"],
    lower_boundary_conditions_particles=["none", "periodic"],
    upper_boundary_conditions_particles=["absorbing", "periodic"],
    moving_window_velocity=[0.0, 0.0],
    warpx_max_grid_size=64,
)

solver = picmi.ElectromagneticSolver(grid=grid, cfl=1.0)

##########################
# diagnostics
##########################

field_diag1 = picmi.FieldDiagnostic(
    name="diag1",
    grid=grid,
    period=diagnostic_intervals,
    data_list=["Er", "Ez", "Bt", "Jr", "Jz", "part_per_cell"],
)

part_diag1 = picmi.ParticleDiagnostic(
    name="diag1",
    period=diagnostic_intervals,
    species=[electrons],
    data_list=["weighting", "momentum"],
)

##########################
# simulation setup
##########################

sim = picmi.Simulation(
    solver=solver,
    max_steps=40,
    verbose=1,
    warpx_current_deposition_algo="esirkepov",
    warpx_field_gathering_algo="energy-conserving",
    warpx_particle_pusher_algo="boris",
    warpx_use_filter=0,
)

sim.add_species(
    electrons,
    layout=picmi.GriddedLayout(n_macroparticle_per_cell=[2, 16, 2], grid=grid),
)
sim.add_species(
    protons, layout=picmi.GriddedLayout(n_macroparticle_per_cell=[2, 16, 2], grid=grid)
)

sim.add_diagnostic(field_diag1)
sim.add_diagnostic(part_diag1)

##########################
# simulation run
##########################

# write_inputs will create an inputs file that can be used to run
# with the compiled version.
# sim.write_input_file(file_name='inputsrz_from_PICMI')

# Alternatively, sim.step will run WarpX, controlling it from Python
sim.step()


# Below is WarpX specific code to check the results.


def calcEr(z, r, k0, w0, wp, t, epsilons):
    """
    Return the radial electric field as an array
    of the same length as z and r, in the half-plane theta=0
    """
    Er_array = (
        epsilons[0]
        * constants.m_e
        * constants.c
        / constants.q_e
        * 2
        * r
        / w0**2
        * np.exp(-(r**2) / w0**2)
        * np.sin(k0 * z)
        * np.sin(wp * t)
        - epsilons[1]
        * constants.m_e
        * constants.c
        / constants.q_e
        * 2
        / w0
        * np.exp(-(r**2) / w0**2)
        * np.sin(k0 * z)
        * np.sin(wp * t)
        + epsilons[1]
        * constants.m_e
        * constants.c
        / constants.q_e
        * 4
        * r**2
        / w0**3
        * np.exp(-(r**2) / w0**2)
        * np.sin(k0 * z)
        * np.sin(wp * t)
        - epsilons[2]
        * constants.m_e
        * constants.c
        / constants.q_e
        * 8
        * r
        / w0**2
        * np.exp(-(r**2) / w0**2)
        * np.sin(k0 * z)
        * np.sin(wp * t)
        + epsilons[2]
        * constants.m_e
        * constants.c
        / constants.q_e
        * 8
        * r**3
        / w0**4
        * np.exp(-(r**2) / w0**2)
        * np.sin(k0 * z)
        * np.sin(wp * t)
    )
    return Er_array


def calcEz(z, r, k0, w0, wp, t, epsilons):
    """
    Return the longitudinal electric field as an array
    of the same length as z and r, in the half-plane theta=0
    """
    Ez_array = (
        -epsilons[0]
        * constants.m_e
        * constants.c
        / constants.q_e
        * k0
        * np.exp(-(r**2) / w0**2)
        * np.cos(k0 * z)
        * np.sin(wp * t)
        - epsilons[1]
        * constants.m_e
        * constants.c
        / constants.q_e
        * k0
        * 2
        * r
        / w0
        * np.exp(-(r**2) / w0**2)
        * np.cos(k0 * z)
        * np.sin(wp * t)
        - epsilons[2]
        * constants.m_e
        * constants.c
        / constants.q_e
        * k0
        * 4
        * r**2
        / w0**2
        * np.exp(-(r**2) / w0**2)
        * np.cos(k0 * z)
        * np.sin(wp * t)
    )
    return Ez_array


# Current time of the simulation
t0 = sim.extension.warpx.gett_new(0)

# Get the raw field data. Note that these are the real and imaginary
# parts of the fields for each azimuthal mode.
Er_sim_wrap = sim.fields.get("Efield_aux", dir="r", level=0)
Ez_sim_wrap = sim.fields.get("Efield_aux", dir="z", level=0)
Er_sim_modes = Er_sim_wrap[...]
Ez_sim_modes = Ez_sim_wrap[...]

rr_Er = Er_sim_wrap.mesh("r")
zz_Er = Er_sim_wrap.mesh("z")
rr_Ez = Ez_sim_wrap.mesh("r")
zz_Ez = Ez_sim_wrap.mesh("z")

rr_Er = rr_Er[:, np.newaxis] * np.ones(zz_Er.shape[0])[np.newaxis, :]
zz_Er = zz_Er[np.newaxis, :] * np.ones(rr_Er.shape[0])[:, np.newaxis]
rr_Ez = rr_Ez[:, np.newaxis] * np.ones(zz_Ez.shape[0])[np.newaxis, :]
zz_Ez = zz_Ez[np.newaxis, :] * np.ones(rr_Ez.shape[0])[:, np.newaxis]

# Sum the real components to get the field along x-axis (theta = 0)
Er_sim = Er_sim_modes[:, :, 0] + np.sum(Er_sim_modes[:, :, 1::2], axis=2)
Ez_sim = Ez_sim_modes[:, :, 0] + np.sum(Ez_sim_modes[:, :, 1::2], axis=2)

# The analytical solutions
Er_th = calcEr(zz_Er, rr_Er, k0, w0, wp, t0, [epsilon0, epsilon1, epsilon2])
Ez_th = calcEz(zz_Ez, rr_Ez, k0, w0, wp, t0, [epsilon0, epsilon1, epsilon2])

max_error_Er = abs(Er_sim - Er_th).max() / abs(Er_th).max()
max_error_Ez = abs(Ez_sim - Ez_th).max() / abs(Ez_th).max()
print("Max error Er %e" % max_error_Er)
print("Max error Ez %e" % max_error_Ez)

# Plot the last field from the loop (Er at iteration 40)
fig, ax = plt.subplots(3)
im = ax[0].imshow(Er_sim, aspect="auto", origin="lower")
fig.colorbar(im, ax=ax[0], orientation="vertical")
ax[0].set_title("Er, last iteration (simulation)")
ax[1].imshow(Er_th, aspect="auto", origin="lower")
fig.colorbar(im, ax=ax[1], orientation="vertical")
ax[1].set_title("Er, last iteration (theory)")
im = ax[2].imshow((Er_sim - Er_th) / abs(Er_th).max(), aspect="auto", origin="lower")
fig.colorbar(im, ax=ax[2], orientation="vertical")
ax[2].set_title("Er, last iteration (difference)")
plt.savefig("langmuir_multi_rz_multimode_analysis_Er.png")

fig, ax = plt.subplots(3)
im = ax[0].imshow(Ez_sim, aspect="auto", origin="lower")
fig.colorbar(im, ax=ax[0], orientation="vertical")
ax[0].set_title("Ez, last iteration (simulation)")
ax[1].imshow(Ez_th, aspect="auto", origin="lower")
fig.colorbar(im, ax=ax[1], orientation="vertical")
ax[1].set_title("Ez, last iteration (theory)")
im = ax[2].imshow((Ez_sim - Ez_th) / abs(Ez_th).max(), aspect="auto", origin="lower")
fig.colorbar(im, ax=ax[2], orientation="vertical")
ax[2].set_title("Ez, last iteration (difference)")
plt.savefig("langmuir_multi_rz_multimode_analysis_Ez.png")

assert max(max_error_Er, max_error_Ez) < 0.02
```

### Executable: Input File

This example can be run as WarpX **executable** using an input file: `warpx.rz inputs_base_rz`
Check out `Examples/Tests/langmuir/inputs_test_rz_langmuir_multi` for additional input parameters.

```none
# Parameters for the plasma wave
my_constants.max_step = 80
my_constants.epsilon = 0.01
my_constants.n0 = 2.e24  # electron density, #/m^3
my_constants.wp = sqrt(n0*q_e**2/(epsilon0*m_e))  # plasma frequency
my_constants.kp = wp/clight  # plasma wavenumber
my_constants.k0 = 2.*pi/20.e-6  # longitudianl perturbation wavenumber
my_constants.w0 = 5.e-6  # transverse perturbation length
# Note: kp is calculated in SI for a density of 2e24
# k0 is calculated so as to have 2 periods within the 40e-6 wide box.

# Maximum number of time steps
max_step = max_step

# number of grid points
amr.n_cell =   64  128

# Maximum allowable size of each subdomain in the problem domain;
#    this is used to decompose the domain for parallel calculations.
amr.max_grid_size = 64

# Maximum level in hierarchy (for now must be 0, i.e., one level in total)
amr.max_level = 0

# Geometry
geometry.dims = RZ
geometry.prob_lo     =   0.e-6   -20.e-6    # physical domain
geometry.prob_hi     =  20.e-6    20.e-6
boundary.field_lo = none periodic
boundary.field_hi = none periodic
boundary.particle_lo = none periodic
boundary.particle_hi = absorbing periodic

warpx.serialize_initial_conditions = 1

# Verbosity
warpx.verbose = 1

# Algorithms
algo.field_gathering = energy-conserving
algo.current_deposition = esirkepov
warpx.use_filter = 0

# Order of particle shape factors
algo.particle_shape = 1

# CFL
warpx.cfl = 1.0

# Having this turned on makes for a more sensitive test
warpx.do_dive_cleaning = 1

# Particles
particles.species_names = electrons ions

electrons.charge = -q_e
electrons.mass = m_e
electrons.injection_style = "NUniformPerCell"
electrons.num_particles_per_cell_each_dim = 2 2 2
electrons.xmin =   0.e-6
electrons.xmax =  18.e-6
electrons.zmin = -20.e-6
electrons.zmax = +20.e-6

electrons.profile = constant
electrons.density = n0   # number of electrons per m^3
electrons.momentum_distribution_type = parse_momentum_function
electrons.momentum_function_ux(x,y,z) = "epsilon/kp*2*x/w0**2*exp(-(x**2+y**2)/w0**2)*sin(k0*z)"
electrons.momentum_function_uy(x,y,z) = "epsilon/kp*2*y/w0**2*exp(-(x**2+y**2)/w0**2)*sin(k0*z)"
electrons.momentum_function_uz(x,y,z) = "-epsilon/kp*k0*exp(-(x**2+y**2)/w0**2)*cos(k0*z)"


ions.charge = q_e
ions.mass = m_p
ions.injection_style = "NUniformPerCell"
ions.num_particles_per_cell_each_dim = 2 2 2
ions.xmin =   0.e-6
ions.xmax =  18.e-6
ions.zmin = -20.e-6
ions.zmax = +20.e-6

ions.profile = constant
ions.density = n0   # number of ions per m^3
ions.momentum_distribution_type = at_rest

# Diagnostics
warpx.synchronize_velocity_for_diagnostics = 1
diagnostics.diags_names = diag1 diag_parser_filter diag_uniform_filter diag_random_filter
diag1.intervals = max_step/2
diag1.diag_type = Full
diag1.fields_to_plot = jr jz Er Ez Bt

## diag_parser_filter is a diag used to test the particle filter function.
diag_parser_filter.intervals = max_step:max_step:
diag_parser_filter.diag_type = Full
diag_parser_filter.species = electrons
diag_parser_filter.electrons.plot_filter_function(t,x,y,z,ux,uy,uz) = "(uy-uz < 0) *
                                                                 (sqrt(x**2+y**2)<10e-6) * (z > 0)"

## diag_uniform_filter is a diag used to test the particle uniform filter.
diag_uniform_filter.intervals = max_step:max_step:
diag_uniform_filter.diag_type = Full
diag_uniform_filter.species = electrons
diag_uniform_filter.electrons.uniform_stride = 3

## diag_random_filter is a diag used to test the particle random filter.
diag_random_filter.intervals = max_step:max_step:
diag_random_filter.diag_type = Full
diag_random_filter.species = electrons
diag_random_filter.electrons.random_fraction = 0.66
```

### 1D

### Python: Script

#### NOTE
TODO: This input file should be created, like the `inputs_test_1d_langmuir_multi` file.

### Executable: Input File

This example can be run as WarpX **executable** using an input file: `warpx.1d inputs_test_1d_langmuir_multi`

```none
# Maximum number of time steps
max_step = 80

# number of grid points
amr.n_cell =  128

# Maximum allowable size of each subdomain in the problem domain;
#    this is used to decompose the domain for parallel calculations.
amr.max_grid_size = 64

# Maximum level in hierarchy (for now must be 0, i.e., one level in total)
amr.max_level = 0

# Geometry
geometry.dims = 1
geometry.prob_lo     = -20.e-6    # physical domain
geometry.prob_hi     =  20.e-6

# Boundary condition
boundary.field_lo = periodic
boundary.field_hi = periodic

warpx.serialize_initial_conditions = 1

# Verbosity
warpx.verbose = 1

# Algorithms
algo.field_gathering = energy-conserving
algo.current_deposition = esirkepov
warpx.use_filter = 0

# Order of particle shape factors
algo.particle_shape = 1

# CFL
warpx.cfl = 0.8

# Parameters for the plasma wave
my_constants.epsilon = 0.01
my_constants.n0 = 2.e24  # electron and positron densities, #/m^3
my_constants.wp = sqrt(2.*n0*q_e**2/(epsilon0*m_e))  # plasma frequency
my_constants.kp = wp/clight  # plasma wavenumber
my_constants.k = 2.*pi/20.e-6  # perturbation wavenumber
# Note: kp is calculated in SI for a density of 4e24 (i.e. 2e24 electrons + 2e24 positrons)
# k is calculated so as to have 2 periods within the 40e-6 wide box.

# Particles
particles.species_names = electrons positrons

electrons.charge = -q_e
electrons.mass = m_e
electrons.injection_style = "NUniformPerCell"
electrons.num_particles_per_cell_each_dim = 2
electrons.zmin = -20.e-6
electrons.zmax = 20.e-6

electrons.profile = constant
electrons.density = n0   # number of electrons per m^3
electrons.momentum_distribution_type = parse_momentum_function
electrons.momentum_function_ux(x,y,z) = "epsilon * k/kp * sin(k*x) * cos(k*y) * cos(k*z)"
electrons.momentum_function_uy(x,y,z) = "epsilon * k/kp * cos(k*x) * sin(k*y) * cos(k*z)"
electrons.momentum_function_uz(x,y,z) = "epsilon * k/kp * cos(k*x) * cos(k*y) * sin(k*z)"

positrons.charge = q_e
positrons.mass = m_e
positrons.injection_style = "NUniformPerCell"
positrons.num_particles_per_cell_each_dim = 2
positrons.zmin = -20.e-6
positrons.zmax = 20.e-6

positrons.profile = constant
positrons.density = n0   # number of positrons per m^3
positrons.momentum_distribution_type = parse_momentum_function
positrons.momentum_function_ux(x,y,z) = "-epsilon * k/kp * sin(k*x) * cos(k*y) * cos(k*z)"
positrons.momentum_function_uy(x,y,z) = "-epsilon * k/kp * cos(k*x) * sin(k*y) * cos(k*z)"
positrons.momentum_function_uz(x,y,z) = "-epsilon * k/kp * cos(k*x) * cos(k*y) * sin(k*z)"

# Diagnostics
diagnostics.diags_names = diag1 openpmd
diag1.intervals = 40
diag1.diag_type = Full
diag1.fields_to_plot = Bx By Bz Ex Ey Ez jx jy jz rho divE
diag1.electrons.variables = z w ux uy uz
diag1.positrons.variables = z w ux uy uz

openpmd.intervals = 40
openpmd.diag_type = Full
openpmd.format = openpmd
```

## Analyze

We run the following script to analyze correctness:

### 3D

### Script `analysis_3d.py`

```python3
#!/usr/bin/env python3

# Copyright 2019-2022 Jean-Luc Vay, Maxence Thevenet, Remi Lehe, Axel Huebl
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL
#
# This is a script that analyses the simulation results from
# the script `inputs.multi.rt`. This simulates a 3D periodic plasma wave.
# The electric field in the simulation is given (in theory) by:
# $$ E_x = \epsilon \,\frac{m_e c^2 k_x}{q_e}\sin(k_x x)\cos(k_y y)\cos(k_z z)\sin( \omega_p t)$$
# $$ E_y = \epsilon \,\frac{m_e c^2 k_y}{q_e}\cos(k_x x)\sin(k_y y)\cos(k_z z)\sin( \omega_p t)$$
# $$ E_z = \epsilon \,\frac{m_e c^2 k_z}{q_e}\cos(k_x x)\cos(k_y y)\sin(k_z z)\sin( \omega_p t)$$
import os
import re
import sys

import matplotlib.pyplot as plt
import yt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

yt.funcs.mylog.setLevel(50)

import numpy as np
from analysis_utils import check_charge_conservation
from openpmd_viewer import OpenPMDTimeSeries
from scipy.constants import c, e, epsilon_0, m_e

# test name
test_name = os.path.split(os.getcwd())[1]

# this will be the name of the plot file
fn = sys.argv[1]

# Parse test name and check if div(E)/div(B) cleaning (warpx.do_div<e,b>_cleaning=1) is used
div_cleaning = True if re.search("div_cleaning", test_name) else False

# Parameters (these parameters must match the parameters in `inputs.multi.rt`)
epsilon = 0.01
n = 4.0e24
n_osc_x = 2
n_osc_y = 2
n_osc_z = 2
lo = [-20.0e-6, -20.0e-6, -20.0e-6]
hi = [20.0e-6, 20.0e-6, 20.0e-6]
Ncell = [64, 64, 64]

# Wave vector of the wave
kx = 2.0 * np.pi * n_osc_x / (hi[0] - lo[0])
ky = 2.0 * np.pi * n_osc_y / (hi[1] - lo[1])
kz = 2.0 * np.pi * n_osc_z / (hi[2] - lo[2])
# Plasma frequency
wp = np.sqrt((n * e**2) / (m_e * epsilon_0))

k = {"Ex": kx, "Ey": ky, "Ez": kz}
cos = {"Ex": (0, 1, 1), "Ey": (1, 0, 1), "Ez": (1, 1, 0)}


def get_contribution_at_positions(is_cos, k, idim, u):
    if is_cos[idim] == 1:
        return np.cos(k * u)
    else:
        return np.sin(k * u)


def get_contribution(is_cos, k, idim):
    du = (hi[idim] - lo[idim]) / Ncell[idim]
    u = lo[idim] + du * (0.5 + np.arange(Ncell[idim]))
    return get_contribution_at_positions(is_cos, k, idim, u)


def get_theoretical_field(field, t):
    amplitude = epsilon * (m_e * c**2 * k[field]) / e * np.sin(wp * t)
    cos_flag = cos[field]
    x_contribution = get_contribution(cos_flag, kx, 0)
    y_contribution = get_contribution(cos_flag, ky, 1)
    z_contribution = get_contribution(cos_flag, kz, 2)

    E = (
        amplitude
        * x_contribution[:, np.newaxis, np.newaxis]
        * y_contribution[np.newaxis, :, np.newaxis]
        * z_contribution[np.newaxis, np.newaxis, :]
    )

    return E


def get_theoretical_field_at_positions(field, t, x, y, z):
    amplitude = epsilon * (m_e * c**2 * k[field]) / e * np.sin(wp * t)
    cos_flag = cos[field]
    x_contribution = get_contribution_at_positions(cos_flag, kx, 0, x)
    y_contribution = get_contribution_at_positions(cos_flag, ky, 1, y)
    z_contribution = get_contribution_at_positions(cos_flag, kz, 2, z)

    E = amplitude * x_contribution * y_contribution * z_contribution

    return E


# Read the file
ds = yt.load(fn)

# Check that the particle selective output worked:
species = "electrons"
print("ds.field_list", ds.field_list)
for field in ["particle_weight", "particle_momentum_x"]:
    print("assert that this is in ds.field_list", (species, field))
    assert (species, field) in ds.field_list
for field in ["particle_momentum_y", "particle_momentum_z"]:
    print("assert that this is NOT in ds.field_list", (species, field))
    assert (species, field) not in ds.field_list
species = "positrons"
for field in ["particle_momentum_x", "particle_momentum_y"]:
    print("assert that this is NOT in ds.field_list", (species, field))
    assert (species, field) not in ds.field_list

t0 = ds.current_time.to_value()
data = ds.covering_grid(
    level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions
)
edge = np.array(
    [
        (ds.domain_left_edge[2]).item(),
        (ds.domain_right_edge[2]).item(),
        (ds.domain_left_edge[0]).item(),
        (ds.domain_right_edge[0]).item(),
    ]
)

# Check the validity of the fields
error_rel = 0
for field in ["Ex", "Ey", "Ez"]:
    E_sim = data[("mesh", field)].to_ndarray()
    E_th = get_theoretical_field(field, t0)
    max_error = abs(E_sim - E_th).max() / abs(E_th).max()
    print("%s: Max error: %.2e" % (field, max_error))
    error_rel = max(error_rel, max_error)

ts = OpenPMDTimeSeries("./diags/openpmd")
x, y, z = ts.get_particle(["x", "y", "z"], species="electrons", iteration=40)
for field in ["Ex", "Ey", "Ez"]:
    E_sim_particles = ts.get_particle(
        [field.lower()], species="electrons", iteration=40
    )
    E_th_particles = get_theoretical_field_at_positions(field, t0, x, y, z)
    max_error = abs(E_sim_particles - E_th_particles).max() / abs(E_th_particles).max()
    print("%s: Max error at particles: %.2e" % (field, max_error))
    error_rel = max(error_rel, max_error)

# Plot the last field from the loop (Ez at iteration 40)
fig, (ax1, ax2) = plt.subplots(1, 2, dpi=100)
# First plot (slice at y=0)
E_plot = E_sim[:, Ncell[1] // 2 + 1, :]
vmin = E_plot.min()
vmax = E_plot.max()
cax1 = make_axes_locatable(ax1).append_axes("right", size="5%", pad="5%")
im1 = ax1.imshow(E_plot, origin="lower", extent=edge, vmin=vmin, vmax=vmax)
cb1 = fig.colorbar(im1, cax=cax1)
ax1.set_xlabel(r"$z$")
ax1.set_ylabel(r"$x$")
ax1.set_title(r"$E_z$ (sim)")
# Second plot (slice at y=0)
E_plot = E_th[:, Ncell[1] // 2 + 1, :]
vmin = E_plot.min()
vmax = E_plot.max()
cax2 = make_axes_locatable(ax2).append_axes("right", size="5%", pad="5%")
im2 = ax2.imshow(E_plot, origin="lower", extent=edge, vmin=vmin, vmax=vmax)
cb2 = fig.colorbar(im2, cax=cax2)
ax2.set_xlabel(r"$z$")
ax2.set_ylabel(r"$x$")
ax2.set_title(r"$E_z$ (theory)")
# Save figure
fig.tight_layout()
fig.savefig("Langmuir_multi_analysis.png", dpi=200)

tolerance_rel = 5e-2

print("error_rel    : " + str(error_rel))
print("tolerance_rel: " + str(tolerance_rel))

assert error_rel < tolerance_rel

# Additional check on charge conservation for certain cases
# (e.g., Esirkepov or Vay deposition, current correction)
check_charge_conservation(data)

if div_cleaning:
    ds_old = yt.load("diags/diag1000038")
    ds_mid = yt.load("diags/diag1000039")
    ds_new = yt.load(fn)  # this is the last plotfile

    ad_old = ds_old.covering_grid(
        level=0, left_edge=ds_old.domain_left_edge, dims=ds_old.domain_dimensions
    )
    ad_mid = ds_mid.covering_grid(
        level=0, left_edge=ds_mid.domain_left_edge, dims=ds_mid.domain_dimensions
    )
    ad_new = ds_new.covering_grid(
        level=0, left_edge=ds_new.domain_left_edge, dims=ds_new.domain_dimensions
    )

    rho = ad_mid["rho"].v.squeeze()
    divE = ad_mid["divE"].v.squeeze()
    F_old = ad_old["F"].v.squeeze()
    F_new = ad_new["F"].v.squeeze()

    # Check max norm of error on dF/dt = div(E) - rho/epsilon_0
    # (the time interval between the old and new data is 2*dt)
    dt = 1.203645751e-15
    x = F_new - F_old
    y = (divE - rho / epsilon_0) * 2 * dt
    error_rel = np.amax(np.abs(x - y)) / np.amax(np.abs(y))
    tolerance = 1e-2
    print("Check div(E) cleaning:")
    print("error_rel = {}".format(error_rel))
    print("tolerance = {}".format(tolerance))
    assert error_rel < tolerance
```

### 2D

### Script `analysis_2d.py`

```python3
#!/usr/bin/env python3

# Copyright 2019 Jean-Luc Vay, Maxence Thevenet, Remi Lehe
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL
#
# This is a script that analyses the simulation results from
# the script `inputs.multi.rt`. This simulates a 3D periodic plasma wave.
# The electric field in the simulation is given (in theory) by:
# $$ E_x = \epsilon \,\frac{m_e c^2 k_x}{q_e}\sin(k_x x)\cos(k_y y)\cos(k_z z)\sin( \omega_p t)$$
# $$ E_y = \epsilon \,\frac{m_e c^2 k_y}{q_e}\cos(k_x x)\sin(k_y y)\cos(k_z z)\sin( \omega_p t)$$
# $$ E_z = \epsilon \,\frac{m_e c^2 k_z}{q_e}\cos(k_x x)\cos(k_y y)\sin(k_z z)\sin( \omega_p t)$$
import os
import re
import sys

import matplotlib.pyplot as plt
import yt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

yt.funcs.mylog.setLevel(50)

import numpy as np
from analysis_utils import check_charge_conservation
from scipy.constants import c, e, epsilon_0, m_e

# test name
test_name = os.path.split(os.getcwd())[1]

# this will be the name of the plot file
fn = sys.argv[1]

# Parse test name and check if particle_shape = 4 is used
particle_shape_4 = True if re.search("particle_shape_4", test_name) else False

# Parameters (must match the parameters in the inputs)
# FIXME read these parameters from warpx_used_inputs
epsilon = 0.01
n = 4.0e24
n_osc_x = 2
n_osc_z = 2
xmin = -20e-6
xmax = 20.0e-6
Nx = 128
zmin = -20e-6
zmax = 20.0e-6
Nz = 128

# Wave vector of the wave
kx = 2.0 * np.pi * n_osc_x / (xmax - xmin)
kz = 2.0 * np.pi * n_osc_z / (zmax - zmin)
# Plasma frequency
wp = np.sqrt((n * e**2) / (m_e * epsilon_0))

k = {"Ex": kx, "Ez": kz}
cos = {"Ex": (0, 1, 1), "Ez": (1, 1, 0)}


def get_contribution(is_cos, k):
    du = (xmax - xmin) / Nx
    u = xmin + du * (0.5 + np.arange(Nx))
    if is_cos == 1:
        return np.cos(k * u)
    else:
        return np.sin(k * u)


def get_theoretical_field(field, t):
    amplitude = epsilon * (m_e * c**2 * k[field]) / e * np.sin(wp * t)
    cos_flag = cos[field]
    x_contribution = get_contribution(cos_flag[0], kx)
    z_contribution = get_contribution(cos_flag[2], kz)

    E = amplitude * x_contribution[:, np.newaxis] * z_contribution[np.newaxis, :]

    return E


# Read the file
ds = yt.load(fn)
t0 = ds.current_time.to_value()
data = ds.covering_grid(
    level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions
)
edge = np.array(
    [
        (ds.domain_left_edge[1]).item(),
        (ds.domain_right_edge[1]).item(),
        (ds.domain_left_edge[0]).item(),
        (ds.domain_right_edge[0]).item(),
    ]
)

# Check the validity of the fields
error_rel = 0
for field in ["Ex", "Ez"]:
    E_sim = data[("mesh", field)].to_ndarray()[:, :, 0]
    E_th = get_theoretical_field(field, t0)
    max_error = abs(E_sim - E_th).max() / abs(E_th).max()
    print("%s: Max error: %.2e" % (field, max_error))
    error_rel = max(error_rel, max_error)

# Plot the last field from the loop (Ez at iteration 40)
fig, (ax1, ax2) = plt.subplots(1, 2, dpi=100)
# First plot
vmin = E_sim.min()
vmax = E_sim.max()
cax1 = make_axes_locatable(ax1).append_axes("right", size="5%", pad="5%")
im1 = ax1.imshow(E_sim, origin="lower", extent=edge, vmin=vmin, vmax=vmax)
cb1 = fig.colorbar(im1, cax=cax1)
ax1.set_xlabel(r"$z$")
ax1.set_ylabel(r"$x$")
ax1.set_title(r"$E_z$ (sim)")
# Second plot
vmin = E_th.min()
vmax = E_th.max()
cax2 = make_axes_locatable(ax2).append_axes("right", size="5%", pad="5%")
im2 = ax2.imshow(E_th, origin="lower", extent=edge, vmin=vmin, vmax=vmax)
cb2 = fig.colorbar(im2, cax=cax2)
ax2.set_xlabel(r"$z$")
ax2.set_ylabel(r"$x$")
ax2.set_title(r"$E_z$ (theory)")
# Save figure
fig.tight_layout()
fig.savefig("Langmuir_multi_2d_analysis.png", dpi=200)

if particle_shape_4:
    # lower fidelity, due to smoothing
    tolerance_rel = 0.07
else:
    tolerance_rel = 0.0503

print("error_rel    : " + str(error_rel))
print("tolerance_rel: " + str(tolerance_rel))

assert error_rel < tolerance_rel

# Additional check on charge conservation for certain cases
# (e.g., Esirkepov or Vay deposition, current correction)
check_charge_conservation(data)
```

### RZ

### Script `analysis_rz.py`

```python3
#!/usr/bin/env python3

# Copyright 2019 David Grote, Maxence Thevenet
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL
#
# This is a script that analyses the simulation results from
# the script `inputs.multi.rz.rt`. This simulates a RZ periodic plasma wave.
# The electric field in the simulation is given (in theory) by:
# $$ E_r = -\partial_r \phi = \epsilon \,\frac{mc^2}{e}\frac{2\,r}{w_0^2} \exp\left(-\frac{r^2}{w_0^2}\right) \sin(k_0 z) \sin(\omega_p t)
# $$ E_z = -\partial_z \phi = - \epsilon \,\frac{mc^2}{e} k_0 \exp\left(-\frac{r^2}{w_0^2}\right) \cos(k_0 z) \sin(\omega_p t)
# Unrelated to the Langmuir waves, we also test the plotfile particle filter function in this
# analysis script.
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yt

yt.funcs.mylog.setLevel(50)

import numpy as np
import post_processing_utils
from analysis_utils import check_charge_conservation
from scipy.constants import c, e, epsilon_0, m_e

# this will be the name of the plot file
fn = sys.argv[1]

# test name
test_name = os.path.split(os.getcwd())[1]

# Parameters (these parameters must match the parameters in `inputs.multi.rz.rt`)
epsilon = 0.01
n = 2.0e24
w0 = 5.0e-6
n_osc_z = 2
rmin = 0e-6
rmax = 20.0e-6
Nr = 64
zmin = -20e-6
zmax = 20.0e-6
Nz = 128

# Wave vector of the wave
k0 = 2.0 * np.pi * n_osc_z / (zmax - zmin)
# Plasma frequency
wp = np.sqrt((n * e**2) / (m_e * epsilon_0))
kp = wp / c


def Er(z, r, epsilon, k0, w0, wp, t):
    """
    Return the radial electric field as an array
    of the same length as z and r, in the half-plane theta=0
    """
    Er_array = (
        epsilon
        * m_e
        * c**2
        / e
        * 2
        * r
        / w0**2
        * np.exp(-(r**2) / w0**2)
        * np.sin(k0 * z)
        * np.sin(wp * t)
    )
    return Er_array


def Ez(z, r, epsilon, k0, w0, wp, t):
    """
    Return the longitudinal electric field as an array
    of the same length as z and r, in the half-plane theta=0
    """
    Ez_array = (
        -epsilon
        * m_e
        * c**2
        / e
        * k0
        * np.exp(-(r**2) / w0**2)
        * np.cos(k0 * z)
        * np.sin(wp * t)
    )
    return Ez_array


# Read the file
ds = yt.load(fn)
t0 = ds.current_time.to_value()
data = ds.covering_grid(
    level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions
)

# Get cell centered coordinates
dr = (rmax - rmin) / Nr
dz = (zmax - zmin) / Nz
coords = np.indices([Nr, Nz], "d")
rr = rmin + (coords[0] + 0.5) * dr
zz = zmin + (coords[1] + 0.5) * dz

# Check the validity of the fields
overall_max_error = 0
Er_sim = data[("boxlib", "Er")].to_ndarray()[:, :, 0]
Er_th = Er(zz, rr, epsilon, k0, w0, wp, t0)
max_error = abs(Er_sim - Er_th).max() / abs(Er_th).max()
print("Er: Max error: %.2e" % (max_error))
overall_max_error = max(overall_max_error, max_error)

Ez_sim = data[("boxlib", "Ez")].to_ndarray()[:, :, 0]
Ez_th = Ez(zz, rr, epsilon, k0, w0, wp, t0)
max_error = abs(Ez_sim - Ez_th).max() / abs(Ez_th).max()
print("Ez: Max error: %.2e" % (max_error))
overall_max_error = max(overall_max_error, max_error)

# Plot the last field from the loop (Ez at iteration 40)
plt.subplot2grid((1, 2), (0, 0))
plt.imshow(Ez_sim)
plt.colorbar()
plt.title("Ez, last iteration\n(simulation)")
plt.subplot2grid((1, 2), (0, 1))
plt.imshow(Ez_th)
plt.colorbar()
plt.title("Ez, last iteration\n(theory)")
plt.tight_layout()
plt.savefig(test_name + "_analysis.png")

error_rel = overall_max_error

tolerance_rel = 0.12

print("error_rel    : " + str(error_rel))
print("tolerance_rel: " + str(tolerance_rel))

assert error_rel < tolerance_rel

# Additional check on charge conservation for certain cases
# (e.g., Esirkepov or Vay deposition, current correction)
check_charge_conservation(data)

## In the final past of the test, we verify that the diagnostic particle filter function works as
## expected in RZ geometry. For this, we only use the last simulation timestep.

dim = "rz"
species_name = "electrons"

# if test_name equals test_rz_langmuir_multi_psatd or test_rz_langmuir_multi_psatd_current_correction,
# we skip the check of the particle momentum along 'x', since it is not included in the output
skip_component = None
if test_name in [
    "test_rz_langmuir_multi_psatd",
    "test_rz_langmuir_multi_psatd_current_correction",
]:
    skip_component = "particle_momentum_x"

parser_filter_fn = "diags/diag_parser_filter000080"
parser_filter_expression = "(py-pz < 0) * (r<10e-6) * (z > 0)"
post_processing_utils.check_particle_filter(
    fn,
    parser_filter_fn,
    parser_filter_expression,
    dim,
    species_name,
    skip_component,
)

uniform_filter_fn = "diags/diag_uniform_filter000080"
uniform_filter_expression = "ids%3 == 0"
post_processing_utils.check_particle_filter(
    fn,
    uniform_filter_fn,
    uniform_filter_expression,
    dim,
    species_name,
    skip_component,
)

random_filter_fn = "diags/diag_random_filter000080"
random_fraction = 0.66
post_processing_utils.check_random_filter(
    fn,
    random_filter_fn,
    random_fraction,
    dim,
    species_name,
    skip_component,
)
```

### 1D

### Script `analysis_1d.py`

```python3
#!/usr/bin/env python3

# Copyright 2019-2022 Jean-Luc Vay, Maxence Thevenet, Remi Lehe, Prabhat Kumar, Axel Huebl
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL
#
# This is a script that analyses the simulation results from
# the script `inputs.multi.rt`. This simulates a 1D periodic plasma wave.
# The electric field in the simulation is given (in theory) by:
# $$ E_z = \epsilon \,\frac{m_e c^2 k_z}{q_e}\sin(k_z z)\sin( \omega_p t)$$
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yt

yt.funcs.mylog.setLevel(50)

import numpy as np
from analysis_utils import check_charge_conservation
from scipy.constants import c, e, epsilon_0, m_e

# test name
test_name = os.path.split(os.getcwd())[1]

# this will be the name of the plot file
fn = sys.argv[1]

# Parameters (these parameters must match the parameters in `inputs.multi.rt`)
epsilon = 0.01
n = 4.0e24
n_osc_z = 2
zmin = -20e-6
zmax = 20.0e-6
Nz = 128

# Wave vector of the wave
kz = 2.0 * np.pi * n_osc_z / (zmax - zmin)
# Plasma frequency
wp = np.sqrt((n * e**2) / (m_e * epsilon_0))

k = {"Ez": kz}
cos = {"Ez": (1, 1, 0)}


def get_contribution(is_cos, k):
    du = (zmax - zmin) / Nz
    u = zmin + du * (0.5 + np.arange(Nz))
    if is_cos == 1:
        return np.cos(k * u)
    else:
        return np.sin(k * u)


def get_theoretical_field(field, t):
    amplitude = epsilon * (m_e * c**2 * k[field]) / e * np.sin(wp * t)
    cos_flag = cos[field]
    z_contribution = get_contribution(cos_flag[2], kz)

    E = amplitude * z_contribution

    return E


# Read the file
ds = yt.load(fn)
t0 = ds.current_time.to_value()
data = ds.covering_grid(
    level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions
)
# Check the validity of the fields
error_rel = 0
for field in ["Ez"]:
    E_sim = data[("mesh", field)].to_ndarray()[:, 0, 0]
    E_th = get_theoretical_field(field, t0)
    max_error = abs(E_sim - E_th).max() / abs(E_th).max()
    print("%s: Max error: %.2e" % (field, max_error))
    error_rel = max(error_rel, max_error)

# Plot the last field from the loop (Ez at iteration 80)
plt.subplot2grid((1, 2), (0, 0))
plt.plot(E_sim)
# plt.colorbar()
plt.title("Ez, last iteration\n(simulation)")
plt.subplot2grid((1, 2), (0, 1))
plt.plot(E_th)
# plt.colorbar()
plt.title("Ez, last iteration\n(theory)")
plt.tight_layout()
plt.savefig("langmuir_multi_1d_analysis.png")

tolerance_rel = 0.05

print("error_rel    : " + str(error_rel))
print("tolerance_rel: " + str(tolerance_rel))

assert error_rel < tolerance_rel

# Additional check on charge conservation for certain cases
# (e.g., Esirkepov or Vay deposition, current correction)
check_charge_conservation(data)
```

## Visualize

#### NOTE
This section is TODO.
