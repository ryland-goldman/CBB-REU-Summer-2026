<a id="examples-plasma-mirror"></a>

# Plasma-Mirror

This example shows how to model a plasma mirror, using a planar target of solid density [[7](../../examples.md#id19), [8](../../examples.md#id20)].

Although laser-solid interaction modeling requires full 3D modeling for adequate description of the dynamics at play, this example models a 2D example.
2D modeling provide a qualitative overview of the dynamics, but mostly saves computational costs since the plasma frequency (and Debye length) of the surface plasma determines the resolution need in laser-solid interaction modeling.

#### NOTE
TODO: The Python (PICMI) input file needs to be created.

## Run

This example can be run **either** as:

* **Python** script: (*TODO*) or
* WarpX **executable** using an input file: `warpx.2d inputs_2d`

For [MPI-parallel](https://www.mpi-forum.org) runs, prefix these lines with `mpiexec -n 4 ...` or `srun -n 4 ...`, depending on the system.

### Python: Script

#### NOTE
TODO: This input file should be created following the `inputs_2d` file.

### Executable: Input File

```none
#################################
####### GENERAL PARAMETERS ######
#################################
max_step = 20 #1000
amr.n_cell = 256 128 #1024 512
amr.max_grid_size = 128
amr.blocking_factor = 32
amr.max_level = 0
geometry.dims = 2
geometry.prob_lo = -100.e-6   0.     # physical domain
geometry.prob_hi =  100.e-6   100.e-6
warpx.verbose = 1
warpx.serialize_initial_conditions = 1

#################################
####### Boundary condition ######
#################################
boundary.field_lo = pml pml
boundary.field_hi = pml pml

#################################
############ NUMERICS ###########
#################################
my_constants.zc    = 20.e-6
my_constants.zp    = 20.05545177444479562e-6
my_constants.lgrad = .08e-6
my_constants.nc    = 1.74e27
my_constants.zp2   = 24.e-6
my_constants.zc2   = 24.05545177444479562e-6
warpx.cfl = 1.0
warpx.use_filter = 1
algo.load_balance_intervals = 66

# Order of particle shape factors
algo.particle_shape = 3

#################################
############ PLASMA #############
#################################
particles.species_names = electrons ions

electrons.charge = -q_e
electrons.mass = m_e
electrons.injection_style = NUniformPerCell
electrons.num_particles_per_cell_each_dim = 2 2
electrons.momentum_distribution_type = "gaussian"
electrons.ux_th = .01
electrons.uz_th = .01
electrons.zmin = "zc-lgrad*log(400)"
electrons.zmax = 25.47931e-6
electrons.profile = parse_density_function
electrons.density_function(x,y,z) = "if(z<zp, nc*exp((z-zc)/lgrad), if(z<=zp2, 2.*nc, nc*exp(-(z-zc2)/lgrad)))"

ions.charge = q_e
ions.mass = m_p
ions.injection_style = NUniformPerCell
ions.num_particles_per_cell_each_dim = 2 2
ions.momentum_distribution_type = "at_rest"
ions.zmin = 19.520e-6
ions.zmax = 25.47931e-6
ions.profile = parse_density_function
ions.density_function(x,y,z) = "if(z<zp, nc*exp((z-zc)/lgrad), if(z<=zp2, 2.*nc, nc*exp(-(z-zc2)/lgrad)))"

#################################
############# LASER #############
#################################
lasers.names        = laser1
laser1.position     = 0. 0. 5.e-6 # This point is on the laser plane
laser1.direction    = 0. 0. 1.     # The plane normal direction
laser1.polarization = 1. 0. 0.     # The main polarization vector
laser1.e_max        = 4.e12        # Maximum amplitude of the laser field (in V/m)
laser1.wavelength = 0.8e-6         # The wavelength of the laser (in meters)
laser1.profile      = Gaussian
laser1.profile_waist = 5.e-6      # The waist of the laser (in meters)
laser1.profile_duration = 15.e-15  # The duration of the laser (in seconds)
laser1.profile_t_peak = 25.e-15    # The time at which the laser reaches its peak (in seconds)
laser1.profile_focal_distance = 15.e-6  # Focal distance from the antenna (in meters)

# Diagnostics
diagnostics.diags_names = diag1
diag1.intervals = 10
diag1.diag_type = Full
```

## Analyze

#### NOTE
This section is TODO.

## Visualize

#### NOTE
This section is TODO.
