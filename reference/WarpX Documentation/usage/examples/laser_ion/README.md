<a id="examples-laser-ion"></a>

# Laser-Ion Acceleration with a Planar Target

This example shows how to model laser-ion acceleration with planar targets of solid density [[4](../../examples.md#id17), [5](../../examples.md#id18), [6](../../examples.md#id16)].
The acceleration mechanism in this scenario depends on target parameters.

Although laser-ion acceleration requires full 3D modeling for adequate description of the acceleration dynamics, especially the acceleration field lengths and decay times, this example models a 2D example.
2D modeling can often hint at a qualitative overview of the dynamics, but mostly saves computational costs since the plasma frequency (and Debye length) of the plasma determines the resolution need in laser-solid interaction modeling.

#### NOTE
The resolution of this 2D case is extremely low by default.
This includes spatial and temporal resolution, but also the number of macro-particles per cell representing the target density for proper phase space sampling.
You will need a computing cluster for adequate resolution of the target density, see comments in the input file.

## Run

This example can be run **either** as:

* **Python** script: `mpiexec -n 2 python3 inputs_test_2d_laser_ion_acc_picmi.py` or
* WarpX **executable** using an input file: `mpiexec -n 2 warpx.2d inputs_test_2d_laser_ion_acc`

### Python: Script

```python3
#!/usr/bin/env python3

from pywarpx import picmi

# Physical constants
c = picmi.constants.c
q_e = picmi.constants.q_e

# We only run 100 steps for tests
# Disable `max_step` below to run until the physical `stop_time`.
max_step = 100
# time-scale with highly kinetic dynamics
stop_time = 0.2e-12

# proper resolution for 30 n_c (dx<=3.33nm) incl. acc. length
# (>=6x V100)
# --> choose larger `max_grid_size` and `blocking_factor` for 1 to 8 grids per GPU accordingly
# nx = 7488
# nz = 14720

# Number of cells
nx = 384
nz = 512

# Domain decomposition (deactivate `warpx_numprocs` in `picmi.Simulation` for this to take effect)
max_grid_size = 64
blocking_factor = 32

# Physical domain
xmin = -7.5e-06
xmax = 7.5e-06
zmin = -5.0e-06
zmax = 25.0e-06

# Create grid
grid = picmi.Cartesian2DGrid(
    number_of_cells=[nx, nz],
    lower_bound=[xmin, zmin],
    upper_bound=[xmax, zmax],
    lower_boundary_conditions=["open", "open"],
    upper_boundary_conditions=["open", "open"],
    lower_boundary_conditions_particles=["absorbing", "absorbing"],
    upper_boundary_conditions_particles=["absorbing", "absorbing"],
    warpx_max_grid_size=max_grid_size,
    warpx_blocking_factor=blocking_factor,
)

# Particles: plasma parameters
# critical plasma density
nc = 1.742e27  # [m^-3]  1.11485e21 * 1.e6 / 0.8**2
# number density: "fully ionized" electron density as reference
#   [material 1] cryogenic H2
n0 = 30.0  # [n_c]
#   [material 2] liquid crystal
# n0 = 192
#   [material 3] PMMA
# n0 = 230
#   [material 4] Copper (ion density: 8.49e28/m^3; times ionization level)
# n0 = 1400
plasma_density = n0 * nc
preplasma_L = 0.05e-6  # [m] scale length (>0)
preplasma_Lcut = 2.0e-6  # [m] hard cutoff from surface
plasma_r0 = 2.5e-6  # [m] radius or half-thickness
plasma_eps_z = 0.05e-6  # [m] small offset in z to make zmin, zmax interval larger than 2*(r0 + Lcut)
plasma_creation_limit_z = (
    plasma_r0 + preplasma_Lcut + plasma_eps_z
)  # [m] upper limit in z for particle creation

plasma_xmin = None
plasma_ymin = None
plasma_zmin = -plasma_creation_limit_z
plasma_xmax = None
plasma_ymax = None
plasma_zmax = plasma_creation_limit_z

density_expression_str = f"{plasma_density}*((abs(z)<={plasma_r0}) + (abs(z)<{plasma_r0}+{preplasma_Lcut}) * (abs(z)>{plasma_r0}) * exp(-(abs(z)-{plasma_r0})/{preplasma_L}))"

slab_with_ramp_dist_hydrogen = picmi.AnalyticDistribution(
    density_expression=density_expression_str,
    lower_bound=[plasma_xmin, plasma_ymin, plasma_zmin],
    upper_bound=[plasma_xmax, plasma_ymax, plasma_zmax],
)

# thermal velocity spread for electrons in gamma*beta
ux_th = 0.01
uz_th = 0.01

slab_with_ramp_dist_electrons = picmi.AnalyticDistribution(
    density_expression=density_expression_str,
    lower_bound=[plasma_xmin, plasma_ymin, plasma_zmin],
    upper_bound=[plasma_xmax, plasma_ymax, plasma_zmax],
    # if `momentum_expressions` and `momentum_spread_expressions` are unset,
    # a Gaussian momentum distribution is assumed given that `rms_velocity` has any non-zero elements
    rms_velocity=[c * ux_th, 0.0, c * uz_th],  # thermal velocity spread in m/s
)

electrons = picmi.Species(
    particle_type="electron",
    name="electrons",
    initial_distribution=slab_with_ramp_dist_electrons,
)

hydrogen = picmi.Species(
    particle_type="proton",
    name="hydrogen",
    initial_distribution=slab_with_ramp_dist_hydrogen,
    warpx_add_real_attributes={"orig_x": "x", "orig_z": "z"},
)

# Laser
# e_max = a0 * 3.211e12 / lambda_0[mu]
#   a0 = 16, lambda_0 = 0.8mu -> e_max = 64.22 TV/m
e_max = 64.22e12
position_z = -4.0e-06
profile_t_peak = 50.0e-15
profile_focal_distance = 4.0e-06
laser = picmi.GaussianLaser(
    wavelength=0.8e-06,
    waist=4.0e-06,
    duration=30.0e-15,
    focal_position=[0, 0, profile_focal_distance + position_z],
    centroid_position=[0, 0, position_z - c * profile_t_peak],
    propagation_direction=[0, 0, 1],
    polarization_direction=[1, 0, 0],
    E0=e_max,
    fill_in=False,
)
laser_antenna = picmi.LaserAntenna(
    position=[0.0, 0.0, position_z], normal_vector=[0, 0, 1]
)

# Electromagnetic solver
solver = picmi.ElectromagneticSolver(
    grid=grid,
    method="Yee",
    cfl=0.999,
    divE_cleaning=0,
    # warpx_pml_ncell=10
)

# Diagnostics
particle_diag = picmi.ParticleDiagnostic(
    name="diagInst",
    period=100,
    warpx_format="openpmd",
    warpx_openpmd_backend="h5",
    # demonstration of a spatial and momentum filter
    warpx_plot_filter_function="(uz>=0) * (x<1.0e-6) * (x>-1.0e-6)",
)
# reduce resolution of output fields
coarsening_ratio = [4, 4]
ncell_field = []
for ncell_comp, cr in zip([nx, nz], coarsening_ratio):
    ncell_field.append(int(ncell_comp / cr))
field_diag = picmi.FieldDiagnostic(
    name="diagInst",
    grid=grid,
    period=100,
    number_of_cells=ncell_field,
    data_list=["B", "E", "J", "rho", "rho_electrons", "rho_hydrogen"],
    warpx_format="openpmd",
    warpx_openpmd_backend="h5",
)

field_time_avg_diag = picmi.TimeAveragedFieldDiagnostic(
    name="diagTimeAvg",
    grid=grid,
    period=100,
    number_of_cells=ncell_field,
    data_list=["B", "E", "J", "rho", "rho_electrons", "rho_hydrogen"],
    warpx_format="openpmd",
    warpx_openpmd_backend="h5",
    warpx_time_average_mode="dynamic_start",
    warpx_average_period_time=2.67e-15,
)

particle_fw_diag = picmi.ParticleDiagnostic(
    name="openPMDfw",
    period=100,
    warpx_format="openpmd",
    warpx_openpmd_backend="h5",
    warpx_plot_filter_function="(uz>=0) * (x<1.0e-6) * (x>-1.0e-6)",
)

particle_bw_diag = picmi.ParticleDiagnostic(
    name="openPMDbw",
    period=100,
    warpx_format="openpmd",
    warpx_openpmd_backend="h5",
    warpx_plot_filter_function="(uz<0)",
)

# histograms with 2.0 degree acceptance angle in fw direction
# 2 deg * pi / 180 : 0.03490658503 rad
# half-angle +/-   : 0.017453292515 rad
histuH_rdiag = picmi.ReducedDiagnostic(
    diag_type="ParticleHistogram",
    name="histuH",
    species=hydrogen,
    bin_number=1000,
    bin_min=0.0,
    bin_max=0.474,  # 100 MeV protons
    histogram_function="u2=ux*ux+uy*uy+uz*uz; if(u2>0, sqrt(u2), 0.0)",
    filter_function="u2=ux*ux+uy*uy+uz*uz; if(u2>0, abs(acos(uz / sqrt(u2))) <= 0.017453, 0)",
)

histue_rdiag = picmi.ReducedDiagnostic(
    diag_type="ParticleHistogram",
    name="histue",
    species=electrons,
    bin_number=1000,
    bin_min=0.0,
    bin_max=197.0,  # 100 MeV electrons
    histogram_function="u2=ux*ux+uy*uy+uz*uz; if(u2>0, sqrt(u2), 0.0)",
    filter_function="u2=ux*ux+uy*uy+uz*uz; if(u2>0, abs(acos(uz / sqrt(u2))) <= 0.017453, 0)",
)

# just a test entry to make sure that the histogram filter is purely optional:
# this one just records uz of all hydrogen ions, independent of their pointing
histuzAll_rdiag = picmi.ReducedDiagnostic(
    diag_type="ParticleHistogram",
    name="histuzAll",
    species=hydrogen,
    bin_number=1000,
    bin_min=-0.474,
    bin_max=0.474,
    histogram_function="uz",
)

field_probe_z_rdiag = picmi.ReducedDiagnostic(
    diag_type="FieldProbe",
    name="FieldProbe_Z",
    integrate=0,
    probe_geometry="Line",
    x_probe=0.0,
    z_probe=-5.0e-6,
    x1_probe=0.0,
    z1_probe=25.0e-6,
    resolution=3712,
)

field_probe_scat_point_rdiag = picmi.ReducedDiagnostic(
    diag_type="FieldProbe",
    name="FieldProbe_ScatPoint",
    integrate=0,
    probe_geometry="Point",
    x_probe=0.0,
    z_probe=15.0e-6,
)

field_probe_scat_line_rdiag = picmi.ReducedDiagnostic(
    diag_type="FieldProbe",
    name="FieldProbe_ScatLine",
    integrate=1,
    probe_geometry="Line",
    x_probe=-2.5e-6,
    z_probe=15.0e-6,
    x1_probe=2.5e-6,
    z1_probe=15e-6,
    resolution=201,
)

load_balance_costs_rdiag = picmi.ReducedDiagnostic(
    diag_type="LoadBalanceCosts",
    name="LBC",
)

# Set up simulation
sim = picmi.Simulation(
    solver=solver,
    max_time=stop_time,  # need to remove `max_step` to run this far
    verbose=1,
    particle_shape="cubic",
    warpx_numprocs=[1, 2],  # deactivate `numprocs` for dynamic load balancing
    warpx_use_filter=1,
    warpx_reduced_diags_intervals=100,
    warpx_load_balance_intervals=100,
    warpx_load_balance_costs_update="heuristic",
)

# Add plasma electrons
sim.add_species(
    electrons,
    layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=[2, 2]),
    # for more realistic simulations, try to avoid that macro-particles represent more than 1 n_c
    # layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=[4,8])
)

# Add hydrogen ions
sim.add_species(
    hydrogen,
    layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=[2, 2]),
    # for more realistic simulations, try to avoid that macro-particles represent more than 1 n_c
    # layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=[4,8])
)

# Add laser
sim.add_laser(laser, injection_method=laser_antenna)

# Add full diagnostics
sim.add_diagnostic(particle_diag)
sim.add_diagnostic(field_diag)
sim.add_diagnostic(field_time_avg_diag)
sim.add_diagnostic(particle_fw_diag)
sim.add_diagnostic(particle_bw_diag)
# Add reduced diagnostics
sim.add_diagnostic(histuH_rdiag)
sim.add_diagnostic(histue_rdiag)
sim.add_diagnostic(histuzAll_rdiag)
sim.add_diagnostic(field_probe_z_rdiag)
sim.add_diagnostic(field_probe_scat_point_rdiag)
sim.add_diagnostic(field_probe_scat_line_rdiag)
sim.add_diagnostic(load_balance_costs_rdiag)
# TODO: make ParticleHistogram2D available

# Write input file that can be used to run with the compiled version
sim.write_input_file(file_name="inputs_2d_picmi")

# Initialize inputs and WarpX instance
sim.initialize_inputs()
sim.initialize_warpx()

# Advance simulation until last time step
sim.step(max_step)
```

### Executable: Input File

```none
#################################
# Domain, Resolution & Numerics
#

# We only run 100 steps for tests
# Disable `max_step` below to run until the physical `stop_time`.
max_step = 100
# time-scale with highly kinetic dynamics
stop_time = 0.2e-12            # [s]
# time-scale for converged ion energy
#   notes: - effective acc. time depends on laser pulse
#          - ions will start to leave the box
#stop_time = 1.0e-12           # [s]

# quick tests at ultra-low res. (for CI, and local computer)
amr.n_cell = 384 512

# proper resolution for 10 n_c excl. acc. length
# (>=1x V100)
#amr.n_cell = 2688 3712

# proper resolution for 30 n_c (dx<=3.33nm) incl. acc. length
# (>=6x V100)
#amr.n_cell = 7488 14720

# simulation box, no MR
#   note: increase z (space & cells) for converging ion energy
amr.max_level = 0
geometry.dims = 2
geometry.prob_lo = -7.5e-6 -5.e-6
geometry.prob_hi =  7.5e-6 25.e-6

# Boundary condition
boundary.field_lo = pml pml
boundary.field_hi = pml pml

# Order of particle shape factors
algo.particle_shape = 3

# improved plasma stability for 2D with very low initial target temperature
# when using Esirkepov current deposition with energy-conserving field gather
interpolation.galerkin_scheme = 0

# numerical tuning
warpx.cfl = 0.999
warpx.use_filter = 1          # bilinear current/charge filter


#################################
# Performance Tuning
#
# simple tuning:
#   the numprocs product must be equal to the number of MPI ranks and splits
#   the domain on the coarsest level equally into grids;
#   slicing in the 2nd dimension is preferred for ideal performance
warpx.numprocs = 1 2   # 2 MPI ranks
#warpx.numprocs = 1 4  # 4 MPI ranks

# detail tuning instead of warpx.numprocs:
#   It is important to have enough cells in a block & grid, otherwise
#   performance will suffer.
#   Use larger values for GPUs, try to fill a GPU well with memory and place
#   few large grids on each device (you can go as low as 1 large grid / device
#   if you do not need load balancing).
#   Slicing in the 2nd dimension is preferred for ideal performance
#amr.blocking_factor = 64
#amr.max_grid_size_x = 2688
#amr.max_grid_size_y = 128  # this is confusingly named and means z in 2D

# load balancing
#   The grid & block parameters above are needed for load balancing:
#   an average of ~10 grids per MPI rank (and device) are a good granularity
#   to allow efficient load-balancing as the simulation evolves
algo.load_balance_intervals = 100
algo.load_balance_costs_update = Heuristic

# particle bin-sorting on GPU (ideal defaults not investigated in 2D)
#   Try larger values than the defaults below and report back! :)
#warpx.sort_intervals = 4    # default on CPU: -1 (off); on GPU: 4
#warpx.sort_bin_size = 1 1 1


#################################
# Target Profile
#

#   definitions for target extent and pre-plasma
my_constants.L    = 0.05e-6            # [m] scale length (>0)
my_constants.Lcut = 2.0e-6             # [m] hard cutoff from surface
my_constants.r0 = 2.5e-6               # [m] radius or half-thickness
my_constants.eps_z = 0.05e-6           # [m] small offset in z to make zmin, zmax interval larger than 2*(r0 + Lcut)
my_constants.zmax = r0 + Lcut + eps_z  # [m] upper limit in z for particle creation

particles.species_names = electrons hydrogen

# particle species
hydrogen.species_type = hydrogen
hydrogen.injection_style = NUniformPerCell
hydrogen.num_particles_per_cell_each_dim = 2 2
# for more realistic simulations, try to avoid that macro-particles represent more than 1 n_c
#hydrogen.num_particles_per_cell_each_dim = 4 8
hydrogen.momentum_distribution_type = at_rest
# minimum and maximum z position between which particles are initialized
# --> should be set for dense targets limit memory consumption during initialization
hydrogen.zmin = -zmax
hydrogen.zmax = zmax
hydrogen.profile = parse_density_function
hydrogen.addRealAttributes = orig_x orig_z
hydrogen.attribute.orig_x(x,y,z,ux,uy,uz,t) = "x"
hydrogen.attribute.orig_z(x,y,z,ux,uy,uz,t) = "z"

electrons.species_type = electron
electrons.injection_style = NUniformPerCell
electrons.num_particles_per_cell_each_dim = 2 2
# for more realistic simulations, try to avoid that macro-particles represent more than 1 n_c
#electrons.num_particles_per_cell_each_dim = 4 8
electrons.momentum_distribution_type = "gaussian"
electrons.ux_th = .01
electrons.uz_th = .01
# minimum and maximum z position between which particles are initialized
# --> should be set for dense targets limit memory consumption during initialization
electrons.zmin = -zmax
electrons.zmax = zmax

# ionization physics (field ionization/ADK)
#   [i1] none (fully pre-ionized):
electrons.profile = parse_density_function
#   [i2] field ionization (ADK):
#hydrogen.do_field_ionization = 1
#hydrogen.physical_element = H
#hydrogen.ionization_initial_level = 0
#hydrogen.ionization_product_species = electrons
#electrons.profile = constant
#electrons.density = 0.0

# collisional physics (binary MC model after Nanbu/Perez)
#collisions.collision_names = c_eH c_ee c_HH
#c_eH.species = electrons hydrogen
#c_ee.species = electrons electrons
#c_HH.species = hydrogen hydrogen
#c_eH.CoulombLog = 15.9
#c_ee.CoulombLog = 15.9
#c_HH.CoulombLog = 15.9

# number density: "fully ionized" electron density as reference
#   [material 1] cryogenic H2
my_constants.nc    = 1.742e27  # [m^-3]  1.11485e21 * 1.e6 / 0.8**2
my_constants.n0    = 30.0      # [n_c]
#   [material 2] liquid crystal
#my_constants.n0    = 192
#   [material 3] PMMA
#my_constants.n0    = 230
#   [material 4] Copper (ion density: 8.49e28/m^3; times ionization level)
#my_constants.n0    = 1400

# density profiles (target extent, pre-plasma and cutoffs defined above particle species list)

# [target 1] flat foil (thickness = 2*r0)
electrons.density_function(x,y,z) = "nc*n0*(
    if(abs(z)<=r0, 1.0, if(abs(z)<r0+Lcut, exp((-abs(z)+r0)/L), 0.0)) )"
hydrogen.density_function(x,y,z) = "nc*n0*(
    if(abs(z)<=r0, 1.0, if(abs(z)<r0+Lcut, exp((-abs(z)+r0)/L), 0.0)) )"

# [target 2] cylinder
#electrons.density_function(x,y,z) = "nc*n0*(
#    ((x*x+z*z)<=(r0*r0)) +
#    (sqrt(x*x+z*z)>r0)*(sqrt(x*x+z*z)<r0+Lcut)*exp( (-sqrt(x*x+z*z)+r0)/L ) )"
#hydrogen.density_function(x,y,z) = "nc*n0*(
#    ((x*x+z*z)<=(r0*r0)) +
#    (sqrt(x*x+z*z)>r0)*(sqrt(x*x+z*z)<r0+Lcut)*exp( (-sqrt(x*x+z*z)+r0)/L ) )"

# [target 3] sphere
#electrons.density_function(x,y,z) = "nc*n0*(
#    ((x*x+y*y+z*z)<=(r0*r0)) +
#    (sqrt(x*x+y*y+z*z)>r0)*(sqrt(x*x+y*y+z*z)<r0+Lcut)*exp( (-sqrt(x*x+y*y+z*z)+r0)/L ) )"
#hydrogen.density_function(x,y,z) = "nc*n0*(
#    ((x*x+y*y+z*z)<=(r0*r0)) +
#    (sqrt(x*x+y*y+z*z)>r0)*(sqrt(x*x+y*y+z*z)<r0+Lcut)*exp( (-sqrt(x*x+y*y+z*z)+r0)/L ) )"


#################################
# Laser Pulse Profile
#
lasers.names        = laser1
laser1.position     = 0. 0. -4.0e-6     # point the laser plane (antenna)
laser1.direction    = 0. 0. 1.          # the plane's (antenna's) normal direction
laser1.polarization = 1. 0. 0.          # the main polarization vector
laser1.a0           = 16.0              # maximum amplitude of the laser field [V/m]
laser1.wavelength   = 0.8e-6            # central wavelength of the laser pulse [m]
laser1.profile      = Gaussian
laser1.profile_waist = 4.e-6            # beam waist (E(w_0)=E_0/e) [m]
laser1.profile_duration = 30.e-15       # pulse length (E(tau)=E_0/e; tau=tau_E=FWHM_I/1.17741) [s]
laser1.profile_t_peak = 50.e-15         # time until peak intensity reached at the laser plane [s]
laser1.profile_focal_distance = 4.0e-6  # focal distance from the antenna [m]

# e_max = a0 * 3.211e12 / lambda_0[mu]
#   a0 = 16, lambda_0 = 0.8mu -> e_max = 64.22 TV/m


#################################
# Diagnostics
#
diagnostics.diags_names = diagInst diagTimeAvg openPMDfw openPMDbw

# instantaneous field and particle diagnostic
diagInst.intervals = 100,96:100  # second interval only for CI testing the time-averaged diags
diagInst.diag_type = Full
diagInst.fields_to_plot = Ex Ey Ez Bx By Bz jx jy jz rho rho_electrons rho_hydrogen
# reduce resolution of output fields
diagInst.coarsening_ratio = 4 4
# demonstration of a spatial and momentum filter
diagInst.electrons.plot_filter_function(t,x,y,z,ux,uy,uz) = (uz>=0) * (x<1.0e-6) * (x>-1.0e-6)
diagInst.hydrogen.plot_filter_function(t,x,y,z,ux,uy,uz) = (uz>=0) * (x<1.0e-6) * (x>-1.0e-6)
diagInst.format = openpmd
diagInst.openpmd_backend = h5

# time-averaged particle and field diagnostic
diagTimeAvg.intervals = 100
diagTimeAvg.diag_type = TimeAveraged
diagTimeAvg.time_average_mode = dynamic_start
#diagTimeAvg.average_period_time = 2.67e-15  # period of 800 nm light waves
diagTimeAvg.average_period_steps = 5  # use only either `time` or `steps`
diagTimeAvg.write_species = 0
diagTimeAvg.fields_to_plot = Ex Ey Ez Bx By Bz jx jy jz rho rho_electrons rho_hydrogen
# reduce resolution of output fields
diagTimeAvg.coarsening_ratio = 4 4
diagTimeAvg.format = openpmd
diagTimeAvg.openpmd_backend = h5

openPMDfw.intervals = 100
openPMDfw.diag_type = Full
openPMDfw.fields_to_plot = Ex Ey Ez Bx By Bz jx jy jz rho rho_electrons rho_hydrogen
# reduce resolution of output fields
openPMDfw.coarsening_ratio = 4 4
openPMDfw.format = openpmd
openPMDfw.openpmd_backend = h5
# demonstration of a spatial and momentum filter
openPMDfw.electrons.plot_filter_function(t,x,y,z,ux,uy,uz) = (uz>=0) * (x<1.0e-6) * (x>-1.0e-6)
openPMDfw.hydrogen.plot_filter_function(t,x,y,z,ux,uy,uz) = (uz>=0) * (x<1.0e-6) * (x>-1.0e-6)

openPMDbw.intervals = 100
openPMDbw.diag_type = Full
openPMDbw.fields_to_plot = rho_hydrogen
# reduce resolution of output fields
openPMDbw.coarsening_ratio = 4 4
openPMDbw.format = openpmd
openPMDbw.openpmd_backend = h5
# demonstration of a momentum filter
openPMDbw.electrons.plot_filter_function(t,x,y,z,ux,uy,uz) = (uz<0)
openPMDbw.hydrogen.plot_filter_function(t,x,y,z,ux,uy,uz) = (uz<0)


#################################
# Reduced Diagnostics
#

# histograms with 2.0 degree acceptance angle in fw direction
# 2 deg * pi / 180 : 0.03490658503 rad
# half-angle +/-   : 0.017453292515 rad
warpx.reduced_diags_names                   = histuH histue histuzAll FieldProbe_Z FieldProbe_ScatPoint FieldProbe_ScatLine LBC PhaseSpaceIons PhaseSpaceElectrons

histuH.type                                 = ParticleHistogram
histuH.intervals                            = 100
histuH.species                              = hydrogen
histuH.bin_number                           = 1000
histuH.bin_min                              =  0.0
histuH.bin_max                              =  0.474  # 100 MeV protons
histuH.histogram_function(t,x,y,z,ux,uy,uz) = "u2=ux*ux+uy*uy+uz*uz; if(u2>0, sqrt(u2), 0.0)"
histuH.filter_function(t,x,y,z,ux,uy,uz) = "u2=ux*ux+uy*uy+uz*uz; if(u2>0, abs(acos(uz / sqrt(u2))) <= 0.017453, 0)"

histue.type                                 = ParticleHistogram
histue.intervals                            = 100
histue.species                              = electrons
histue.bin_number                           = 1000
histue.bin_min                              = 0.0
histue.bin_max                              = 197  # 100 MeV electrons
histue.histogram_function(t,x,y,z,ux,uy,uz) = "u2=ux*ux+uy*uy+uz*uz; if(u2>0, sqrt(u2), 0.0)"
histue.filter_function(t,x,y,z,ux,uy,uz) = "u2=ux*ux+uy*uy+uz*uz; if(u2>0, abs(acos(uz / sqrt(u2))) <= 0.017453, 0)"

# just a test entry to make sure that the histogram filter is purely optional:
# this one just records uz of all hydrogen ions, independent of their pointing
histuzAll.type                                 = ParticleHistogram
histuzAll.intervals                            = 100
histuzAll.species                              = hydrogen
histuzAll.bin_number                           = 1000
histuzAll.bin_min                              = -0.474
histuzAll.bin_max                              =  0.474
histuzAll.histogram_function(t,x,y,z,ux,uy,uz) = "uz"

FieldProbe_Z.type = FieldProbe
FieldProbe_Z.intervals = 100
FieldProbe_Z.integrate = 0
FieldProbe_Z.probe_geometry = Line
FieldProbe_Z.x_probe = 0.0
FieldProbe_Z.z_probe = -5.0e-6
FieldProbe_Z.x1_probe = 0.0
FieldProbe_Z.z1_probe = 25.0e-6
FieldProbe_Z.resolution = 3712

FieldProbe_ScatPoint.type = FieldProbe
FieldProbe_ScatPoint.intervals = 1
FieldProbe_ScatPoint.integrate = 0
FieldProbe_ScatPoint.probe_geometry = Point
FieldProbe_ScatPoint.x_probe = 0.0
FieldProbe_ScatPoint.z_probe = 15e-6

FieldProbe_ScatLine.type = FieldProbe
FieldProbe_ScatLine.intervals = 100
FieldProbe_ScatLine.integrate = 1
FieldProbe_ScatLine.probe_geometry = Line
FieldProbe_ScatLine.x_probe = -2.5e-6
FieldProbe_ScatLine.z_probe = 15e-6
FieldProbe_ScatLine.x1_probe = 2.5e-6
FieldProbe_ScatLine.z1_probe = 15e-6
FieldProbe_ScatLine.resolution = 201

# check computational load per box
LBC.type = LoadBalanceCosts
LBC.intervals = 100

PhaseSpaceIons.type                                 = ParticleHistogram2D
PhaseSpaceIons.intervals                            = 100
PhaseSpaceIons.species                              = hydrogen
PhaseSpaceIons.bin_number_abs                       = 1000
PhaseSpaceIons.bin_number_ord                       = 1000
PhaseSpaceIons.bin_min_abs                          = -5.e-6
PhaseSpaceIons.bin_max_abs                          = 25.e-6
PhaseSpaceIons.bin_min_ord                          = -0.474
PhaseSpaceIons.bin_max_ord                          = 0.474
PhaseSpaceIons.histogram_function_abs(t,x,y,z,ux,uy,uz,w) = "z"
PhaseSpaceIons.histogram_function_ord(t,x,y,z,ux,uy,uz,w) = "uz"
PhaseSpaceIons.value_function(t,x,y,z,ux,uy,uz,w) = "w"
# PhaseSpaceIons.filter_function(t,x,y,z,ux,uy,uz,w) = "u2=ux*ux+uy*uy+uz*uz; if(u2>0, abs(acos(uz / sqrt(u2))) <= 0.017453, 0)"

PhaseSpaceElectrons.type                                 = ParticleHistogram2D
PhaseSpaceElectrons.intervals                            = 100
PhaseSpaceElectrons.species                              = electrons
PhaseSpaceElectrons.bin_number_abs                       = 1000
PhaseSpaceElectrons.bin_number_ord                       = 1000
PhaseSpaceElectrons.bin_min_abs                          = -5.e-6
PhaseSpaceElectrons.bin_max_abs                          = 25.e-6
PhaseSpaceElectrons.bin_min_ord                          = -197
PhaseSpaceElectrons.bin_max_ord                          = 197
PhaseSpaceElectrons.histogram_function_abs(t,x,y,z,ux,uy,uz,w) = "z"
PhaseSpaceElectrons.histogram_function_ord(t,x,y,z,ux,uy,uz,w) = "uz"
PhaseSpaceElectrons.value_function(t,x,y,z,ux,uy,uz,w) = "w"
PhaseSpaceElectrons.filter_function(t,x,y,z,ux,uy,uz,w) = "sqrt(x*x+y*y) < 1e-6"

#################################
# Physical Background
#
# This example is modeled after a target similar to the hydrogen jet here:
#   [1] https://doi.org/10.1038/s41598-017-10589-3
#   [2] https://arxiv.org/abs/1903.06428
#
authors = "Axel Huebl <axelhuebl@lbl.gov>"
```

## Analyze

<a id="fig-tnsa-ps-electrons-pinhole"></a>
![Longitudinal phase space of forward-moving electrons in a 2 degree opening angle.](https://user-images.githubusercontent.com/5416860/295003882-c755fd47-4bb3-4439-9319-c48214cbaafd.png)

<a id="fig-tnsa-ps-protons-pinhole"></a>
![Longitudinal phase space of forward-moving protons in a 2 degree opening angle.](https://user-images.githubusercontent.com/5416860/295003988-dea3dfb7-0d55-4616-b32d-061fb429f9ac.png)

Time-resolved phase electron space analysis as in [Fig. 3](#fig-tnsa-ps-electrons-pinhole) gives information about, e.g., how laser energy is locally converted into electron kinetic energy.
Later in time, ion phase spaces like [Fig. 4](#fig-tnsa-ps-protons-pinhole) can reveal where accelerated ion populations originate.

### Script `analysis_histogram_2D.py`

```python3
#!/usr/bin/env python3

# This script displays a 2D histogram.

import argparse

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from openpmd_viewer import OpenPMDTimeSeries

parser = argparse.ArgumentParser(
    description="Process a 2D histogram name and an integer."
)
parser.add_argument("hist2D", help="Folder name of the reduced diagnostic.")
parser.add_argument(
    "iter",
    help="Iteration number of the simulation that is plotted. Enter a number from the list of iterations or 'All' if you want all plots.",
)
args = parser.parse_args()

path = "diags/reducedfiles/" + args.hist2D

ts = OpenPMDTimeSeries(path)

it = ts.iterations
data, info = ts.get_field(field="data", iteration=0, plot=True)
print("The available iterations of the simulation are:", it)
print("The axes of the histogram are (0: ordinate ; 1: abscissa):", info.axes)
print("The data shape is:", data.shape)

# Add the simulation time to the title once this information
# is available in the "info" FieldMetaInformation object.
if args.iter == "All":
    for it_idx, i in enumerate(it):
        plt.figure()
        data, info = ts.get_field(field="data", iteration=i, plot=False)
        abscissa_name = info.axes[1]  # This might be 'z' or something else
        abscissa_values = getattr(info, abscissa_name, None)
        ordinate_name = info.axes[0]  # This might be 'z' or something else
        ordinate_values = getattr(info, ordinate_name, None)

        plt.pcolormesh(
            abscissa_values / 1e-6,
            ordinate_values,
            data,
            norm=colors.LogNorm(),
            rasterized=True,
        )
        plt.title(args.hist2D + f" Time: {ts.t[it_idx]:.2e} s  (Iteration: {i:d})")
        plt.xlabel(info.axes[1] + r" ($\mu$m)")
        plt.ylabel(info.axes[0] + r" ($m_\mathrm{species} c$)")
        plt.colorbar()
        plt.tight_layout()
        plt.savefig("Histogram_2D_" + args.hist2D + "_iteration_" + str(i) + ".png")
else:
    i = int(args.iter)
    it_idx = np.where(i == it)[0][0]
    plt.figure()
    data, info = ts.get_field(field="data", iteration=i, plot=False)
    abscissa_name = info.axes[1]  # This might be 'z' or something else
    abscissa_values = getattr(info, abscissa_name, None)
    ordinate_name = info.axes[0]  # This might be 'z' or something else
    ordinate_values = getattr(info, ordinate_name, None)

    plt.pcolormesh(
        abscissa_values / 1e-6,
        ordinate_values,
        data,
        norm=colors.LogNorm(),
        rasterized=True,
    )
    plt.title(args.hist2D + f" Time: {ts.t[it_idx]:.2e} s  (Iteration: {i:d})")
    plt.xlabel(info.axes[1] + r" ($\mu$m)")
    plt.ylabel(info.axes[0] + r" ($m_\mathrm{species} c$)")
    plt.colorbar()
    plt.tight_layout()
    plt.savefig("Histogram_2D_" + args.hist2D + "_iteration_" + str(i) + ".png")
```

## Visualize

#### NOTE
The following images for densities and electromagnetic fields were created with a run on 64 NVidia A100 GPUs featuring a total number of cells of `nx = 8192` and `nz = 16384`, as well as 64 particles per cell per species.

<a id="fig-tnsa-densities"></a>
![Particle densities for electrons (top), protons (middle), and electrons again in logarithmic scale (bottom).](https://user-images.githubusercontent.com/5416860/296338802-8059c39c-0be8-4e4d-b41b-f976b626bd7f.png)

Particle density output illustrates the evolution of the target in time and space.
Logarithmic scales can help to identify where the target becomes transparent for the laser pulse (bottom panel in [Fig. 5](#fig-tnsa-densities) ).

<a id="fig-tnsa-fields"></a>
![Electromagnetic field visualization for E_x (top), B_y (middle), and E_z (bottom).](https://user-images.githubusercontent.com/5416860/296338609-a49eee7f-6793-4b55-92f1-0b887e6437ab.png)

Electromagnetic field output shows where the laser field is strongest at a given point in time, and where accelerating fields build up [Fig. 6](#fig-tnsa-fields).

### Script `plot_2d.py`

```python3
#!/usr/bin/env python3

# Copyright 2023 The WarpX Community
#
# This file is part of WarpX.
#
# Authors: Marco Garten
# License: BSD-3-Clause-LBNL
#
# This script plots the densities and fields of a 2D laser-ion acceleration simulation.


import argparse
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.constants as sc
from matplotlib.colors import TwoSlopeNorm
from openpmd_viewer import OpenPMDTimeSeries

plt.rcParams.update({"font.size": 16})


def create_analysis_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def visualize_density_iteration(ts, iteration, out_dir):
    """
    Visualize densities and fields of a single iteration.

    :param ts: OpenPMDTimeSeries
    :param iteration: Output iteration (simulation timestep)
    :param out_dir: Directory for PNG output
    :return:
    """
    # Physics parameters
    lambda_L = 800e-9  # Laser wavelength in meters
    omega_L = 2 * np.pi * sc.c / lambda_L  # Laser frequency in seconds
    n_c = (
        sc.m_e * sc.epsilon_0 * omega_L**2 / sc.elementary_charge**2
    )  # Critical plasma density in meters^(-3)
    micron = 1e-6

    # Simulation parameters
    n_e0 = 30
    n_max = 2 * n_e0
    nr = 1  # Number to decrease resolution

    # Data fetching
    it = iteration
    ii = np.where(ts.iterations == it)[0][0]

    time = ts.t[ii]
    rho_e, rho_e_info = ts.get_field(field="rho_electrons", iteration=it)
    rho_d, rho_d_info = ts.get_field(field="rho_hydrogen", iteration=it)

    # Rescale to critical density
    rho_e = rho_e / (sc.elementary_charge * n_c)
    rho_d = rho_d / (sc.elementary_charge * n_c)

    # Axes setup
    fig, axs = plt.subplots(3, 1, figsize=(5, 8))
    xax, zax = rho_e_info.x, rho_e_info.z

    # Plotting
    # Electron density
    im0 = axs[0].pcolormesh(
        zax[::nr] / micron,
        xax[::nr] / micron,
        -rho_e.T[::nr, ::nr],
        vmin=0,
        vmax=n_max,
        cmap="Reds",
        rasterized=True,
    )
    plt.colorbar(im0, ax=axs[0], label=r"$n_\mathrm{\,e}\ (n_\mathrm{c})$")

    # Hydrogen density
    im1 = axs[1].pcolormesh(
        zax[::nr] / micron,
        xax[::nr] / micron,
        rho_d.T[::nr, ::nr],
        vmin=0,
        vmax=n_max,
        cmap="Blues",
        rasterized=True,
    )
    plt.colorbar(im1, ax=axs[1], label=r"$n_\mathrm{\,H}\ (n_\mathrm{c})$")

    # Masked electron density
    divnorm = TwoSlopeNorm(vmin=-7.0, vcenter=0.0, vmax=2)
    masked_data = np.ma.masked_where(rho_e.T == 0, rho_e.T)
    my_cmap = plt.cm.PiYG_r.copy()
    my_cmap.set_bad(color="black")
    im2 = axs[2].pcolormesh(
        zax[::nr] / micron,
        xax[::nr] / micron,
        np.log(-masked_data[::nr, ::nr]),
        norm=divnorm,
        cmap=my_cmap,
        rasterized=True,
    )
    plt.colorbar(
        im2,
        ax=axs[2],
        ticks=[-6, -3, 0, 1, 2],
        extend="both",
        label=r"$\log n_\mathrm{\,e}\ (n_\mathrm{c})$",
    )

    # Axis labels and title
    for ax in axs:
        ax.set_aspect(1.0)
        ax.set_ylabel(r"$x$ ($\mu$m)")
    for ax in axs[:-1]:
        ax.set_xticklabels([])
    axs[2].set_xlabel(r"$z$ ($\mu$m)")
    fig.suptitle(f"Iteration: {it}, Time: {time / 1e-15:.1f} fs")

    plt.tight_layout()

    plt.savefig(f"{out_dir}/densities_{it:06d}.png")


def visualize_field_iteration(ts, iteration, out_dir):
    # Additional parameters
    nr = 1  # Number to decrease resolution
    micron = 1e-6

    # Data fetching
    it = iteration
    ii = np.where(ts.iterations == it)[0][0]
    time = ts.t[ii]

    Ex, Ex_info = ts.get_field(field="E", coord="x", iteration=it)
    Exmax = np.max(np.abs([np.min(Ex), np.max(Ex)]))
    By, By_info = ts.get_field(field="B", coord="y", iteration=it)
    Bymax = np.max(np.abs([np.min(By), np.max(By)]))
    Ez, Ez_info = ts.get_field(field="E", coord="z", iteration=it)
    Ezmax = np.max(np.abs([np.min(Ez), np.max(Ez)]))

    # Axes setup
    fig, axs = plt.subplots(3, 1, figsize=(5, 8))
    xax, zax = Ex_info.x, Ex_info.z

    # Plotting
    im0 = axs[0].pcolormesh(
        zax[::nr] / micron,
        xax[::nr] / micron,
        Ex.T[::nr, ::nr],
        vmin=-Exmax,
        vmax=Exmax,
        cmap="RdBu",
        rasterized=True,
    )

    plt.colorbar(im0, ax=axs[00], label=r"$E_x$ (V/m)")

    im1 = axs[1].pcolormesh(
        zax[::nr] / micron,
        xax[::nr] / micron,
        By.T[::nr, ::nr],
        vmin=-Bymax,
        vmax=Bymax,
        cmap="RdBu",
        rasterized=True,
    )
    plt.colorbar(im1, ax=axs[1], label=r"$B_y$ (T)")

    im2 = axs[2].pcolormesh(
        zax[::nr] / micron,
        xax[::nr] / micron,
        Ez.T[::nr, ::nr],
        vmin=-Ezmax,
        vmax=Ezmax,
        cmap="RdBu",
        rasterized=True,
    )
    plt.colorbar(im2, ax=axs[2], label=r"$E_z$ (V/m)")

    # Axis labels and title
    for ax in axs:
        ax.set_aspect(1.0)
        ax.set_ylabel(r"$x$ ($\mu$m)")
    for ax in axs[:-1]:
        ax.set_xticklabels([])
    axs[2].set_xlabel(r"$z$ ($\mu$m)")
    fig.suptitle(f"Iteration: {it}, Time: {time / 1e-15:.1f} fs")

    plt.tight_layout()

    plt.savefig(f"{out_dir}/fields_{it:06d}.png")


def visualize_particle_histogram_iteration(
    diag_name="histuH", species="hydrogen", iteration=1000, out_dir="./analysis"
):
    it = iteration

    if species == "hydrogen":
        # proton rest energy in eV
        mc2 = sc.m_p / sc.electron_volt * sc.c**2
    elif species == "electron":
        mc2 = sc.m_e / sc.electron_volt * sc.c**2
    else:
        raise NotImplementedError(
            "The only implemented presets for this analysis script are `electron` or `hydrogen`."
        )

    fs = 1.0e-15
    MeV = 1.0e6

    df = pd.read_csv(f"./diags/reducedfiles/{diag_name}.txt", delimiter=r"\s+")
    # the columns look like this:
    #     #[0]step() [1]time(s) [2]bin1=0.000220() [3]bin2=0.000660() [4]bin3=0.001100()

    # matches words, strings surrounded by " ' ", dots, minus signs and e for scientific notation in numbers
    nested_list = [re.findall(r"[\w'\.]+", col) for col in df.columns]

    index = pd.MultiIndex.from_tuples(
        nested_list, names=("column#", "name", "bin value")
    )

    df.columns = index

    steps = df.values[:, 0].astype(int)
    ii = np.where(steps == it)[0][0]
    time = df.values[:, 1]
    data = df.values[:, 2:]
    edge_vals = np.array([float(row[2]) for row in df.columns[2:]])
    edges_MeV = (np.sqrt(edge_vals**2 + 1) - 1) * mc2 / MeV

    time_fs = time / fs

    fig, ax = plt.subplots(1, 1)

    ax.plot(edges_MeV, data[ii, :])
    ax.set_yscale("log")
    ax.set_ylabel(r"d$N$/d$\mathcal{E}$ (arb. u.)")
    ax.set_xlabel(r"$\mathcal{E}$ (MeV)")

    fig.suptitle(f"{species} - Iteration: {it}, Time: {time_fs[ii]:.1f} fs")

    plt.tight_layout()
    plt.savefig(f"./{out_dir}/{diag_name}_{it:06d}.png")


if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(
        description="Visualize Laser-Ion Accelerator Densities and Fields"
    )
    parser.add_argument(
        "-d",
        "--diag_dir",
        type=str,
        default="./diags/diagInst",
        help="Directory containing density and field diagnostics",
    )
    parser.add_argument(
        "-i",
        "--iteration",
        type=int,
        default=None,
        help="Specific iteration to visualize",
    )
    parser.add_argument(
        "-hn",
        "--histogram_name",
        type=str,
        default="histuH",
        help="Name of histogram diagnostic to visualize",
    )
    parser.add_argument(
        "-hs",
        "--histogram_species",
        type=str,
        default="hydrogen",
        help="Particle species in the visualized histogram diagnostic",
    )
    args = parser.parse_args()

    # Create analysis directory
    analysis_dir = "analysis"
    create_analysis_dir(analysis_dir)

    # Loading the time series
    ts = OpenPMDTimeSeries(args.diag_dir)

    if args.iteration is not None:
        visualize_density_iteration(ts, args.iteration, analysis_dir)
        visualize_field_iteration(ts, args.iteration, analysis_dir)
        visualize_particle_histogram_iteration(
            args.histogram_name, args.histogram_species, args.iteration, analysis_dir
        )
    else:
        for it in ts.iterations:
            visualize_density_iteration(ts, it, analysis_dir)
            visualize_field_iteration(ts, it, analysis_dir)
            visualize_particle_histogram_iteration(
                args.histogram_name, args.histogram_species, it, analysis_dir
            )
```
