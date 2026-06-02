<a id="usage-picmi"></a>

<a id="usage-picmi-run"></a>

# Inputs: PICMI Python Script

This documents on how to use WarpX as a Python script (e.g., `python3 PICMI_script.py`).

WarpX uses the [PICMI standard](https://github.com/picmi-standard/picmi) for its Python input files.
Complete example input files can be found in [the examples section](examples.md#usage-examples).

In the input file, instances of classes are created defining the various aspects of the simulation.
A variable of type [`pywarpx.picmi.Simulation`](#pywarpx.picmi.Simulation) is the central object to which all other options are passed, defining the simulation time, field solver, registered species, etc.

Once the simulation is fully configured, it can be used in one of two modes.
**Interactive** use is the most common and can be [extended with custom runtime functionality](workflows/python_extend.md#usage-python-extend):

### Interactive

[`step()`](#pywarpx.picmi.Simulation.step): run directly from Python

### Preprocessor

[`write_input_file()`](#pywarpx.picmi.Simulation.write_input_file): create an [inputs file for a WarpX executable](parameters.md#running-cpp-parameters)

When run directly from Python, one can also extend WarpX with further custom user logic.
See the [detailed workflow page](workflows/python_extend.md#usage-python-extend) on how to extend WarpX from Python.

<a id="usage-picmi-parameters"></a>

## Simulation and Grid Setup

### *class* pywarpx.picmi.Simulation(solver=None, time_step_size=None, max_steps=None, max_time=None, verbose=None, particle_shape='linear', gamma_boost=None, load_balancing=None, \*\*kw)

Creates a Simulation object

* **Parameters:**
  * **solver** (*field solver instance*) – This is the field solver to be used in the simulation.
    It should be an instance of field solver classes.
  * **time_step_size** (*float*) – Absolute time step size of the simulation [s].
    Needed if the CFL is not specified elsewhere.
  * **max_steps** (*integer*) – Maximum number of time steps.
    Specify either this, or max_time, or use the step function directly.
  * **max_time** (*float*) – Maximum physical time to run the simulation [s].
    Specify either this, or max_steps, or use the step function directly.
  * **verbose** (*integer* *,* *optional*) – Verbosity flag. A larger integer results in more verbose output
  * **particle_shape** ( *{'NGP'* *,*  *'linear'* *,*  *'quadratic'* *,*  *'cubic'}*) – Default particle shape for species added to this simulation
  * **gamma_boost** – Lorentz factor of the boosted simulation frame.
    Note that all input values should be in the lab frame.

See [Input Parameters](https://warpx.readthedocs.io/en/latest/usage/parameters.html) for more information.

* **Parameters:**
  * **warpx_evolve_scheme** (*solver scheme instance* *,* *optional*) – Which evolve scheme to use
  * **warpx_current_deposition_algo** ( *{'direct'* *,*  *'esirkepov'* *,* *and 'vay'}* *,* *optional*) – Current deposition algorithm. The default depends on conditions.
  * **warpx_charge_deposition_algo** ( *{'standard'}* *,* *optional*) – Charge deposition algorithm.
  * **warpx_field_gathering_algo** ( *{'energy-conserving'* *,*  *'momentum-conserving'}* *,* *optional*) – Field gathering algorithm. The default depends on conditions.
  * **warpx_particle_pusher_algo** ( *{'boris'* *,*  *'vay'* *,*  *'higuera'}* *,* *default='boris'*) – Particle pushing algorithm.
  * **warpx_use_filter** (*bool* *,* *optional*) – Whether to use filtering. The default depends on the conditions.
  * **warpx_grid_type** ( *{'collocated'* *,*  *'staggered'* *,*  *'hybrid'}* *,* *default='staggered'*) – Whether to use a collocated grid (all fields defined at the cell nodes),
    a staggered grid (fields defined on a Yee grid), or a hybrid grid
    (fields and currents are interpolated back and forth between a staggered grid
    and a collocated grid, must be used with momentum-conserving field gathering algorithm).
  * **warpx_do_current_centering** (*bool* *,* *optional*) – If true, the current is deposited on a nodal grid and then centered
    to a staggered grid (Yee grid), using finite-order interpolation.
    Default: warpx.do_current_centering=0 with collocated or staggered grids,
    warpx.do_current_centering=1 with hybrid grids.
  * **warpx_field_centering_nox/noy/noz** (*integer* *,* *optional*) – The order of interpolation used with staggered or hybrid grids (`warpx_grid_type=staggered`
    or `warpx_grid_type=hybrid`) and momentum-conserving field gathering
    (`warpx_field_gathering_algo=momentum-conserving`) to interpolate the
    electric and magnetic fields from the cell centers to the cell nodes,
    before gathering the fields from the cell nodes to the particle positions.
    Default: `warpx_field_centering_no<x,y,z>=2` with staggered grids,
    `warpx_field_centering_no<x,y,z>=8` with hybrid grids (typically necessary
    to ensure stability in boosted-frame simulations of relativistic plasmas and beams).
  * **warpx_current_centering_nox/noy/noz** (*integer* *,* *optional*) – The order of interpolation used with hybrid grids (`warpx_grid_type=hybrid`)
    to interpolate the currents from the cell nodes to the cell centers when
    `warpx_do_current_centering=1`, before pushing the Maxwell fields on staggered grids.
    Default: `warpx_current_centering_no<x,y,z>=8` with hybrid grids (typically necessary
    to ensure stability in boosted-frame simulations of relativistic plasmas and beams).
  * **warpx_serialize_initial_conditions** (*bool* *,* *default=False*) – Controls the random numbers used for initialization.
    This parameter should only be used for testing and continuous integration.
  * **warpx_random_seed** (*string* *or* *int* *,* *optional*) – (See documentation)
  * **warpx_do_dynamic_scheduling** (*bool* *,* *default=True*) – Whether to do dynamic scheduling with OpenMP
  * **warpx_roundrobin_sfc** (*bool* *,* *default=False*) – Whether to use the RRSFC strategy for making DistributionMapping
  * **warpx_load_balance_intervals** (*string* *,* *default='0'*) – The intervals for doing load balancing
  * **warpx_load_balance_efficiency_ratio_threshold** (*float* *,* *default=1.1*) – (See documentation)
  * **warpx_load_balance_with_sfc** (*bool* *,* *default=0*) – (See documentation)
  * **warpx_load_balance_knapsack_factor** (*float* *,* *default=1.24*) – (See documentation)
  * **warpx_load_balance_costs_update** ( *{'heuristic'* *or*  *'timers'}* *,* *optional*) – (See documentation)
  * **warpx_costs_heuristic_particles_wt** (*float* *,* *optional*) – (See documentation)
  * **warpx_costs_heuristic_cells_wt** (*float* *,* *optional*) – (See documentation)
  * **warpx_use_fdtd_nci_corr** (*bool* *,* *optional*) – Whether to use the NCI correction when using the FDTD solver
  * **warpx_amr_check_input** (*bool* *,* *optional*) – Whether AMReX should perform checks on the input
    (primarily related to the max grid size and blocking factors)
  * **warpx_amr_restart** (*string* *,* *optional*) – The name of the restart to use
  * **warpx_amrex_the_arena_is_managed** (*bool* *,* *optional*) – Whether to use managed memory in the AMReX Arena
  * **warpx_amrex_the_arena_init_size** (*long int* *,* *optional*) – The amount of memory in bytes to allocate in the Arena.
  * **warpx_amrex_use_gpu_aware_mpi** (*bool* *,* *optional*) – Whether to use GPU-aware MPI communications
  * **warpx_do_device_synchronize** (*bool* *,* *optional*) – Whether to synchronize GPU threads at ends of profiling regions.
    Note that if this is set to False, the TinyProfiler table can be
    misleading.
  * **warpx_zmax_plasma_to_compute_max_step** (*float* *,* *optional*) – Sets the simulation run time based on the maximum z value
  * **warpx_compute_max_step_from_btd** (*bool* *,* *default=0*) – If specified, automatically calculates the number of iterations
    required in the boosted frame for all back-transformed diagnostics
    to be completed.
  * **warpx_collisions** (*collision instance* *,* *optional*) – The collision instance specifying the particle collisions
  * **warpx_collisions_split_momentum_push** (*bool* *,* *default=1*) – If true, collisions are performed in the middle of the momentum push,
    which is split into two substeps.
    This improves energy conservation, as demonstrated in
    (Vay et al., Phys. Rev. E 111, 2025).
    This is only implemented for the explicit evolve scheme
    and is not available for the implicit evolve schemes.
  * **warpx_embedded_boundary** (*embedded boundary instance* *,* *optional*) – 
  * **warpx_break_signals** (*list* *of* *strings*) – Signals on which to break
  * **warpx_checkpoint_signals** (*list* *of* *strings*) – Signals on which to write out a checkpoint
  * **warpx_synchronize_velocity** (*bool* *,* *default=False*) – Flags whether the particle velocities are synchronized in time with
    the positions in the diagnostics. When False, the particles are
    one half step behind the positions (except for the final diagnostic).
  * **warpx_numprocs** (*list* *of* *ints* *(**1 in 1D* *,* *2 in 2D* *,* *3 in 3D* *)*) – Domain decomposition on the coarsest level.
    The domain will be chopped into the exact number of pieces in each dimension as specified by this parameter.
    [https://warpx.readthedocs.io/en/latest/usage/parameters.html#distribution-across-mpi-ranks-and-parallelization](https://warpx.readthedocs.io/en/latest/usage/parameters.html#distribution-across-mpi-ranks-and-parallelization)
    [https://warpx.readthedocs.io/en/latest/usage/domain_decomposition.html#simple-method](https://warpx.readthedocs.io/en/latest/usage/domain_decomposition.html#simple-method)
  * **warpx_sort_intervals** (*string* *,* *optional* *(**defaults: -1 on CPU; 4 on GPU* *)*) – Using the Intervals parser syntax, this string defines the timesteps at which particles are sorted. If <=0, do not sort particles.
    It is turned on on GPUs for performance reasons (to improve memory locality).
  * **warpx_sort_particles_for_deposition** (*bool* *,* *optional* *(**default: true for the CUDA backend* *,* *otherwise false* *)*) – This option controls the type of sorting used if particle sorting is turned on, i.e. if sort_intervals is not <=0.
    If true, particles will be sorted by cell to optimize deposition with many particles per cell, in the order x -> y -> z -> ppc.
    If false, particles will be sorted by bin, using the sort_bin_size parameter below, in the order ppc -> x -> y -> z.
    true is recommended for best performance on NVIDIA GPUs, especially if there are many particles per cell.
  * **warpx_sort_idx_type** (*list* *of* *int* *,* *optional* *(**default: 0 0 0* *)*) – 

    This controls the type of grid used to sort the particles when sort_particles_for_deposition is true.
    Possible values are:
    * idx_type = {0, 0, 0}: Sort particles to a cell centered grid,
    * idx_type = {1, 1, 1}: Sort particles to a node centered grid,
    * idx_type = {2, 2, 2}: Compromise between a cell and node centered grid.

    In 2D (XZ and RZ), only the first two elements are read. In 1D, only the first element is read.
  * **warpx_sort_bin_size** (*list* *of* *int* *,* *optional* *(**default 1 1 1* *)*) – If sort_intervals is activated and sort_particles_for_deposition is false, particles are sorted in bins of sort_bin_size cells.
    In 2D, only the first two elements are read.
  * **warpx_used_inputs_file** (*string* *,* *optional*) – The name of the text file that the used input parameters is written to,
  * **warpx_reduced_diags_path** (*string* *,* *optional*) – Sets the default path for reduced diagnostic output files
  * **warpx_reduced_diags_extension** (*string* *,* *optional*) – Sets the default extension for reduced diagnostic output files
  * **warpx_reduced_diags_intervals** (*string* *,* *optional*) – Sets the default intervals for reduced diagnostic output files
  * **warpx_reduced_diags_separator** (*string* *,* *optional*) – Sets the default separator for reduced diagnostic output files
  * **warpx_reduced_diags_precision** (*integer* *,* *optional*) – Sets the default precision for reduced diagnostic output files

#### add_applied_field(applied_field)

Add an applied field

* **Parameters:**
  **applied_field** (*applied field instance*) – One of the applied field instance.
  Specifies the properties of the applied field.

#### add_diagnostic(diagnostic)

Add a diagnostic

* **Parameters:**
  **diagnostic** (*diagnostic instance*) – One of the diagnostic instances.

#### add_interaction(interaction)

Add an interaction

* **Parameters:**
  **interaction** (*interaction instance*) – One of the interaction objects.

#### add_laser(laser, injection_method)

Add a laser pulses that to be injected in the simulation

* **Parameters:**
  * **laser_profile** (*laser instance*) – One of laser profile instances.
    Specifies the **physical** properties of the laser pulse
    (e.g. spatial and temporal profile, wavelength, amplitude, etc.).
  * **injection_method** (*laser injection instance* *,* *optional*) – Specifies how the laser is injected (numerically) into the simulation
    (e.g. through a laser antenna, or directly added to the mesh).
    This argument describes an **algorithm**, not a physical object.
    It is up to each code to define the default method
    of injection, if the user does not provide injection_method.

#### add_species(species, layout, initialize_self_field=None)

Add species to be used in the simulation

* **Parameters:**
  * **species** (*species instance*) – An instance of one of the PICMI species objects.
    Defines species to be added from the *physical* point of view
    (e.g. charge, mass, initial distribution of particles).
  * **layout** (*layout instance*) – An instance of one of the PICMI particle layout objects.
    Defines how particles are added into the simulation, from the *numerical* point of view.
  * **initialize_self_field** (*bool* *,* *optional*) – Whether the initial space-charge fields of this species
    is calculated and added to the simulation

#### step(nsteps=None, mpi_comm=None)

Run the simulation for nsteps timesteps

* **Parameters:**
  **nsteps** (*integer* *,* *default=1*) – The number of timesteps

#### write_input_file(file_name='inputs')

Write the parameters of the simulation, as defined in the PICMI input,
into a code-specific input file.

This can be used for codes that are not Python-driven (e.g. compiled,
pure C++ or Fortran codes) and expect a text input in a given format.

* **Parameters:**
  **file_name** (*string*) – The path to the file that will be created

### *class* pywarpx.picmi.Cartesian3DGrid(number_of_cells=None, lower_bound=None, upper_bound=None, lower_boundary_conditions=None, upper_boundary_conditions=None, nx=None, ny=None, nz=None, xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None, bc_xmin=None, bc_xmax=None, bc_ymin=None, bc_ymax=None, bc_zmin=None, bc_zmax=None, moving_window_velocity=None, refined_regions=[], lower_bound_particles=None, upper_bound_particles=None, xmin_particles=None, xmax_particles=None, ymin_particles=None, ymax_particles=None, zmin_particles=None, zmax_particles=None, lower_boundary_conditions_particles=None, upper_boundary_conditions_particles=None, bc_xmin_particles=None, bc_xmax_particles=None, bc_ymin_particles=None, bc_ymax_particles=None, bc_zmin_particles=None, bc_zmax_particles=None, guard_cells=None, pml_cells=None, \*\*kw)

Three dimensional Cartesian grid
Parameters can be specified either as vectors or separately.
(If both are specified, the vector is used.)

* **Parameters:**
  * **number_of_cells** (*vector* *of* *integers*) – Number of cells along each axis (number of nodes is number_of_cells+1)
  * **lower_bound** (*vector* *of* *floats*) – Position of the node at the lower bound [m]
  * **upper_bound** (*vector* *of* *floats*) – Position of the node at the upper bound [m]
  * **lower_boundary_conditions** (*vector* *of* *strings*) – Conditions at lower boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann
  * **upper_boundary_conditions** (*vector* *of* *strings*) – Conditions at upper boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann
  * **nx** (*integer*) – Number of cells along X (number of nodes=nx+1)
  * **ny** (*integer*) – Number of cells along Y (number of nodes=ny+1)
  * **nz** (*integer*) – Number of cells along Z (number of nodes=nz+1)
  * **xmin** (*float*) – Position of first node along X [m]
  * **xmax** (*float*) – Position of last node along X [m]
  * **ymin** (*float*) – Position of first node along Y [m]
  * **ymax** (*float*) – Position of last node along Y [m]
  * **zmin** (*float*) – Position of first node along Z [m]
  * **zmax** (*float*) – Position of last node along Z [m]
  * **bc_xmin** (*string*) – Boundary condition at min X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann
  * **bc_xmax** (*string*) – Boundary condition at max X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann
  * **bc_ymin** (*string*) – Boundary condition at min Y: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann
  * **bc_ymax** (*string*) – Boundary condition at max Y: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann
  * **bc_zmin** (*string*) – Boundary condition at min Z: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann
  * **bc_zmax** (*string*) – Boundary condition at max Z: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann
  * **moving_window_velocity** (*vector* *of* *floats* *,* *optional*) – Moving frame velocity [m/s]
  * **refined_regions** (*list* *of* *lists* *,* *optional*) – List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor],
    with level being the refinement level, with 1 being the first level of refinement, 2 being the second etc,
    lo and hi being vectors of length 3 specifying the extent of the region,
    and refinement_factor defaulting to [2,2,2] (relative to next lower level)
  * **lower_bound_particles** (*vector* *of* *floats* *,* *optional*) – Position of particle lower bound [m]
  * **upper_bound_particles** (*vector* *of* *floats* *,* *optional*) – Position of particle upper bound [m]
  * **xmin_particles** (*float* *,* *optional*) – Position of min particle boundary along X [m]
  * **xmax_particles** (*float* *,* *optional*) – Position of max particle boundary along X [m]
  * **ymin_particles** (*float* *,* *optional*) – Position of min particle boundary along Y [m]
  * **ymax_particles** (*float* *,* *optional*) – Position of max particle boundary along Y [m]
  * **float** (*zmin_particles*) – Position of min particle boundary along Z [m]
  * **optional** – Position of min particle boundary along Z [m]
  * **zmax_particles** (*float* *,* *optional*) – Position of max particle boundary along Z [m]
  * **lower_boundary_conditions_particles** (*vector* *of* *strings* *,* *optional*) – Conditions at lower boundaries for particles, periodic, absorbing, reflect or thermal
  * **upper_boundary_conditions_particles** (*vector* *of* *strings* *,* *optional*) – Conditions at upper boundaries for particles, periodic, absorbing, reflect or thermal
  * **bc_xmin_particles** (*string* *,* *optional*) – Boundary condition at min X for particles: One of periodic, absorbing, reflect, thermal
  * **bc_xmax_particles** (*string* *,* *optional*) – Boundary condition at max X for particles: One of periodic, absorbing, reflect, thermal
  * **bc_ymin_particles** (*string* *,* *optional*) – Boundary condition at min Y for particles: One of periodic, absorbing, reflect, thermal
  * **bc_ymax_particles** (*string* *,* *optional*) – Boundary condition at max Y for particles: One of periodic, absorbing, reflect, thermal
  * **bc_zmin_particles** (*string* *,* *optional*) – Boundary condition at min Z for particles: One of periodic, absorbing, reflect, thermal
  * **bc_zmax_particles** (*string* *,* *optional*) – Boundary condition at max Z for particles: One of periodic, absorbing, reflect, thermal
  * **guard_cells** (*vector* *of* *integers* *,* *optional*) – Number of guard cells used along each direction
  * **pml_cells** (*vector* *of* *integers* *,* *optional*) – Number of Perfectly Matched Layer (PML) cells along each direction

### References

absorbing_silver_mueller: A local absorbing boundary condition that works best under normal incidence angle.
Based on the Silver-Mueller Radiation Condition, e.g., in

* A. K. Belhora and L. Pichon, “Maybe Efficient Absorbing Boundary Conditions for the Finite Element Solution of 3D Scattering Problems,” 1995,
  [https://doi.org/10.1109/20.376322](https://doi.org/10.1109/20.376322)
* B Engquist and A. Majdat, “Absorbing boundary conditions for numerical simulation of waves,” 1977,
  [https://doi.org/10.1073/pnas.74.5.1765](https://doi.org/10.1073/pnas.74.5.1765)
* R. Lehe, “Electromagnetic wave propagation in Particle-In-Cell codes,” 2016,
  US Particle Accelerator School (USPAS) Summer Session, Self-Consistent Simulations of Beam and Plasma Systems
  [https://people.nscl.msu.edu/~lund/uspas/scs_2016/lec_adv/A1b_EM_Waves.pdf](https://people.nscl.msu.edu/~lund/uspas/scs_2016/lec_adv/A1b_EM_Waves.pdf)
  > Implementation specific documentation

See [Input Parameters](https://warpx.readthedocs.io/en/latest/usage/parameters.html) for more information.

* **Parameters:**
  * **warpx_max_grid_size** (*integer* *,* *default=32*) – Maximum block size in either direction
  * **warpx_max_grid_size_x** (*integer* *,* *optional*) – Maximum block size in x direction
  * **warpx_max_grid_size_y** (*integer* *,* *optional*) – Maximum block size in z direction
  * **warpx_max_grid_size_z** (*integer* *,* *optional*) – Maximum block size in z direction
  * **warpx_blocking_factor** (*integer* *,* *optional*) – Blocking factor (which controls the block size)
  * **warpx_blocking_factor_x** (*integer* *,* *optional*) – Blocking factor (which controls the block size) in the x direction
  * **warpx_blocking_factor_y** (*integer* *,* *optional*) – Blocking factor (which controls the block size) in the z direction
  * **warpx_blocking_factor_z** (*integer* *,* *optional*) – Blocking factor (which controls the block size) in the z direction
  * **warpx_potential_lo_x** (*float* *,* *default=0.*) – Electrostatic potential on the lower x boundary
  * **warpx_potential_hi_x** (*float* *,* *default=0.*) – Electrostatic potential on the upper x boundary
  * **warpx_potential_lo_y** (*float* *,* *default=0.*) – Electrostatic potential on the lower z boundary
  * **warpx_potential_hi_y** (*float* *,* *default=0.*) – Electrostatic potential on the upper z boundary
  * **warpx_potential_lo_z** (*float* *,* *default=0.*) – Electrostatic potential on the lower z boundary
  * **warpx_potential_hi_z** (*float* *,* *default=0.*) – Electrostatic potential on the upper z boundary
  * **warpx_start_moving_window_step** (*int* *,* *default=0*) – The timestep at which the moving window starts
  * **warpx_end_moving_window_step** (*int* *,* *default=-1*) – The timestep at which the moving window ends. If -1, the moving window
    will continue until the end of the simulation.
  * **warpx_boundary_u_th** (*dict* *,* *default=None*) – If a thermal boundary is used for particles, this dictionary should
    specify the thermal speed for each species in the form {<species>: u_th}.
    Note: u_th = sqrt(T\*q_e/mass)/clight with T in eV.

### *class* pywarpx.picmi.Cartesian2DGrid(number_of_cells=None, lower_bound=None, upper_bound=None, lower_boundary_conditions=None, upper_boundary_conditions=None, nx=None, ny=None, xmin=None, xmax=None, ymin=None, ymax=None, bc_xmin=None, bc_xmax=None, bc_ymin=None, bc_ymax=None, moving_window_velocity=None, refined_regions=[], lower_bound_particles=None, upper_bound_particles=None, xmin_particles=None, xmax_particles=None, ymin_particles=None, ymax_particles=None, lower_boundary_conditions_particles=None, upper_boundary_conditions_particles=None, bc_xmin_particles=None, bc_xmax_particles=None, bc_ymin_particles=None, bc_ymax_particles=None, guard_cells=None, pml_cells=None, \*\*kw)

Two dimensional Cartesian grid
Parameters can be specified either as vectors or separately.
(If both are specified, the vector is used.)

* **Parameters:**
  * **number_of_cells** (*vector* *of* *integers*) – Number of cells along each axis (number of nodes is number_of_cells+1)
  * **lower_bound** (*vector* *of* *floats*) – Position of the node at the lower bound [m]
  * **upper_bound** (*vector* *of* *floats*) – Position of the node at the upper bound [m]
  * **lower_boundary_conditions** (*vector* *of* *strings*) – Conditions at lower boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann
  * **upper_boundary_conditions** (*vector* *of* *strings*) – Conditions at upper boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann
  * **nx** (*integer*) – Number of cells along X (number of nodes=nx+1)
  * **ny** (*integer*) – Number of cells along Y (number of nodes=ny+1)
  * **xmin** (*float*) – Position of first node along X [m]
  * **xmax** (*float*) – Position of last node along X [m]
  * **ymin** (*float*) – Position of first node along Y [m]
  * **ymax** (*float*) – Position of last node along Y [m]
  * **bc_xmin** (*vector* *of* *strings*) – Boundary condition at min X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann
  * **bc_xmax** (*vector* *of* *strings*) – Boundary condition at max X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann
  * **bc_ymin** (*vector* *of* *strings*) – Boundary condition at min Y: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann
  * **bc_ymax** (*vector* *of* *strings*) – Boundary condition at max Y: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann
  * **moving_window_velocity** (*vector* *of* *floats* *,* *optional*) – Moving frame velocity [m/s]
  * **refined_regions** (*list* *of* *lists* *,* *optional*) – List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor],
    with level being the refinement level, with 1 being the first level of refinement, 2 being the second etc,
    lo and hi being vectors of length 2 specifying the extent of the region,
    and refinement_factor defaulting to [2,2] (relative to next lower level)
  * **lower_bound_particles** (*vector* *of* *floats* *,* *optional*) – Position of particle lower bound [m]
  * **upper_bound_particles** (*vector* *of* *floats* *,* *optional*) – Position of particle upper bound [m]
  * **xmin_particles** (*float* *,* *optional*) – Position of min particle boundary along X [m]
  * **xmax_particles** (*float* *,* *optional*) – Position of max particle boundary along X [m]
  * **ymin_particles** (*float* *,* *optional*) – Position of min particle boundary along Y [m]
  * **ymax_particles** (*float* *,* *optional*) – Position of max particle boundary along Y [m]
  * **lower_boundary_conditions_particles** (*vector* *of* *strings* *,* *optional*) – Conditions at lower boundaries for particles, periodic, absorbing, reflect or thermal
  * **upper_boundary_conditions_particles** (*vector* *of* *strings* *,* *optional*) – Conditions at upper boundaries for particles, periodic, absorbing, reflect or thermal
  * **bc_xmin_particles** (*string* *,* *optional*) – Boundary condition at min X for particles: One of periodic, absorbing, reflect, thermal
  * **bc_xmax_particles** (*string* *,* *optional*) – Boundary condition at max X for particles: One of periodic, absorbing, reflect, thermal
  * **bc_ymin_particles** (*string* *,* *optional*) – Boundary condition at min Y for particles: One of periodic, absorbing, reflect, thermal
  * **bc_ymax_particles** (*string* *,* *optional*) – Boundary condition at max Y for particles: One of periodic, absorbing, reflect, thermal
  * **guard_cells** (*vector* *of* *integers* *,* *optional*) – Number of guard cells used along each direction
  * **pml_cells** (*vector* *of* *integers* *,* *optional*) – Number of Perfectly Matched Layer (PML) cells along each direction

### References

absorbing_silver_mueller: A local absorbing boundary condition that works best under normal incidence angle.
Based on the Silver-Mueller Radiation Condition, e.g., in

* A. K. Belhora and L. Pichon, “Maybe Efficient Absorbing Boundary Conditions for the Finite Element Solution of 3D Scattering Problems,” 1995,
  [https://doi.org/10.1109/20.376322](https://doi.org/10.1109/20.376322)
* B Engquist and A. Majdat, “Absorbing boundary conditions for numerical simulation of waves,” 1977,
  [https://doi.org/10.1073/pnas.74.5.1765](https://doi.org/10.1073/pnas.74.5.1765)
* R. Lehe, “Electromagnetic wave propagation in Particle-In-Cell codes,” 2016,
  US Particle Accelerator School (USPAS) Summer Session, Self-Consistent Simulations of Beam and Plasma Systems
  [https://people.nscl.msu.edu/~lund/uspas/scs_2016/lec_adv/A1b_EM_Waves.pdf](https://people.nscl.msu.edu/~lund/uspas/scs_2016/lec_adv/A1b_EM_Waves.pdf)
  > Implementation specific documentation

See [Input Parameters](https://warpx.readthedocs.io/en/latest/usage/parameters.html) for more information.

* **Parameters:**
  * **warpx_max_grid_size** (*integer* *,* *default=32*) – Maximum block size in either direction
  * **warpx_max_grid_size_x** (*integer* *,* *optional*) – Maximum block size in x direction
  * **warpx_max_grid_size_y** (*integer* *,* *optional*) – Maximum block size in z direction
  * **warpx_blocking_factor** (*integer* *,* *optional*) – Blocking factor (which controls the block size)
  * **warpx_blocking_factor_x** (*integer* *,* *optional*) – Blocking factor (which controls the block size) in the x direction
  * **warpx_blocking_factor_y** (*integer* *,* *optional*) – Blocking factor (which controls the block size) in the z direction
  * **warpx_potential_lo_x** (*float* *,* *default=0.*) – Electrostatic potential on the lower x boundary
  * **warpx_potential_hi_x** (*float* *,* *default=0.*) – Electrostatic potential on the upper x boundary
  * **warpx_potential_lo_z** (*float* *,* *default=0.*) – Electrostatic potential on the lower z boundary
  * **warpx_potential_hi_z** (*float* *,* *default=0.*) – Electrostatic potential on the upper z boundary
  * **warpx_start_moving_window_step** (*int* *,* *default=0*) – The timestep at which the moving window starts
  * **warpx_end_moving_window_step** (*int* *,* *default=-1*) – The timestep at which the moving window ends. If -1, the moving window
    will continue until the end of the simulation.
  * **warpx_boundary_u_th** (*dict* *,* *default=None*) – If a thermal boundary is used for particles, this dictionary should
    specify the thermal speed for each species in the form {<species>: u_th}.
    Note: u_th = sqrt(T\*q_e/mass)/clight with T in eV.

### *class* pywarpx.picmi.Cartesian1DGrid(number_of_cells=None, lower_bound=None, upper_bound=None, lower_boundary_conditions=None, upper_boundary_conditions=None, nx=None, xmin=None, xmax=None, bc_xmin=None, bc_xmax=None, moving_window_velocity=None, refined_regions=[], lower_bound_particles=None, upper_bound_particles=None, xmin_particles=None, xmax_particles=None, lower_boundary_conditions_particles=None, upper_boundary_conditions_particles=None, bc_xmin_particles=None, bc_xmax_particles=None, guard_cells=None, pml_cells=None, \*\*kw)

One-dimensional Cartesian grid
Parameters can be specified either as vectors or separately.
(If both are specified, the vector is used.)

* **Parameters:**
  * **number_of_cells** (*vector* *of* *integers*) – Number of cells along each axis (number of nodes is number_of_cells+1)
  * **lower_bound** (*vector* *of* *floats*) – Position of the node at the lower bound [m]
  * **upper_bound** (*vector* *of* *floats*) – Position of the node at the upper bound [m]
  * **lower_boundary_conditions** (*vector* *of* *strings*) – Conditions at lower boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann
  * **upper_boundary_conditions** (*vector* *of* *strings*) – Conditions at upper boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann
  * **nx** (*integer*) – Number of cells along X (number of nodes=nx+1)
  * **xmin** (*float*) – Position of first node along X [m]
  * **xmax** (*float*) – Position of last node along X [m]
  * **bc_xmin** (*vector* *of* *strings*) – Boundary condition at min X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann
  * **bc_xmax** (*vector* *of* *strings*) – Boundary condition at max X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann
  * **moving_window_velocity** (*vector* *of* *floats* *,* *optional*) – Moving frame velocity [m/s]
  * **refined_regions** (*list* *of* *lists* *,* *optional*) – List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor],
    with level being the refinement level, with 1 being the first level of refinement, 2 being the second etc,
    lo and hi being vectors of length 2 specifying the extent of the region,
    and refinement_factor defaulting to [2,2] (relative to next lower level)
  * **lower_bound_particles** (*vector* *of* *floats* *,* *optional*) – Position of particle lower bound [m]
  * **upper_bound_particles** (*vector* *of* *floats* *,* *optional*) – Position of particle upper bound [m]
  * **xmin_particles** (*float* *,* *optional*) – Position of min particle boundary along X [m]
  * **xmax_particles** (*float* *,* *optional*) – Position of max particle boundary along X [m]
  * **lower_boundary_conditions_particles** (*vector* *of* *strings* *,* *optional*) – Conditions at lower boundaries for particles, periodic, absorbing, reflect or thermal
  * **upper_boundary_conditions_particles** (*vector* *of* *strings* *,* *optional*) – Conditions at upper boundaries for particles, periodic, absorbing, reflect or thermal
  * **bc_xmin_particles** (*string* *,* *optional*) – Boundary condition at min X for particles: One of periodic, absorbing, reflect, thermal
  * **bc_xmax_particles** (*string* *,* *optional*) – Boundary condition at max X for particles: One of periodic, absorbing, reflect, thermal
  * **guard_cells** (*vector* *of* *integers* *,* *optional*) – Number of guard cells used along each direction
  * **pml_cells** (*vector* *of* *integers* *,* *optional*) – Number of Perfectly Matched Layer (PML) cells along each direction

### References

absorbing_silver_mueller: A local absorbing boundary condition that works best under normal incidence angle.
Based on the Silver-Mueller Radiation Condition, e.g., in

* A. K. Belhora and L. Pichon, “Maybe Efficient Absorbing Boundary Conditions for the Finite Element Solution of 3D Scattering Problems,” 1995,
  [https://doi.org/10.1109/20.376322](https://doi.org/10.1109/20.376322)
* B Engquist and A. Majdat, “Absorbing boundary conditions for numerical simulation of waves,” 1977,
  [https://doi.org/10.1073/pnas.74.5.1765](https://doi.org/10.1073/pnas.74.5.1765)
* R. Lehe, “Electromagnetic wave propagation in Particle-In-Cell codes,” 2016,
  US Particle Accelerator School (USPAS) Summer Session, Self-Consistent Simulations of Beam and Plasma Systems
  [https://people.nscl.msu.edu/~lund/uspas/scs_2016/lec_adv/A1b_EM_Waves.pdf](https://people.nscl.msu.edu/~lund/uspas/scs_2016/lec_adv/A1b_EM_Waves.pdf)
  > Implementation specific documentation

See [Input Parameters](https://warpx.readthedocs.io/en/latest/usage/parameters.html) for more information.

* **Parameters:**
  * **warpx_max_grid_size** (*integer* *,* *default=32*) – Maximum block size in either direction
  * **warpx_max_grid_size_x** (*integer* *,* *optional*) – Maximum block size in longitudinal direction
  * **warpx_blocking_factor** (*integer* *,* *optional*) – Blocking factor (which controls the block size)
  * **warpx_blocking_factor_x** (*integer* *,* *optional*) – Blocking factor (which controls the block size) in the longitudinal direction
  * **warpx_potential_lo_z** (*float* *,* *default=0.*) – Electrostatic potential on the lower longitudinal boundary
  * **warpx_potential_hi_z** (*float* *,* *default=0.*) – Electrostatic potential on the upper longitudinal boundary
  * **warpx_start_moving_window_step** (*int* *,* *default=0*) – The timestep at which the moving window starts
  * **warpx_end_moving_window_step** (*int* *,* *default=-1*) – The timestep at which the moving window ends. If -1, the moving window
    will continue until the end of the simulation.
  * **warpx_boundary_u_th** (*dict* *,* *default=None*) – If a thermal boundary is used for particles, this dictionary should
    specify the thermal speed for each species in the form {<species>: u_th}.
    Note: u_th = sqrt(T\*q_e/mass)/clight with T in eV.

### *class* pywarpx.picmi.CylindricalGrid(number_of_cells=None, lower_bound=None, upper_bound=None, lower_boundary_conditions=None, upper_boundary_conditions=None, nr=None, nz=None, n_azimuthal_modes=None, rmin=None, rmax=None, zmin=None, zmax=None, bc_rmin=None, bc_rmax=None, bc_zmin=None, bc_zmax=None, moving_window_velocity=None, refined_regions=[], lower_bound_particles=None, upper_bound_particles=None, rmin_particles=None, rmax_particles=None, zmin_particles=None, zmax_particles=None, lower_boundary_conditions_particles=None, upper_boundary_conditions_particles=None, bc_rmin_particles=None, bc_rmax_particles=None, bc_zmin_particles=None, bc_zmax_particles=None, guard_cells=None, pml_cells=None, \*\*kw)

Axisymmetric, cylindrical grid
Parameters can be specified either as vectors or separately.
(If both are specified, the vector is used.)

* **Parameters:**
  * **number_of_cells** (*vector* *of* *integers*) – Number of cells along each axis (number of nodes is number_of_cells+1)
  * **lower_bound** (*vector* *of* *floats*) – Position of the node at the lower bound [m]
  * **upper_bound** (*vector* *of* *floats*) – Position of the node at the upper bound [m]
  * **lower_boundary_conditions** (*vector* *of* *strings*) – Conditions at lower boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann
  * **upper_boundary_conditions** (*vector* *of* *strings*) – Conditions at upper boundaries, periodic, open, dirichlet, absorbing_silver_mueller, or neumann
  * **nr** (*integer*) – Number of cells along R (number of nodes=nr+1)
  * **nz** (*integer*) – Number of cells along Z (number of nodes=nz+1)
  * **n_azimuthal_modes** (*integer*) – Number of azimuthal modes
  * **rmin** (*float*) – Position of first node along R [m]
  * **rmax** (*float*) – Position of last node along R [m]
  * **zmin** (*float*) – Position of first node along Z [m]
  * **zmax** (*float*) – Position of last node along Z [m]
  * **bc_rmin** (*vector* *of* *strings*) – Boundary condition at min R: One of open, dirichlet, absorbing_silver_mueller, or neumann
  * **bc_rmax** (*vector* *of* *strings*) – Boundary condition at max R: One of open, dirichlet, absorbing_silver_mueller, or neumann
  * **bc_zmin** (*vector* *of* *strings*) – Boundary condition at min Z: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann
  * **bc_zmax** (*vector* *of* *strings*) – Boundary condition at max Z: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann
  * **moving_window_velocity** (*vector* *of* *floats* *,* *optional*) – Moving frame velocity [m/s]
  * **refined_regions** (*list* *of* *lists* *,* *optional*) – List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor],
    with level being the refinement level, with 1 being the first level of refinement, 2 being the second etc,
    lo and hi being vectors of length 2 specifying the extent of the region,
    and refinement_factor defaulting to [2,2] (relative to next lower level)
  * **lower_bound_particles** (*vector* *of* *floats* *,* *optional*) – Position of particle lower bound [m]
  * **upper_bound_particles** (*vector* *of* *floats* *,* *optional*) – Position of particle upper bound [m]
  * **rmin_particles** (*float* *,* *optional*) – Position of min particle boundary along R [m]
  * **rmax_particles** (*float* *,* *optional*) – Position of max particle boundary along R [m]
  * **zmin_particles** (*float* *,* *optional*) – Position of min particle boundary along Z [m]
  * **zmax_particles** (*float* *,* *optional*) – Position of max particle boundary along Z [m]
  * **lower_boundary_conditions_particles** (*vector* *of* *strings* *,* *optional*) – Conditions at lower boundaries for particles, periodic, absorbing, reflect or thermal
  * **upper_boundary_conditions_particles** (*vector* *of* *strings* *,* *optional*) – Conditions at upper boundaries for particles, periodic, absorbing, reflect or thermal
  * **bc_rmin_particles** (*string* *,* *optional*) – Boundary condition at min R for particles: One of periodic, absorbing, reflect, thermal
  * **bc_rmax_particles** (*string* *,* *optional*) – Boundary condition at max R for particles: One of periodic, absorbing, reflect, thermal
  * **bc_zmin_particles** (*string* *,* *optional*) – Boundary condition at min Z for particles: One of periodic, absorbing, reflect, thermal
  * **bc_zmax_particles** (*string* *,* *optional*) – Boundary condition at max Z for particles: One of periodic, absorbing, reflect, thermal
  * **guard_cells** (*vector* *of* *integers* *,* *optional*) – Number of guard cells used along each direction
  * **pml_cells** (*vector* *of* *integers* *,* *optional*) – Number of Perfectly Matched Layer (PML) cells along each direction

### References

absorbing_silver_mueller: A local absorbing boundary condition that works best under normal incidence angle.
Based on the Silver-Mueller Radiation Condition, e.g., in

* A. K. Belhora and L. Pichon, “Maybe Efficient Absorbing Boundary Conditions for the Finite Element Solution of 3D Scattering Problems,” 1995,
  [https://doi.org/10.1109/20.376322](https://doi.org/10.1109/20.376322)
* B Engquist and A. Majdat, “Absorbing boundary conditions for numerical simulation of waves,” 1977,
  [https://doi.org/10.1073/pnas.74.5.1765](https://doi.org/10.1073/pnas.74.5.1765)
* R. Lehe, “Electromagnetic wave propagation in Particle-In-Cell codes,” 2016,
  US Particle Accelerator School (USPAS) Summer Session, Self-Consistent Simulations of Beam and Plasma Systems
  [https://people.nscl.msu.edu/~lund/uspas/scs_2016/lec_adv/A1b_EM_Waves.pdf](https://people.nscl.msu.edu/~lund/uspas/scs_2016/lec_adv/A1b_EM_Waves.pdf)
  > Implementation specific documentation

This assumes that WarpX was compiled with USE_RZ = TRUE

See [Input Parameters](https://warpx.readthedocs.io/en/latest/usage/parameters.html) for more information.

* **Parameters:**
  * **warpx_max_grid_size** (*integer* *,* *default=32*) – Maximum block size in either direction
  * **warpx_max_grid_size_x** (*integer* *,* *optional*) – Maximum block size in radial direction
  * **warpx_max_grid_size_y** (*integer* *,* *optional*) – Maximum block size in longitudinal direction
  * **warpx_blocking_factor** (*integer* *,* *optional*) – Blocking factor (which controls the block size)
  * **warpx_blocking_factor_x** (*integer* *,* *optional*) – Blocking factor (which controls the block size) in the radial direction
  * **warpx_blocking_factor_y** (*integer* *,* *optional*) – Blocking factor (which controls the block size) in the longitudinal direction
  * **warpx_potential_lo_r** (*float* *,* *default=0.*) – Electrostatic potential on the lower radial boundary
  * **warpx_potential_hi_r** (*float* *,* *default=0.*) – Electrostatic potential on the upper radial boundary
  * **warpx_potential_lo_z** (*float* *,* *default=0.*) – Electrostatic potential on the lower longitudinal boundary
  * **warpx_potential_hi_z** (*float* *,* *default=0.*) – Electrostatic potential on the upper longitudinal boundary
  * **warpx_reflect_all_velocities** (*bool default=False*) – Whether the sign of all of the particle velocities are changed upon
    reflection on a boundary, or only the velocity normal to the surface
  * **warpx_start_moving_window_step** (*int* *,* *default=0*) – The timestep at which the moving window starts
  * **warpx_end_moving_window_step** (*int* *,* *default=-1*) – The timestep at which the moving window ends. If -1, the moving window
    will continue until the end of the simulation.
  * **warpx_boundary_u_th** (*dict* *,* *default=None*) – If a thermal boundary is used for particles, this dictionary should
    specify the thermal speed for each species in the form {<species>: u_th}.
    Note: u_th = sqrt(T\*q_e/mass)/clight with T in eV.

### *class* pywarpx.picmi.EmbeddedBoundary(implicit_function=None, stl_file=None, stl_scale=None, stl_center=None, stl_reverse_normal=False, potential=None, cover_multiple_cuts=None, \*\*kw)

Custom class to handle set up of embedded boundaries specific to WarpX.
If embedded boundary initialization is added to picmistandard this can be
changed to inherit that functionality. The geometry can be specified either as
an implicit function or as an STL file (ASCII or binary). In the latter case the
geometry specified in the STL file can be scaled, translated and inverted.

* **Parameters:**
  * **implicit_function** (*string*) – Analytic expression describing the embedded boundary
  * **stl_file** (*string*) – STL file path (string), file contains the embedded boundary geometry
  * **stl_scale** (*float*) – Factor by which the STL geometry is scaled
  * **stl_center** (*vector* *of* *floats*) – Vector by which the STL geometry is translated (in meters)
  * **stl_reverse_normal** (*bool*) – If True inverts the orientation of the STL geometry
  * **potential** (*string* *,* *default=0.*) – Analytic expression defining the potential. Can only be specified
    when the solver is electrostatic.
  * **cover_multiple_cuts** (*bool* *,* *default=None*) – Whether to cover cells with multiple cuts.
    (If False, this will raise an error if some cells have multiple cuts)
  * **arguments.** (*Parameters used in the analytic expressions should be given as additional keyword*) – 

Field solvers define the updates of electric and magnetic fields.

### *class* pywarpx.picmi.ElectromagneticSolver(grid, method=None, stencil_order=None, cfl=None, source_smoother=None, field_smoother=None, subcycling=None, galilean_velocity=None, divE_cleaning=None, divB_cleaning=None, pml_divE_cleaning=None, pml_divB_cleaning=None, \*\*kw)

Electromagnetic field solver

* **Parameters:**
  * **grid** (*grid instance*) – Grid object for the diagnostic
  * **method** ( *{'Yee'* *,*  *'CKC'* *,*  *'Lehe'* *,*  *'PSTD'* *,*  *'PSATD'* *,*  *'GPSTD'* *,*  *'DS'* *,*  *'ECT'}*) – 

    The advance method use to solve Maxwell’s equations. The default method is code dependent.
    - ’Yee’: standard solver using the staggered Yee grid ([https://doi.org/10.1109/TAP.1966.1138693](https://doi.org/10.1109/TAP.1966.1138693))
    - ’CKC’: solver with the extended Cole-Karkkainen-Cowan stencil with better dispersion properties
      ([https://doi.org/10.1103/PhysRevSTAB.16.041303](https://doi.org/10.1103/PhysRevSTAB.16.041303))
    - ’Lehe’: CKC-style solver with modified dispersion ([https://doi.org/10.1103/PhysRevSTAB.16.021301](https://doi.org/10.1103/PhysRevSTAB.16.021301))
    - ’PSTD’: Spectral solver with finite difference in time domain, e.g., Q. H. Liu, Letters 15 (3) (1997) 158–165
    - ’PSATD’: Spectral solver with analytic in time domain ([https://doi.org/10.1016/j.jcp.2013.03.010](https://doi.org/10.1016/j.jcp.2013.03.010))
    - ’DS’: Directional Splitting after Yasuhiko Sentoku ([https://doi.org/10.1140/epjd/e2014-50162-y](https://doi.org/10.1140/epjd/e2014-50162-y))
    - ’ECT’: Enlarged Cell Technique solver, allowing internal conductors ([https://doi.org/10.1109/APS.2005.1551259](https://doi.org/10.1109/APS.2005.1551259))
  * **stencil_order** (*vector* *of* *integers*) – Order of stencil for each axis (-1=infinite)
  * **cfl** (*float* *,* *optional*) – Fraction of the Courant-Friedrich-Lewy criteria [1]
  * **source_smoother** (*smoother instance* *,* *optional*) – Smoother object to apply to the sources
  * **field_smoother** (*smoother instance* *,* *optional*) – Smoother object to apply to the fields
  * **subcycling** (*integer* *,* *optional*) – Level of subcycling for the GPSTD solver
  * **galilean_velocity** (*vector* *of* *floats* *,* *optional*) – Velocity of Galilean reference frame [m/s]
  * **divE_cleaning** (*bool* *,* *optional*) – Solver uses div(E) cleaning if True
  * **divB_cleaning** (*bool* *,* *optional*) – Solver uses div(B) cleaning if True
  * **pml_divE_cleaning** (*bool* *,* *optional*) – Solver uses div(E) cleaning in the PML if True
  * **pml_divB_cleaning** – Solver uses div(B) cleaning in the PML if True

See [Input Parameters](https://warpx.readthedocs.io/en/latest/usage/parameters.html) for more information.

* **Parameters:**
  * **warpx_pml_ncell** (*integer* *,* *optional*) – The depth of the PML, in number of cells
  * **warpx_periodic_single_box_fft** (*bool* *,* *default=False*) – Whether to do the spectral solver FFTs assuming a single
    simulation block
  * **warpx_current_correction** (*bool* *,* *default=True*) – Whether to do the current correction for the spectral solver.
    See documentation for exceptions to the default value.
  * **warpx_psatd_update_with_rho** (*bool* *,* *optional*) – Whether to update with the actual rho for the spectral solver
  * **warpx_psatd_do_time_averaging** (*bool* *,* *optional*) – Whether to do the time averaging for the spectral solver
  * **warpx_psatd_JRhom** (*str*) – This determines whether the PSATD JRhom algorithm is used.
    The parameter is a string composed by two characters and one digit.
    The first character represents the time dependency of J within the
    time step over which the electromagnetic fields are evolved, e.g.,
    “C” for constant in time, “L” for linear in time, “Q” for quadratic
    in time.
    The second character represents the time dependency of rho within the
    time step over which the electromagnetic fields are evolved, following
    the same naming convention as for J.
    The last digit is an integer that represents the number of subintervals
    used in the JRhom algorithm.
    Examples: “CL1” (equivalent to the standard PSATD PIC algorithm),
    “CL2”, “LL4”, etc.
    By default, the string is empty and the JRhom algorithm is not used.
  * **warpx_do_pml_in_domain** (*bool* *,* *default=False*) – Whether to do the PML boundaries within the domain (versus
    in the guard cells)
  * **warpx_pml_has_particles** (*bool* *,* *default=False*) – Whether to allow particles in the PML region
  * **warpx_do_pml_j_damping** (*bool* *,* *default=False*) – Whether to do damping of J in the PML

### *class* pywarpx.picmi.ElectrostaticSolver(grid, method=None, required_precision=None, maximum_iterations=None, \*\*kw)

Electrostatic field solver

* **Parameters:**
  * **grid** (*grid instance*) – Grid object for the diagnostic
  * **method** (*string*) – One of ‘FFT’, or ‘Multigrid’
  * **required_precision** (*float* *,* *optional*) – Level of precision required for iterative solvers
  * **maximum_iterations** – Maximum number of iterations for iterative solvers

See [Input Parameters](https://warpx.readthedocs.io/en/latest/usage/parameters.html) for more information.

The standard PICMI parameters required_precision and maximum_iterations control the
MLMG Poisson solver convergence for the labframe electrostatic solvers. When warpx_magnetostatic=True,
these parameters are used as defaults for the magnetostatic solver but can be overridden
with the explicit warpx_magnetostatic_\* parameters.

* **Parameters:**
  * **warpx_relativistic** (*bool* *,* *default=False*) – Whether to use the relativistic solver or lab frame solver
  * **warpx_absolute_tolerance** (*float* *,* *default=0.*) – Absolute tolerance on the labframe electrostatic solver
  * **warpx_self_fields_verbosity** (*integer* *,* *default=2*) – Level of verbosity for the labframe electrostatic solver
  * **warpx_magnetostatic** (*bool* *,* *default=False*) – Whether to also solve for self-consistent magnetic fields from currents.
  * **warpx_magnetostatic_required_precision** (*float* *,* *optional*) – Relative precision for the magnetostatic solver. If not specified,
    defaults to the value of required_precision.
  * **warpx_magnetostatic_absolute_tolerance** (*float* *,* *optional*) – Absolute tolerance for the magnetostatic solver. If not specified,
    defaults to the value of warpx_absolute_tolerance.
  * **warpx_magnetostatic_max_iters** (*integer* *,* *optional*) – Maximum iterations for the magnetostatic solver. If not specified,
    defaults to the value of maximum_iterations.
  * **warpx_magnetostatic_verbosity** (*integer* *,* *optional*) – Verbosity level for the magnetostatic solver. If not specified,
    defaults to the value of warpx_self_fields_verbosity.
  * **warpx_effective_potential** (*bool* *,* *default=False*) – Whether to use the effective potential Poisson solver (EP-PIC)
  * **warpx_effective_potential_factor** (*float* *,* *default=4*) – If the effective potential Poisson solver is used, this sets the value
    of C_EP (the method is marginally stable at C_EP = 1)
  * **warpx_effective_potential_time_filter_param** (*float* *,* *default=0.1*) – Time filtering parameter used to filter sigma in the effective
    potential scheme. sigma is updated using:
    sigma^n = warpx_effective_potential_time_filter_param \* sigma^n + (1 - warpx_effective_potential_time_filter_param) \* sigma^n-1
  * **warpx_effective_potential_density_floor** (*float* *,* *default=0*) – If given, this value will be used as the minimum density during the
    local calculation of sigma.
  * **warpx_dt_update_interval** (*integer* *,* *optional* *(**default = -1* *)*) – How frequently the timestep is updated. Adaptive timestepping is disabled when this is <= 0.
  * **warpx_cfl** (*float* *,* *optional*) – Fraction of the CFL condition for particle velocity vs grid size, used to set the timestep when warpx_dt_update_interval > 0.
  * **warpx_max_dt** (*float* *,* *optional*) – The maximum allowable timestep when warpx_dt_update_interval > 0.

### *class* pywarpx.picmi.HybridPICSolver(grid, Te=None, n0=None, gamma=None, n_floor=None, plasma_resistivity=None, plasma_hyper_resistivity=None, substeps=None, use_rkf45=None, substep_rtol=None, substep_atol=None, substep_safety=None, substep_max_growth=None, max_substep_attempts=None, holmstrom_vacuum_region=None, Jx_external_function=None, Jy_external_function=None, Jz_external_function=None, A_external=None, do_external_diva_cleaning=None, \*\*kw)

Hybrid-PIC solver based on Ohm’s law.
See [Theory Section](https://warpx.readthedocs.io/en/latest/theory/kinetic_fluid_hybrid_model.html) for more information.

* **Parameters:**
  * **Te** (*float*) – Electron temperature in eV.
  * **n0** (*float*) – Reference plasma density in m^-3.
  * **gamma** (*float* *,* *default=5/3*) – Exponent in calculation of electron pressure.
  * **n_floor** (*float* *,* *optional*) – Minimum density used in Ohm’s law calculation.
  * **plasma_resistivity** (*float* *or* *str*) – Value or expression to use for the plasma resistivity in Ohm\*m.
    Can be a constant value or an expression depending on `rho` (charge density),
    `J` (current density magnitude), and `t` (simulation time).
  * **plasma_hyper_resistivity** (*float* *or* *str*) – Value or expression to use for the plasma hyper-resistivity in Ohm\*m^3.
    Can be a constant value or an expression depending on `rho` (charge density)
    and `B` (magnetic field magnitude).
  * **substeps** (*int* *,* *default=10*) – Total number of substeps used to advance the B-field over one full
    timestep (split evenly between the two half-steps, so `substeps/2`
    RK4 steps are taken per half-step, each of duration
    `dt / substeps`). Must be divisible by 2; if not, the value is
    automatically rounded up to the next even number.
    When `use_rkf45=True`, this is instead used only as the initial
    substep count estimate for the adaptive solver.
  * **use_rkf45** (*bool* *,* *default=False*) – If True, use the adaptive Runge-Kutta-Fehlberg 4(5) (RKF45)
    integrator (Fehlberg 1969, NASA Technical Report R-315,
    [https://ntrs.nasa.gov/citations/19690021375](https://ntrs.nasa.gov/citations/19690021375)) for the B-field substep
    advance, with step-size control governed by `substep_rtol` and
    `substep_atol`. If False, use the fixed-step classical RK4
    integrator with `substeps` total substeps per timestep.
  * **substep_rtol** (*float* *,* *default=1e-4*) – Relative tolerance for the RKF45 adaptive step-size control.
    Only used when `use_rkf45=True`.
  * **substep_atol** (*float* *,* *default=1e-8*) – Absolute tolerance for the RKF45 adaptive step-size control.
    Only used when `use_rkf45=True`.
  * **substep_safety** (*float* *,* *default=0.9*) – Safety factor applied to the step-size adjustment formula.
    Only used when `use_rkf45=True`.
  * **substep_max_growth** (*float* *,* *default=5.0*) – Maximum factor by which the substep size may grow after an accepted
    step. Only used when `use_rkf45=True`.
  * **max_substep_attempts** (*int* *,* *default=250*) – Maximum number of substep attempts (accepted + rejected combined) per
    half-step before the simulation aborts. Only used when
    `use_rkf45=True`.
  * **holmstrom_vacuum_region** (*bool* *,* *default=False*) – Flag to determine handling of vacuum region (where rho < n_floor\*q_e). Setting to True will solve the simplified Generalized Ohm’s Law dropping the Hall and pressure terms in the vacuum region. See [Holmstrom (2013)](https://arxiv.org/abs/1301.0272v1).
    This flag is useful for suppressing vacuum region fluctuations. A large resistivity value must be used when rho <= rho_floor.
  * **Jx/y/z_external_function** (*str*) – Function of space and time specifying external (non-plasma) currents.
  * **A_external** (*dict*) – Function of space and time specifying external (non-plasma) vector potential fields.
    It is expected that a nested dictionary will be passed in for each separate vector potential that may have
    different spatial configuration or time dependence. Each field entry should contain either implicit functions
    with (x,y,z) dependence for ‘Ax_external_function’, ‘Ay_external_function’,
    ‘Az_external_function’, plus ‘A_time_external_function’ with (t) dependence, or
    alternatively ‘load_from_file’: True with a ‘path’ to an OpenPMD file along with
    ‘A_time_external_function’.
  * **do_external_diva_cleaning** (*bool* *(**default=True* *)*) – This flag can be used to disable divA cleaning. This may be necessary when using a non-periodic
    external A with periodic field boundary conditions.

### Notes

**Required Parameters:**

- `Te` must be specified when using the hybrid solver.
- `n0` should be specified if `gamma != 1`.

**Best Practices:**

- *Grid type:* Setting `warpx_grid_type='collocated'` is recommended.
- *Particle shape:* Linear particles (`algo.particle_shape = 1`) are recommended.

**Constraints and Limitations:**

- *Mesh refinement:* Only one level is supported (no AMR). The solver will abort if more than one level is used.
- *RZ geometry:* Only the m=0 azimuthal mode is supported in RZ geometry.
- *External vector potential:* If `A_external` is provided, it must be non-empty.
- *Time-dependent A fields:* When using expressions for external vector potentials, time variation must be specified via `A_time_external_function`, not directly in the `A[x,y,z]_external_function` expressions.

For complete parameter documentation, see the [Input Parameters section](https://warpx.readthedocs.io/en/latest/usage/parameters.html#maxwell-solver-kinetic-fluid-hybrid).

Object that allows smoothing of fields.

### *class* pywarpx.picmi.BinomialSmoother(n_pass=None, compensation=None, stride=None, alpha=None, \*\*kw)

Describes a binomial smoother operator (applied to grids)

* **Parameters:**
  * **n_pass** (*vector* *of* *integers*) – Number of passes along each axis
  * **compensation** (*vector* *of* *booleans* *,* *optional*) – Flags whether to apply comensation along each axis
  * **stride** (*vector* *of* *integers* *,* *optional*) – Stride along each axis
  * **alpha** (*vector* *of* *floats* *,* *optional*) – Smoothing coefficients along each axis

## Evolve Schemes

These define the scheme use to evolve the fields and particles.
An instance of one of these would be passed as the evolve_scheme into the Simulation.

### *class* pywarpx.picmi.ExplicitEvolveScheme

Sets up the explicit evolve scheme

### *class* pywarpx.picmi.ThetaImplicitEMEvolveScheme(nonlinear_solver, theta=None)

Sets up the “theta implicit” electromagnetic evolve scheme

* **Parameters:**
  * **nonlinear_solver** (*nonlinear solver instance*) – The nonlinear solver to use for the iterations
  * **theta** (*float* *,* *optional*) – The “theta” parameter, determining the level of implicitness

### *class* pywarpx.picmi.SemiImplicitEMEvolveScheme(nonlinear_solver)

Sets up the “semi-implicit” electromagnetic evolve scheme

* **Parameters:**
  **nonlinear_solver** (*nonlinear solver instance*) – The nonlinear solver to use for the iterations

There are several support classes use to specify components of the evolve schemes

### *class* pywarpx.picmi.PicardNonlinearSolver(verbose=None, require_convergence=None, max_iterations=None, relative_tolerance=None, absolute_tolerance=None, diagnostic_file=None, diagnostic_interval=None)

Sets up the iterative Picard nonlinear solver for the implicit evolve scheme

* **Parameters:**
  * **verbose** (*bool* *,* *default=True*) – Whether there is verbose output from the solver
  * **require_convergence** (*bool* *,* *default True*) – Whether convergence is required. If True and convergence is not obtained, the code will exit.
  * **max_iterations** (*integer* *,* *default=100*) – Maximum number of iterations
  * **relative_tolerance** (*float* *,* *default=1.e-6*) – Relative tolerance of the convergence
  * **absolute_tolerance** (*float* *,* *default=0.*) – Absoluate tolerence of the convergence
  * **diagnostic_file** (*string* *,* *optional*) – File name where solver diagnostics are written
  * **diagnostic_interval** (*string* *,* *optional*) – The intervals for writing out solver diagnostics to the diagnostic file

### *class* pywarpx.picmi.NewtonNonlinearSolver(verbose=None, linear_solver=None, require_convergence=None, max_iterations=None, relative_tolerance=None, absolute_tolerance=None, diagnostic_file=None, diagnostic_interval=None, max_particle_iterations=None, particle_tolerance=None, particle_suborbits=None, print_unconverged_particle_detail=None, use_mass_matrices_jacobian=None, skip_particle_picard_init=None, use_mass_matrices_pc=None, mass_matrices_pc_width=None, pc_type=None)

Sets up the iterative Newton nonlinear solver for the implicit evolve scheme

* **Parameters:**
  * **verbose** (*bool* *,* *default=True*) – Whether there is verbose output from the solver
  * **linear_solver** (*linear solver instance* *,* *optional*) – Specifies input arguments to the linear solver
  * **require_convergence** (*bool* *,* *default True*) – Whether convergence is required. If True and convergence is not obtained, the code will exit.
  * **max_iterations** (*integer* *,* *default=100*) – Maximum number of iterations
  * **relative_tolerance** (*float* *,* *default=1.e-6*) – Relative tolerance of the convergence
  * **absolute_tolerance** (*float* *,* *default=0.*) – Absoluate tolerence of the convergence
  * **diagnostic_file** (*string* *,* *optional*) – File name where solver diagnostics are written
  * **diagnostic_interval** (*string* *,* *optional*) – The intervals for writing out solver diagnostics to the diagnostic file
  * **max_particle_iterations** (*integer* *,* *optional*) – The maximum number of particle iterations
  * **particle_tolerance** (*float* *,* *optional*) – The tolerance of parrticle quantities for convergence
  * **particle_suborbits** (*bool* *,* *optional*) – Whether to use particle suborbits during the solve
  * **print_unconverged_particle_detail** (*bool* *,* *optional*) – Whether to print the details of unconverged particles during suborbits
  * **use_mass_matrices_jacobian** (*bool* *,* *optional*) – Whether to use mass-matrices during the linear stage of PS-JFNK
  * **skip_particle_picard_init** (*bool* *,* *optional*) – When use_mass_matrices_jacobian is True, whether to skip the particle picard iteration on the initial Newton step
  * **use_mass_matrices_pc** (*bool* *,* *optional*) – Whether to capture the plasma response in the preconditioner
  * **mass_matrices_pc_width** (*int* *,* *optional*) – When use_mass_matrices_pc is True, the width of the preconditioner mass matrices
  * **pc_type** (*preconditioner instance* *,* *optional*) – The preconditioner type, An instance of either CurlCurlMLMGPreconditioner, JacobiPreconditioner, or PETScPreconditioner

### *class* pywarpx.picmi.GMRESLinearSolver(verbose_int=None, restart_length=None, absolute_tolerance=None, relative_tolerance=None, max_iterations=None)

Sets up the iterative GMRES linear solver for the implicit Newton nonlinear solver

* **Parameters:**
  * **verbose_int** (*integer* *,* *default=2*) – Level of verbosity of output
  * **restart_length** (*integer* *,* *default=30*) – How often to restart the GMRES iterations
  * **max_iterations** (*integer* *,* *default=1000*) – Maximum number of iterations
  * **relative_tolerance** (*float* *,* *default=1.e-4*) – Relative tolerance of the convergence
  * **absolute_tolerance** (*float* *,* *default=0.*) – Absoluate tolerence of the convergence

### *class* pywarpx.picmi.PETScKSPLinearSolver

Sets up the petsc_ksp linear solver for the implicit Newton nonlinear solver

### *class* pywarpx.picmi.CurlCurlMLMGPreconditioner(verbose, bottom_verbose, agglomeration, consolidation, max_iter, max_coarsening_level, relative_tolerance, absolute_tolerance)

Sets up the curl-curl multigrid preconditioner used during the nonlinear solver

* **Parameters:**
  * **verbose** (*bool* *,* *optional*) – Whether there is verbose output from the solver
  * **bottom_verbose** (*bool* *,* *optional*) – Whether there is verbose output from the bottom solver
  * **agglomeration** (*bool* *,* *optional*) – 
  * **consolidation** (*bool* *,* *optional*) – 
  * **max_iter** (*int* *,* *optional*) – Maximum number of iterations
  * **max_coarsening_level** (*int* *,* *optional*) – Maximum coarsening level
  * **relative_tolerance** (*float* *,* *optional*) – Relative tolerance of the convergence
  * **absolute_tolerance** (*float* *,* *optional*) – Absoluate tolerence of the convergence

### *class* pywarpx.picmi.JacobiPreconditioner(verbose, max_iter, relative_tolerance, absolute_tolerance)

Sets up the point Jacobi preconditioner used during the nonlinear solver

* **Parameters:**
  * **verbose** (*bool* *,* *optional*) – Whether there is verbose output from the solver
  * **max_iter** (*int* *,* *optional*) – Maximum number of iterations
  * **relative_tolerance** (*float* *,* *optional*) – Relative tolerance of the convergence
  * **absolute_tolerance** (*float* *,* *optional*) – Absoluate tolerence of the convergence

### *class* pywarpx.picmi.PETScPreconditioner(type, asm_overlap, sub_type, ilu_factor_levels, hypre_type, euclid_factor_levels)

Sets up the PETSc preconditioner used during the nonlinear solver

* **Parameters:**
  * **type** (*string* *,* *optional*) – PETSc solver type, one of “lu”, “asm”, or “hypre”
  * **asm_overlap** (*int* *,* *optional*) – Parameter for type is “asm”
  * **sub_type** (*string* *,* *optional*) – When type is “asm”, one of “ilu” or “lu”, defailt “ilu”
  * **ilu_factor_levels** (*int* *,* *optional*) – When type is “asm”, and sub_type is “ilu”
  * **hypre_type** (*string* *,* *optional*) – When type is “hypre”, default “euclid”
  * **euclid_factor_levels** (*string* *,* *optional*) – When type is “hypre” and hypre_type is “euclid”

## Constants

For convenience, the PICMI interface defines the following constants,
which can be used directly inside any PICMI script. The values are in SI units.

- `picmi.constants.c`: The speed of light in vacuum.
- `picmi.constants.ep0`: The vacuum permittivity $\epsilon_0$
- `picmi.constants.mu0`: The vacuum permeability $\mu_0$
- `picmi.constants.q_e`: The elementary charge (absolute value of the charge of an electron).
- `picmi.constants.m_e`: The electron mass
- `picmi.constants.m_p`: The proton mass

## Applied fields

Instances of the classes below need to be passed to the method add_applied_field of the Simulation class.

### *class* pywarpx.picmi.AnalyticInitialField(\*\*kw)

Describes an analytic applied field

The expressions should be in terms of the position and time, written as ‘x’, ‘y’, ‘z’, ‘t’.
Parameters can be used in the expression with the values given as additional keyword arguments.
Expressions should be relative to the lab frame.

* **Parameters:**
  * **Ex_expression** (*string* *,* *optional*) – Analytic expression describing Ex field [V/m]
  * **Ey_expression** (*string* *,* *optional*) – Analytic expression describing Ey field [V/m]
  * **Ez_expression** (*string* *,* *optional*) – Analytic expression describing Ez field [V/m]
  * **Bx_expression** (*string* *,* *optional*) – Analytic expression describing Bx field [T]
  * **By_expression** (*string* *,* *optional*) – Analytic expression describing By field [T]
  * **Bz_expression** (*string* *,* *optional*) – Analytic expression describing Bz field [T]
  * **lower_bound** (*vector* *,* *optional*) – Lower bound of the region where the field is applied [m].
  * **upper_bound** – Upper bound of the region where the field is applied [m]

Field Initializer that takes an implicit function to be loaded as an initial E/B field.

* **Parameters:**
  * **warpx_do_initial_div_cleaning** (*bool* *,* *default=True*) – Flag that controls whether or not to execute the Projection based B-field divergence cleaner.
  * **warpx_projection_div_cleaner_atol** (*float*) – Controls the absolute tolerance used in the divergence cleaner solve.
  * **warpx_projection_div_cleaner_rtol** (*float*) – Controls the relative tolerance used in the divergence cleaner solve.

### *class* pywarpx.picmi.ConstantAppliedField(Ex=None, Ey=None, Ez=None, Bx=None, By=None, Bz=None, lower_bound=[None, None, None], upper_bound=[None, None, None], \*\*kw)

Describes a constant applied field

* **Parameters:**
  * **Ex** (*float* *,* *default=0.*) – Constant Ex field [V/m]
  * **Ey** (*float* *,* *default=0.*) – Constant Ey field [V/m]
  * **Ez** (*float* *,* *default=0.*) – Constant Ez field [V/m]
  * **Bx** (*float* *,* *default=0.*) – Constant Bx field [T]
  * **By** (*float* *,* *default=0.*) – Constant By field [T]
  * **Bz** (*float* *,* *default=0.*) – Constant Bz field [T]
  * **lower_bound** (*vector* *,* *optional*) – Lower bound of the region where the field is applied [m].
  * **upper_bound** (*vector* *,* *optional*) – Upper bound of the region where the field is applied [m]

### *class* pywarpx.picmi.AnalyticAppliedField(Ex_expression=None, Ey_expression=None, Ez_expression=None, Bx_expression=None, By_expression=None, Bz_expression=None, lower_bound=[None, None, None], upper_bound=[None, None, None], \*\*kw)

Describes an analytic applied field

The expressions should be in terms of the position and time, written as ‘x’, ‘y’, ‘z’, ‘t’.
Parameters can be used in the expression with the values given as additional keyword arguments.
Expressions should be relative to the lab frame.

* **Parameters:**
  * **Ex_expression** (*string* *,* *optional*) – Analytic expression describing Ex field [V/m]
  * **Ey_expression** (*string* *,* *optional*) – Analytic expression describing Ey field [V/m]
  * **Ez_expression** (*string* *,* *optional*) – Analytic expression describing Ez field [V/m]
  * **Bx_expression** (*string* *,* *optional*) – Analytic expression describing Bx field [T]
  * **By_expression** (*string* *,* *optional*) – Analytic expression describing By field [T]
  * **Bz_expression** (*string* *,* *optional*) – Analytic expression describing Bz field [T]
  * **lower_bound** (*vector* *,* *optional*) – Lower bound of the region where the field is applied [m].
  * **upper_bound** (*vector* *,* *optional*) – Upper bound of the region where the field is applied [m]

### *class* pywarpx.picmi.LoadInitialField(read_fields_from_path, load_B=True, load_E=True, \*\*kw)

The data read in is used to initialize the E and B fields on the grid at the start of the simulation.
The expected format is the file is OpenPMD with axes (x,y,z) in Cartesian, or (r,z) in Cylindrical geometry.

* **Parameters:**
  * **read_fields_from_path** (*string*) – Path to file with field data
  * **load_B** (*bool* *,* *default=True*) – If False, do not load magnetic field
  * **load_E** – If False, do not load electric field

Field Initializer that loads the initial field from a file.

* **Parameters:**
  * **warpx_do_initial_div_cleaning** (*bool* *,* *default=True*) – Flag that controls whether or not to execute the Projection based B-field divergence cleaner.
  * **warpx_projection_div_cleaner_atol** (*float*) – Controls the absolute tolerance used in the divergence cleaner solve.
  * **warpx_projection_div_cleaner_rtol** (*float*) – Controls the relative tolerance used in the divergence cleaner solve.

### *class* pywarpx.picmi.PlasmaLens(period, starts, lengths, strengths_E=None, strengths_B=None, \*\*kw)

Custom class to setup a plasma lens lattice.
The applied fields are dependent only on the transverse position.

* **Parameters:**
  * **period** (*float*) – Periodicity of the lattice (in lab frame, in meters)
  * **starts** (*list* *of* *floats*) – The start of each lens relative to the periodic repeat
  * **lengths** (*list* *of* *floats*) – The length of each lens
  * **strengths_E=None** (*list* *of* *floats* *,* *default = 0.*) – The electric field strength of each lens
  * **strengths_B=None** (*list* *of* *floats* *,* *default = 0.*) – The magnetic field strength of each lens

The field that is applied depends on the transverse position of the particle, (x,y)

- Ex = x\*strengths_E
- Ey = y\*strengths_E
- Bx = +y\*strengths_B
- By = -x\*strengths_B

### *class* pywarpx.picmi.Mirror(x_front_location=None, y_front_location=None, z_front_location=None, depth=None, number_of_cells=None, \*\*kw)

Describes a perfectly reflecting mirror, where the E and B fields are zeroed
out in a plane of finite thickness.

* **Parameters:**
  * **x_front_location** (*float* *,* *optional* *(**see comment below* *)*) – Location in x of the front of the nirror [m]
  * **y_front_location** (*float* *,* *optional* *(**see comment below* *)*) – Location in y of the front of the nirror [m]
  * **z_front_location** (*float* *,* *optional* *(**see comment below* *)*) – Location in z of the front of the nirror [m]
  * **depth** (*float* *,* *optional* *(**see comment below* *)*) – Depth of the mirror [m]
  * **number_of_cells** (*integer* *,* *optional* *(**see comment below* *)*) – Minimum numer of cells zeroed out

Only one of the [x,y,z]_front_location should be specified. The mirror will be set
perpendicular to the respective direction and infinite in the others.
The depth of the mirror will be the maximum of the specified depth and number_of_cells,
or the code’s default value if neither are specified.

## Diagnostics

### *class* pywarpx.picmi.ParticleDiagnostic(period, species=None, data_list=None, write_dir=None, step_min=None, step_max=None, parallelio=None, name=None, \*\*kw)

Defines the particle diagnostics in the simulation frame

* **Parameters:**
  * **period** (*integer*) – Period of time steps that the diagnostic is performed
  * **species** (*species instance* *or* *list* *of* *species instances* *,* *optional*) – Species to write out. If not specified, all species are written.
    Note that the name attribute must be defined for the species.
  * **data_list** (*list* *of* *strings* *,* *optional*) – The data to be written out. Possible values ‘position’, ‘momentum’, ‘weighting’.
    Defaults to the output list of the implementing code.
  * **write_dir** (*string* *,* *optional*) – Directory where data is to be written
  * **step_min** (*integer* *,* *default=0*) – Minimum step at which diagnostics could be written
  * **step_max** (*integer* *,* *default=unbounded*) – Maximum step at which diagnostics could be written
  * **parallelio** (*bool* *,* *optional*) – If set to True, particle diagnostics are dumped in parallel
  * **name** – Sets the base name for the diagnostic output files

See [Input Parameters](https://warpx.readthedocs.io/en/latest/usage/parameters.html) for more information.

* **Parameters:**
  * **warpx_format** ( *{plotfile* *,* *checkpoint* *,* *openpmd* *,* *ascent* *,* *sensei}* *,* *optional*) – Diagnostic file format
  * **warpx_openpmd_backend** ( *{bp* *,* *h5* *,* *json}* *,* *optional*) – Openpmd backend file format
  * **warpx_openpmd_encoding** ( *'v'* *(**variable based* *)* *,*  *'f'* *(**file based* *) or*  *'g'* *(**group based* *)* *,* *optional*) – Only read if `<diag_name>.format = openpmd`. openPMD file output encoding.
    File based: one file per timestep (slower), group/variable based: one file for all steps (faster)).
    Variable based is an experimental feature with ADIOS2. Default: ‘f’.
  * **warpx_file_prefix** (*string* *,* *optional*) – Prefix on the diagnostic file name
  * **warpx_file_min_digits** (*integer* *,* *optional*) – Minimum number of digits for the time step number in the file name
  * **warpx_random_fraction** (*float* *or* *dict* *,* *optional*) – Random fraction of particles to include in the diagnostic. If a float
    is given the same fraction will be used for all species, if a dictionary
    is given the keys should be species with the value specifying the random
    fraction for that species.
  * **warpx_uniform_stride** (*integer* *or* *dict* *,* *optional*) – Stride to down select to the particles to include in the diagnostic.
    If an integer is given the same stride will be used for all species, if
    a dictionary is given the keys should be species with the value
    specifying the stride for that species.
  * **warpx_dump_last_timestep** (*bool* *,* *optional*) – If true, the last timestep is dumped regardless of the diagnostic period/intervals.
  * **warpx_plot_filter_function** (*string* *,* *optional*) – Analytic expression to down select the particles to in the diagnostic
  * **warpx_verbose** (*int* *,* *optional*) – Verbosity level to use for printing diagnostic output information.

### *class* pywarpx.picmi.FieldDiagnostic(grid, period, data_list=None, write_dir=None, step_min=None, step_max=None, number_of_cells=None, lower_bound=None, upper_bound=None, parallelio=None, name=None, \*\*kw)

Defines the electromagnetic field diagnostics in the simulation frame

* **Parameters:**
  * **grid** (*grid instance*) – Grid object for the diagnostic
  * **period** (*integer*) – Period of time steps that the diagnostic is performed
  * **data_list** (*list* *of* *strings* *,* *optional*) – List of quantities to write out. Possible values ‘rho’, ‘E’, ‘B’, ‘J’, ‘Ex’ etc.
    Defaults to the output list of the implementing code.
  * **write_dir** (*string* *,* *optional*) – Directory where data is to be written
  * **step_min** (*integer* *,* *default=0*) – Minimum step at which diagnostics could be written
  * **step_max** (*integer* *,* *default=unbounded*) – Maximum step at which diagnostics could be written
  * **number_of_cells** (*vector* *of* *integers* *,* *optional*) – Number of cells in each dimension.
    If not given, will be obtained from grid.
  * **lower_bound** (*vector* *of* *floats* *,* *optional*) – Lower corner of diagnostics box in each direction.
    If not given, will be obtained from grid.
  * **upper_bound** (*vector* *of* *floats* *,* *optional*) – Higher corner of diagnostics box in each direction.
    If not given, will be obtained from grid.
  * **parallelio** (*bool* *,* *optional*) – If set to True, field diagnostics are dumped in parallel
  * **name** – Sets the base name for the diagnostic output files

See [Input Parameters](https://warpx.readthedocs.io/en/latest/usage/parameters.html) for more information.

* **Parameters:**
  * **warpx_plot_raw_fields** (*bool* *,* *optional*) – Flag whether to dump the raw fields
  * **warpx_plot_raw_fields_guards** (*bool* *,* *optional*) – Flag whether the raw fields should include the guard cells
  * **warpx_format** ( *{plotfile* *,* *checkpoint* *,* *openpmd* *,* *ascent* *,* *sensei}* *,* *optional*) – Diagnostic file format
  * **warpx_openpmd_backend** ( *{bp* *,* *h5* *,* *json}* *,* *optional*) – Openpmd backend file format
  * **warpx_openpmd_encoding** ( *'v'* *(**variable based* *)* *,*  *'f'* *(**file based* *) or*  *'g'* *(**group based* *)* *,* *optional*) – Only read if `<diag_name>.format = openpmd`. openPMD file output encoding.
    File based: one file per timestep (slower), group/variable based: one file for all steps (faster)).
    Variable based is an experimental feature with ADIOS2. Default: ‘f’.
  * **warpx_file_prefix** (*string* *,* *optional*) – Prefix on the diagnostic file name
  * **warpx_file_min_digits** (*integer* *,* *optional*) – Minimum number of digits for the time step number in the file name
  * **warpx_dump_rz_modes** (*bool* *,* *optional*) – Flag whether to dump the data for all RZ modes
  * **warpx_dump_last_timestep** (*bool* *,* *optional*) – If true, the last timestep is dumped regardless of the diagnostic period/intervals.
  * **warpx_particle_fields_to_plot** (*list* *of* *ParticleFieldDiagnostics*) – List of ParticleFieldDiagnostic classes to install in the simulation. Error
    checking is handled in the class itself.
  * **warpx_particle_fields_species** (*list* *of* *strings* *,* *optional*) – Species for which to calculate particle_fields_to_plot functions. Fields will
    be calculated separately for each specified species. If not passed, default is
    all of the available particle species.
  * **warpx_verbose** (*int* *,* *optional*) – Verbosity level to use for printing diagnostic output information.

### *class* pywarpx.picmi.TimeAveragedFieldDiagnostic(grid, period, data_list=None, write_dir=None, step_min=None, step_max=None, number_of_cells=None, lower_bound=None, upper_bound=None, parallelio=None, name=None, \*\*kw)

Defines the electromagnetic field diagnostics in the simulation frame

* **Parameters:**
  * **grid** (*grid instance*) – Grid object for the diagnostic
  * **period** (*integer*) – Period of time steps that the diagnostic is performed
  * **data_list** (*list* *of* *strings* *,* *optional*) – List of quantities to write out. Possible values ‘rho’, ‘E’, ‘B’, ‘J’, ‘Ex’ etc.
    Defaults to the output list of the implementing code.
  * **write_dir** (*string* *,* *optional*) – Directory where data is to be written
  * **step_min** (*integer* *,* *default=0*) – Minimum step at which diagnostics could be written
  * **step_max** (*integer* *,* *default=unbounded*) – Maximum step at which diagnostics could be written
  * **number_of_cells** (*vector* *of* *integers* *,* *optional*) – Number of cells in each dimension.
    If not given, will be obtained from grid.
  * **lower_bound** (*vector* *of* *floats* *,* *optional*) – Lower corner of diagnostics box in each direction.
    If not given, will be obtained from grid.
  * **upper_bound** (*vector* *of* *floats* *,* *optional*) – Higher corner of diagnostics box in each direction.
    If not given, will be obtained from grid.
  * **parallelio** (*bool* *,* *optional*) – If set to True, field diagnostics are dumped in parallel
  * **name** – Sets the base name for the diagnostic output files

See [Input Parameters](https://warpx.readthedocs.io/en/latest/usage/parameters.html) for more information.

* **Parameters:**
  * **warpx_plot_raw_fields** (*bool* *,* *optional*) – Flag whether to dump the raw fields
  * **warpx_plot_raw_fields_guards** (*bool* *,* *optional*) – Flag whether the raw fields should include the guard cells
  * **warpx_format** ( *{plotfile* *,* *checkpoint* *,* *openpmd* *,* *ascent* *,* *sensei}* *,* *optional*) – Diagnostic file format
  * **warpx_openpmd_backend** ( *{bp* *,* *h5* *,* *json}* *,* *optional*) – Openpmd backend file format
  * **warpx_openpmd_encoding** ( *'v'* *(**variable based* *)* *,*  *'f'* *(**file based* *) or*  *'g'* *(**group based* *)* *,* *optional*) – Only read if `<diag_name>.format = openpmd`. openPMD file output encoding.
    File based: one file per timestep (slower), group/variable based: one file for all steps (faster)).
    Variable based is an experimental feature with ADIOS2. Default: ‘f’.
  * **warpx_file_prefix** (*string* *,* *optional*) – Prefix on the diagnostic file name
  * **warpx_file_min_digits** (*integer* *,* *optional*) – Minimum number of digits for the time step number in the file name
  * **warpx_dump_rz_modes** (*bool* *,* *optional*) – Flag whether to dump the data for all RZ modes
  * **warpx_dump_last_timestep** (*bool* *,* *optional*) – If true, the last timestep is dumped regardless of the diagnostic period/intervals.
  * **warpx_particle_fields_to_plot** (*list* *of* *ParticleFieldDiagnostics*) – List of ParticleFieldDiagnostic classes to install in the simulation. Error
    checking is handled in the class itself.
  * **warpx_particle_fields_species** (*list* *of* *strings* *,* *optional*) – Species for which to calculate particle_fields_to_plot functions. Fields will
    be calculated separately for each specified species. If not passed, default is
    all of the available particle species.
  * **warpx_verbose** – Verbosity level to use for printing diagnostic output information.

See [Input Parameters](https://warpx.readthedocs.io/en/latest/usage/parameters.html) for more information.

* **Parameters:**
  * **warpx_time_average_mode** (*str*) – 

    Type of time averaging diagnostic
    Supported values include `"none"`, `"fixed_start"`, and `"dynamic_start"`
    > * `"none"` for no averaging (instantaneous fields)
    > * `"fixed_start"` for a diagnostic that averages all fields between the current output step and a fixed point in time
    > * `"dynamic_start"` for a constant averaging period and output at different points in time (non-overlapping)
  * **warpx_average_period_steps** (*int* *,* *optional*) – Configures the number of time steps in an averaging period.
    Set this only in the `"dynamic_start"` mode and only if `warpx_average_period_time` has not already been set.
    Will be ignored in the `"fixed_start"` mode (with warning).
  * **warpx_average_period_time** (*float* *,* *optional*) – Configures the time (SI units) in an averaging period.
    Set this only in the `"dynamic_start"` mode and only if `average_period_steps` has not already been set.
    Will be ignored in the `"fixed_start"` mode (with warning).
  * **warpx_average_start_steps** (*int* *,* *optional*) – Configures the time step at which time-averaging begins.
    Set this only in the `"fixed_start"` mode.
    Will be ignored in the `"dynamic_start"` mode (with warning).

### pywarpx.picmi.ElectrostaticFieldDiagnostic

alias of [`FieldDiagnostic`](#pywarpx.picmi.FieldDiagnostic)

### *class* pywarpx.picmi.Checkpoint(period=1, write_dir=None, name=None, \*\*kw)

Sets up checkpointing of the simulation, allowing for later restarts

See [Input Parameters](https://warpx.readthedocs.io/en/latest/usage/parameters.html) for more information.

* **Parameters:**
  * **warpx_file_prefix** (*string*) – The prefix to the checkpoint directory names
  * **warpx_file_min_digits** (*integer*) – Minimum number of digits for the time step number in the checkpoint
    directory name.
  * **warpx_verbose** (*int* *,* *optional*) – Verbosity level to use for printing diagnostic output information.

### *class* pywarpx.picmi.ReducedDiagnostic(diag_type, name=None, period=None, path=None, extension=None, separator=None, \*\*kw)

Sets up a reduced diagnostic in the simulation.

See [Input Parameters](https://warpx.readthedocs.io/en/latest/usage/parameters.html#reduced-diagnostics)
for more information.

* **Parameters:**
  * **diag_type** (*string*) – The type of reduced diagnostic. See the link above for all the different
    types of reduced diagnostics available.
  * **name** (*string*) – The name of this diagnostic which will also be the name of the data
    file written to disk.
  * **period** (*integer*) – The simulation step interval at which to output this diagnostic.
  * **path** (*string*) – The file path in which the diagnostic file should be written.
  * **extension** (*string*) – The file extension used for the diagnostic output.
  * **separator** (*string*) – The separator between row values in the output file.
  * **species** (*species instance*) – The name of the species for which to calculate the diagnostic, required for
    diagnostic types ‘BeamRelevant’, ‘ParticleHistogram’, ‘ParticleHistogram2D’, and ‘ParticleExtrema’
  * **bin_number** (*integer*) – For diagnostic type ‘ParticleHistogram’, the number of bins used for the histogram
  * **bin_max** (*float*) – For diagnostic type ‘ParticleHistogram’, the maximum value of the bins
  * **bin_min** (*float*) – For diagnostic type ‘ParticleHistogram’, the minimum value of the bins
  * **normalization** ( *{'unity_particle_weight'* *,*  *'max_to_unity'* *,*  *'area_to_unity'}* *,* *optional*) – For diagnostic type ‘ParticleHistogram’, normalization method of the histogram.
  * **histogram_function** (*string*) – For diagnostic type ‘ParticleHistogram’, the function evaluated to produce the histogram data
  * **filter_function** (*string* *,* *optional*) – For diagnostic types ‘ParticleHistogram’ and ‘ParticleHistogram2D’, the function to filter whether particles are included in the histogram
  * **bin_max_abs** (*float*) – For diagnostic type ‘ParticleHistogram2D’, the maximum value of the bins for the abscissa axis.
  * **bin_max_ord** (*float*) – For diagnostic type ‘ParticleHistogram2D’, the maximum value of the bins for the ordinate axis.
  * **bin_min_abs** (*float*) – For diagnostic type ‘ParticleHistogram2D’, the minimum value of the bins for the abscissa axis.
  * **bin_min_ord** (*float*) – For diagnostic type ‘ParticleHistogram2D’, the minimum value of the bins for the ordinate axis.
  * **bin_number_abs** (*integer*) – For diagnostic type ‘ParticleHistogram2D’, the number of bins used for the histogram for the abscissa axis.
  * **bin_number_ord** (*integer*) – For diagnostic type ‘ParticleHistogram2D’, the number of bins used for the histogram for the ordinate axis.
  * **histogram_function_abs** (*string*) – For diagnostic type ‘ParticleHistogram2D’, the histogram function for the abscissa axis.
  * **histogram_function_ord** (*string*) – For diagnostic type ‘ParticleHistogram2D’, the histogram function for the ordinate axis.
  * **value_function** (*string* *,* *optional*) – For diagnostic type ‘ParticleHistogram2D’, the expression for the weight used to calculate the histogram.
  * **reduced_function** (*string*) – For diagnostic type ‘FieldReduction’, the function of the fields to evaluate
  * **weighting_function** (*string* *,* *optional*) – For diagnostic type ‘ChargeOnEB’, the function to weight contributions to the total charge
  * **reduction_type** ( *{'Maximum'* *,*  *'Minimum'* *, or*  *'Integral'}*) – For diagnostic type ‘FieldReduction’, the type of reduction
  * **probe_geometry** ( *{'Point'* *,*  *'Line'* *,*  *'Plane'}* *,* *default='Point'*) – For diagnostic type ‘FieldProbe’, the geometry of the probe
  * **integrate** (*bool* *,* *default=false*) – For diagnostic type ‘FieldProbe’, whether the field is integrated
  * **do_moving_window_FP** (*bool* *,* *default=False*) – For diagnostic type ‘FieldProbe’, whether the moving window is followed
  * **x_probe** (*floats*) – For diagnostic type ‘FieldProbe’, a probe location. For ‘Point’, the location of the point. For ‘Line’, the start of the
    line. For ‘Plane’, the center of the square detector.
  * **y_probe** (*floats*) – For diagnostic type ‘FieldProbe’, a probe location. For ‘Point’, the location of the point. For ‘Line’, the start of the
    line. For ‘Plane’, the center of the square detector.
  * **z_probe** (*floats*) – For diagnostic type ‘FieldProbe’, a probe location. For ‘Point’, the location of the point. For ‘Line’, the start of the
    line. For ‘Plane’, the center of the square detector.
  * **interp_order** (*integer*) – For diagnostic type ‘FieldProbe’, the interpolation order for ‘Line’ and ‘Plane’
  * **resolution** (*integer*) – For diagnostic type ‘FieldProbe’, the number of points along the ‘Line’ or along each edge of the square ‘Plane’
  * **x1_probe** (*floats*) – For diagnostic type ‘FieldProbe’, the end point for ‘Line’
  * **y1_probe** (*floats*) – For diagnostic type ‘FieldProbe’, the end point for ‘Line’
  * **z1_probe** (*floats*) – For diagnostic type ‘FieldProbe’, the end point for ‘Line’
  * **detector_radius** (*float*) – For diagnostic type ‘FieldProbe’, the detector “radius” (half edge length) of the ‘Plane’
  * **target_normal_x** (*floats*) – For diagnostic type ‘FieldProbe’, the normal vector to the ‘Plane’. Only applicable in 3D
  * **target_normal_y** (*floats*) – For diagnostic type ‘FieldProbe’, the normal vector to the ‘Plane’. Only applicable in 3D
  * **target_normal_z** (*floats*) – For diagnostic type ‘FieldProbe’, the normal vector to the ‘Plane’. Only applicable in 3D
  * **target_up_x** (*floats*) – For diagnostic type ‘FieldProbe’, the vector specifying up in the ‘Plane’
  * **target_up_y** (*floats*) – For diagnostic type ‘FieldProbe’, the vector specifying up in the ‘Plane’
  * **target_up_z** (*floats*) – For diagnostic type ‘FieldProbe’, the vector specifying up in the ‘Plane’

Lab-frame diagnostics diagnostics are used when running boosted-frame simulations.

### *class* pywarpx.picmi.LabFrameParticleDiagnostic(grid, num_snapshots, dt_snapshots, data_list=None, time_start=0.0, species=None, write_dir=None, parallelio=None, name=None, \*\*kw)

Defines the particle diagnostics in the lab frame

* **Parameters:**
  * **grid** (*grid instance*) – Grid object for the diagnostic
  * **num_snapshots** (*integer*) – Number of lab frame snapshots to make
  * **dt_snapshots** (*float*) – Time between each snapshot in lab frame
  * **species** (*species instance* *or* *list* *of* *species instances* *,* *optional*) – Species to write out. If not specified, all species are written.
    Note that the name attribute must be defined for the species.
  * **data_list** (*list* *of* *strings* *,* *optional*) – The data to be written out. Possible values ‘position’, ‘momentum’, ‘weighting’.
    Defaults to the output list of the implementing code.
  * **time_start** (*float* *,* *default=0*) – Time for the first snapshot in lab frame
  * **write_dir** (*string* *,* *optional*) – Directory where data is to be written
  * **parallelio** (*bool* *,* *optional*) – If set to True, particle diagnostics are dumped in parallel
  * **name** – Sets the base name for the diagnostic output files

See [Input Parameters](https://warpx.readthedocs.io/en/latest/usage/parameters.html#backtransformed-diagnostics)
for more information.

* **Parameters:**
  * **warpx_format** (*string* *,* *optional*) – Passed to <diagnostic name>.format
  * **warpx_openpmd_backend** (*string* *,* *optional*) – Passed to <diagnostic name>.openpmd_backend
  * **warpx_openpmd_encoding** ( *'f'* *(**file based* *) or*  *'g'* *(**group based* *)* *,* *optional*) – Only read if `<diag_name>.format = openpmd`. openPMD file output encoding.
    File based: one file per timestep (slower), group/variable based: one file for all steps (faster)).
    Default: ‘f’.
  * **warpx_file_prefix** (*string* *,* *optional*) – Passed to <diagnostic name>.file_prefix
  * **warpx_intervals** (*integer* *or* *string*) – Selects the snapshots to be made, instead of using “num_snapshots” which
    makes all snapshots. “num_snapshots” is ignored.
  * **warpx_file_min_digits** (*integer* *,* *optional*) – Passed to <diagnostic name>.file_min_digits
  * **warpx_buffer_size** (*integer* *,* *optional*) – Passed to <diagnostic name>.buffer_size
  * **warpx_verbose** (*int* *,* *optional*) – Verbosity level to use for printing diagnostic output information.

### *class* pywarpx.picmi.LabFrameFieldDiagnostic(grid, num_snapshots, dt_snapshots, data_list=None, z_subsampling=1, time_start=0.0, write_dir=None, parallelio=None, name=None, \*\*kw)

Defines the electromagnetic field diagnostics in the lab frame

* **Parameters:**
  * **grid** (*grid instance*) – Grid object for the diagnostic
  * **num_snapshots** (*integer*) – Number of lab frame snapshots to make
  * **dt_snapshots** (*float*) – Time between each snapshot in lab frame
  * **data_list** (*list* *of* *strings* *,* *optional*) – List of quantities to write out. Possible values ‘rho’, ‘E’, ‘B’, ‘J’, ‘Ex’ etc.
    Defaults to the output list of the implementing code.
  * **z_subsampling** (*integer* *,* *default=1*) – A factor which is applied on the resolution of the lab frame reconstruction
  * **time_start** (*float* *,* *default=0*) – Time for the first snapshot in lab frame
  * **write_dir** (*string* *,* *optional*) – Directory where data is to be written
  * **parallelio** (*bool* *,* *optional*) – If set to True, field diagnostics are dumped in parallel
  * **name** – Sets the base name for the diagnostic output files

See [Input Parameters](https://warpx.readthedocs.io/en/latest/usage/parameters.html#backtransformed-diagnostics)
for more information.

* **Parameters:**
  * **warpx_format** (*string* *,* *optional*) – Passed to <diagnostic name>.format
  * **warpx_openpmd_backend** (*string* *,* *optional*) – Passed to <diagnostic name>.openpmd_backend
  * **warpx_openpmd_encoding** ( *'f'* *(**file based* *) or*  *'g'* *(**group based* *)* *,* *optional*) – Only read if `<diag_name>.format = openpmd`. openPMD file output encoding.
    File based: one file per timestep (slower), group/variable based: one file for all steps (faster)).
    Default: ‘f’.
  * **warpx_file_prefix** (*string* *,* *optional*) – Passed to <diagnostic name>.file_prefix
  * **warpx_intervals** (*integer* *or* *string*) – Selects the snapshots to be made, instead of using “num_snapshots” which
    makes all snapshots. “num_snapshots” is ignored.
  * **warpx_file_min_digits** (*integer* *,* *optional*) – Passed to <diagnostic name>.file_min_digits
  * **warpx_buffer_size** (*integer* *,* *optional*) – Passed to <diagnostic name>.buffer_size
  * **warpx_lower_bound** (*vector* *of* *floats* *,* *optional*) – Passed to <diagnostic name>.lower_bound
  * **warpx_upper_bound** (*vector* *of* *floats* *,* *optional*) – Passed to <diagnostic name>.upper_bound
  * **warpx_verbose** (*int* *,* *optional*) – Verbosity level to use for printing diagnostic output information.

## Particles

Species objects are a collection of particles with similar properties.
For instance, background plasma electrons, background plasma ions and an externally injected beam could each be their own particle species.

### *class* pywarpx.picmi.Species(particle_type=None, name=None, charge_state=None, charge=None, mass=None, initial_distribution=None, particle_shape=None, density_scale=None, method=None, \*\*kw)

Sets up the species to be simulated.
The species charge and mass can be specified by setting the particle type or by setting them directly.
If the particle type is specified, the charge or mass can be set to override the value from the type.

* **Parameters:**
  * **particle_type** (*string* *,* *optional*) – A string specifying an elementary particle, atom, or other, as defined in
    the openPMD 2 species type extension, openPMD-standard/EXT_SpeciesType.md
  * **name** (*string* *,* *optional*) – Name of the species
  * **method** ( *{'Boris'* *,*  *'Vay'* *,*  *'Higuera-Cary'* *,*  *'Li'* *,*  *'free-streaming'* *,*  *'LLRK4'}*) – 

    The particle advance method to use. Code-specific method can be specified using ‘other:<method>’. The default is code
    dependent.
    - ’Boris’: Standard “leap-frog” Boris advance
    - ’Vay’:
    - ’Higuera-Cary’:
    - ’Li’ :
    - ’free-streaming’: Advance with no fields
    - ’LLRK4’: Landau-Lifschitz radiation reaction formula with RK-4)
  * **charge_state** (*float* *,* *optional*) – Charge state of the species (applies only to atoms) [1]
  * **charge** (*float* *,* *optional*) – Particle charge, required when type is not specified, otherwise determined from type [C]
  * **mass** (*float* *,* *optional*) – Particle mass, required when type is not specified, otherwise determined from type [kg]
  * **initial_distribution** (*distribution instance*) – The initial distribution loaded at t=0. Must be one of the standard distributions implemented.
  * **density_scale** (*float* *,* *optional*) – A scale factor on the density given by the initial_distribution
  * **particle_shape** – Particle shape used for deposition and gather.
    If not specified, the value from the Simulation object will be used.
    Other values maybe specified that are code dependent.

See [Input Parameters](https://warpx.readthedocs.io/en/latest/usage/parameters.html) for more information.

* **Parameters:**
  * **warpx_boost_adjust_transverse_positions** (*bool* *,* *default=False*) – Whether to adjust transverse positions when apply the boost
    to the simulation frame
  * **warpx_self_fields_required_precision** (*float* *,* *default=1.e-11*) – Relative precision on the electrostatic solver
    (when using the relativistic solver)
  * **warpx_self_fields_absolute_tolerance** (*float* *,* *default=0.*) – Absolute precision on the electrostatic solver
    (when using the relativistic solver)
  * **warpx_self_fields_max_iters** (*integer* *,* *default=200*) – Maximum number of iterations for the electrostatic
    solver for the species
  * **warpx_self_fields_verbosity** (*integer* *,* *default=2*) – Level of verbosity for the electrostatic solver
  * **warpx_save_previous_position** (*bool* *,* *default=False*) – Whether to save the old particle positions
  * **warpx_do_not_deposit** (*bool* *,* *default=False*) – Whether or not to deposit the charge and current density for
    for this species
  * **warpx_do_not_push** (*bool* *,* *default=False*) – Whether or not to push this species
  * **warpx_do_not_gather** (*bool* *,* *default=False*) – Whether or not to gather the fields from grids for this species
  * **warpx_radial_numpercell_power** (*float* *,* *default=0.*) – With cylindrical geometry, specifies the radial power of the number of particles per cell
  * **warpx_random_theta** (*bool* *,* *default=True*) – Whether or not to add random angle to the particles in theta
    when in RZ mode.
  * **warpx_reflection_model_xlo** (*string* *,* *default='0.'*) – Expression (in terms of the velocity “v”) specifying the probability
    that the particle will reflect on the lower x boundary
  * **warpx_reflection_model_xhi** (*string* *,* *default='0.'*) – Expression (in terms of the velocity “v”) specifying the probability
    that the particle will reflect on the upper x boundary
  * **warpx_reflection_model_ylo** (*string* *,* *default='0.'*) – Expression (in terms of the velocity “v”) specifying the probability
    that the particle will reflect on the lower y boundary
  * **warpx_reflection_model_yhi** (*string* *,* *default='0.'*) – Expression (in terms of the velocity “v”) specifying the probability
    that the particle will reflect on the upper y boundary
  * **warpx_reflection_model_zlo** (*string* *,* *default='0.'*) – Expression (in terms of the velocity “v”) specifying the probability
    that the particle will reflect on the lower z boundary
  * **warpx_reflection_model_zhi** (*string* *,* *default='0.'*) – Expression (in terms of the velocity “v”) specifying the probability
    that the particle will reflect on the upper z boundary
  * **warpx_save_particles_at_xlo** (*bool* *,* *default=False*) – Whether to save particles lost at the lower x boundary
  * **warpx_save_particles_at_xhi** (*bool* *,* *default=False*) – Whether to save particles lost at the upper x boundary
  * **warpx_save_particles_at_ylo** (*bool* *,* *default=False*) – Whether to save particles lost at the lower y boundary
  * **warpx_save_particles_at_yhi** (*bool* *,* *default=False*) – Whether to save particles lost at the upper y boundary
  * **warpx_save_particles_at_zlo** (*bool* *,* *default=False*) – Whether to save particles lost at the lower z boundary
  * **warpx_save_particles_at_zhi** (*bool* *,* *default=False*) – Whether to save particles lost at the upper z boundary
  * **warpx_save_particles_at_eb** (*bool* *,* *default=False*) – Whether to save particles lost at the embedded boundary
  * **warpx_do_resampling** (*bool* *,* *default=False*) – Whether particles will be resampled
  * **warpx_resampling_min_ppc** (*int* *,* *default=1*) – Cells with fewer particles than this number will be
    skipped during resampling.
  * **warpx_resampling_algorithm_target_weight** (*float*) – Weight that the product particles from resampling will not exceed.
  * **warpx_resampling_trigger_intervals** (*bool* *,* *default=0*) – Timesteps at which to resample
  * **warpx_resampling_trigger_max_avg_ppc** (*int* *,* *default=infinity*) – Resampling will be done when the average number of
    particles per cell exceeds this number
  * **warpx_resampling_algorithm** (*str* *,* *default="leveling_thinning"*) – Resampling algorithm to use.
  * **warpx_resampling_algorithm_velocity_grid_type** (*str* *,* *default="spherical"*) – Type of grid to use when clustering particles in velocity space. Only
    applicable with the velocity_coincidence_thinning algorithm.
  * **warpx_resampling_algorithm_delta_ur** (*float*) – Size of velocity window used for clustering particles during grid-based
    merging, with velocity_grid_type == “spherical”.
  * **warpx_resampling_algorithm_n_theta** (*int*) – Number of bins to use in theta when clustering particle velocities
    during grid-based merging, with velocity_grid_type == “spherical”.
  * **warpx_resampling_algorithm_n_phi** (*int*) – Number of bins to use in phi when clustering particle velocities
    during grid-based merging, with velocity_grid_type == “spherical”.
  * **warpx_resampling_algorithm_delta_u** (*array* *of* *floats* *or* *float*) – Size of velocity window used in ux, uy and uz for clustering particles
    during grid-based merging, with velocity_grid_type == “cartesian”. If
    a single number is given the same du value will be used in all three
    directions.
  * **warpx_add_int_attributes** (*dict*) – Dictionary of extra integer particle attributes initialized from an
    expression that is a function of the variables (x, y, z, ux, uy, uz, t).
  * **warpx_add_real_attributes** (*dict*) – Dictionary of extra real particle attributes initialized from an
    expression that is a function of the variables (x, y, z, ux, uy, uz, t).
  * **warpx_do_temperature_deposition** (*bool* *,* *default=False*) – This flag is set per species to do another pass to deposit temperature
    on each timestep if required. Currently only works with Ohm’s Law Hybrid Solver.

### *class* pywarpx.picmi.MultiSpecies(particle_types=None, names=None, charge_states=None, charges=None, masses=None, proportions=None, initial_distribution=None, particle_shape=None, \*\*kw)

INCOMPLETE: proportions argument is not implemented
Multiple species that are initialized with the same distribution.
Each parameter can be list, giving a value for each species, or a single value which is given to all species.
The species charge and mass can be specified by setting the particle type or by setting them directly.
If the particle type is specified, the charge or mass can be set to override the value from the type.

* **Parameters:**
  * **particle_types** (*list* *of* *strings* *,* *optional*) – A string specifying an elementary particle, atom, or other, as defined in
    the openPMD 2 species type extension, openPMD-standard/EXT_SpeciesType.md
  * **names** (*list* *of* *strings* *,* *optional*) – Names of the species
  * **charge_states** (*list* *of* *floats* *,* *optional*) – Charge states of the species (applies only to atoms)
  * **charges** (*list* *of* *floats* *,* *optional*) – Particle charges, required when type is not specified, otherwise determined from type [C]
  * **masses** (*list* *of* *floats* *,* *optional*) – Particle masses, required when type is not specified, otherwise determined from type [kg]
  * **proportions** (*list* *of* *floats* *,* *optional*) – Proportions of the initial distribution made up by each species
  * **initial_distribution** (*distribution instance*) – Initial particle distribution, applied to all species
  * **particle_shape** ( *{'NGP'* *,*  *'linear'* *,*  *'quadratic'* *,*  *'cubic'}*) – Particle shape used for deposition and gather.
    If not specified, the value from the Simulation object will be used.
    Other values maybe specified that are code dependent.

Particle distributions can be used for to initialize particles in a particle species.

### *class* pywarpx.picmi.GaussianBunchDistribution(n_physical_particles, rms_bunch_size, rms_velocity=[0.0, 0.0, 0.0], centroid_position=[0.0, 0.0, 0.0], centroid_velocity=[0.0, 0.0, 0.0], velocity_divergence=[0.0, 0.0, 0.0], \*\*kw)

Describes a Gaussian distribution of particles

* **Parameters:**
  * **n_physical_particles** (*integer*) – Number of physical particles in the bunch
  * **rms_bunch_size** (*vector* *of* *length 3* *of* *floats*) – RMS bunch size at t=0 [m]
  * **rms_velocity** (*vector* *of* *length 3* *of* *floats* *,* *default=* *[**0.* *,**0.* *,**0.* *]*) – RMS velocity spread at t=0 [m/s]
  * **centroid_position** (*vector* *of* *length 3* *of* *floats* *,* *default=* *[**0.* *,**0.* *,**0.* *]*) – Position of the bunch centroid at t=0 [m]
  * **centroid_velocity** (*vector* *of* *length 3* *of* *floats* *,* *default=* *[**0.* *,**0.* *,**0.* *]*) – Velocity (gamma\*V) of the bunch centroid at t=0 [m/s]
  * **velocity_divergence** (*vector* *of* *length 3* *of* *floats* *,* *default=* *[**0.* *,**0.* *,**0.* *]*) – Expansion rate of the bunch at t=0 [m/s/m]

### *class* pywarpx.picmi.UniformDistribution(density, lower_bound=[None, None, None], upper_bound=[None, None, None], rms_velocity=[0.0, 0.0, 0.0], directed_velocity=[0.0, 0.0, 0.0], fill_in=None, \*\*kw)

Describes a uniform density distribution of particles

* **Parameters:**
  * **density** (*float*) – Physical number density [m^-3]
  * **lower_bound** (*vector* *of* *length 3* *of* *floats* *,* *optional*) – Lower bound of the distribution [m]
  * **upper_bound** (*vector* *of* *length 3* *of* *floats* *,* *optional*) – Upper bound of the distribution [m]
  * **rms_velocity** (*vector* *of* *length 3* *of* *floats* *,* *default=* *[**0.* *,**0.* *,**0.* *]*) – Thermal velocity spread [m/s]
  * **directed_velocity** (*vector* *of* *length 3* *of* *floats* *,* *default=* *[**0.* *,**0.* *,**0.* *]*) – Directed, average, proper velocity [m/s]
  * **fill_in** (*bool* *,* *optional*) – Flags whether to fill in the empty spaced opened up when the grid moves

### *class* pywarpx.picmi.AnalyticDistribution(density_expression, momentum_expressions=[None, None, None], momentum_spread_expressions=[None, None, None], lower_bound=[None, None, None], upper_bound=[None, None, None], rms_velocity=[0.0, 0.0, 0.0], directed_velocity=[0.0, 0.0, 0.0], fill_in=None, \*\*kw)

Describes a plasma with density following a provided analytic expression

* **Parameters:**
  * **density_expression** (*string*) – Analytic expression describing physical number density (string) [m^-3].
    Expression should be in terms of the position, written as ‘x’, ‘y’, and ‘z’.
    Parameters can be used in the expression with the values given as keyword arguments.
  * **momentum_expressions** (*list* *of* *strings*) – Analytic expressions describing the gamma\*velocity for each axis [m/s].
    Expressions should be in terms of the position, written as ‘x’, ‘y’, and ‘z’.
    Parameters can be used in the expression with the values given as keyword arguments.
    For any axis not supplied (set to None), directed_velocity will be used.
  * **momentum_spread_expressions** (*list* *of* *strings*) – Analytic expressions describing the gamma\*velocity Gaussian thermal spread sigma for each axis [m/s].
    Expressions should be in terms of the position, written as ‘x’, ‘y’, and ‘z’.
    Parameters can be used in the expression with the values given as keyword arguments.
    For any axis not supplied (set to None), zero will be used.
  * **lower_bound** (*vector* *of* *length 3* *of* *floats* *,* *optional*) – Lower bound of the distribution [m]
  * **upper_bound** (*vector* *of* *length 3* *of* *floats* *,* *optional*) – Upper bound of the distribution [m]
  * **rms_velocity** (*vector* *of* *length 3* *of* *floats* *,* *detault=* *[**0.* *,**0.* *,**0.* *]*) – Thermal velocity spread [m/s]
  * **directed_velocity** (*vector* *of* *length 3* *of* *floats* *,* *detault=* *[**0.* *,**0.* *,**0.* *]*) – Directed, average, proper velocity [m/s]
  * **fill_in** (*bool* *,* *optional*) – Flags whether to fill in the empty spaced opened up when the grid moves

This example will create a distribution where the density is n0 below rmax and zero elsewhere.:

```cpp
.. code-block:: python
```

> dist = AnalyticDistribution(density_expression=’((x\*\*2+y\*\*2)<rmax\*\*2)\*n0’,
> : > rmax = 1.,
>   > n0 = 1.e20,
>   > …)
>   <br/>
>   Implementation specific documentation
* **Parameters:**
  * **warpx_density_min** (*float*) – Minimum plasma density. No particle is injected where the density is
    below this value.
  * **warpx_density_max** (*float*) – Maximum plasma density. The density at each point is the minimum between
    the value given in the profile, and density_max.
  * **warpx_momentum_spread_expressions** (*list* *of* *string*) – Analytic expressions describing the gamma\*velocity spread for each axis [m/s].
    Expressions should be in terms of the position, written as ‘x’, ‘y’, and ‘z’.
    Parameters can be used in the expression with the values given as keyword arguments.
    For any axis not supplied (set to None), zero will be used.

### *class* pywarpx.picmi.UniformFluxDistribution(flux, flux_normal_axis, surface_flux_position, flux_direction, lower_bound=[None, None, None], upper_bound=[None, None, None], rms_velocity=[0.0, 0.0, 0.0], directed_velocity=[0.0, 0.0, 0.0], flux_tmin=None, flux_tmax=None, gaussian_flux_momentum_distribution=None, \*\*kw)

Describes a flux of particles emitted from a plane

* **Parameters:**
  * **flux** (*string*) – Analytic expression describing flux of particles [m^-2.s^-1]
    Expression should be in terms of the position and time, written as ‘x’, ‘y’, ‘z’, and ‘t’.
  * **flux_normal_axis** (*string*) – x, y, or z for 3D, x or z for 2D, or r, t, or z in RZ geometry
  * **surface_flux_position** (*double*) – location of the injection plane [m] along the direction
    specified by flux_normal_axis
  * **flux_direction** (*int*) – Direction of the flux relative to the plane: -1 or +1
  * **lower_bound** (*vector* *of* *floats* *,* *optional*) – Lower bound of the distribution [m]
  * **upper_bound** (*vector* *of* *floats* *,* *optional*) – Upper bound of the distribution [m]
  * **rms_velocity** (*vector* *of* *floats* *,* *default=* *[**0.* *,**0.* *,**0.* *]*) – Thermal velocity spread [m/s]
  * **directed_velocity** (*vector* *of* *floats* *,* *default=* *[**0.* *,**0.* *,**0.* *]*) – Directed, average, proper velocity [m/s]
  * **flux_tmin** (*float* *,* *optional*) – Time at which the flux injection will be turned on.
  * **flux_tmax** (*float* *,* *optional*) – Time at which the flux injection will be turned off.
  * **gaussian_flux_momentum_distribution** – If True, the momentum distribution is v\*Gaussian,
    in the direction normal to the plane. Otherwise,
    the momentum distribution is simply Gaussian.
  * **warpx_inject_from_embedded_boundary** (*bool*) – When true, the flux is injected from the embedded boundaries instead
    of a plane.

### *class* pywarpx.picmi.AnalyticFluxDistribution(flux, flux_normal_axis, surface_flux_position, flux_direction, lower_bound=[None, None, None], upper_bound=[None, None, None], rms_velocity=[0.0, 0.0, 0.0], directed_velocity=[0.0, 0.0, 0.0], flux_tmin=None, flux_tmax=None, gaussian_flux_momentum_distribution=None, \*\*kw)

Describes a flux of particles emitted from a plane

* **Parameters:**
  * **flux** (*string*) – Analytic expression describing flux of particles [m^-2.s^-1]
    Expression should be in terms of the position and time, written as ‘x’, ‘y’, ‘z’, and ‘t’.
  * **flux_normal_axis** (*string*) – x, y, or z for 3D, x or z for 2D, or r, t, or z in RZ geometry
  * **surface_flux_position** (*double*) – location of the injection plane [m] along the direction
    specified by flux_normal_axis
  * **flux_direction** (*int*) – Direction of the flux relative to the plane: -1 or +1
  * **lower_bound** (*vector* *of* *floats* *,* *optional*) – Lower bound of the distribution [m]
  * **upper_bound** (*vector* *of* *floats* *,* *optional*) – Upper bound of the distribution [m]
  * **rms_velocity** (*vector* *of* *floats* *,* *default=* *[**0.* *,**0.* *,**0.* *]*) – Thermal velocity spread [m/s]
  * **directed_velocity** (*vector* *of* *floats* *,* *default=* *[**0.* *,**0.* *,**0.* *]*) – Directed, average, proper velocity [m/s]
  * **flux_tmin** (*float* *,* *optional*) – Time at which the flux injection will be turned on.
  * **flux_tmax** (*float* *,* *optional*) – Time at which the flux injection will be turned off.
  * **gaussian_flux_momentum_distribution** – If True, the momentum distribution is v\*Gaussian,
    in the direction normal to the plane. Otherwise,
    the momentum distribution is simply Gaussian.
  * **warpx_inject_from_embedded_boundary** (*bool*) – When true, the flux is injected from the embedded boundaries instead
    of a plane.

### *class* pywarpx.picmi.ParticleListDistribution(x=0.0, y=0.0, z=0.0, ux=0.0, uy=0.0, uz=0.0, weight=0.0, \*\*kw)

Load particles at the specified positions and velocities

* **Parameters:**
  * **x** (*float* *,* *default=0.*) – List of x positions of the particles [m]
  * **y** (*float* *,* *default=0.*) – List of y positions of the particles [m]
  * **z** (*float* *,* *default=0.*) – List of z positions of the particles [m]
  * **ux** (*float* *,* *default=0.*) – List of ux positions of the particles (ux = gamma\*vx) [m/s]
  * **uy** (*float* *,* *default=0.*) – List of uy positions of the particles (uy = gamma\*vy) [m/s]
  * **uz** (*float* *,* *default=0.*) – List of uz positions of the particles (uz = gamma\*vz) [m/s]
  * **weight** (*float*) – Particle weight or list of weights, number of real particles per simulation particle

### *class* pywarpx.picmi.FromFileDistribution(file_path, \*\*kw)

Load particles from an openPMD file.

The openPMD file must contain the attributes position, momentum, weighting.

Particle layouts determine how to microscopically place macro particles in a grid cell.

### *class* pywarpx.picmi.GriddedLayout(n_macroparticle_per_cell, grid=None, \*\*kw)

Specifies a gridded layout of particles

* **Parameters:**
  * **n_macroparticle_per_cell** (*vector* *of* *integers*) – Number of particles per cell along each axis
  * **grid** (*grid instance* *,* *optional*) – Grid object specifying the grid to follow.
    If not specified, the underlying grid of the code is used.

### *class* pywarpx.picmi.PseudoRandomLayout(n_macroparticles=None, n_macroparticles_per_cell=None, seed=None, grid=None, \*\*kw)

Specifies a pseudo-random layout of the particles

* **Parameters:**
  * **n_macroparticles** (*integer*) – Total number of macroparticles to load.
    Either this argument or n_macroparticles_per_cell should be supplied.
  * **n_macroparticles_per_cell** (*integer*) – Number of macroparticles to load per cell.
    Either this argument or n_macroparticles should be supplied.
  * **seed** (*integer* *,* *optional*) – Pseudo-random number generator seed
  * **grid** (*grid instance* *,* *optional*) – Grid object specifying the grid to follow for n_macroparticles_per_cell.
    If not specified, the underlying grid of the code is used.

Other operations related to particles:

### *class* pywarpx.picmi.CoulombCollisions(name, species, CoulombLog=None, ndt_supercycle=None, ndt_subcycle=None, \*\*kw)

Custom class to handle setup of binary Coulomb collisions in WarpX. If
collision initialization is added to picmistandard this can be changed to
inherit that functionality.

* **Parameters:**
  * **name** (*string*) – Name of instance (used in the inputs file)
  * **species** (*list* *of* *species instances*) – The species involved in the collision. Must be of length 2.
  * **CoulombLog** (*float* *,* *optional*) – Value of the Coulomb log to use in the collision cross section.
    If not supplied, it is calculated from the local conditions.
  * **ndt_supercycle** (*integer* *,* *optional*) – Run collision once every ndt_supercycle PIC time steps
    (dt_collision = ndt_supercycle \* dt_PIC). Must be >= 1.
    Mutually exclusive with ndt_subcycle. Default is 1.
  * **ndt_subcycle** (*integer* *,* *optional*) – Run collision ndt_subcycle times per PIC time step
    (dt_collision = dt_PIC / ndt_subcycle). Must be >= 1.
    Mutually exclusive with ndt_supercycle.

### *class* pywarpx.picmi.DSMCCollisions(name, species, scattering_processes, product_species=None, ndt_supercycle=None, ndt_subcycle=None, \*\*kw)

Custom class to handle setup of DSMC collisions in WarpX. If collision
initialization is added to picmistandard this can be changed to inherit
that functionality.

* **Parameters:**
  * **name** (*string*) – Name of instance (used in the inputs file)
  * **species** (*species instance*) – The species involved in the collision
  * **scattering_processes** (*dictionary*) – The scattering process to use and any needed information
  * **product_species** (*list*) – The species produced by collision processes (currently both
    ionization and charge-exchange require defining the product species).
  * **ndt_supercycle** (*integer* *,* *optional*) – Run collision once every ndt_supercycle PIC time steps
    (dt_collision = ndt_supercycle \* dt_PIC). Must be >= 1.
    Mutually exclusive with ndt_subcycle. Default is 1.
  * **ndt_subcycle** (*integer* *,* *optional*) – Run collision ndt_subcycle times per PIC time step
    (dt_collision = dt_PIC / ndt_subcycle). Must be >= 1.
    Mutually exclusive with ndt_supercycle.

### *class* pywarpx.picmi.MCCCollisions(name, species, background_density, background_temperature, scattering_processes, background_mass=None, max_background_density=None, ndt_supercycle=None, ndt_subcycle=None, \*\*kw)

Custom class to handle setup of MCC collisions in WarpX. If collision
initialization is added to picmistandard this can be changed to inherit
that functionality.

* **Parameters:**
  * **name** (*string*) – Name of instance (used in the inputs file)
  * **species** (*species instance*) – The species involved in the collision
  * **background_density** (*float* *or* *string*) – The density of the background. An string expression as a function of (x, y, z, t) can be used.
  * **background_temperature** (*float* *or* *string*) – The temperature of the background. An string expression as a function of (x, y, z, t) can be used.
  * **scattering_processes** (*dictionary*) – The scattering process to use and any needed information
  * **background_mass** (*float* *,* *optional*) – The mass of the background particle. If not supplied, the default depends
    on the type of scattering process.
  * **max_background_density** (*float*) – The maximum background density. When the background_density is an expression, this must also
    be specified.
  * **ndt_supercycle** (*integer* *,* *optional*) – Run collision once every ndt_supercycle PIC time steps
    (dt_collision = ndt_supercycle \* dt_PIC). Must be >= 1.
    Mutually exclusive with ndt_subcycle. Default is 1.
  * **ndt_subcycle** (*integer* *,* *optional*) – Run collision ndt_subcycle times per PIC time step
    (dt_collision = dt_PIC / ndt_subcycle). Must be >= 1.
    Mutually exclusive with ndt_supercycle.

### *class* pywarpx.picmi.FieldIonization(model, ionized_species, product_species, \*\*kw)

Field ionization on an ion species

* **Parameters:**
  * **model** (*string*) – Ionization model, e.g. “ADK”
  * **ionized_species** (*species instance*) – Species that is ionized
  * **product_species** – Species in which ionized electrons are stored.

WarpX only has ADK ionization model implemented.

## Laser Pulses

Laser profiles can be used to initialize laser pulses in the simulation.

### *class* pywarpx.picmi.GaussianLaser(wavelength, waist, duration, propagation_direction, polarization_direction, focal_position, centroid_position, a0=None, E0=None, phi0=None, zeta=None, beta=None, phi2=None, name=None, fill_in=True, \*\*kw)

Specifies a Gaussian laser distribution.

More precisely, the electric field **near the focal plane** is given by:

$$
E(\boldsymbol{x},t) = a_0\times E_0\,
\exp\left( -\frac{r^2}{w_0^2} - \frac{(z-z_0-ct)^2}{c^2\tau^2} \right)
\cos[ k_0( z - z_0 - ct ) - \phi_{cep} ]
$$

where $k_0 = 2\pi/\lambda_0$ is the wavevector and where
$E_0 = m_e c^2 k_0 / q_e$ is the field amplitude for $a_0=1$.

#### NOTE
The additional terms that arise **far from the focal plane**
(Gouy phase, wavefront curvature, …) are not included in the above
formula for simplicity, but are of course taken into account by
the code, when initializing the laser pulse away from the focal plane.

* **Parameters:**
  * **wavelength** (*float*) – Laser wavelength [m], defined as $\lambda_0$ in the above formula
  * **waist** (*float*) – Waist of the Gaussian pulse at focus [m], defined as $w_0$ in the above formula
  * **duration** (*float*) – Duration of the Gaussian pulse [s], defined as $\tau$ in the above formula
  * **propagation_direction** (*unit vector* *of* *length 3* *of* *floats*) – Direction of propagation [1]
  * **polarization_direction** (*unit vector* *of* *length 3* *of* *floats*) – Direction of polarization [1]
  * **focal_position** (*vector* *of* *length 3* *of* *floats*) – Position of the laser focus [m]
  * **centroid_position** (*vector* *of* *length 3* *of* *floats*) – Position of the laser centroid at time 0 [m]
  * **a0** (*float*) – Normalized vector potential at focus
    Specify either a0 or E0 (E0 takes precedence).
  * **E0** (*float*) – Maximum amplitude of the laser field [V/m]
    Specify either a0 or E0 (E0 takes precedence).
  * **phi0** (*float*) – Carrier envelope phase (CEP) [rad]
  * **zeta** (*float*) – Spatial chirp at focus (in the lab frame) [m.s]
  * **beta** (*float*) – Angular dispersion at focus (in the lab frame) [rad.s]
  * **phi2** (*float*) – Temporal chirp at focus (in the lab frame) [s^2]
  * **fill_in** (*bool* *,* *default=True*) – Flags whether to fill in the empty spaced opened up when the grid moves
  * **name** (*string* *,* *optional*) – Optional name of the laser

### *class* pywarpx.picmi.AnalyticLaser(field_expression, wavelength, propagation_direction, polarization_direction, amax=None, Emax=None, name=None, fill_in=True, \*\*kw)

Specifies a laser with an analytically described distribution

* **Parameters:**
  * **name=None** (*string* *,* *optional*) – Optional name of the laser
  * **field_expression** (*string*) – Analytic expression describing the electric field of the laser [V/m]
    Expression should be in terms of the position, ‘X’, ‘Y’, in the plane orthogonal
    to the propagation direction, and ‘t’ the time. The expression should describe
    the full field, including the oscillitory component.
    Parameters can be used in the expression with the values given as keyword arguments.
  * **wavelength** (*float*) – Laser wavelength.
    This should be built into the expression, but some codes require a specified value for numerical purposes.
  * **propagation_direction** (*unit vector* *of* *length 3* *of* *floats*) – Direction of propagation [1]
  * **polarization_direction** (*unit vector* *of* *length 3* *of* *floats*) – Direction of polarization [1]
  * **amax** (*float* *,* *optional*) – Maximum normalized vector potential.
    Specify either amax or Emax (Emax takes precedence).
    This should be built into the expression, but some codes require a specified value for numerical purposes.
  * **Emax** (*float* *,* *optional*) – Maximum amplitude of the laser field [V/m].
    Specify either amax or Emax (Emax takes precedence).
    This should be built into the expression, but some codes require a specified value for numerical purposes.
  * **fill_in** (*bool* *,* *default=True*) – Flags whether to fill in the empty spaced opened up when the grid moves

Laser injectors control where to initialize laser pulses on the simulation grid.

### *class* pywarpx.picmi.LaserAntenna(position, normal_vector=None, \*\*kw)

Specifies the laser antenna injection method

* **Parameters:**
  * **position** (*vector* *of* *strings*) – Position of antenna launching the laser [m]
  * **normal_vector** (*vector* *of* *strings* *,* *optional*) – Vector normal to antenna plane, defaults to the laser direction
    of propagation [1]
