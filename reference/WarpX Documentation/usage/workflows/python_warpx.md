# Accessing global WarpX functionalities (e.g., extract timestep)

An important object is `sim.extension.warpx`, which is the Python equivalent to the
C++ `WarpX` simulation class and gives access to global functionalities:

### *class* WarpX

#### getistep(lev: int)

Get the current step on mesh-refinement level `lev`.

#### gett_new(lev: int)

Get the current physical time on mesh-refinement level `lev`.

#### getdt(lev: int)

Get the current physical time step size on mesh-refinement level `lev`.

#### multi_particle_container()

#### get_particle_boundary_buffer()

#### set_potential_on_domain_boundary(potential_[lo/hi]_[x/y/z]: str)

The potential on the domain boundaries can be modified when using the electrostatic solver.
This function updates the strings and function parsers which set the domain
boundary potentials during the Poisson solve.

#### set_potential_on_eb(potential: str)

The embedded boundary (EB) conditions can be modified when using the electrostatic solver.
This set the EB potential string and updates the function parser.

#### evolve(numsteps=-1)

Evolve the simulation the specified number of steps.

#### step(numsteps=-1)

An alias to the evolve method.

#### finalize(finalize_mpi=1)

Call finalize for WarpX and AMReX. Registered to run at program exit.
