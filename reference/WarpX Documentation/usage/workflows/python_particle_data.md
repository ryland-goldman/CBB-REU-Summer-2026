# Accessing particles data

## Selecting a given particle species

The simulation’s particles are accessed through the `sim.particles` attribute, where `sim` is obtained as shown
in the [How to run a simulation with Python extensions](python_extend.md#usage-python-extend-run-simulation) section. Specific particle species (e.g. electrons, ions, photons, etc.)
are selected with the `sim.particles.get` method, as shown in the example below.

```python
# Preparation: set up the sim object
#   sim = picmi.Simulation(...)
#   ...

# Extract the electrons particle species
electrons = sim.particles.get("electrons")
```

The function `sim.particles.get` returns an object of type `WarpXParticleContainer`, from which the
data of individual particles can be accessed or modified as described below.

## Accessing/modifying the underlying particle data

There are several ways to access and modify the particle data, i.e. the positions (`'x'`, `'y'`, `'z'` in 3D Cartesian geometry),
normalized momenta (`'ux'`, `'uy'`, `'uz'`), the particle weights (`'w'`), and the unique IDs (`'idcpu'`) of each individual particles.

#### NOTE
In geometries other than 3D Cartesian, particle positions are defined by different variables.
For example, in RZ geometry, positions are accessed as `'r'` and `'z'` (with an additional `'theta'` attribute),
while in RCYLINDER geometry, only `'r'` is available (with `'theta'`), and in RSPHERE geometry, `'r'` is available
(with `'theta'` and `'phi'`). See [Dimensionality](../../developers/dimensionality.md#developers-dimensionality) for a complete table of position attributes
available in each geometry.

The different methods below differ in their user-friendliness, flexibility and performance overhead.
(For more information, see the [pyamrex documentation](http://pyamrex.readthedocs.io/en/latest/usage/compute.html#particles).)

### Global access through pandas DataFrame (read-only)

The method `to_df` of the `WarpXParticleContainer` object returns a
[pandas DataFrame](https://pandas.pydata.org/docs/user_guide/dsintro.html#dataframe) containing the particle data.
More specifically, the keys of the DataFrame are the particle attributes (e.g., `'ux'`, `'w'`, `'idcpu'`),
and the corresponding arrays have one element per particle, and contain the particles of that species across all
boxes and tiles (on the current MPI rank) and across all mesh refinement levels.

#### WARNING
The data in the DataFrame is a copy of the particle data, and therefore modifying it will not modify
the actual particle data in the simulation.

#### NOTE
The method `to_df` is very convenient because it automatically concatenates all particles across boxes and tiles,
and across all mesh refinement levels. However, this implies significant performance overheads, as it incurs copies
and CPU-GPU data transfers. This method is thus mostly meant for debugging and visualization purposes,
and not for performance-critical operations.

```python
# Preparation: set up the simulation
#   sim = picmi.Simulation(...)
#   ...

# Extract the electrons particle species
electrons = sim.particles.get("electrons")

# local particles (returns only particles on the current MPI rank)
df = electrons.to_df(local=True)  # this is a copy!
print('Available attributes: ', df.columns)
print('Number of particles: ', len(df))

# print position x (one element per particle)
print('Position x: ', df['x'])

# Warning: because `df` is a copy, modifying it will
# not modify the actual particle data
df['x'] += 0.1 # This does not modify the actual particle data
```

### Explicit loop over boxes/tiles

This method provides similar capabilities to the pandas DataFrame approach, but adapts to
the data structure of the particles in WarpX (i.e. particles are organized per box/tile and
mesh refinement level). Unlike the pandas DataFrame approach, this avoids unneeded copies
and CPU-GPU data transfers. As a result, this method offers significantly higher performance,
especially for large-scale parallel simulations and GPU-accelerated runs. The data is accessed
by explicitly looping over mesh-refinement levels and individual grid blocks (boxes), giving
direct access to the underlying particle data arrays for each local block.

```python
# Preparation: set up the simulation
#   sim = picmi.Simulation(...)
#   xp, _ = load_cupy()
#   ...

# Extract the electrons particle species
electrons = sim.particles.get("electrons")

# iterate over boxes/tiles on level 0
for pti in electrons.iterator(level=0):

    # print position x (one element per particle)
    print('Position x: ', pti['x'])

    # increment position by a random value
    # using numpy/cupy syntax
    pti["x"][:] += xp.random.random( len(pti['x'][:]) )
```

In the above example, `xp` represents either the `numpy` or `cupy` package.
See [Writing portable Python code that can be executed on CPU and GPU](python_portable.md#usage-python-portable) for more details on writing portable Python code.

## Adding new particles

New particles can be added to a given species by using the method `add_particles` of the `WarpXParticleContainer` object,
using the following syntax:

### pywarpx.extensions.WarpXParticleContainer.add_particles(self, x=None, y=None, z=None, ux=None, uy=None, uz=None, w=None, unique_particles=True, \*\*kwargs)

A function for adding particles to the WarpX simulation.

* **Parameters:**
  * **species_name** (*str*) – The type of species for which particles will be added
  * **x** (*arrays* *or* *scalars*) – The particle positions (m) (default = 0.)
  * **y** (*arrays* *or* *scalars*) – The particle positions (m) (default = 0.)
  * **z** (*arrays* *or* *scalars*) – The particle positions (m) (default = 0.)
  * **ux** (*arrays* *or* *scalars*) – The particle proper velocities (m/s) (default = 0.)
  * **uy** (*arrays* *or* *scalars*) – The particle proper velocities (m/s) (default = 0.)
  * **uz** (*arrays* *or* *scalars*) – The particle proper velocities (m/s) (default = 0.)
  * **w** (*array* *or* *scalars*) – Particle weights (default = 0.)
  * **unique_particles** (*bool*) – True means the added particles are duplicated by each process;
    False means the number of added particles is independent of
    the number of processes (default = True)
  * **kwargs** (*dict*) – Containing an entry for all the extra particle attribute arrays. If
    an attribute is not given it will be set to 0.

### See this function used in a full example

```python
#!/usr/bin/env python3
# --- Input file for particle-boundary interaction testing in RZ.
# --- This input is a simple case of reflection
# --- of one electron on the surface of a sphere.

import numpy as np
import scipy.constants as scc

from pywarpx import callbacks, particle_containers, picmi
from pywarpx.LoadThirdParty import load_cupy

##########################
# numerics parameters
##########################

dt = 1.0e-11

# --- Nb time steps

max_steps = 23
diagnostic_interval = 1

# --- grid

nr = 64
nz = 64

rmin = 0.0
rmax = 2
zmin = -2
zmax = 2

##########################
# numerics components
##########################

grid = picmi.CylindricalGrid(
    number_of_cells=[nr, nz],
    n_azimuthal_modes=1,
    lower_bound=[rmin, zmin],
    upper_bound=[rmax, zmax],
    lower_boundary_conditions=["none", "dirichlet"],
    upper_boundary_conditions=["dirichlet", "dirichlet"],
    lower_boundary_conditions_particles=["none", "reflecting"],
    upper_boundary_conditions_particles=["absorbing", "reflecting"],
)


solver = picmi.ElectrostaticSolver(
    grid=grid, method="Multigrid", warpx_absolute_tolerance=1e-7
)

embedded_boundary = picmi.EmbeddedBoundary(
    implicit_function="-(x**2+y**2+z**2-radius**2)", radius=0.2
)

##########################
# physics components
##########################

# one particle
e_dist = picmi.ParticleListDistribution(
    x=0.0, y=0.0, z=-0.25, ux=0.5e10, uy=0.0, uz=1.0e10, weight=1
)

electrons = picmi.Species(
    particle_type="electron",
    name="electrons",
    initial_distribution=e_dist,
    warpx_save_particles_at_eb=1,
)

##########################
# diagnostics
##########################

field_diag = picmi.FieldDiagnostic(
    name="diag1",
    grid=grid,
    period=diagnostic_interval,
    data_list=["Er", "Ez", "phi", "rho", "rho_electrons"],
    warpx_format="openpmd",
)

part_diag = picmi.ParticleDiagnostic(
    name="diag1",
    period=diagnostic_interval,
    species=[electrons],
    warpx_format="openpmd",
)

##########################
# simulation setup
##########################

sim = picmi.Simulation(
    solver=solver,
    time_step_size=dt,
    max_steps=max_steps,
    warpx_embedded_boundary=embedded_boundary,
    warpx_amrex_the_arena_is_managed=1,
)

sim.add_species(
    electrons,
    layout=picmi.GriddedLayout(n_macroparticle_per_cell=[10, 1, 1], grid=grid),
)
sim.add_diagnostic(part_diag)
sim.add_diagnostic(field_diag)

sim.initialize_inputs()
sim.initialize_warpx()

##########################
# python particle data access
##########################
xp, _ = load_cupy()


def concat(list_of_arrays):
    if len(list_of_arrays) == 0:
        # Return a 1d array of size 0
        return xp.empty(0)
    else:
        return xp.concatenate(list_of_arrays)


def to_numpy(arr):
    if hasattr(arr, "get"):
        return arr.get()
    else:
        return arr


def mirror_reflection():
    buffer = particle_containers.ParticleBoundaryBufferWrapper()  # boundary buffer

    # STEP 1: extract the different parameters of the boundary buffer (normal, time, position)
    lev = 0  # level 0 (no mesh refinement here)
    delta_t = concat(
        buffer.get_particle_scraped_this_step(
            "electrons", "eb", "deltaTimeScraped", lev
        )
    )
    r = concat(buffer.get_particle_scraped_this_step("electrons", "eb", "r", lev))
    theta = concat(
        buffer.get_particle_scraped_this_step("electrons", "eb", "theta", lev)
    )

    z = concat(buffer.get_particle_scraped_this_step("electrons", "eb", "z", lev))
    x = r * xp.cos(theta)  # from RZ coordinates to 3D coordinates
    y = r * xp.sin(theta)
    ux = concat(buffer.get_particle_scraped_this_step("electrons", "eb", "ux", lev))
    uy = concat(buffer.get_particle_scraped_this_step("electrons", "eb", "uy", lev))
    uz = concat(buffer.get_particle_scraped_this_step("electrons", "eb", "uz", lev))
    w = concat(buffer.get_particle_scraped_this_step("electrons", "eb", "w", lev))
    nx = concat(buffer.get_particle_scraped_this_step("electrons", "eb", "nx", lev))
    ny = concat(buffer.get_particle_scraped_this_step("electrons", "eb", "ny", lev))
    nz = concat(buffer.get_particle_scraped_this_step("electrons", "eb", "nz", lev))

    # STEP 2: use these parameters to inject particle from the same position in the plasma
    electrons = sim.particles.get("electrons")  # general particle container

    ####this part is specific to the case of simple reflection.
    un = ux * nx + uy * ny + uz * nz
    ux_reflect = -2 * un * nx + ux  # for a "mirror reflection" u(sym)=-2(u.n)n+u
    uy_reflect = -2 * un * ny + uy
    uz_reflect = -2 * un * nz + uz

    x = to_numpy(x)
    y = to_numpy(y)
    z = to_numpy(z)
    w = to_numpy(w)
    delta_t = to_numpy(delta_t)
    ux_reflect = to_numpy(ux_reflect)
    uy_reflect = to_numpy(uy_reflect)
    uz_reflect = to_numpy(uz_reflect)

    inv_c2 = 1.0 / (scc.c**2)
    inv_gamma = 1.0 / np.sqrt(
        1.0 + (ux_reflect**2 + uy_reflect**2 + uz_reflect**2) * inv_c2
    )
    dt_remaining = dt - delta_t

    electrons.add_particles(
        x=x + dt_remaining * ux_reflect * inv_gamma,
        y=y + dt_remaining * uy_reflect * inv_gamma,
        z=z + dt_remaining * uz_reflect * inv_gamma,
        ux=ux_reflect,
        uy=uy_reflect,
        uz=uz_reflect,
        w=w,
    )  # adds the particle in the general particle container at the next step
    #### Can be modified depending on the model of interaction.


callbacks.installafterstep(
    mirror_reflection
)  # mirror_reflection is called at the next step
# using the new particle container modified at the last step

##########################
# simulation run
##########################

sim.step(max_steps)  # the whole process is done "max_steps" times
```
