# Accessing fields data

## Selecting a given field

The simulation’s fields are accessed through the `sim.fields` attribute, where `sim` is obtained as shown
in the [How to run a simulation with Python extensions](python_extend.md#usage-python-extend-run-simulation) section. Specific fields (e.g. electric field, charge density, etc.)
are selected with the `sim.fields.get` method, as shown in the example below.

```python
# Preparation: set up the sim object
#   sim = picmi.Simulation(...)
#   ...

# Extract the Ex field, at level 0 of mesh refinement
Ex = sim.fields.get("Efield_fp", dir="x", level=0)
```

The available field names (e.g. `"Efield_fp"`, `"rho_fp"`, etc.) are listed in the [Field Names](../../developers/fields.md#developers-fields-names) section.
The function `sim.fields.get` returns a [pyamrex](https://pyamrex.readthedocs.io/en/latest/index.html) object of
type [MultiFab](https://pyamrex.readthedocs.io/en/latest/usage/api.html#amrex.space3d.MultiFab), whose field data can be accessed or modified as described below.

## Accessing/modifying the underlying field data

Several ways to access and modify the field data (i.e., the values of the fields on the grid points) are available.
These different methods differ in their user-friendliness, flexibility and performance overhead.

### Pre-defined AMReX methods

The AMReX library defines many functions that can operate on field data, and many of them are accessible
from Python via the [pyamrex](https://pyamrex.readthedocs.io/en/latest/index.html) library. For a
full list of these methods, see the [pyamrex API documentation](https://pyamrex.readthedocs.io/en/latest/usage/api.html#amrex.space3d.MultiFab).

Examples include:

- Finding the maximum value of the field over the entire domain: `Ex.max()`
- Scaling the field by a factor of 2: `Ex.mult(2.)`
- Adding two fields together: `Ex.saxpy(...)` (see [this link](https://pyamrex.readthedocs.io/en/latest/usage/api.html#amrex.space3d.MultiFab.saxpy))

#### NOTE
These methods generally have high performance and low overhead, including when using GPUs and multi-node parallelization, but are limited to the existing functions provided by AMReX.

### See some of these methods used in a full example

```python
#!/usr/bin/env python3
#
# --- Input file for spacecraft charging testing in RZ.
# --- This input defines a conducting sphere (spacecraft) immersed in a thermal
# --- plasma with the same given initial conditions as in the article:
# --- (*) J. Deca, G. Lapenta, R. Marchand, S. Markidis;
# ---     Spacecraft charging analysis with the implicit particle-in-cell code iPic3D.
# ---     Part III. A. pages 3-4
# ---     Phys. Plasmas 1 October 2013; 20 (10): 102902. https://doi.org/10.1063/1.4826951.
# --- The conducting sphere starts with an initial potential of 1V and will interact with
# --- the surrounding plasma, initially static. The charging of the spacecraft - by accumulation
# --- of electrons - leads to a decrease of the potential on the surface over the time
# --- until reaching an equilibrium floating potential of ~144.5 V (*).

import numpy as np
import scipy.constants as scc
from mpi4py import MPI as mpi

from pywarpx import picmi
from pywarpx.callbacks import installafterEsolve, installafterInitEsolve
from pywarpx.particle_containers import ParticleBoundaryBufferWrapper


# Utilities
class SpaceChargeFieldCorrector(object):
    """
    Class used by the callback functions to calculate the
    correct field around the spacecraft, at each timestep
    (taking into account the charge that has been collected on the spacecraft)
    """

    def __init__(self):
        self.saved_first_iteration_fields = False
        self.spacecraft_potential = 1.0  # Initial voltage: 1V
        self.spacecraft_capacitance = None

    def correct_space_charge_fields(self, q=None):
        """
        Function that will be called at each iteration,
        after each electrostatic solve in WarpX
        """
        assert self.saved_first_iteration_fields

        # Compute the charge that WarpX thinks there is on the spacecraft
        # from phi and rho after the Poisson solver
        q_v = compute_virtual_charge_on_spacecraft()
        if q is None:
            q = compute_actual_charge_on_spacecraft()

        # Correct fields so as to recover the actual charge
        warpx = sim.extension.warpx
        Er = sim.fields.get("Efield_fp", dir="r", level=0)
        normalized_Er = sim.fields.get("normalized_Er", level=0)
        Er.saxpy(q - q_v, normalized_Er, 0, 0, 1, 0)
        Ez = sim.fields.get("Efield_fp", dir="z", level=0)
        normalized_Ez = sim.fields.get("normalized_Ez", level=0)
        Ez.saxpy(q - q_v, normalized_Ez, 0, 0, 1, 0)
        phi = sim.fields.get("phi_fp", level=0)
        normalized_phi = sim.fields.get("normalized_phi", level=0)
        phi.saxpy(q - q_v, normalized_phi, 0, 0, 1, 0)

        self.spacecraft_potential += (q - q_v) * self.spacecraft_capacitance
        warpx.set_potential_on_eb("%f" % self.spacecraft_potential)
        print("Setting potential to %f" % self.spacecraft_potential)

        # Confirm that the charge on the spacecraft is now correct
        compute_virtual_charge_on_spacecraft()

    def save_normalized_vacuum_Efields(
        self,
    ):
        # Compute the charge that WarpX thinks there is on the spacecraft
        # from phi and rho after the Poisson solver
        q_v = compute_virtual_charge_on_spacecraft()
        self.spacecraft_capacitance = 1.0 / q_v  # the potential was set to 1V

        phi = sim.fields.get("phi_fp", level=0)
        Er = sim.fields.get("Efield_fp", dir="r", level=0)
        Ez = sim.fields.get("Efield_fp", dir="z", level=0)
        # Allocate the fields `normalized_Er`, `normalized_Ez`, and `normalized_phi
        # in WarpX's multifab register. This allows to get these fields at later
        # iterations with sim.fields.get( ... ).
        # These new fields are automatically redistributed when doing load balancing.
        normalized_Er = sim.fields.alloc_init(
            name="normalized_Er",
            level=0,
            ba=Er.box_array(),
            dm=Er.dm(),
            ncomp=Er.n_comp,
            ngrow=Er.n_grow_vect,
            initial_value=0.0,
            redistribute=True,
            redistribute_on_remake=True,
        )

        normalized_Ez = sim.fields.alloc_init(
            name="normalized_Ez",
            level=0,
            ba=Ez.box_array(),
            dm=Ez.dm(),
            ncomp=Ez.n_comp,
            ngrow=Ez.n_grow_vect,
            initial_value=0.0,
            redistribute=True,
            redistribute_on_remake=True,
        )

        normalized_phi = sim.fields.alloc_init(
            name="normalized_phi",
            level=0,
            ba=phi.box_array(),
            dm=phi.dm(),
            ncomp=phi.n_comp,
            ngrow=phi.n_grow_vect,
            initial_value=0.0,
            redistribute=True,
            redistribute_on_remake=True,
        )

        # Record fields

        normalized_Er.copymf(Er, 0, 0, 1, Er.n_grow_vect)
        normalized_Er.mult(1 / q_v, 0)
        normalized_Ez.copymf(Ez, 0, 0, 1, Ez.n_grow_vect)
        normalized_Ez.mult(1 / q_v, 0)
        normalized_phi.copymf(phi, 0, 0, 1, phi.n_grow_vect)
        normalized_phi.mult(1 / q_v, 0)

        self.saved_first_iteration_fields = True
        self.correct_space_charge_fields(q=0)


def compute_virtual_charge_on_spacecraft():
    """
    Given that we asked WarpX to solve the Poisson
    equation with phi=1 on the spacecraft and phi=0
    on the boundary of the domain, compute the charge
    that WarpX thinks there should be on the spacecraft.
    """
    warpx = sim.extension.warpx
    rho = sim.fields.get("rho_fp", level=0)
    phi = sim.fields.get("phi_fp", level=0)

    dr, dz = warpx.Geom(lev=0).data().CellSize()

    # Compute integral of grad phi over surfaces of the domain
    nr = phi.shape[0]
    r = np.linspace(rmin, rmax, nr, endpoint=False) + (rmax - rmin) / (
        2 * nr
    )  # shift of the r points because the derivaties are calculated in the middle
    face_z0 = (
        2 * np.pi * 1.0 / dz * ((phi[:, 0] - phi[:, 1]) * r).sum() * dr
    )  # here I am assuming that phi is a numpy array that can handle elementwise mult
    face_zend = 2 * np.pi * 1.0 / dz * ((phi[:, -1] - phi[:, -2]) * r).sum() * dr
    face_rend = 2 * np.pi * 1.0 / dr * ((phi[-1, :] - phi[-2, :]) * rmax).sum() * dz
    grad_phi_integral = face_z0 + face_zend + face_rend

    # Compute integral of rho over volume of the domain
    # (i.e. total charge of the plasma particles)
    rho_integral = (
        (rho[1 : nr - 1, 1 : nz - 1] * r[1 : nr - 1, np.newaxis]).sum()
        * 2
        * np.pi
        * dr
        * dz
    )

    # Compute charge of the spacecraft, based on Gauss theorem
    q_spacecraft = -rho_integral - scc.epsilon_0 * grad_phi_integral
    print("Virtual charge on the spacecraft: %e" % q_spacecraft)
    return q_spacecraft


def compute_actual_charge_on_spacecraft():
    """
    Compute the actual charge on the spacecraft,
    by counting how many electrons and protons
    were collected by the WarpX embedded boundary (EB)
    """
    charge = {"electrons": -scc.e, "protons": scc.e}
    q_spacecraft = 0
    particle_buffer = ParticleBoundaryBufferWrapper()
    for species in charge.keys():
        weights = particle_buffer.get_particle_boundary_buffer(species, "eb", "w", 0)
        sum_weights_over_tiles = sum([w.sum() for w in weights])

        # Reduce across all MPI ranks
        ntot = float(mpi.COMM_WORLD.allreduce(sum_weights_over_tiles, op=mpi.SUM))
        print("Total number of %s collected on spacecraft: %e" % (species, ntot))
        q_spacecraft += ntot * charge[species]

    print("Actual charge on the spacecraft: %e" % q_spacecraft)
    return q_spacecraft


##########################
# numerics parameters
##########################

dt = 1.27e-8

# --- Nb time steps
max_steps = 1000
diagnostic_interval = 10

# --- grid
nr = 40
nz = 80

rmin = 0.0
rmax = 3
zmin = -3
zmax = 3

number_per_cell = 5
number_per_cell_each_dim = [10, 1, 1]


##########################
# physics components
##########################

n = 7.0e9  # plasma density #particles/m^3
Te = 85  # Electron temp in eV
Ti = 0.05 * Te  # Ion temp in eV
qe = picmi.constants.q_e  # elementary charge
m_e = picmi.constants.m_e  # electron mass
m_i = 1836.0 * m_e  # mass of ion
v_eth = (qe * Te / m_e) ** 0.5
v_pth = (qe * Ti / m_i) ** 0.5

# nothing to change in the distribution function?
e_dist = picmi.UniformDistribution(density=n, rms_velocity=[v_eth, v_eth, v_eth])
e_dist2 = picmi.UniformFluxDistribution(
    flux=n * v_eth / (2 * np.pi) ** 0.5,  # Flux for Gaussian with vmean=0
    surface_flux_position=3,
    flux_direction=-1,
    flux_normal_axis="r",
    gaussian_flux_momentum_distribution=True,
    rms_velocity=[v_eth, v_eth, v_eth],
)
electrons = picmi.Species(
    particle_type="electron",
    name="electrons",
    initial_distribution=[e_dist, e_dist2],
    warpx_save_particles_at_eb=1,
)

p_dist = picmi.UniformDistribution(density=n, rms_velocity=[v_pth, v_pth, v_pth])
p_dist2 = picmi.UniformFluxDistribution(
    flux=n * v_pth / (2 * np.pi) ** 0.5,  # Flux for Gaussian with vmean=0
    surface_flux_position=3,
    flux_direction=-1,
    flux_normal_axis="r",
    gaussian_flux_momentum_distribution=True,
    rms_velocity=[v_pth, v_pth, v_pth],
)
protons = picmi.Species(
    particle_type="proton",
    name="protons",
    initial_distribution=[p_dist, p_dist2],
    warpx_save_particles_at_eb=1,
)


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
    implicit_function="-(x**2+y**2+z**2-radius**2)",
    potential=1.0,  # arbitrary value ; this will be corrected by a callback function
    radius=0.3277,
)


##########################
# diagnostics
##########################

field_diag = picmi.FieldDiagnostic(
    name="diag1",
    grid=grid,
    period=diagnostic_interval,
    data_list=["Er", "Ez", "phi", "rho", "rho_electrons", "rho_protons"],
    warpx_format="openpmd",
)

part_diag = picmi.ParticleDiagnostic(
    name="diag1",
    period=diagnostic_interval,
    species=[electrons, protons],
    warpx_format="openpmd",
)

part_scraping_boundary_diag = picmi.ParticleBoundaryScrapingDiagnostic(
    name="diag2",
    period=-1,  # only at the end, because we also use the buffers in this test
    species=[electrons, protons],
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
    warpx_random_seed=1,
)

layout1 = picmi.GriddedLayout(
    n_macroparticle_per_cell=number_per_cell_each_dim, grid=grid
)
layout2 = picmi.PseudoRandomLayout(n_macroparticles_per_cell=number_per_cell, grid=grid)
sim.add_species(electrons, layout=[layout1, layout2])

sim.add_species(protons, layout=[layout1, layout2])

sim.add_diagnostic(field_diag)
sim.add_diagnostic(part_diag)
sim.add_diagnostic(part_scraping_boundary_diag)

##########################
# simulation run
##########################

spc = SpaceChargeFieldCorrector()

installafterInitEsolve(spc.save_normalized_vacuum_Efields)
installafterEsolve(spc.correct_space_charge_fields)

sim.step(max_steps)
```

### Numpy-like global indexing

The field data in a `MultiFab` object can also be accessed via global indexing.
Using standard array indexing with square brackets, the data can be accessed using indices that are relative to the full domain (across the `MultiFab` and across processors).
When the data is fetched the result is a `numpy` array that contains a copy of the data, and when using multiple processors is broadcast to all processors (and is a global operation).

#### WARNING
Global indexing is convenient and user-friendly, but has significant performance overheads,
since it potentially incurs MPI communications and CPU-GPU copies under the hood.
This method is thus mostly meant for debugging and visualization purposes,
and not for performance-critical operations.

For indices within the domain, values from valid cells are always returned. The ghost cells at the exterior of the domain are
accessed using imaginary numbers, with negative values accessing the lower ghost cells, and positive the upper ghost cells.
This example will return the `Bz` field at all valid interior points along `x` at the specified `y` and `z` indices.

```python
Bz = sim.fields.get("Bfield_fp", dir=2, level=0)
Bz_along_x = Bz[:,5,6]
```

The same global indexing can be done to set values. This example will set the values over a range in `y` and `z` at the
specified `x`. The data will be scattered appropriately to the underlying FABs. Setting values is a local operation.

```python
Jy = sim.fields.get("current_fp", dir=1, level=0)
Jy[5,6:20,8:30] = 7.
```

In the example below, 7 is added to all of the values along `x`, including both valid and ghost cells (specified by using the empty tuple,
`()`), the first ghost cell at the lower boundary in `y`, and the last valid cell and first upper ghost cell in `z`.
Note that the `+=` will be a global operation.

```python
Jx = sim.fields.get("current_fp", dir=0, level=0)
Jx[(),-1j,-1:2j] += 7.
```

Instead of setting values with a scalar value, you can also set values using an array. The array shape must match the selected region.
The array must be either a `numpy` array (if WarpX is run on CPU) or a `cupy` array (if WarpX is run on GPU).
To write portable code that works on both CPU and GPU, it is recommended to use the load_cupy package.
See [Writing portable Python code that can be executed on CPU and GPU](python_portable.md#usage-python-portable) for more details on writing portable Python code.

```python
from pywarpx.LoadThirdParty import load_cupy
xp = load_cupy()

Jy = sim.fields.get("current_fp", dir=1, level=0)
# Create random values with shape matching the selected region
random_values = xp.random.random((14, 22))  # shape: (20-6, 30-8)
Jy[5,6:20,8:30] = random_values
```

To fetch the data from all of the valid cells of all dimensions, the ellipsis can be used, `Jx[...]`.
Similarly, to fetch all of the data including valid cells and ghost cells, use an empty tuple, `Jx[()]`.
The code does error checking to ensure that the specified indices are within the bounds of the global domain.

Finally, the `mesh` method returns the physical coordinates of the mesh along a specified direction,
with appropriate centering based on the field’s staggering. This is useful for plotting,
analysis, or when you need to know the physical positions corresponding to field values.

```python
Ex = sim.fields.get("Efield_fp", dir="x", level=0)
x_coords = Ex.mesh("x")
y_coords = Ex.mesh("y")
z_coords = Ex.mesh("z")
```

The method accepts a direction string (`"x"`, `"y"`, `"z"` in 3D; `"r"`, `"z"` in RZ geometry)
and an optional `include_ghosts` parameter (default `False`) to include ghost cell coordinates.
The returned array contains the physical coordinates of the mesh points along the specified direction,
properly accounting for the field’s cell-centered or face-centered staggering.

### Explicit loop over boxes

This method provides similar capabilities to the numpy-like global indexing approach, but operates
only on local data within each MPI rank. Unlike global indexing, which may involve MPI communications
and CPU-GPU data transfers under the hood, this approach performs all operations locally on each processor.
As a result, this method offers significantly higher performance, especially for large-scale parallel simulations
and GPU-accelerated runs. The data is accessed by explicitly looping over mesh-refinement levels and
individual grid blocks (boxes), giving you direct access to the underlying `numpy` or `cupy` arrays for each local block.

The example below accesses the $Ex(x,y,z)$ field at level 0 after every time step and sets all of the values to `42`.
This shows how to loop over levels and grid blocks.

```python
from pywarpx import picmi
from pywarpx.callbacks import callfromafterstep

# Preparation: set up the simulation
#   sim = picmi.Simulation(...)
#   ...

# Extract the Ex field, at level 0 of mesh refinement
Ex = sim.fields.get("Efield_fp", dir="x", level=0)

# compute on Ex
# iterate over mesh-refinement levels
for lev in range(warpx.finest_level + 1):
    # grow (aka guard/ghost/halo) regions
    ngv = Ex.n_grow_vect

    # get every local block of the field
    for mfi in Ex:
        # global index space box, including guards
        bx = mfi.tilebox().grow(ngv)
        print(bx)  # note: global index space of this block

    # numpy/cupy representation of the field data, including
    # the guard/ghost region
    Ex_arr = Ex.array(mfi).to_xp()

    # notes on indexing in Ex:
    # - numpy/cupy use locally zero-based indexing
    # - layout is F_CONTIGUOUS by default, just like AMReX

    # notes:
    # Only the next lines are the "HOT LOOP" of the computation.
    # For efficiency, we use array operation for speed.
    Ex_arr[()] = 42.0
```

For further details on how to [access GPU data](https://pyamrex.readthedocs.io/en/latest/usage/zerocopy.html) or compute on `Ex`, please see the [pyAMReX documentation](https://pyamrex.readthedocs.io/en/latest/usage/compute.html#fields).

## Defining a new custom field

For some use cases, it is sometimes needed to create new custom fields (in addition to the [existing fields in WarpX](../../developers/fields.md#developers-fields-names)).
New `MultiFab` objects can be created at the Python level. Using this method, the new `MultiFab` will be handled in the same way as WarpX’s internal `MultiFab`.
For example, their data will be automatically redistributed during load balancing (when the flags are set as shown in the example).

In the example below, a new `MultiFab` is created with the same properties as `Ex`.

```python
Ex = sim.fields.get("Efield_fp", dir=0, level=0)
normalized_Ex = sim.fields.alloc_init(name="normalized_Ex",
                                      dir=0,
                                      level=0,
                                      ba=Ex.box_array(),
                                      dm=Ex.dm(),
                                      ncomp=Ex.n_comp,
                                      ngrow=Ex.n_grow_vect,
                                      initial_value=0.,
                                      redistribute=True,
                                      redistribute_on_remake=True)
```

### See this function used in a full example

```python
#!/usr/bin/env python3
#
# --- Input file for spacecraft charging testing in RZ.
# --- This input defines a conducting sphere (spacecraft) immersed in a thermal
# --- plasma with the same given initial conditions as in the article:
# --- (*) J. Deca, G. Lapenta, R. Marchand, S. Markidis;
# ---     Spacecraft charging analysis with the implicit particle-in-cell code iPic3D.
# ---     Part III. A. pages 3-4
# ---     Phys. Plasmas 1 October 2013; 20 (10): 102902. https://doi.org/10.1063/1.4826951.
# --- The conducting sphere starts with an initial potential of 1V and will interact with
# --- the surrounding plasma, initially static. The charging of the spacecraft - by accumulation
# --- of electrons - leads to a decrease of the potential on the surface over the time
# --- until reaching an equilibrium floating potential of ~144.5 V (*).

import numpy as np
import scipy.constants as scc
from mpi4py import MPI as mpi

from pywarpx import picmi
from pywarpx.callbacks import installafterEsolve, installafterInitEsolve
from pywarpx.particle_containers import ParticleBoundaryBufferWrapper


# Utilities
class SpaceChargeFieldCorrector(object):
    """
    Class used by the callback functions to calculate the
    correct field around the spacecraft, at each timestep
    (taking into account the charge that has been collected on the spacecraft)
    """

    def __init__(self):
        self.saved_first_iteration_fields = False
        self.spacecraft_potential = 1.0  # Initial voltage: 1V
        self.spacecraft_capacitance = None

    def correct_space_charge_fields(self, q=None):
        """
        Function that will be called at each iteration,
        after each electrostatic solve in WarpX
        """
        assert self.saved_first_iteration_fields

        # Compute the charge that WarpX thinks there is on the spacecraft
        # from phi and rho after the Poisson solver
        q_v = compute_virtual_charge_on_spacecraft()
        if q is None:
            q = compute_actual_charge_on_spacecraft()

        # Correct fields so as to recover the actual charge
        warpx = sim.extension.warpx
        Er = sim.fields.get("Efield_fp", dir="r", level=0)
        normalized_Er = sim.fields.get("normalized_Er", level=0)
        Er.saxpy(q - q_v, normalized_Er, 0, 0, 1, 0)
        Ez = sim.fields.get("Efield_fp", dir="z", level=0)
        normalized_Ez = sim.fields.get("normalized_Ez", level=0)
        Ez.saxpy(q - q_v, normalized_Ez, 0, 0, 1, 0)
        phi = sim.fields.get("phi_fp", level=0)
        normalized_phi = sim.fields.get("normalized_phi", level=0)
        phi.saxpy(q - q_v, normalized_phi, 0, 0, 1, 0)

        self.spacecraft_potential += (q - q_v) * self.spacecraft_capacitance
        warpx.set_potential_on_eb("%f" % self.spacecraft_potential)
        print("Setting potential to %f" % self.spacecraft_potential)

        # Confirm that the charge on the spacecraft is now correct
        compute_virtual_charge_on_spacecraft()

    def save_normalized_vacuum_Efields(
        self,
    ):
        # Compute the charge that WarpX thinks there is on the spacecraft
        # from phi and rho after the Poisson solver
        q_v = compute_virtual_charge_on_spacecraft()
        self.spacecraft_capacitance = 1.0 / q_v  # the potential was set to 1V

        phi = sim.fields.get("phi_fp", level=0)
        Er = sim.fields.get("Efield_fp", dir="r", level=0)
        Ez = sim.fields.get("Efield_fp", dir="z", level=0)
        # Allocate the fields `normalized_Er`, `normalized_Ez`, and `normalized_phi
        # in WarpX's multifab register. This allows to get these fields at later
        # iterations with sim.fields.get( ... ).
        # These new fields are automatically redistributed when doing load balancing.
        normalized_Er = sim.fields.alloc_init(
            name="normalized_Er",
            level=0,
            ba=Er.box_array(),
            dm=Er.dm(),
            ncomp=Er.n_comp,
            ngrow=Er.n_grow_vect,
            initial_value=0.0,
            redistribute=True,
            redistribute_on_remake=True,
        )

        normalized_Ez = sim.fields.alloc_init(
            name="normalized_Ez",
            level=0,
            ba=Ez.box_array(),
            dm=Ez.dm(),
            ncomp=Ez.n_comp,
            ngrow=Ez.n_grow_vect,
            initial_value=0.0,
            redistribute=True,
            redistribute_on_remake=True,
        )

        normalized_phi = sim.fields.alloc_init(
            name="normalized_phi",
            level=0,
            ba=phi.box_array(),
            dm=phi.dm(),
            ncomp=phi.n_comp,
            ngrow=phi.n_grow_vect,
            initial_value=0.0,
            redistribute=True,
            redistribute_on_remake=True,
        )

        # Record fields

        normalized_Er.copymf(Er, 0, 0, 1, Er.n_grow_vect)
        normalized_Er.mult(1 / q_v, 0)
        normalized_Ez.copymf(Ez, 0, 0, 1, Ez.n_grow_vect)
        normalized_Ez.mult(1 / q_v, 0)
        normalized_phi.copymf(phi, 0, 0, 1, phi.n_grow_vect)
        normalized_phi.mult(1 / q_v, 0)

        self.saved_first_iteration_fields = True
        self.correct_space_charge_fields(q=0)


def compute_virtual_charge_on_spacecraft():
    """
    Given that we asked WarpX to solve the Poisson
    equation with phi=1 on the spacecraft and phi=0
    on the boundary of the domain, compute the charge
    that WarpX thinks there should be on the spacecraft.
    """
    warpx = sim.extension.warpx
    rho = sim.fields.get("rho_fp", level=0)
    phi = sim.fields.get("phi_fp", level=0)

    dr, dz = warpx.Geom(lev=0).data().CellSize()

    # Compute integral of grad phi over surfaces of the domain
    nr = phi.shape[0]
    r = np.linspace(rmin, rmax, nr, endpoint=False) + (rmax - rmin) / (
        2 * nr
    )  # shift of the r points because the derivaties are calculated in the middle
    face_z0 = (
        2 * np.pi * 1.0 / dz * ((phi[:, 0] - phi[:, 1]) * r).sum() * dr
    )  # here I am assuming that phi is a numpy array that can handle elementwise mult
    face_zend = 2 * np.pi * 1.0 / dz * ((phi[:, -1] - phi[:, -2]) * r).sum() * dr
    face_rend = 2 * np.pi * 1.0 / dr * ((phi[-1, :] - phi[-2, :]) * rmax).sum() * dz
    grad_phi_integral = face_z0 + face_zend + face_rend

    # Compute integral of rho over volume of the domain
    # (i.e. total charge of the plasma particles)
    rho_integral = (
        (rho[1 : nr - 1, 1 : nz - 1] * r[1 : nr - 1, np.newaxis]).sum()
        * 2
        * np.pi
        * dr
        * dz
    )

    # Compute charge of the spacecraft, based on Gauss theorem
    q_spacecraft = -rho_integral - scc.epsilon_0 * grad_phi_integral
    print("Virtual charge on the spacecraft: %e" % q_spacecraft)
    return q_spacecraft


def compute_actual_charge_on_spacecraft():
    """
    Compute the actual charge on the spacecraft,
    by counting how many electrons and protons
    were collected by the WarpX embedded boundary (EB)
    """
    charge = {"electrons": -scc.e, "protons": scc.e}
    q_spacecraft = 0
    particle_buffer = ParticleBoundaryBufferWrapper()
    for species in charge.keys():
        weights = particle_buffer.get_particle_boundary_buffer(species, "eb", "w", 0)
        sum_weights_over_tiles = sum([w.sum() for w in weights])

        # Reduce across all MPI ranks
        ntot = float(mpi.COMM_WORLD.allreduce(sum_weights_over_tiles, op=mpi.SUM))
        print("Total number of %s collected on spacecraft: %e" % (species, ntot))
        q_spacecraft += ntot * charge[species]

    print("Actual charge on the spacecraft: %e" % q_spacecraft)
    return q_spacecraft


##########################
# numerics parameters
##########################

dt = 1.27e-8

# --- Nb time steps
max_steps = 1000
diagnostic_interval = 10

# --- grid
nr = 40
nz = 80

rmin = 0.0
rmax = 3
zmin = -3
zmax = 3

number_per_cell = 5
number_per_cell_each_dim = [10, 1, 1]


##########################
# physics components
##########################

n = 7.0e9  # plasma density #particles/m^3
Te = 85  # Electron temp in eV
Ti = 0.05 * Te  # Ion temp in eV
qe = picmi.constants.q_e  # elementary charge
m_e = picmi.constants.m_e  # electron mass
m_i = 1836.0 * m_e  # mass of ion
v_eth = (qe * Te / m_e) ** 0.5
v_pth = (qe * Ti / m_i) ** 0.5

# nothing to change in the distribution function?
e_dist = picmi.UniformDistribution(density=n, rms_velocity=[v_eth, v_eth, v_eth])
e_dist2 = picmi.UniformFluxDistribution(
    flux=n * v_eth / (2 * np.pi) ** 0.5,  # Flux for Gaussian with vmean=0
    surface_flux_position=3,
    flux_direction=-1,
    flux_normal_axis="r",
    gaussian_flux_momentum_distribution=True,
    rms_velocity=[v_eth, v_eth, v_eth],
)
electrons = picmi.Species(
    particle_type="electron",
    name="electrons",
    initial_distribution=[e_dist, e_dist2],
    warpx_save_particles_at_eb=1,
)

p_dist = picmi.UniformDistribution(density=n, rms_velocity=[v_pth, v_pth, v_pth])
p_dist2 = picmi.UniformFluxDistribution(
    flux=n * v_pth / (2 * np.pi) ** 0.5,  # Flux for Gaussian with vmean=0
    surface_flux_position=3,
    flux_direction=-1,
    flux_normal_axis="r",
    gaussian_flux_momentum_distribution=True,
    rms_velocity=[v_pth, v_pth, v_pth],
)
protons = picmi.Species(
    particle_type="proton",
    name="protons",
    initial_distribution=[p_dist, p_dist2],
    warpx_save_particles_at_eb=1,
)


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
    implicit_function="-(x**2+y**2+z**2-radius**2)",
    potential=1.0,  # arbitrary value ; this will be corrected by a callback function
    radius=0.3277,
)


##########################
# diagnostics
##########################

field_diag = picmi.FieldDiagnostic(
    name="diag1",
    grid=grid,
    period=diagnostic_interval,
    data_list=["Er", "Ez", "phi", "rho", "rho_electrons", "rho_protons"],
    warpx_format="openpmd",
)

part_diag = picmi.ParticleDiagnostic(
    name="diag1",
    period=diagnostic_interval,
    species=[electrons, protons],
    warpx_format="openpmd",
)

part_scraping_boundary_diag = picmi.ParticleBoundaryScrapingDiagnostic(
    name="diag2",
    period=-1,  # only at the end, because we also use the buffers in this test
    species=[electrons, protons],
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
    warpx_random_seed=1,
)

layout1 = picmi.GriddedLayout(
    n_macroparticle_per_cell=number_per_cell_each_dim, grid=grid
)
layout2 = picmi.PseudoRandomLayout(n_macroparticles_per_cell=number_per_cell, grid=grid)
sim.add_species(electrons, layout=[layout1, layout2])

sim.add_species(protons, layout=[layout1, layout2])

sim.add_diagnostic(field_diag)
sim.add_diagnostic(part_diag)
sim.add_diagnostic(part_scraping_boundary_diag)

##########################
# simulation run
##########################

spc = SpaceChargeFieldCorrector()

installafterInitEsolve(spc.save_normalized_vacuum_Efields)
installafterEsolve(spc.correct_space_charge_fields)

sim.step(max_steps)
```
