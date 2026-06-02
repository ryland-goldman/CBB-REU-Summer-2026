<a id="usage-examples"></a>

# Examples

This section allows you to **download input files** that correspond to different physical situations.

We provide two kinds of inputs:

* PICMI python input files, [with parameters described here](https://picmi-standard.github.io).
* AMReX `inputs` files, [with parameters described here](parameters.md#running-cpp-parameters),

For a complete list of all example input files, also have a look at our [Examples/](https://github.com/BLAST-WarpX/warpx/tree/development/Examples) directory.
It contains folders and subfolders with self-describing names that you can try.
All these input files are automatically tested, so they should always be up-to-date.

## Plasma-Based Acceleration

* [Laser-Wakefield Acceleration of Electrons](examples/lwfa/README.md)
* [Beam-Driven Wakefield Acceleration of Electrons](examples/pwfa/README.md)
* [In-Depth: PWFA](pwfa.md)

## Laser-Plasma Interaction

* [Laser-Ion Acceleration with a Planar Target](examples/laser_ion/README.md)
* [Plasma-Mirror](examples/plasma_mirror/README.md)

## Particle Accelerator & Beam Physics

* [Gaussian Beam](examples/gaussian_beam/README.md)
* [Beam-beam collision](examples/beam_beam_collision/README.md)
* [Free-electron laser](examples/free_electron_laser/README.md)
* [Ion-Beam Extraction from a Plasma Source](examples/ion_beam_extraction/README.md)
* [Thomson Parabola Spectrometer](examples/thomson_parabola_spectrometer/README.md)

## High Energy Astrophysical Plasma Physics

* [Ohm Solver: Magnetic Reconnection](examples/ohm_solver_magnetic_reconnection/README.md)

## Fundamental Plasma Physics

* [Langmuir Waves](examples/langmuir/README.md)
* [Capacitive Discharge](examples/capacitive_discharge/README.md)
* [Pierce Diode at the Child–Langmuir Limit](examples/pierce_diode/README.md)

<a id="examples-hybrid-model"></a>

### Kinetic-fluid Hybrid Models

WarpX includes a reduced plasma model in which electrons are treated as a massless
fluid while ions are kinetically evolved, and Ohm’s law is used to calculate
the electric field. This model is appropriate for problems in which ion kinetics
dominate (ion cyclotron waves, for instance). See the
[theory section](../theory/models_algorithms/kinetic_fluid_hybrid_model.md#theory-kinetic-fluid-hybrid-model) for more details. Several
examples and benchmarks of this kinetic-fluid hybrid model are provided below.
A few of the examples are replications of the verification tests described in
Muñoz *et al.* [[1](#id10)]. The hybrid-PIC model was added to WarpX in
[PR #3665](https://github.com/BLAST-WarpX/warpx/pull/3665) - the figures in the
examples below were generated at that time.

* [Ohm solver: Electromagnetic modes](examples/ohm_solver_em_modes/README.md)
* [Ohm solver: Cylindrical normal modes](examples/ohm_solver_em_modes/README.md#ohm-solver-cylindrical-normal-modes)
* [Ohm solver: Ion Beam R Instability](examples/ohm_solver_ion_beam_instability/README.md)
* [Ohm solver: Ion Landau Damping](examples/ohm_solver_ion_Landau_damping/README.md)

## High-Performance Computing and Numerics

The following examples are commonly used to study the performance of WarpX, e.g., for computing efficiency, scalability, and I/O patterns.
While all prior examples are used for such studies as well, the examples here need less explanation on the physics, less-detail tuning on load balancing, and often simply scale (weak or strong) by changing the number of cells, AMReX block size and number of compute units.

* [Uniform Plasma](examples/uniform_plasma/README.md)

## Manipulating fields via Python

#### NOTE
TODO: The section needs to be sorted into either science cases (above) or later sections ([workflows and Python API details](workflows/python_extend.md#usage-python-extend)).

An example of using Python to access the simulation charge density, solve the Poisson equation (using `superLU`) and write the resulting electrostatic potential back to the simulation is given in the input file below. This example uses the `fields.py` module included in the `pywarpx` library.

* [`Direct Poisson solver example`](../../../Examples/Physics_applications/capacitive_discharge/inputs_test_2d_background_mcc_picmi.py)

An example of initializing the fields by accessing their data through Python, advancing the simulation for a chosen number of time steps, and plotting the fields again through Python. The simulation runs with 128 regular cells, 8 guard cells, and 10 PML cells, in each direction. Moreover, it uses div(E) and div(B) cleaning both in the regular grid and in the PML and initializes all available electromagnetic fields (E,B,F,G) identically.

* [`Unit pulse with PML`](../../../Examples/Tests/python_wrappers/inputs_test_2d_python_wrappers_picmi.py)

## Many Further Examples, Demos and Tests

* [Field Ionization](examples/field_ionization/README.md)

WarpX runs over 200 integration tests on a variety of modeling cases, which validate and demonstrate its functionality.
Please see the [Examples/Tests/](https://github.com/BLAST-WarpX/warpx/tree/development/Examples/Tests) directory for many more examples.

## Example References
