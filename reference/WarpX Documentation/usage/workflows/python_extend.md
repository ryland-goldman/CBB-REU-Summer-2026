<a id="usage-python-extend"></a>

# Extend a Simulation with Python

## Overview

WarpX’s Python bindings let you integrate Python code directly into a WarpX simulation.
Through this interface, you can **access and modify simulation data** – such as particle properties, field values – as the simulation runs.
This versatility opens the door to a wide range of workflows, including:

> - **Adding a custom physics module** (for instance, a specific collision model) that may not yet be available in WarpX’s C++ implementation, and that can be quickly implemented in Python.
> - **Coupling WarpX with another simulation tool** that has a Python interface, enabling both codes to operate on the same particle or field data.
> - **Incorporating AI-based surrogate models** built in Python (e.g., with PyTorch or TensorFlow) to emulate complex physical processes.

If your custom Python code uses high-performance, GPU-accelerated libraries – such as [cupy](https://cupy.dev/), [pytorch](https://pytorch.org/),
or [numba](https://numba.pydata.org/) – the extra computations are unlikely to significantly impact simulation speed.
Note that WarpX’s Python bindings provide direct access to particle and field data without creating copies, resulting in very low overhead.

<a id="usage-python-extend-run-simulation"></a>

## How to run a simulation with Python extensions

- **Install WarpX with support for the Python interface**: for instance, if you [compile WarpX from source](../../install/cmake.md#install-build-code), this involves using `-DWarpX_PYTHON=ON`.
- **Write a Python script that extends the simulation**: this can be done starting from a simulation defined either with a [parameter list](../parameters.md#running-cpp-parameters) or with the [PICMI Python interface](../python.md#usage-picmi).
  The Python script typically contains [callback functions](#usage-python-extend-callbacks) that [access/modify](#usage-python-extend-data-access) the simulation data (see the sections below for more details).

### Parameter List

When starting from a [parameter list](../parameters.md#running-cpp-parameters), write a Python script that loads the parameter list file using the `load_inputs_file()` method:

```python3
from pywarpx import warpx

sim = warpx
sim.load_inputs_file("./inputs_test_3d_laser_acceleration")

# register callbacks ...

# advance simulation until the last time step
sim.step()
```

### Full Example

```python3
#!/usr/bin/env python3
#
# Starting from an inputs file, define a WarpX simulation
# and extend it with Python logic.

from pywarpx import warpx
from pywarpx.callbacks import callfromafterstep

sim = warpx
sim.load_inputs_file("./inputs_test_3d_laser_acceleration")


# Optional: Define callbacks, e.g., after every step
@callfromafterstep
def my_simple_callback():
    """This simple callback uses particle container and MultiFab objects,
    https://warpx.readthedocs.io/en/latest/usage/workflows/python_extend.html#particles
    and
    https://warpx.readthedocs.io/en/latest/usage/workflows/python_extend.html#fields
    """
    print("  my_simple_callback")

    # electrons: access (and potentially manipulate)
    electrons = sim.particles.get("electrons")
    print(f"    {electrons}")

    # electric field: access (and potentially manipulate)
    Ex = sim.fields.get("Efield_fp", dir="x", level=0)
    print(f"    {Ex}")


@callfromafterstep
def my_advanced_callback():
    """This callback dives deeper using pyAMReX methods and data containers directly.
    https://pyamrex.readthedocs.io/en/latest/usage/compute.html
    """
    print("  my_advanced_callback")

    # the pyAMReX module
    amr = sim.extension.amr
    amr.Print(f"    {amr.ParallelDescriptor.NProcs()} MPI process(es) active")

    # electrons: access (and potentially manipulate)
    electrons = sim.particles.get("electrons")
    print(f"    {electrons}")

    # electric field: access (and potentially manipulate)
    Ex_mf = sim.fields.get("Efield_fp", dir="x", level=0)
    print(f"    {Ex_mf}")


# Advance simulation until the last time step
sim.step()
```

### PICMI

When starting from a [PICMI Python script](../python.md#usage-picmi), simply add the Python code that extends the simulation to this script, before the call to [`step()`](../python.md#pywarpx.picmi.Simulation.step).

```python3
# Preparation: set up the simulation
#   sim = picmi.Simulation(...)
#   ...

# register callbacks ...

sim.step(nsteps=1000)
```

- **Then, run the simulation by executing the Python script**: for instance using `mpirun` or `srun` on an HPC system.

```bash
mpirun -np <n_ranks> python <python_script>
```

<a id="usage-python-extend-callbacks"></a>

## Callback Functions

Installing [callback functions](https://en.wikipedia.org/wiki/Callback_(computer_programming)) will execute a given Python function at a
specific location in the WarpX simulation loop. The syntax to use in order to define callback functions is described in the links below.

* [Callback Locations](python_callbacks.md)
* [`installcallback()`](python_callbacks.md#pywarpx.callbacks.installcallback)
* [`isinstalled()`](python_callbacks.md#pywarpx.callbacks.isinstalled)
* [`uninstallcallback()`](python_callbacks.md#pywarpx.callbacks.uninstallcallback)

<a id="usage-python-extend-data-access"></a>

## Accessing simulation data through Python

While the simulation is running, the Python code (e.g. the code in the callback functions) will have read and write access the WarpX simulation data.
The specific Python syntax to access this data is described in the following sections.

* [Accessing fields data](python_field_data.md)
* [Accessing particles data](python_particle_data.md)
* [Accessing the particles that hit the boundaries](python_particle_boundary_data.md)
* [Accessing global WarpX functionalities (e.g., extract timestep)](python_warpx.md)
* [Writing portable Python code that can be executed on CPU and GPU](python_portable.md)
