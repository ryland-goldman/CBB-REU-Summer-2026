# Callback Locations

These are the functions which allow installing user created functions so that
they are called at various places along the time step.

The following three functions allow the user to install, uninstall and verify
the different call back types.

* [`installcallback()`](#pywarpx.callbacks.installcallback): Installs a function to be called at that specified time
* [`uninstallcallback()`](#pywarpx.callbacks.uninstallcallback): Uninstalls the function (so it won’t be called anymore)
* [`isinstalled()`](#pywarpx.callbacks.isinstalled): Checks if the function is installed

These functions all take a callback location name (string) and function or
instance method as an argument. Note that if an instance method is used, an
extra reference to the method’s object is saved.

Functions can be called at the following times:

* `loadExternalFields`: during `WarpX::LoadExternalFields` to write `B/Efield_fp_external` values
* `beforeInitEsolve`: before the initial solve for the E fields (i.e. before the PIC loop starts)
* `afterInitEsolve`: after the initial solve for the E fields (i.e. before the PIC loop starts)
* `afterinit`: immediately after the init is complete
* `beforeEsolve`: before the solve for E fields (not called during init E solve, use beforeInitEsolve to apply to first solve)
* `poissonsolver`: In place of the computePhi call but only in an electrostatic simulation
* `afterEsolve`: after the solve for E fields (not called after init E solve, use afterInitEsolve to apply to first solve)
* `afterBpush`: after the B field advance for electromagnetic solvers
* `afterEpush`: after the E field advance for electromagnetic solvers
* `beforedeposition`: before the particle deposition (for charge and/or current)
* `afterdeposition`: after particle deposition (for charge and/or current)
* `beforestep`: before the time step
* `afterstep`: after the time step
* `afterdiagnostics`: after diagnostic output
* `oncheckpointsignal`: on a checkpoint signal
* `onbreaksignal`: on a break signal. These callbacks will be the last ones executed before the simulation ends.
* `particlescraper`: before particle boundary conditions are applied
* `particleloader`: at the time that the standard particle loader is called
* `particleinjection`: called when particle injection happens, after the position
  advance and before deposition is called, allowing a user
  defined particle distribution to be injected each time step

Example that calls the Python function `myplots` after each step:

```python3
from pywarpx.callbacks import installcallback

def myplots():
    # do something here

installcallback('afterstep', myplots)

# run simulation
sim.step(nsteps=100)
```

The install can also be done using a [Python decorator](https://docs.python.org/3/glossary.html#term-decorator), which has the prefix `callfrom`.
To use a decorator, the syntax is as follows. This will install the function `myplots` to be called after each step.
The above example is quivalent to the following:

```python3
from pywarpx.callbacks import callfromafterstep

@callfromafterstep
def myplots():
    # do something here

# run simulation
sim.step(nsteps=100)
```

### pywarpx.callbacks.installcallback(name, f)

Installs a function to be called at that specified time.

Adds a function to the list of functions called by this callback.

### pywarpx.callbacks.isinstalled(name, f)

Checks if a function is installed for this callback.

### pywarpx.callbacks.uninstallcallback(name, f)

Uninstalls the function (so it won’t be called anymore).

Removes the function from the list of functions called by this callback.
