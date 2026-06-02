## 2. Programs to Prepare Input Files

The following programs help prepare input files for IMPACT-Z simulations:

**`Engscan.py`**: does a single cavity phase scan.

**`RFcoef.f90`**: prepares Fourier expansion coefficients of RF or solenoid field to be used in the simulation (with Lorentz nonlinear integrator).

**`phaseOptZ.py`**: sets up the RF driven phase of the cavities with the user specified design phases.

---

### Notes

1. The current version of the code is for serial single processor computer with Fortran90 compiler. To run the code on a parallel computer with MPI, the user has to comment out the line `"use mpistub"` in `Contrl/Input.f90`, `DataStruct/Data.f90`, `DataStruct/Pgrid.f90`, `DataStruct/PhysConst.f90`, and `Func/Timer.f90`. The user also has to remove the `mpif.h` file under the `Appl`, `Control`, `DataStruct`, and `Func` directories. The user also has to modify the Makefile to remove the `mpistub.o` inside the file and to use the appropriate parallel Fortran90 compiler such as `mpif90`.

2. The subroutines in `FFT.f90`: `realft`, `four1`, and `sinft`, can be replaced with functions from the Numerical Recipe or some equivalent 1D FFT functions.

3. To compile the code on Windows PC without `"make"` function, one needs to move all `*.f90` files to one directory and use a compiler (`g95` or `gfortran`) to compile the code in a single line.

4. There is a `ImpactTv2.pdf` document that contains a more detailed physics description.
