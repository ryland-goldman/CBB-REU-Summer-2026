# 4 Pre Processing

There are three preprocessing programs.

## RFcoefext

The program **RFcoefext** takes discrete field data as provided on axis in the file `rfdata.in`, extends these data by symmetry to make a new set of data describing the fields on an extended domain, and generates the Fourier coefficients for these extended fields based on a Fourier series expansion. These are stored in the file `rfdatax` for IMPACT-T simulation use. (Here, "x" in `rfdatax` has to be replaced by the corresponding "number" in the `ImpactT.in` input line for that element.) This file is used when only half of the injector is given. It contains both the real part and the imaginary part of the complex Fourier coefficients.

> *Note: the user should use another output file `rfdata.out` to check whether the input Fourier mode number is sufficient or not by comparing the plot of col. 2 vs. col. 1 with those from `rfdata.in`. Basically, `rfdata.out` contains the reconstructed and extended RF field and its derivatives using the Fourier coefficients.*

## RFcoefcls

The program **RFcoefcls** takes the discrete field data, generates the Fourier coefficients of these data directly and stores them in `rfdataxxx`. It also generates a shifted field data in `rfdata.tmp` based on these coefficients. This file will be used to regenerate Fourier coefficients for the shifted field data, which will be used to model a traveling wave field using the summation of two standing wave fields (one field is the shift of the other field). The `rfdataxx` will be used in beam line elements such as SolRF. The output file `rfdata.out` is the field reconstructed using those Fourier coefficients. This will be used as a good check of the accuracy of the Fourier expansion approximation.

## ImpactTphase.py

A Python scripting program, **ImpactTphase.py**, scans through the initial driven phase of a single RF cavity in the `ImpactT.in` file. It can also be used to scan through other parameters. The user needs to go into the code to specify the exact location of the parameter inside the file. Some other small changes (e.g. name of executable file) might also be needed. For the phase scan, the output of this scripting file is `Engout` with column one the driven phase, column five the final energy.

## PhaseOpt.py

Another python scripting program, **PhaseOpt.py**, can be used to automatically find the driven phase for all RF cavities in the `ImpactT.in` file. Here, the user needs to specify the design phase of each cavity in the `ImpactT.in` file. Minor modification of the program is needed to specify the line range of the lattice inside the `ImpactT.in` file. The single on-axis particle can be specified using the waterbag distribution with 0 rms sizes and transverse offsets except that the longitudinal offset corresponds to `-v0 * t_laser / 2`, where `v0` can be found from the initial beam kinetic energy "Bkenergy", `t_laser` is the total laser pulse length (in seconds).

## chicaneImpt.f90

There is a F90 program, **chicaneImpt.f90**, which will help generate four input files that will be used in the Impact-T simulation of a chicane. The user is encouraged to read this small program and to set the inputs appropriately.
