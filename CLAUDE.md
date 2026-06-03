# CLAUDE.md

This file provides guidance to Claude Code when working with code in this repository.

## Overview

The code in this folder is meant to contain simulations for a Research Experience for Undergraduates (REU) program at the Cornell Center for Bright Beams (CBB) and Cornell Laboratory for Accelerator ScienceS and Education (CLASSE). The project focus is to build a beam simulation for the Cornell High Energy Synchotron Source (CHESS). 

Please use the CBB conda environment by running `conda activate CBB` (Miniforge is installed at ~/miniforge3; if conda isn't on your PATH yet, first run source ~/miniforge3/bin/activate).

## Reference Materials

Before coding, read any relevant documentation (listed below) or papers (listed in `reference/Papers/README.md`) to better understand the task.

### Simulation Codes

| Tool | Location | Purpose |
|------|----------|---------|
| **IMPACT-T** | `reference/Impact-T Documentation/README.md` | 3D relativistic particle tracking with space charge, wakefields, and CSR. Parallel implementation, used in photoinjector design. |
| **IMPACT-Z** | `reference/Impact-Z Documentation/README.md` | 3D parallel PIC code for intense beams through drifts, quadrupoles, solenoids, bending magnets, multipoles, and RF cavities. |
| **GPT** | `reference/GPT Documentation/README.md` | General Particle Tracer — 3D charged particle dynamics including space charge. Uses GDF file format for I/O. |
| **WarpX** | `reference/WarpX Documentation/README.md` | Massively parallel PIC code (EM and electrostatic). Supports GPU backends (CUDA/HIP/SYCL), adaptive mesh refinement, Python interface via `pywarpx`/PICMI. |
| **G4beamline** | `reference/G4beamline Documentation/README.md` | Geant4-based beamline simulation — command-driven input file, full physics lists, virtual detectors, NTuples, and 3D visualization. |
| **BMAD** | `reference/BMAD Documentation/README.md` | Fortran90 subroutine library for reading MAD-format lattice files, computing Twiss parameters, and tracking particles. Developed at Cornell (CESR/CLASSE). Supports Taylor maps, Runge-Kutta, symplectic integrators, and PTC interface. MAD-X User Manual at `reference/BMAD Documentation/MAD-X User Manual/README.md`. |
| **Linac Sim GUI** | `reference/Linac Simulation Documentation/README.md` | Adam Bartnik's CESR Linac simulation GUI (Java). Chains a custom 1D cathode code → GPT (space charge, cylindrical symmetry) → BMAD (high-energy, 3D). Includes fieldmaps for the thermionic gun, prebunchers, solenoid lenses, and SLAC-design linac cavities. |
| **LUME-Impact** | `reference/lume-impact Documentation/README.md` | Python interface for IMPACT-T and IMPACT-Z. Provides `Impact` and `ImpactZ` classes for input configuration, execution, output parsing, and plotting. Integrates with openPMD-beamphysics and BMAD. |

### Key Concepts

- **PIC (Particle-In-Cell)**: The computational method used by IMPACT-Z and WarpX — particles tracked on a mesh, fields solved on grid.
- **Space charge**: Coulomb self-repulsion of the beam, dominant at low energy and high current; all four codes model it.
- **CSR (Coherent Synchrotron Radiation)**: Wakefield from relativistic bunches in bending magnets; modeled in IMPACT-T.
- **GDF**: GPT's native binary data format; convert to/from ASCII with `GDF2A`/`ASCI2GDF`.
- **IMPACT-Z output files**: Named `fort.18`, `fort.24`–`fort.30`, `fort.32` — see `reference/Impact-Z Documentation/output_files/` for column definitions.

### Adding New Papers

When saving a new paper to `reference/Papers/`, add a summary entry to `reference/Papers/README.md` following the existing format (title, file, author, abstract summary).
