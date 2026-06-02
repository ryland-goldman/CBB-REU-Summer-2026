# IMPACT-Z User Manual

**Version 2.7** — Ji Qiang, Lawrence Berkeley National Laboratory

IMPACT-Z is a 3D parallel/serial Particle-In-Cell (PIC) code based on multi-layer object-oriented design. It can treat intense beams propagating through drifts, magnetic quadrupoles, magnetic solenoids, bending magnets, multipoles, and RF cavities using either a map integrator or nonlinear Lorentz integrator.

This markdown-based manual was automatically generated from a PDF version and may have issues with transcription. If there are formatting errors, please see [ImpactZusermanual.pdf](ImpactZusermanual.pdf).

---

## Introduction

| Section | File |
|---------|------|
| 1. Introduction | [introduction.md](introduction/introduction.md) |
| 2. Programs to prepare input files | [programs_to_prepare_input_files.md](introduction/programs_to_prepare_input_files.md) |

---

## Input Files

The main IMPACT-Z input file is called `ImpactZ.in`. If field data is needed, the file `rfdataN.in` is required (N is the field ID set in the element input information). For initial distribution, you can either choose distributions within the code (Gaussian, Waterbag, KV, etc.) or read from `particle.in`.

| Section | File |
|---------|------|
| 3. Input Files overview | [1.1_impactz_in.md](input_files/1.1_impactz_in.md) |
| 1.1. ImpactZ.in | [1.1_impactz_in.md](input_files/1.1_impactz_in.md) |
| 1.1.1. Beam Section | [1.1.1_beam_section.md](input_files/1.1.1_beam_section.md) |
| 1.1.2. Lattice Section | [1.1.2_lattice_section.md](input_files/1.1.2_lattice_section.md) |
| 1.2. rfdataN.in | [1.2_rfdatan_in.md](input_files/1.2_rfdatan_in.md) |
| 1.3. particle.in | [1.3_particle_in.md](input_files/1.3_particle_in.md) |

---

## Run IMPACT-Z

| Section | File |
|---------|------|
| 2.1. Run IMPACT-Z | [2.1_run_impactz.md](run_impactz/2.1_run_impactz.md) |
| 2.2. Console | [2.2_console.md](run_impactz/2.2_console.md) |
| 2.3. Stop IMPACT-Z | [2.3_stop_impactz.md](run_impactz/2.3_stop_impactz.md) |

---

## Output Files

There are several output files produced by IMPACT-Z: `fort.18`, `fort.24`, `fort.25`, `fort.26`, `fort.27`, `fort.28`, `fort.29`, `fort.30`, and `fort.32`.

| Section | File |
|---------|------|
| 5.1. fort.18: reference particle information | [5.1_fort18.md](output_files/5.1_fort18.md) |
| 5.2. fort.24/25/26: RMS size information | [5.2_fort24_fort25_fort26.md](output_files/5.2_fort24_fort25_fort26.md) |
| 5.3. fort.27: maximum amplitude information | [5.3_fort27.md](output_files/5.3_fort27.md) |
| 5.4. fort.28: load balance and loss diagnostic | [5.4_fort28.md](output_files/5.4_fort28.md) |
| 5.5. fort.29: cubic root of 3rd moments | [5.5_fort29.md](output_files/5.5_fort29.md) |
| 5.6. fort.30: square root of 4th moments | [5.6_fort30.md](output_files/5.6_fort30.md) |
| 5.7. fort.32: number of particles per charge state | [5.7_fort32.md](output_files/5.7_fort32.md) |
| 5.8. Distribution unit | [5.8_distribution_unit.md](output_files/5.8_distribution_unit.md) |
