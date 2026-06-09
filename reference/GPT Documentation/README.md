# GPT User Manual

The General Particle Tracer (GPT) code is a well-established simulation platform for the study of charged particle dynamics in electromagnetic fields. The code is completely 3D, including the space-charge model. GPT can be conveniently customized without compromising its ease of use, accuracy or simulation speed, because of its modern implementation.

The first chapter of this manual describes the GPT kernel in detail. The second chapter is a tutorial introduction to GPT where simplified but practical examples are discussed. Although we recommend reading the first chapter first, some may prefer to start directly with the tutorial. Chapter three gives an overview of all pre- and postprocessing commands, many of which are also covered in the tutorial. Chapter four describes all the built-in GPT elements and keywords that can be used in the inputfile.

Customizing GPT is the subject of a separate manual, the GPT Programmer's Reference. This markdown-based manual was automatically generated from a PDF version and may have issues with transcription. If there are formatting errors, please see [UserManual.pdf](UserManual.pdf).

> **Running GPT in this repo:** the `gpt` executable is **not installed locally**. To run a GPT deck (e.g. to regenerate a `fieldmaps/*.gdf` map), use the `gpt_remote` wrapper (`/usr/local/bin/gpt_remote`), which mirrors the working directory to a remote host over SSH, runs the real `gpt` there with your exact arguments, and syncs results back. It is a drop-in for the `gpt` CLI: `gpt_remote -o out.gdf input.in`. See the root [`README.md`](../../README.md) ("Running GPT (remote wrapper)") for details and the `GPT_REMOTE` / `GPT_REMOTE_BIN` / `GPT_SSH_TIMEOUT` overrides.

---

## Preface

Particle tracing is a very powerful tool in the design of charged-particle accelerators and beam lines. However, the computers used to be too slow to trace the number of particles needed to obtain reliable statistics in a reasonable amount of computing time. Therefore the particle tracing method was often abandoned and especially matrix/optical methods became popular. These alternative methods however are limited because they fail to simulate an actual beam. Most of these methods yield inaccurate results or cannot be applied at all when space-charge effects are dominant and/or paraxial approximations can not be used.

Currently it is possible to trace millions of particles through complex electromagnetic fields. Because of the enormous advances in computer technology all 3D effects and space-charge forces can be taken into account. This offers much more accurate results compared to matrix and paraxial codes. The General Particle Tracer (GPT) code offers all of these capabilities in a user-friendly and easy-to-customize package.

In order to be able to read this manual, basic knowledge of electromagnetism theory and the special theory of relativity is required. Experience with accelerator or beam-line design is useful, but is not necessary. For information about how to write custom GPT code, we refer to the GPT Custom Elements documentation.

---

## Introduction

The GPT package consists of the GPT executable, which performs the actual calculations, and an extensive set of pre- and post-processing tools including data analysis tools and a graphical user interface.

| Section | File |
|---------|------|
| Overview | [overview.md](introduction/overview.md) |
| The GPT executable | [the_gpt_executable.md](introduction/the_gpt_executable.md) |
| Pre- and post-processing | [pre_and_post_processing.md](introduction/pre_and_post_processing.md) |

---

## Chapter 1 — The GPT Code

This chapter describes the GPT kernel in detail: the inputfile format, how to run GPT, coordinate systems, space-charge models, initial particle distributions, the equations of motion, the Runge-Kutta integrator, output mechanisms, the GDF file format, collector design, and automatic solving.

| Section | File |
|---------|------|
| 1.1 The GPT inputfile | [1.1_the_gpt_inputfile.md](chapter_1_the_gpt_code/1.1_the_gpt_inputfile.md) |
| 1.2 Running GPT | [1.2_running_gpt.md](chapter_1_the_gpt_code/1.2_running_gpt.md) |
| 1.3 Error messages | [1.3_error_messages.md](chapter_1_the_gpt_code/1.3_error_messages.md) |
| 1.4 Coordinate Systems | [1.4_coordinate_systems.md](chapter_1_the_gpt_code/1.4_coordinate_systems.md) |
| 1.5 Space charge | [1.5_space_charge.md](chapter_1_the_gpt_code/1.5_space_charge.md) |
| 1.6 Initial particle distribution | [1.6_initial_particle_distribution.md](chapter_1_the_gpt_code/1.6_initial_particle_distribution.md) |
| 1.7 Equations of motion | [1.7_equations_of_motion.md](chapter_1_the_gpt_code/1.7_equations_of_motion.md) |
| 1.8 Runge-Kutta | [1.8_runge_kutta.md](chapter_1_the_gpt_code/1.8_runge_kutta.md) |
| 1.9 Output | [1.9_output.md](chapter_1_the_gpt_code/1.9_output.md) |
| 1.10 GDF | [1.10_gdf.md](chapter_1_the_gpt_code/1.10_gdf.md) |
| 1.11 Collector design | [1.11_collector_design.md](chapter_1_the_gpt_code/1.11_collector_design.md) |
| 1.12 Automatic solving | [1.12_automatic_solving.md](chapter_1_the_gpt_code/1.12_automatic_solving.md) |

---

## Chapter 2 — Tutorial

This chapter gives a tutorial introduction to the GPT package. All the files in the tutorial can be found in the GPT distribution in the `tutorial` directory. Within that directory, a subdirectory or folder exists for every section in the tutorial containing the described files.

| Section | File |
|---------|------|
| 2.1 Quadrupole focusing | [2.1_quadrupole_focusing.md](chapter_2_tutorial/2.1_quadrupole_focusing.md) |
| 2.2 The GPTwin user interface | [2.2_the_gptwin_user_interface.md](chapter_2_tutorial/2.2_the_gptwin_user_interface.md) |
| 2.3 Scanning the beam energy | [2.3_scanning_the_beam_energy.md](chapter_2_tutorial/2.3_scanning_the_beam_energy.md) |
| 2.4 Automatic solving | [2.4_automatic_solving.md](chapter_2_tutorial/2.4_automatic_solving.md) |
| 2.5 Electrostatic accelerator | [2.5_electrostatic_accelerator.md](chapter_2_tutorial/2.5_electrostatic_accelerator.md) |
| 2.6 Accuracy | [2.6_accuracy.md](chapter_2_tutorial/2.6_accuracy.md) |
| 2.7 Magnetic mirror | [2.7_magnetic_mirror.md](chapter_2_tutorial/2.7_magnetic_mirror.md) |
| 2.8 Element Coordinate System | [2.8_element_coordinate_system.md](chapter_2_tutorial/2.8_element_coordinate_system.md) |
| 2.9 Initial particle distribution | [2.9_initial_particle_distribution.md](chapter_2_tutorial/2.9_initial_particle_distribution.md) |
| 2.10 Photo-cathode, starting particles as function of time | [2.10_photo_cathode.md](chapter_2_tutorial/2.10_photo_cathode.md) |
| 2.11 Collector design | [2.11_collector_design.md](chapter_2_tutorial/2.11_collector_design.md) |

---

## Chapter 3 — GDF

GPT Datafile Format, GDF, is a file format developed for programs producing large amounts of numerical data. It is the native output format of the General Particle Tracer, GPT, allowing it to automatically scan parameters and perform efficient data processing and visualization. GDF-files are machine independent and can be transferred between PC's and Linux machines.

This chapter describes the basics of a GDF file and the following utility programs:

| Program | Description |
|---------|-------------|
| ASCI2GDF | Converts an ASCII file to a GDF file. |
| ASTRA2GDF | Convert output from the ASTRA program to input for GPT. |
| FISH2GDF | Converts SUPERFISH output to a GDF file. |
| FISHFILE | Calculates "unrolled" scatter statistics or creates a GPT inputfile based on a Superfish file. |
| GDF2A | Writes the contents of a GDF file in ASCII to the terminal or a file. The file can be processed by an ASCII-based editor or imported into a spreadsheet program. |
| GDFFLATTEN | Removes one level from the hierarchical structure of a GDF file by combining the results. |
| GDF2DXF | Converts a GDF file to a point, line or vector drawing for further customization in a 2D or 3D drawing package. |
| GDF2GDF | Combines several GDF files into a single file, optionally sorting the top-level group. |
| GDF2HIS | Calculates histograms from the contents of a GDF file. The output is in GDF format and can be plotted or converted to ASCII. |
| GDF2SDDS | Converts a GDF file to an SDDS file. |
| GDF2VTP | Convert a GDF file to ParaView format. |
| GDFA | Removes one level of hierarchy in a GDF file by performing user-specified calculations. The calculations typically produce averages and standard deviations of arrays, reducing the total amount of data considerably. |
| GDFMGO | Multi-objective genetic optimizer for GPT. |
| GDFSELECT | Select a specific or a range of groups from a GDF file. |
| GDFSOLVE | Multi-dimensional root finder and optimizer for GPT. |
| GDFTRANS | Transpose a GDF file to calculate trajectories as function of any parameter. |
| MR | Scan one or more parameters and concatenate the output in a single hierarchical file. |
| RAW2GDF | Converts a raw binary datafile to a GDF file. |

| Section | File |
|---------|------|
| 3.1 Basics | [3.1_basics.md](chapter_3_gdf/3.1_basics.md) |
| 3.2 ASCI2GDF | [3.2_asci2gdf.md](chapter_3_gdf/3.2_asci2gdf.md) |
| 3.3 ASTRA2GDF | [3.3_astra2gdf.md](chapter_3_gdf/3.3_astra2gdf.md) |
| 3.4 FISH2GDF | [3.4_fish2gdf.md](chapter_3_gdf/3.4_fish2gdf.md) |
| 3.5 FISHFILE | [3.5_fishfile.md](chapter_3_gdf/3.5_fishfile.md) |
| 3.6 GDF2A | [3.6_gdf2a.md](chapter_3_gdf/3.6_gdf2a.md) |
| 3.7 GDF2DXF | [3.7_gdf2dxf.md](chapter_3_gdf/3.7_gdf2dxf.md) |
| 3.8 GDF2GDF | [3.8_gdf2gdf.md](chapter_3_gdf/3.8_gdf2gdf.md) |
| 3.9 GDF2HIS | [3.9_gdf2his.md](chapter_3_gdf/3.9_gdf2his.md) |
| 3.10 GDF2SDDS | [3.10_gdf2sdds.md](chapter_3_gdf/3.10_gdf2sdds.md) |
| 3.11 GDF2VTP | [3.11_gdf2vtp.md](chapter_3_gdf/3.11_gdf2vtp.md) |
| 3.12 GDFA | [3.12_gdfa.md](chapter_3_gdf/3.12_gdfa.md) |
| 3.13 GDFFLATTEN | [3.13_gdfflatten.md](chapter_3_gdf/3.13_gdfflatten.md) |
| 3.14 GDFMGO | [3.14_gdfmgo.md](chapter_3_gdf/3.14_gdfmgo.md) |
| 3.15 GDFSELECT | [3.15_gdfselect.md](chapter_3_gdf/3.15_gdfselect.md) |
| 3.16 GDFSOLVE | [3.16_gdfsolve.md](chapter_3_gdf/3.16_gdfsolve.md) |
| 3.17 GDFTRANS | [3.17_gdftrans.md](chapter_3_gdf/3.17_gdftrans.md) |
| 3.18 MR | [3.18_mr.md](chapter_3_gdf/3.18_mr.md) |
| 3.19 RAW2GDF | [3.19_raw2gdf.md](chapter_3_gdf/3.19_raw2gdf.md) |

---

## Chapter 4 — GPT Reference

This chapter serves as a reference of the simulation variables, keywords and elements. It is not recommended for first reading. The inputfile syntax, operator precedence and predefined constants are provided in section 4.1.

If necessary, you can also write your own element or modify an existing one. How this is done is explained in the GPT Programmer's Reference.

| Section | File |
|---------|------|
| 4.1 Inputfile syntax | [4.1_inputfile_syntax.md](chapter_4_gpt_reference/4.1_inputfile_syntax.md) |
| 4.2 Initial particle distribution | [4.2_initial_particle_distribution.md](chapter_4_gpt_reference/4.2_initial_particle_distribution.md) |
| 4.3 Timestep and output control | [4.3_timestep_and_output_control.md](chapter_4_gpt_reference/4.3_timestep_and_output_control.md) |
| 4.4 Accelerating structures | [4.4_accelerating_structures.md](chapter_4_gpt_reference/4.4_accelerating_structures.md) |
| 4.5 Static electric fields | [4.5_static_electric_fields.md](chapter_4_gpt_reference/4.5_static_electric_fields.md) |
| 4.6 Static magnetic fields | [4.6_static_magnetic_fields.md](chapter_4_gpt_reference/4.6_static_magnetic_fields.md) |
| 4.7 Field maps | [4.7_field_maps.md](chapter_4_gpt_reference/4.7_field_maps.md) |
| 4.8 Spacecharge | [4.8_spacecharge.md](chapter_4_gpt_reference/4.8_spacecharge.md) |
| 4.9 Scattering | [4.9_scattering.md](chapter_4_gpt_reference/4.9_scattering.md) |
| 4.10 Remove particles | [4.10_remove_particles.md](chapter_4_gpt_reference/4.10_remove_particles.md) |
| 4.11 Miscellaneous | [4.11_miscellaneous.md](chapter_4_gpt_reference/4.11_miscellaneous.md) |
| 4.12 Free Electron Laser (FEL), Wakefields, and CSR | [4.12_fel_wakefields_csr.md](chapter_4_gpt_reference/4.12_fel_wakefields_csr.md) |
| 4.13 Obsolete elements | [4.13_obsolete_elements.md](chapter_4_gpt_reference/4.13_obsolete_elements.md) |
