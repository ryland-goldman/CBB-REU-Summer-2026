# MAD-X User's Reference Manual

**The MAD-X Program** (Methodical Accelerator Design) — Version 5.05.02

EUROPEAN LABORATORY FOR PARTICLE PHYSICS

*Laurent Deniau (editor), Hans Grote, Ghislain Roy, Frank Schmidt*

---

`MAD-X` is a general-purpose tool for charged-particle optics design and studies in alternating-gradient accelerators and beam lines. It can handle medium size to very large accelerators and solves various problems on such machines.

`MAD-X` is the successor of `MAD-8` and was specifically adapted to the needs of the design of the LHC. The `PTC` library of E. Forest is also embedded in `MAD-X` as an addition to better support small and low energy accelerators. A large part of the present document is based on the `MAD-8` documentation originally written and published by F.C. Iselin.

Comments and corrections may be sent to [`mad@cern.ch`](mailto:mad@cern.ch?subject=[user's guide]).

---

## Table of Contents

- **Chapter 1 — Control**
  - [1.1 Conventions](chapter_1_control/1.1_conventions.md)
  - [1.2 Command Format](chapter_1_control/1.2_command_format.md)
  - [1.3 Program Flow Statements](chapter_1_control/1.3_program_flow_statements.md)
  - [1.4 General Control Statements](chapter_1_control/1.4_general_control_statements.md)
  - [1.5 File Handling Statements](chapter_1_control/1.5_file_handling_statements.md)
  - [1.6 Table Handling Statements](chapter_1_control/1.6_table_handling_statements.md)
  - [1.7 Beam Handling Statements](chapter_1_control/1.7_beam_handling_statements.md)
  - [1.8 Sequence Editor](chapter_1_control/1.8_sequence_editor.md)
- **Chapter 2 — Elements, Beamlines and Sequences**
  - [2.1 Definition of Elements](chapter_2_elements_beamlines_sequences/2.1_definition_of_elements.md)
  - [2.2 Element Types](chapter_2_elements_beamlines_sequences/2.2_element_types.md)
  - [2.3 Range and Class Selection](chapter_2_elements_beamlines_sequences/2.3_range_and_class_selection.md)
  - [2.4 Beam Lines](chapter_2_elements_beamlines_sequences/2.4_beam_lines.md)
  - [2.5 Sequences](chapter_2_elements_beamlines_sequences/2.5_sequences.md)
- **Chapter 3 — Input and Output**
  - [3.1 TFS File Format](chapter_3_input_output/3.1_tfs_file_format.md)
  - [3.2 Conversion to SixTrack](chapter_3_input_output/3.2_conversion_to_sixtrack.md)
  - [3.3 SXF File Format](chapter_3_input_output/3.3_sxf_file_format.md)
  - [3.4 Plotting Data](chapter_3_input_output/3.4_plotting_data.md)
- **Chapter 4 — MAD-X Modules**
  - [4.1 SURVEY](chapter_4_modules/4.1_survey.md)
  - [4.2 Twiss Module](chapter_4_modules/4.2_twiss.md)
  - [4.3 Matching Module](chapter_4_modules/4.3_matching.md)
  - [4.4 EMIT: Equilibrium Emittances](chapter_4_modules/4.4_emit_equilibrium_emittances.md)
  - [4.5 Physical Aperture](chapter_4_modules/4.5_physical_aperture.md)
  - [4.6 Slicing a Sequence into Thin Lenses](chapter_4_modules/4.6_slicing_thin_lenses.md)
  - [4.7 Error Definitions](chapter_4_modules/4.7_error_definitions.md)
  - [4.8 Orbit Correction](chapter_4_modules/4.8_orbit_correction.md)
  - [4.9 SODD: Second Order Detuning and Distortion](chapter_4_modules/4.9_sodd.md)
  - [4.10 Touschek Lifetime and Scattering Rates](chapter_4_modules/4.10_touschek_lifetime.md)
  - [4.11 Intra-Beam Scattering](chapter_4_modules/4.11_intra_beam_scattering.md)
  - [4.12 Particle Tracking](chapter_4_modules/4.12_particle_tracking.md)
- **Chapter 5 — PTC Commands**
  - [5.1 PTC Set-up Parameters](chapter_5_ptc_commands/5.1_setup_parameters.md)
  - [5.2 Thick-Lens Tracking Module](chapter_5_ptc_commands/5.2_thick_lens_tracking.md)
  - [5.3 Ripken Optics Parameters](chapter_5_ptc_commands/5.3_ripken_optics_parameters.md)
  - [5.4 Non-Linear Machine Parameters](chapter_5_ptc_commands/5.4_non_linear_machine_parameters.md)
  - [5.5 MAD-X–PTC Auxiliaries](chapter_5_ptc_commands/5.5_ptc_auxiliaries.md)
- **Chapter 6 — Trailing Material**
  - [6.1 Known Differences to Other Programs](chapter_6_trailing_material/6.1_known_differences.md)
  - [6.2 MAD-X Recipes and Pitfalls](chapter_6_trailing_material/6.2_recipes_and_pitfalls.md)
  - [6.3 Contributors to MAD-X](chapter_6_trailing_material/6.3_contributors.md)
  - [6.4 Change Log](chapter_6_trailing_material/6.4_change_log.md)
