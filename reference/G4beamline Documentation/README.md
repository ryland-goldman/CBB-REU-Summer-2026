# G4beamline User's Guide

G4beamline is a Geant4-based program that simulates the passage of particles through beamlines and other systems. It uses a command-based input file to describe beamline geometry, fields, and physics processes.

**Version:** 3.08  
**Author:** Tom Roberts, Muons, Inc.  
**Date:** August 2022

This markdown-based documentation was automatically converted from [G4beamlineUsersGuide-2.pdf](G4beamlineUsersGuide-2.pdf). If there are formatting errors, please refer to the original PDF.

---

## Chapter 1 — Introduction

| Section | File |
|---------|------|
| 1 Introduction | [introduction.md](chapter_1_introduction/introduction.md) |
| 1.1 License and Distribution | [1.1_license_and_distribution.md](chapter_1_introduction/1.1_license_and_distribution.md) |

---

## Chapter 2 — Quick Start Guide

| Section | File |
|---------|------|
| 2.1 Installation | [2.1_installation.md](chapter_2_quick_start/2.1_installation.md) |
| 2.2 Initial Test | [2.2_initial_test.md](chapter_2_quick_start/2.2_initial_test.md) |
| 2.3 Using the G4beamline Command Line | [2.3_command_line.md](chapter_2_quick_start/2.3_command_line.md) |
| 2.4 Companion Programs for Plotting | [2.4_companion_programs.md](chapter_2_quick_start/2.4_companion_programs.md) |

---

## Chapter 3 — Basic Concepts

| Section | File |
|---------|------|
| 3.1 Particle tracking in simulations | [3.1_particle_tracking.md](chapter_3_basic_concepts/3.1_particle_tracking.md) |
| 3.2 Units | [3.2_units.md](chapter_3_basic_concepts/3.2_units.md) |
| 3.3 Geometry | [3.3_geometry.md](chapter_3_basic_concepts/3.3_geometry.md) |
| 3.4 Coordinates | [3.4_coordinates.md](chapter_3_basic_concepts/3.4_coordinates.md) |
| 3.5 Rotations | [3.5_rotations.md](chapter_3_basic_concepts/3.5_rotations.md) |
| 3.6 Materials | [3.6_materials.md](chapter_3_basic_concepts/3.6_materials.md) |
| 3.7 Electromagnetic Fields | [3.7_electromagnetic_fields.md](chapter_3_basic_concepts/3.7_electromagnetic_fields.md) |
| 3.8 Visualization | [3.8_visualization.md](chapter_3_basic_concepts/3.8_visualization.md) |
| 3.9 Obtaining Results – Virtual Detectors and NTuples | [3.9_virtual_detectors.md](chapter_3_basic_concepts/3.9_virtual_detectors.md) |
| 3.10 Random Number Generator | [3.10_random_number_generator.md](chapter_3_basic_concepts/3.10_random_number_generator.md) |
| 3.11 Tuning the Beamline | [3.11_tuning_beamline.md](chapter_3_basic_concepts/3.11_tuning_beamline.md) |
| 3.12 The Geant4 User Interface | [3.12_geant4_user_interface.md](chapter_3_basic_concepts/3.12_geant4_user_interface.md) |
| 3.13 Event Time Limit | [3.13_event_time_limit.md](chapter_3_basic_concepts/3.13_event_time_limit.md) |
| 3.14 Wall Clock Limit | [3.14_wall_clock_limit.md](chapter_3_basic_concepts/3.14_wall_clock_limit.md) |
| 3.15 Event and Track Numbering (EventID, TrackID) | [3.15_event_track_numbering.md](chapter_3_basic_concepts/3.15_event_track_numbering.md) |
| 3.16 Signals | [3.16_signals.md](chapter_3_basic_concepts/3.16_signals.md) |

---

## Chapter 4 — Important Values that Affect the Validity and Accuracy of Simulations

| Section | File |
|---------|------|
| 4.1 Physics List | [4.1_physics_list.md](chapter_4_important_values/4.1_physics_list.md) |
| 4.2 Tracking Accuracy Parameters | [4.2_tracking_accuracy.md](chapter_4_important_values/4.2_tracking_accuracy.md) |
| 4.3 Electromagnetic Field Map Tolerance | [4.3_field_map_tolerance.md](chapter_4_important_values/4.3_field_map_tolerance.md) |
| 4.4 Secondary Creation Threshold in Physics Processes | [4.4_secondary_creation_threshold.md](chapter_4_important_values/4.4_secondary_creation_threshold.md) |
| 4.5 Track Cuts | [4.5_track_cuts.md](chapter_4_important_values/4.5_track_cuts.md) |
| 4.6 Fringe Fields | [4.6_fringe_fields.md](chapter_4_important_values/4.6_fringe_fields.md) |
| 4.7 Stepping and Hard Edges | [4.7_stepping_hard_edges.md](chapter_4_important_values/4.7_stepping_hard_edges.md) |
| 4.8 Radiative Decays | [4.8_radiative_decays.md](chapter_4_important_values/4.8_radiative_decays.md) |
| 4.9 Miscellaneous | [4.9_miscellaneous.md](chapter_4_important_values/4.9_miscellaneous.md) |

---

## Chapter 5 — Input File Description

| Section | File |
|---------|------|
| 5.1 Expressions | [5.1_expressions.md](chapter_5_input_file/5.1_expressions.md) |
| 5.2 Element Names | [5.2_element_names.md](chapter_5_input_file/5.2_element_names.md) |
| 5.3 Parameters | [5.3_parameters.md](chapter_5_input_file/5.3_parameters.md) |
| 5.4 Pillbox Geometry and Dimensions | [5.4_pillbox_geometry.md](chapter_5_input_file/5.4_pillbox_geometry.md) |
| 5.5 Absorber Geometry and Dimensions | [5.5_absorber_geometry.md](chapter_5_input_file/5.5_absorber_geometry.md) |
| 5.6 Cornerarc Geometry | [5.6_cornerarc_geometry.md](chapter_5_input_file/5.6_cornerarc_geometry.md) |
| 5.7 genericquad Aperture | [5.7_genericquad_aperture.md](chapter_5_input_file/5.7_genericquad_aperture.md) |
| 5.8 G4beamline Commands by Type | [5.8_commands_by_type.md](chapter_5_input_file/5.8_commands_by_type.md) |

---

## Chapter 6 — G4beamline Commands (Alphabetical)

| Section | File |
|---------|------|
| 6.1 absorber | [6.1_absorber.md](chapter_6_commands/6.1_absorber.md) |
| 6.2 beam | [6.2_beam.md](chapter_6_commands/6.2_beam.md) |
| 6.3 beamlossntuple | [6.3_beamlossntuple.md](chapter_6_commands/6.3_beamlossntuple.md) |
| 6.4 boolean | [6.4_boolean.md](chapter_6_commands/6.4_boolean.md) |
| 6.5 box | [6.5_box.md](chapter_6_commands/6.5_box.md) |
| 6.6 bug1021 | [6.6_bug1021.md](chapter_6_commands/6.6_bug1021.md) |
| 6.7 coil | [6.7_coil.md](chapter_6_commands/6.7_coil.md) |
| 6.8 collective | [6.8_collective.md](chapter_6_commands/6.8_collective.md) |
| 6.9 corner | [6.9_corner.md](chapter_6_commands/6.9_corner.md) |
| 6.10 cornerarc | [6.10_cornerarc.md](chapter_6_commands/6.10_cornerarc.md) |
| 6.11 cosmicraybeam | [6.11_cosmicraybeam.md](chapter_6_commands/6.11_cosmicraybeam.md) |
| 6.12 cylinder | [6.12_cylinder.md](chapter_6_commands/6.12_cylinder.md) |
| 6.13 define | [6.13_define.md](chapter_6_commands/6.13_define.md) |
| 6.14 demo | [6.14_demo.md](chapter_6_commands/6.14_demo.md) |
| 6.15 detector | [6.15_detector.md](chapter_6_commands/6.15_detector.md) |
| 6.16 do | [6.16_do.md](chapter_6_commands/6.16_do.md) |
| 6.17 elementdb | [6.17_elementdb.md](chapter_6_commands/6.17_elementdb.md) |
| 6.18 endgroup | [6.18_endgroup.md](chapter_6_commands/6.18_endgroup.md) |
| 6.19 eventcuts | [6.19_eventcuts.md](chapter_6_commands/6.19_eventcuts.md) |
| 6.20 exit | [6.20_exit.md](chapter_6_commands/6.20_exit.md) |
| 6.21 expandworld | [6.21_expandworld.md](chapter_6_commands/6.21_expandworld.md) |
| 6.22 extrusion | [6.22_extrusion.md](chapter_6_commands/6.22_extrusion.md) |
| 6.23 fieldexpr | [6.23_fieldexpr.md](chapter_6_commands/6.23_fieldexpr.md) |
| 6.24 fieldlines | [6.24_fieldlines.md](chapter_6_commands/6.24_fieldlines.md) |
| 6.25 fieldmap | [6.25_fieldmap.md](chapter_6_commands/6.25_fieldmap.md) |
| 6.26 fieldntuple | [6.26_fieldntuple.md](chapter_6_commands/6.26_fieldntuple.md) |
| 6.27 for | [6.27_for.md](chapter_6_commands/6.27_for.md) |
| 6.28 g4ui | [6.28_g4ui.md](chapter_6_commands/6.28_g4ui.md) |
| 6.29 genericbend | [6.29_genericbend.md](chapter_6_commands/6.29_genericbend.md) |
| 6.30 genericquad | [6.30_genericquad.md](chapter_6_commands/6.30_genericquad.md) |
| 6.31 genericsectorbend | [6.31_genericsectorbend.md](chapter_6_commands/6.31_genericsectorbend.md) |
| 6.32 geometry | [6.32_geometry.md](chapter_6_commands/6.32_geometry.md) |
| 6.33 group | [6.33_group.md](chapter_6_commands/6.33_group.md) |
| 6.34 helicaldipole | [6.34_helicaldipole.md](chapter_6_commands/6.34_helicaldipole.md) |
| 6.35 helicalharmonic | [6.35_helicalharmonic.md](chapter_6_commands/6.35_helicalharmonic.md) |
| 6.36 help | [6.36_help.md](chapter_6_commands/6.36_help.md) |
| 6.37 idealsectorbend | [6.37_idealsectorbend.md](chapter_6_commands/6.37_idealsectorbend.md) |
| 6.38 if | [6.38_if.md](chapter_6_commands/6.38_if.md) |
| 6.39 include | [6.39_include.md](chapter_6_commands/6.39_include.md) |
| 6.40 isotropicsource | [6.40_isotropicsource.md](chapter_6_commands/6.40_isotropicsource.md) |
| 6.41 label | [6.41_label.md](chapter_6_commands/6.41_label.md) |
| 6.42 lilens | [6.42_lilens.md](chapter_6_commands/6.42_lilens.md) |
| 6.43 list | [6.43_list.md](chapter_6_commands/6.43_list.md) |
| 6.44 man | [6.44_man.md](chapter_6_commands/6.44_man.md) |
| 6.45 material | [6.45_material.md](chapter_6_commands/6.45_material.md) |
| 6.46 movie | [6.46_movie.md](chapter_6_commands/6.46_movie.md) |
| 6.47 multipole | [6.47_multipole.md](chapter_6_commands/6.47_multipole.md) |
| 6.48 muminuscapturefix | [6.48_muminuscapturefix.md](chapter_6_commands/6.48_muminuscapturefix.md) |
| 6.49 muonium | [6.49_muonium.md](chapter_6_commands/6.49_muonium.md) |
| 6.50 newparticlentuple | [6.50_newparticlentuple.md](chapter_6_commands/6.50_newparticlentuple.md) |
| 6.51 ntuple | [6.51_ntuple.md](chapter_6_commands/6.51_ntuple.md) |
| 6.52 output | [6.52_output.md](chapter_6_commands/6.52_output.md) |
| 6.53 param | [6.53_param.md](chapter_6_commands/6.53_param.md) |
| 6.54 particlecolor | [6.54_particlecolor.md](chapter_6_commands/6.54_particlecolor.md) |
| 6.55 particlefilter | [6.55_particlefilter.md](chapter_6_commands/6.55_particlefilter.md) |
| 6.56 particlesource | [6.56_particlesource.md](chapter_6_commands/6.56_particlesource.md) |
| 6.57 physics | [6.57_physics.md](chapter_6_commands/6.57_physics.md) |
| 6.58 pillbox | [6.58_pillbox.md](chapter_6_commands/6.58_pillbox.md) |
| 6.59 place | [6.59_place.md](chapter_6_commands/6.59_place.md) |
| 6.60 polycone | [6.60_polycone.md](chapter_6_commands/6.60_polycone.md) |
| 6.61 printf | [6.61_printf.md](chapter_6_commands/6.61_printf.md) |
| 6.62 printfield | [6.62_printfield.md](chapter_6_commands/6.62_printfield.md) |
| 6.63 probefield | [6.63_probefield.md](chapter_6_commands/6.63_probefield.md) |
| 6.64 profile | [6.64_profile.md](chapter_6_commands/6.64_profile.md) |
| 6.65 randomseed | [6.65_randomseed.md](chapter_6_commands/6.65_randomseed.md) |
| 6.66 rdecaysource | [6.66_rdecaysource.md](chapter_6_commands/6.66_rdecaysource.md) |
| 6.67 reference | [6.67_reference.md](chapter_6_commands/6.67_reference.md) |
| 6.68 rfdevice | [6.68_rfdevice.md](chapter_6_commands/6.68_rfdevice.md) |
| 6.69 sample | [6.69_sample.md](chapter_6_commands/6.69_sample.md) |
| 6.70 setdecay | [6.70_setdecay.md](chapter_6_commands/6.70_setdecay.md) |
| 6.71 showmaterial | [6.71_showmaterial.md](chapter_6_commands/6.71_showmaterial.md) |
| 6.72 solenoid | [6.72_solenoid.md](chapter_6_commands/6.72_solenoid.md) |
| 6.73 sourceonly | [6.73_sourceonly.md](chapter_6_commands/6.73_sourceonly.md) |
| 6.74 spacecharge | [6.74_spacecharge.md](chapter_6_commands/6.74_spacecharge.md) |
| 6.75 spacechargelw | [6.75_spacechargelw.md](chapter_6_commands/6.75_spacechargelw.md) |
| 6.76 sphere | [6.76_sphere.md](chapter_6_commands/6.76_sphere.md) |
| 6.77 start | [6.77_start.md](chapter_6_commands/6.77_start.md) |
| 6.78 survey | [6.78_survey.md](chapter_6_commands/6.78_survey.md) |
| 6.79 tess | [6.79_tess.md](chapter_6_commands/6.79_tess.md) |
| 6.80 tessellatedsolid | [6.80_tessellatedsolid.md](chapter_6_commands/6.80_tessellatedsolid.md) |
| 6.81 test | [6.81_test.md](chapter_6_commands/6.81_test.md) |
| 6.82 timentuple | [6.82_timentuple.md](chapter_6_commands/6.82_timentuple.md) |
| 6.83 torus | [6.83_torus.md](chapter_6_commands/6.83_torus.md) |
| 6.84 totalenergy | [6.84_totalenergy.md](chapter_6_commands/6.84_totalenergy.md) |
| 6.85 trace | [6.85_trace.md](chapter_6_commands/6.85_trace.md) |
| 6.86 trackcolor | [6.86_trackcolor.md](chapter_6_commands/6.86_trackcolor.md) |
| 6.87 trackcuts | [6.87_trackcuts.md](chapter_6_commands/6.87_trackcuts.md) |
| 6.88 tracker | [6.88_tracker.md](chapter_6_commands/6.88_tracker.md) |
| 6.89 trackermode | [6.89_trackermode.md](chapter_6_commands/6.89_trackermode.md) |
| 6.90 trackerplane | [6.90_trackerplane.md](chapter_6_commands/6.90_trackerplane.md) |
| 6.91 trap | [6.91_trap.md](chapter_6_commands/6.91_trap.md) |
| 6.92 tube | [6.92_tube.md](chapter_6_commands/6.92_tube.md) |
| 6.93 tubs | [6.93_tubs.md](chapter_6_commands/6.93_tubs.md) |
| 6.94 tune | [6.94_tune.md](chapter_6_commands/6.94_tune.md) |
| 6.95 usertrackfilter | [6.95_usertrackfilter.md](chapter_6_commands/6.95_usertrackfilter.md) |
| 6.96 virtualdetector | [6.96_virtualdetector.md](chapter_6_commands/6.96_virtualdetector.md) |
| 6.97 zntuple | [6.97_zntuple.md](chapter_6_commands/6.97_zntuple.md) |

---

## Chapter 7 — Examples

| Section | File |
|---------|------|
| 7.1 Example1 – Simple Tracking and Virtualdetectors | [7.1_example1.md](chapter_7_examples/7.1_example1.md) |
| 7.2 ExampleN02.g4bl – The Geant4 ExampleN02 | [7.2_exampleN02.md](chapter_7_examples/7.2_exampleN02.md) |
| 7.3 FieldLines.g4bl | [7.3_fieldlines.md](chapter_7_examples/7.3_fieldlines.md) |
| 7.4 Study2Cooling.g4bl | [7.4_study2cooling.md](chapter_7_examples/7.4_study2cooling.md) |
| 7.5 MultipleScattering.g4bl | [7.5_multiple_scattering.md](chapter_7_examples/7.5_multiple_scattering.md) |
| 7.6 TungstenTarget.g4bl | [7.6_tungsten_target.md](chapter_7_examples/7.6_tungsten_target.md) |
| 7.7 MICE_StageVI.g4bl | [7.7_mice_stagevi.md](chapter_7_examples/7.7_mice_stagevi.md) |
| 7.8 Idealized_g-2.g4bl | [7.8_idealized_g2.md](chapter_7_examples/7.8_idealized_g2.md) |
| 7.9 SpaceCharge.g4bl | [7.9_spacecharge.md](chapter_7_examples/7.9_spacecharge.md) |
| 7.10 SampleMovie.g4bl | [7.10_sample_movie.md](chapter_7_examples/7.10_sample_movie.md) |
| 7.11 Movie.in | [7.11_movie_in.md](chapter_7_examples/7.11_movie_in.md) |
| 7.12 triplet.sh – tune a quad triplet for point-to-point focus | [7.12_triplet.md](chapter_7_examples/7.12_triplet.md) |
| 7.13 emittancematch.sh – match a quad triplet into a solenoid | [7.13_emittancematch.md](chapter_7_examples/7.13_emittancematch.md) |

---

## Chapter 8 — Tips and Techniques

| Section | File |
|---------|------|
| 8.1 Getting Help on Using G4beamline – Users Forum | [8.1_getting_help.md](chapter_8_tips/8.1_getting_help.md) |
| 8.2 Reporting Bugs in G4beamline | [8.2_reporting_bugs.md](chapter_8_tips/8.2_reporting_bugs.md) |
| 8.3 Requesting New Features in G4beamline | [8.3_requesting_features.md](chapter_8_tips/8.3_requesting_features.md) |
| 8.4 Getting Help on Individual G4beamline Commands | [8.4_command_help.md](chapter_8_tips/8.4_command_help.md) |
| 8.5 In what Directory should I Work? | [8.5_working_directory.md](chapter_8_tips/8.5_working_directory.md) |
| 8.6 Files Common to Multiple Simulations | [8.6_common_files.md](chapter_8_tips/8.6_common_files.md) |
| 8.7 Basic Execution in a Command-Line Environment | [8.7_command_line_execution.md](chapter_8_tips/8.7_command_line_execution.md) |
| 8.8 Basic Execution in a GUI Environment | [8.8_gui_execution.md](chapter_8_tips/8.8_gui_execution.md) |
| 8.9 Putting Shielding into a Simulation | [8.9_shielding.md](chapter_8_tips/8.9_shielding.md) |
| 8.10 How to Debug a Simulation | [8.10_debugging.md](chapter_8_tips/8.10_debugging.md) |
| 8.11 Geant4 Commands | [8.11_geant4_commands.md](chapter_8_tips/8.11_geant4_commands.md) |
| 8.12 Obtaining Plots and Histograms | [8.12_plots_histograms.md](chapter_8_tips/8.12_plots_histograms.md) |
| 8.13 Obtaining Pictures of the System and Events | [8.13_pictures.md](chapter_8_tips/8.13_pictures.md) |
| 8.14 Warning and Error messages – Which Ones can be ignored | [8.14_warnings_errors.md](chapter_8_tips/8.14_warnings_errors.md) |
| 8.15 Secondary Tracks and Particles | [8.15_secondary_tracks.md](chapter_8_tips/8.15_secondary_tracks.md) |
| 8.16 Finding Example Input Files using the XXX command | [8.16_finding_examples.md](chapter_8_tips/8.16_finding_examples.md) |
| 8.17 Parameterizing the Input File | [8.17_parameterizing.md](chapter_8_tips/8.17_parameterizing.md) |
| 8.18 Setting Fields of Magnets | [8.18_setting_fields.md](chapter_8_tips/8.18_setting_fields.md) |
| 8.19 Tuning Bending Magnets | [8.19_tuning_bending.md](chapter_8_tips/8.19_tuning_bending.md) |
| 8.20 Setting the Phase of RF Cavities (pillbox) | [8.20_rf_phase.md](chapter_8_tips/8.20_rf_phase.md) |
| 8.21 Tuning the maxGradient of RF Cavities (pillbox) | [8.21_rf_gradient.md](chapter_8_tips/8.21_rf_gradient.md) |
| 8.22 Multiple Jobs in Parallel | [8.22_parallel_jobs.md](chapter_8_tips/8.22_parallel_jobs.md) |

---
