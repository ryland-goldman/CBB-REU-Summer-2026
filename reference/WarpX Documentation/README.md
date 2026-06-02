# WarpX Documentation

WarpX is an advanced, massively parallel **Particle-In-Cell (PIC)** code built on the [AMReX](https://amrex-codes.github.io/) adaptive mesh refinement framework. It supports electromagnetic and electrostatic simulations across multiple geometries and compute backends, and was awarded the [2022 ACM Gordon Bell Prize](https://www.exascaleproject.org/ecp-supported-collaborative-teams-win-the-2022-acm-gordon-bell-prize-and-special-prize/).

**GitHub:** https://github.com/BLAST-WarpX/warpx

---

## Key Features

- **Field solvers:** Maxwell's equations (FDTD, PSATD), Poisson's equation, Ohm's law / kinetic-fluid hybrid
- **Geometries:** 1D, 2D, 3D Cartesian; cylindrical (RZ); embedded boundaries
- **Compute backends:** Multi-core CPU (OpenMP), NVIDIA/AMD/Intel GPUs (CUDA, HIP, SYCL)
- **Multi-physics:** Ionization, atomic/fusion/collisional physics, QED
- **Advanced numerics:** Explicit and implicit time advance, adaptive mesh refinement, boosted-frame simulations, load balancing
- **Python interface:** Run and extend simulations via `pywarpx` / PICMI; couple to AI/ML frameworks
- **Platforms:** Linux, macOS, Windows; scales to the world's largest supercomputers

---

## Table of Contents

### Installation
| Topic | File |
|-------|------|
| User install (pip / spack / conda) | [install/users.md](install/users.md) |
| CMake build from source | [install/cmake.md](install/cmake.md) |
| HPC systems (Frontier, Perlmutter, LUMI, …) | [install/hpc/](install/hpc/) |
| Batch schedulers (Slurm, PBS, LSF, …) | [install/batch/](install/batch/) |

### Usage
| Topic | File |
|-------|------|
| How to run WarpX | [usage/how_to_run.md](usage/how_to_run.md) |
| All input parameters | [usage/parameters.md](usage/parameters.md) |
| Python (PICMI) interface | [usage/python.md](usage/python.md) |
| Example simulations | [usage/examples.md](usage/examples.md) |
| Plasma wakefield acceleration (PWFA) | [usage/pwfa.md](usage/pwfa.md) |
| FAQ | [usage/faq.md](usage/faq.md) |

### Tutorials & Examples
| Example | File |
|---------|------|
| Laser wakefield acceleration (LWFA) | [usage/examples/lwfa/README.md](usage/examples/lwfa/README.md) |
| Plasma wakefield acceleration (PWFA) | [usage/examples/pwfa/README.md](usage/examples/pwfa/README.md) |
| Langmuir wave | [usage/examples/langmuir/README.md](usage/examples/langmuir/README.md) |
| Laser-ion acceleration | [usage/examples/laser_ion/README.md](usage/examples/laser_ion/README.md) |
| Capacitive discharge | [usage/examples/capacitive_discharge/README.md](usage/examples/capacitive_discharge/README.md) |
| Beam-beam collision | [usage/examples/beam_beam_collision/README.md](usage/examples/beam_beam_collision/README.md) |
| All tutorials index | [tutorials.md](tutorials.md) |

### Workflows
| Topic | File |
|-------|------|
| Extending WarpX with Python | [usage/workflows/python_extend.md](usage/workflows/python_extend.md) |
| Python callbacks | [usage/workflows/python_callbacks.md](usage/workflows/python_callbacks.md) |
| Accessing field data from Python | [usage/workflows/python_field_data.md](usage/workflows/python_field_data.md) |
| Accessing particle data from Python | [usage/workflows/python_particle_data.md](usage/workflows/python_particle_data.md) |
| AI-assisted input design | [usage/workflows/ai_input_design.md](usage/workflows/ai_input_design.md) |
| Optimization with Optimas | [usage/workflows/optimas.md](usage/workflows/optimas.md) |
| Debugging | [usage/workflows/debugging.md](usage/workflows/debugging.md) |
| Domain decomposition | [usage/workflows/domain_decomposition.md](usage/workflows/domain_decomposition.md) |

### Data Analysis
| Topic | File |
|-------|------|
| Output formats (openPMD, HDF5, ADIOS) | [dataanalysis/formats.md](dataanalysis/formats.md) |
| openPMD API | [dataanalysis/openpmdapi.md](dataanalysis/openpmdapi.md) |
| openPMD-viewer | [dataanalysis/openpmdviewer.md](dataanalysis/openpmdviewer.md) |
| yt integration | [dataanalysis/yt.md](dataanalysis/yt.md) |
| 3D visualizations | [dataanalysis/3dvisualizations.md](dataanalysis/3dvisualizations.md) |
| In-situ analysis (Ascent, Catalyst, SENSEI) | [dataanalysis/insitu.md](dataanalysis/insitu.md) |
| ParaView | [dataanalysis/paraview.md](dataanalysis/paraview.md) |
| VisIt | [dataanalysis/visit.md](dataanalysis/visit.md) |

### Theory
| Topic | File |
|-------|------|
| Introduction | [theory/intro.md](theory/intro.md) |
| Electromagnetic PIC | [theory/models_algorithms/electromagnetic_pic.md](theory/models_algorithms/electromagnetic_pic.md) |
| Explicit EM PIC | [theory/models_algorithms/explicit_em_pic.md](theory/models_algorithms/explicit_em_pic.md) |
| Implicit EM PIC | [theory/models_algorithms/implicit_em_pic.md](theory/models_algorithms/implicit_em_pic.md) |
| Electrostatic PIC | [theory/models_algorithms/electrostatic_pic.md](theory/models_algorithms/electrostatic_pic.md) |
| Kinetic-fluid hybrid model | [theory/models_algorithms/kinetic_fluid_hybrid_model.md](theory/models_algorithms/kinetic_fluid_hybrid_model.md) |
| Adaptive mesh refinement | [theory/amr.md](theory/amr.md) |
| Boosted-frame simulations | [theory/boosted_frame.md](theory/boosted_frame.md) |
| Boundary conditions | [theory/boundary_conditions.md](theory/boundary_conditions.md) |
| Collisions | [theory/multiphysics/collisions.md](theory/multiphysics/collisions.md) |
| Ionization | [theory/multiphysics/ionization.md](theory/multiphysics/ionization.md) |
| QED | [theory/multiphysics/qed.md](theory/multiphysics/qed.md) |

### Developer Guide
| Topic | File |
|-------|------|
| Contributing | [developers/contributing.md](developers/contributing.md) |
| Repository organization | [developers/repo_organization.md](developers/repo_organization.md) |
| CMake build (local dev) | [developers/how_to_compile_locally.md](developers/how_to_compile_locally.md) |
| Testing | [developers/how_to_test.md](developers/how_to_test.md) |
| Code architecture: fields | [developers/fields.md](developers/fields.md) |
| Code architecture: particles | [developers/particles.md](developers/particles.md) |
| Code architecture: diagnostics | [developers/diagnostics.md](developers/diagnostics.md) |
| Dimensionality | [developers/dimensionality.md](developers/dimensionality.md) |
| AMReX basics | [developers/amrex_basics.md](developers/amrex_basics.md) |
| Python bindings | [developers/python.md](developers/python.md) |
| Portability (CPU/GPU) | [developers/portability.md](developers/portability.md) |
| Profiling | [developers/how_to_profile.md](developers/how_to_profile.md) |
| Writing docs | [developers/how_to_write_the_docs.md](developers/how_to_write_the_docs.md) |
| Developing with LLMs | [developers/how_to_develop_with_llms.md](developers/how_to_develop_with_llms.md) |
| Developer FAQ | [developers/faq.md](developers/faq.md) |

### Community
| Topic | File |
|-------|------|
| How to acknowledge WarpX | [acknowledge_us.md](acknowledge_us.md) |
| Science highlights | [highlights.md](highlights.md) |
| Glossary | [glossary.md](glossary.md) |
| Code of conduct | [coc.md](coc.md) |
| Governance | [governance.md](governance.md) |
| Acknowledgements | [acknowledgements.md](acknowledgements.md) |

---

## Contact & Support

- **Discussions / Q&A:** https://github.com/BLAST-WarpX/warpx/discussions
- **Bug reports / feature requests:** https://github.com/BLAST-WarpX/warpx/issues
- **Watch / star the repo** to receive updates and support the project
