# openPMD-beamphysics

Python tools for analyzing, converting, and writing particle and field data in the
openPMD standard (beamphysics extension). Provides two core classes used throughout
this repo's IMPACT-T work: **`ParticleGroup`** (particle data) and **`FieldMesh`**
(external field-map data — e.g. `FieldMesh.from_onaxis` for the DC gun field in
`impactt_cathode/`). Source repository: [openPMD-beamphysics on GitHub](https://github.com/ChristopherMayes/openPMD-beamphysics).

The `api/` files give the class signatures, properties, and method lists for
`ParticleGroup` and `FieldMesh` (filled in from the rendered docs site, since the
upstream sources are build-time mkdocstrings directives). The example notebooks below
are the worked-example reference.

---

## Overview

| Section | File |
|---------|------|
| Overview, installation, class summary | [index.md](index.md) |
| Beam-statistics YAML schema (named stat quantities + definitions) | [statistics_schema.md](statistics_schema.md) |

---

## API reference (autodoc stubs)

| Class | File |
|-------|------|
| `ParticleGroup` | [api/particles.md](api/particles.md) |
| `FieldMesh` | [api/fields.md](api/fields.md) |
| Wakefield classes (`ResistiveWallWakefield`, `ImpedanceWakefield`, …) | [api/wakefields.md](api/wakefields.md) |

---

## Examples (notebooks — the practical reference)

### Particles
| Topic | File |
|-------|------|
| Reading particle data (multiple code formats) | [examples/read_examples.ipynb](examples/read_examples.ipynb) |
| Writing particle data | [examples/write_examples.ipynb](examples/write_examples.ipynb) |
| `ParticleGroup` features (stats, slicing, derived quantities) | [examples/particle_examples.ipynb](examples/particle_examples.ipynb) |
| Plotting | [examples/plot_examples.ipynb](examples/plot_examples.ipynb) |
| Units handling | [examples/units.ipynb](examples/units.ipynb) |
| Plot labels | [examples/labels.ipynb](examples/labels.ipynb) |
| Normalized coordinates | [examples/normalized_coordinates.ipynb](examples/normalized_coordinates.ipynb) |
| Bunching factor | [examples/bunching.ipynb](examples/bunching.ipynb) |

### Fields (`FieldMesh`)
| Topic | File |
|-------|------|
| Field-map basics | [examples/fields/field_examples.ipynb](examples/fields/field_examples.ipynb) |
| Format conversion (incl. IMPACT-T `rfdata`) | [examples/fields/field_conversion.ipynb](examples/fields/field_conversion.ipynb) |
| On-axis expansion | [examples/fields/field_expansion.ipynb](examples/fields/field_expansion.ipynb) |
| Field analysis | [examples/fields/field_analysis.ipynb](examples/fields/field_analysis.ipynb) |
| Field tracking | [examples/fields/field_tracking.ipynb](examples/fields/field_tracking.ipynb) |
| Solenoid modeling | [examples/fields/solenoid_modeling.ipynb](examples/fields/solenoid_modeling.ipynb) |
| Corrector modeling | [examples/fields/corrector_modeling.ipynb](examples/fields/corrector_modeling.ipynb) |

### Wakefields & wavefronts
| Topic | File |
|-------|------|
| Resistive-wall wakefield | [examples/wakefields/resistive_wall.ipynb](examples/wakefields/resistive_wall.ipynb) |
| Impedance wakefield | [examples/wakefields/impedance_wakefield.ipynb](examples/wakefields/impedance_wakefield.ipynb) |
| Wavefront basics | [examples/wavefront/wavefront.ipynb](examples/wavefront/wavefront.ipynb) |
| Wavefront drift | [examples/wavefront/advanced_drift.ipynb](examples/wavefront/advanced_drift.ipynb) |
| Experimental PyTorch wavefront | [examples/wavefront/experimental_pytorch_wavefront.ipynb](examples/wavefront/experimental_pytorch_wavefront.ipynb) |

> Sample `examples/data/` (HDF5/field-map binaries) was removed during import — it is
> input data for running the notebooks, not documentation. Re-fetch from the upstream
> repo if a notebook needs to be executed.
