# openPMD-viewer

Tools to load and visualize openPMD file series (typically a time series of field and
particle diagnostics). Used in this repo to read WarpX openPMD output — e.g.
`OpenPMDTimeSeries` / `ts.get_field(...)` in `warpx_cathode/plot_cathode.py`. Usable
both as a **Python API** and as an **interactive Jupyter GUI**. Source repository:
[openPMD-viewer on GitHub](https://github.com/openPMD/openPMD-viewer).

The `api_reference/` files carry the full class/method signatures and parameter docs
(filled in from the rendered ReadTheDocs pages, since the upstream sources are
build-time autodoc directives). The tutorial notebooks are the worked-example reference.

---

## Overview

| Section | File |
|---------|------|
| Overview, installation, usage | [source/index.md](source/index.md) |

---

## API reference

| Class | File |
|-------|------|
| `OpenPMDTimeSeries` — `get_field`, `get_particle`, `iterate`, `slider` | [source/api_reference/generic_interface.md](source/api_reference/generic_interface.md) |
| `ParticleTracker` | [source/api_reference/particle_tracking.md](source/api_reference/particle_tracking.md) |
| `LpaDiagnostics` — energy spread, emittance, current, a0, … | [source/api_reference/lpa_diagnostics.md](source/api_reference/lpa_diagnostics.md) |
| API reference index | [source/api_reference/api_reference.md](source/api_reference/api_reference.md) |

---

## Tutorials (notebooks — the practical reference)

| Topic | File |
|-------|------|
| Introduction to the API | [source/tutorials/1_Introduction-to-the-API.ipynb](source/tutorials/1_Introduction-to-the-API.ipynb) |
| Specific field geometries (thetaMode, slices) | [source/tutorials/2_Specific-field-geometries.ipynb](source/tutorials/2_Specific-field-geometries.ipynb) |
| Introduction to the interactive GUI | [source/tutorials/3_Introduction-to-the-GUI.ipynb](source/tutorials/3_Introduction-to-the-GUI.ipynb) |
| Particle selection | [source/tutorials/4_Particle_selection.ipynb](source/tutorials/4_Particle_selection.ipynb) |
| Laser-plasma (LPA) tools | [source/tutorials/5_Laser-plasma_tools.ipynb](source/tutorials/5_Laser-plasma_tools.ipynb) |
