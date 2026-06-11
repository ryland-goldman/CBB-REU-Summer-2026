# lume-gpt

Python wrapper for **GPT (General Particle Tracer)** within the LUME framework: builds/edits
GPT input decks, runs the `gpt` binary, parses the GDF output into openPMD-beamphysics
`ParticleGroup` objects, and provides autophasing, lattice construction, plotting, and HDF5
archiving. Source repository:
[lume-gpt on GitHub](https://github.com/ColwynGulliford/lume-gpt) (Colwyn Gulliford).
Available on conda-forge (`conda install -c conda-forge lume-gpt`). Requires a separately
licensed **GPT** install — see the [GPT website](http://www.pulsar.nl/gpt/) and the local
`reference/GPT Documentation/`.

This folder is a checkout of the package itself: the code is the reference
([src/gpt/](src/gpt/)), and the [examples/](examples/) notebooks + templates are the
practical worked reference.

> **Naming gotcha:** the distribution is named `lume-gpt` but the *import* is `gpt`
> (`from gpt import GPT`) — `pyproject.toml` declares `name = "gpt"`, packaged from `src/`.

---

## Overview

The central object is `GPT` (in [src/gpt/gpt.py](src/gpt/gpt.py)):

```python
from gpt import GPT
G = GPT(input_file='gpt.in', initial_particles=P)   # P: ParticleGroup (optional)
G.set_variable('sol_1_current', 3.2)                 # deck variables: name=value;
G.run()
G.stat('sigma_x', 'screen')                          # any ParticleGroup stat vs z
G.screen[-1]                                         # final screen as ParticleGroup
G.plot(['sigma_x', 'sigma_y'], x='mean_z')
```

Workflow: `__init__` parses the deck (`parsers.parse_gpt_input_file`) and sets up a temp
working dir (`use_tempdir=True` default; support files — field maps, `setfile` beams — are
referenced absolutely or copied with `copy_support_files=True`). `run()` writes the deck,
invokes the GPT binary, and `load_output()` parses the output GDF into **touts**
(time snapshots) and **screens** (fixed-z planes), each a `ParticleGroup`.

The GPT binary is found via the **`$GPT_BIN`** environment variable (likewise
`$GDF2A_BIN` / `$ASCI2GDF_BIN` for the converters used by `executables.py` and phasing).

Top-level exports (`gpt/__init__.py`): `GPT`, `run_gpt`, `evaluate_gpt`,
`run_gpt_with_distgen`, `evaluate_gpt_with_distgen`.

---

## `GPT` class essentials

| Member | Purpose |
|--------|---------|
| `GPT(input_file, initial_particles=None, gpt_bin='$GPT_BIN', use_tempdir=True, workdir=None, timeout=None, ccs_beg='wcs', spin_tracking=False, n_cpu=1, ...)` | Parse deck + configure run dir. `initial_particles` (a `ParticleGroup`) is written to GDF and wired into the deck's `setfile("beam", ...)` line. |
| `.set_variable(name, val)` / `.set_variables(dict)` | Set deck variables declared as `name=value;` lines. |
| `.run(gpt_verbose=False)` | Execute GPT (timeout + `kill_msgs` watchdog), then parse output. |
| `.tout`, `.screen`, `.particles` | Lists of `ParticleGroup` (time snapshots / z-screens / both). `n_tout`, `n_screen` count them. `tout_ccs` gives touts in the bend-following coordinate system. |
| `.stat(key, data_type='all'|'tout'|'tout_ccs'|'screen')` | Any `ParticleGroup` stat vs output index; special-cases `twiss_*`, field samples `fEx…fBz` (needs `load_fields=True`), spin `spinx/y/z`, `spin_polarization`. `.units(key)` gives units. |
| `.trajectory(pid, data_type='tout')` | Single-particle trajectory arrays by particle ID. |
| `.plot(y=[...], y2=[...], x='mean_z')` | Stat plot with layout (`plot.plot_stats_with_layout`). |
| `.track(particles, s=None)` / `.track1(...)` / `.track1_to_z(...)` / `.track1_in_ccs(...)` | Track a `ParticleGroup` (or one constructed single particle) through the configured deck. |
| `.archive(h5)` / `.load_archive(h5)` / `GPT.from_archive(h5)` | Full input+output round-trip to HDF5 ([src/gpt/archive.py](src/gpt/archive.py)). |
| `.auto_phase()` | Run the marker-comment phasing pass (below) using the centroid of `initial_particles`. |
| `.copy()`, `.fingerprint()` | Deep copy; hash of the configured input for caching. |
| `GPT.from_tao(tao, ...)` | Build a GPT deck from a Bmad/Tao lattice ([src/gpt/gpt_tao.py](src/gpt/gpt_tao.py)). |

Module-level `run_gpt(settings, gpt_input_file, ...)` is the one-call version:
make object → apply `settings` → run.

---

## Autophasing (`gpt_phasing.py`)

RF cavity phasing is driven by **marker comments in the GPT deck** naming the variable to
solve for (the numeric suffix orders the pass):

| Deck marker | Meaning |
|-------------|---------|
| `phasing_amplitude_N = "varname";` | Amplitude variable for cavity group N (restored after phasing). |
| `phasing_on_crest_N = "varname";` | Phase variable solved for max energy gain (on-crest). |
| `phasing_relative_N = "varname";` | Phase offset applied relative to the found crest. |
| `phasing_gamma_N = "varname";` | Expected gamma bookkeeping variable. |

`gpt_phasing(path_to_input_file, gpt_bin='$GPT_BIN', gdf2a_bin='$GDF2A_BIN',
path_to_phasing_dist=None)` replaces the beam with a single centroid particle, brents each
cavity to crest in order, and writes `<deck>.phased.in`. The
[examples/auto_phasing.ipynb](examples/auto_phasing.ipynb) notebook demonstrates; in
`run_gpt_with_distgen(..., auto_phase=True)` it runs automatically with a centroid particle
derived from the distgen distribution (`get_distgen_beam_for_phasing`).

This is the same marker-comment convention used by Adam Bartnik's LinacSim GPT decks (see
`reference/Linac Simulation Documentation/` and the
[examples/templates/cu_injector/gpt.in](examples/templates/cu_injector/gpt.in) deck —
that template *is* the Cornell DC-gun injector front end, with the same
`dcgun_GHV`, `buncher_CTB`, and `solenoid_SLA_L60` field maps this repo's chain rebuilds).

## distgen integration (`gpt_distgen.py`)

`run_gpt_with_distgen(settings, gpt_input_file, distgen_input_file, auto_phase=False, ...)`
runs distgen → writes the GDF beam → runs GPT. `settings` keys are routed automatically:
keys prefixed `distgen:` (e.g. `distgen:n_particle`) update the distgen YAML tree; bare keys
set GPT deck variables (`set_gpt_and_distgen`). `evaluate_gpt_with_distgen` /
`evaluate_gpt` wrap a run and return a merit dict (`merit.default_gpt_merit`: end emittances,
sigmas, energy, charge conservation error) for optimizers (e.g. xopt).

## Programmatic lattices (`lattice.py`, `element.py`, `maps.py`, `bstatic.py`)

An alternative to hand-written decks: build a `Lattice` and write the deck from it.

| Piece | Contents |
|-------|----------|
| `Lattice(name)` | `.add(element, ds, ref_element=, ref_origin='end', element_origin='beg')` s-ordered placement; `.write_gpt_lines(...)` deck emission; indexing by name; bend-aware CCS chaining. |
| `element.py` | Base `Element` (geometry, CCS, placement), `Screen`, `Quad`, `SectorBend`, `Beg`. |
| `maps.py` | GDF field-map elements: `Map1D_E`, `Map1D_B`, `Map1D_TM` (RF), `Map2D_E`, `Map2D_B`, `Map25D_TM`, `Map3D_E`, `Map3D_B`; plus `autophase`/`autophase_track1` (crest a single map), `track_on_axis`, `larmor_angle`, `cavity_voltage`, floor-plan plotting. |
| `bstatic.py` / `estatic.py` | Analytic magnets: `Sectormagnet`, `Rectmagnet`, `Quadrupole`, `QuadF`, `Bzsolenoid`; electrostatic `Erect`. |
| `remove_particles.py` | `Aperture`, `CircularAperture` (GPT `rmax`-style scraping elements). |
| `wien_filter.py` | `WienFilter` element (+ `WienFilter3D` map in `maps.py`). |
| `autoscale.py` | `autophase1`/`autoscale1`: single-particle `ztrack1` passes that phase/scale each auto-element in lattice order. |

Element-by-element notebooks live in [examples/elements/](examples/elements/)
(one per element type, with sample field maps in `examples/elements/fields/`).

## Supporting modules

| Module | Purpose |
|--------|---------|
| [src/gpt/parsers.py](src/gpt/parsers.py) | `parse_gpt_input_file`, `write_gpt_input_file`, support-file path rewriting, GDF readers. |
| [src/gpt/particles.py](src/gpt/particles.py) | Raw GDF dicts → `ParticleGroup`s (`raw_data_to_particle_groups`, `gdf_to_particle_groups`); `identify_species` (mass+charge → name); `particle_stats`. |
| [src/gpt/easygdf.py](src/gpt/easygdf.py) | **Vendored** GDF loader (`load`, `load_initial_distribution`) returning tout/screen dicts. Distinct from the standalone `easygdf` package documented in `reference/easygdf Documentation/` — both exist in this codebase; lume-gpt's requirements pull the standalone one too. |
| [src/gpt/executables.py](src/gpt/executables.py) | `gpt`, `gdf2a`, `asci2gdf` subprocess wrappers; `$GPT_BIN`-style env-var expansion decorator. |
| [src/gpt/evaluate.py](src/gpt/evaluate.py), [src/gpt/merit.py](src/gpt/merit.py) | Merit-function evaluation for optimization loops. |
| [src/gpt/archive.py](src/gpt/archive.py) | openPMD-flavored HDF5 archive read/write. |
| [src/gpt/gpt_tao.py](src/gpt/gpt_tao.py) | Bmad/Tao → GPT deck conversion (fieldmap packing, bends, soft-edge solenoids). |
| [src/gpt/template.py](src/gpt/template.py) | Minimal deck templates (`basic`, `ztrack1`) used by tracking/autoscale. |
| [src/gpt/watcher.py](src/gpt/watcher.py), [src/gpt/fake_gpt.py](src/gpt/fake_gpt.py) | Run watchdog; fake GPT binaries for tests. |

---

## Examples (the practical reference)

| Topic | File |
|-------|------|
| Minimal run + plotting | [examples/simple_gpt_example.ipynb](examples/simple_gpt_example.ipynb), [examples/drift.ipynb](examples/drift.ipynb) |
| Autophasing markers in a real deck | [examples/auto_phasing.ipynb](examples/auto_phasing.ipynb) |
| **Cornell DC-gun injector** (gun + buncher + solenoids) | [examples/cu_injector_example.ipynb](examples/cu_injector_example.ipynb), templates [cu_injector](examples/templates/cu_injector/), [cu_injector_frontend](examples/templates/cu_injector_frontend/), [dcgun](examples/templates/dcgun/) |
| APEX RF gun, SRF injector, TESLA 9-cell | [examples/apex_gun_example.ipynb](examples/apex_gun_example.ipynb), [examples/srf_injector.ipynb](examples/srf_injector.ipynb), [examples/tesla_9cell_cavity_example.ipynb](examples/tesla_9cell_cavity_example.ipynb) |
| Space charge with cathode emission (3D tree) | [examples/spacecharge3dTree_with_cathode.ipynb](examples/spacecharge3dTree_with_cathode.ipynb) |
| GDF → `ParticleGroup` conversion | [examples/gdf_to_particles.ipynb](examples/gdf_to_particles.ipynb) |
| Element/lattice construction (per-element) | [examples/elements/](examples/elements/) |
| Tao/Bmad → GPT FODO | [examples/tao_to_gpt_fodo_example.ipynb](examples/tao_to_gpt_fodo_example.ipynb) |
| Spin tracking | [examples/spin_tracking.ipynb](examples/spin_tracking.ipynb), [examples/gpt_manual_spin_tracking.ipynb](examples/gpt_manual_spin_tracking.ipynb) |
| Parameter sweeps / multiple runs | [examples/multiple_runs.ipynb](examples/multiple_runs.ipynb) |

Template decks (deck + `distgen.in.yaml` + field maps) are under
[examples/templates/](examples/templates/).

---

## Gotchas

- **Import name is `gpt`, not `lume_gpt`** (and the version lookup is `version("gpt")`).
- This checkout imports openPMD-beamphysics as `beamphysics` (`from beamphysics import
  ParticleGroup`) — released lume-gpt versions import `pmd_beamphysics`. Match package
  versions accordingly (same situation in the `GPT_tools` checkout).
- The GPT binary is **not bundled** — set `$GPT_BIN` (and `$GDF2A_BIN`/`$ASCI2GDF_BIN`); a
  GPT license is required.
- Deck variables are only settable if declared as plain `name=value;` lines —
  `set_variable` errors otherwise.
- Touts and screens are different samplings: touts are equal-*time* snapshots (beam spread
  over z), screens are equal-*z* crossings (beam spread over t). Compare against WarpX
  openPMD dumps (which are touts) accordingly.
- `auto_phase` requires `initial_particles` (it phases with the bunch centroid).
- For bend lattices, lab-frame stats come from `tout`; use `tout_ccs` / `s_ccs` for the
  beamline-following frame.
