# GPT_tools

**Adam Bartnik's** analysis/plotting/run-harness layer on top of lume-gpt, distgen, and
openPMD-beamphysics — interactive (ipywidgets) plot GUIs for lume-gpt runs, an extended
`ParticleGroup`, settings-driven run wrappers with auto-phasing, screen postprocessing, an
xopt Pareto-front explorer, and several specialized physics models (THz interaction, image
charge, inverse Compton, tip emission). Source repository:
[GPT_tools on GitHub](https://github.com/AdamCBartnik/GPT_tools) (`pip install` from
source; see [setup.py](setup.py) / [requirements.txt](requirements.txt) — depends on
`lume-gpt` (`gpt`), `distgen`, `openpmd-beamphysics`, `xopt`, `ipywidgets`).

This folder is a checkout of the package itself: the code is the reference
([GPT_tools/](GPT_tools/)); [examples/example.ipynb](examples/example.ipynb) and
[examples/example_population.ipynb](examples/example_population.ipynb) are the worked
demos (with a sample [gpt.in](examples/gpt.in), [distgen.in.yaml](examples/distgen.in.yaml),
and [xopt.in.yaml](examples/xopt.in.yaml)).

Directly relevant to this repo: it is built by the author of LinacSim (see
`reference/Linac Simulation Documentation/`), and its plot conventions (trend vs z,
dist1d/dist2d with stats tables) are the model for comparing this repo's WarpX chain
against GPT runs of the same beamline.

---

## Plotting ([GPT_tools/gpt_plot.py](GPT_tools/gpt_plot.py), [gpt_plot_gui.py](GPT_tools/gpt_plot_gui.py))

All take a lume-gpt `GPT` object (or its data) and auto-format units/labels
([nicer_units.py](GPT_tools/nicer_units.py)):

| Function | Purpose |
|----------|---------|
| `gpt_plot(gpt_data, var1, var2, ...)` | Trend plot of any stat(s) vs `mean_z`/`t`, e.g. `gpt_plot(G, 'mean_z', 'sigma_x')`; supports two y-axes, survivor filtering (`show_survivors_at_z=`), and screen markers. |
| `gpt_plot_dist1d(pmd, var, plot_type='charge', ...)` | 1D histogram of a screen/tout `ParticleGroup` with a stats table (emittances, sigmas, charge). |
| `gpt_plot_dist2d(pmd, var1, var2, plot_type='histogram', ...)` | 2D scatter/histogram with density coloring + stats table. |
| `gpt_plot_trajectory(gpt_data, var1, var2, nlines=, ...)` | Individual particle trajectories through the touts. |
| `gpt_plot_gui(gpt_data)` | **ipywidgets GUI** wrapping all of the above: dropdown variables, screen slider (with "special screen" detection), log/linear, slicing — the quickest way to explore a finished `GPT` run in Jupyter. |

Plot/stat helpers live in [GPT_tools/tools.py](GPT_tools/tools.py)
(`get_screen_data`, `special_screens`, weighted stats, unit scaling) and
[SnappingCursor.py](GPT_tools/SnappingCursor.py) (interactive cursor used by the GUIs).

## Extended ParticleGroup ([GPT_tools/ParticleGroupExtension.py](GPT_tools/ParticleGroupExtension.py))

`ParticleGroupExtension(ParticleGroup)` adds derived properties usable in any stat/plot
call: `core_emit_x/y/4d` (Gaussian-core emittance fit), `slice_emit_x/y/4d` (with
`n_slices`/`slice_key`), `sqrt_norm_emit_4d`, `root_norm_emit_6d`, `action_x/y/4d`,
`transverse_energy`, `ptrans`, `rp`, `r_centered`, `pr_centered`, …

| Helper | Purpose |
|--------|---------|
| `convert_gpt_data(G)` | Wrap every tout/screen of a `GPT` object in `ParticleGroupExtension`. |
| `slice_emit(p_list, key)` | Slice-emittance across a divided bunch. |
| `divide_particles(pg, nbins, key)` | Split a bunch into slices by any variable. |
| `core_emit_calc(x, xp, w)` / `core_emit_calc_4d(...)` | Core-emittance fitting primitives. |
| [emittance_vs_fraction.py](GPT_tools/emittance_vs_fraction.py) `emittance_vs_fraction(pg, var, ...)` | Emittance vs enclosed-fraction curve (minimum-bounding-ellipse based). |

## Run harness ([GPT_tools/GPTExtension.py](GPT_tools/GPTExtension.py))

Settings-dict driven wrappers that chain **distgen → lume-gpt (+ auto-phasing) →
postprocess → merit**, one level above `run_gpt_with_distgen`:

| Function | Purpose |
|----------|---------|
| `run_gpt_with_settings(settings, gpt_input_file=, distgen_input_file=, ...)` | Main entry: distgen beam (or cached), phased GPT run, returns the `GPT` object. |
| `run_gpt_with_particlegroup(settings, ...)` | Same but injecting an existing `ParticleGroup` (chain GPT runs stage-to-stage). |
| `evaluate_run_gpt_with_settings` / `evaluate_run_gpt_with_particlegroup` | Run + return a flat merit dict (extended `default_gpt_merit`: end/screen stats, core/slice emittances) — the xopt evaluation functions. |
| `multithread_gpt_with_settings(...)` | Split the bunch and run GPT on N threads, recombining output (`split_particle_group`, `run_one_thread`). |
| `evaluate_gpt_with_stability` / `evaluate_multirun_gpt_with_stability` / `add_jitter_to_settings` | Jitter/stability studies: re-run with randomized settings offsets and collect spreads. |
| `run_gpt_with_THz` / `run_gpt_with_analytic_THz` | Runs inserting a THz interaction kick (below). |
| `clip_to_charge` / `clip_to_emit` / `radius_including_charge` / `radius_including_emit` | Charge-/emittance-fraction clipping utilities used in merits. |

Settings keys follow the lume-gpt routing convention: `distgen:`-prefixed keys edit the
distgen YAML tree, bare keys set GPT deck variables.

## Screen postprocessing ([GPT_tools/postprocess.py](GPT_tools/postprocess.py))

`postprocess_screen(screen, **params)` dispatches: `take_range(var, min, max)`,
`take_slice(var, index, n_slices)`, `clip_to_charge`, `clip_to_emit`,
`remove_correlation(var1, var2, max_power)` (polynomial decorrelation),
`remove_spinning` (subtracts coherent angular momentum), `kill_zero_weight`,
`include_ids` / `id_of_nearest_N` / `random_N` (subselection),
`add_cylindrical_copies` (azimuthal cloning for axisymmetric statistics),
`keep_only_last_forward_pass` (drops back-bending trajectory segments, in GPTExtension).

## Cathode beams ([GPT_tools/cathode_particlegroup.py](GPT_tools/cathode_particlegroup.py))

`get_cathode_particlegroup(settings, distgen_input_file)` — distgen run → `ParticleGroup`
with the `distgen:` settings applied; `get_coreshield_particlegroup` builds a
core+shield two-population cathode beam. [image_charge.py](GPT_tools/image_charge.py)
implements a cathode image-charge/Schottky emission model (MTE/QE vs extraction field,
`MakeMetalParticleGroup`).

## xopt front explorer ([GPT_tools/front_tools.py](GPT_tools/front_tools.py), [front_gui.py](GPT_tools/front_gui.py))

For CNSGA (xopt) optimizer output populations (`pop_*.csv`/json):
`show_fronts(...)` (Pareto fronts across generations), `get_pop`, `find_settings`
(pick the individual nearest a target objective), `reevaluate_population`,
`clamp_population`, and `front_gui(xopt_file, pop_directory)` — an ipywidgets GUI with
snapping cursor to browse individuals and pull their settings. Demo:
[examples/example_population.ipynb](examples/example_population.ipynb).

## Specialized physics models

| Module | Purpose |
|--------|---------|
| [THz_functions.py](GPT_tools/THz_functions.py) | THz-pulse/beam interaction (lump-element kick, parabolic mirror geometry, analytic vs GPT comparison). |
| [compton.py](GPT_tools/compton.py) | `inverse_compton_scatter(PG_e, ...)` — ICS photon ParticleGroup from an electron bunch + laser. |
| [tip_emission.py](GPT_tools/tip_emission.py) / [conical_tip_emission.py](GPT_tools/conical_tip_emission.py) | Nanotip field-emission distributions. |
| [boundary_element_solver.py](GPT_tools/boundary_element_solver.py) | Axisymmetric electrostatic BEM solver (electrode → on-axis field, e.g. for gun field maps). |

---

## Gotchas

- **Not on PyPI/conda-forge** — install from the GitHub checkout (`pip install .`); the
  heavy dependency set (xopt, ipywidgets, fastnumbers, deap via xopt) is only needed for
  the GUI/front modules; plotting + ParticleGroupExtension work with just
  lume-gpt + openpmd-beamphysics + matplotlib.
- The GUIs (`gpt_plot_gui`, `front_gui`) are **Jupyter ipywidgets** apps, not standalone
  windows.
- `convert_gpt_data` deep-copies every tout/screen — fine for typical runs, slow for
  thousands of dumps.
- Several functions mutate the passed screen unless `make_copy=True` (postprocess family) —
  the default differs per function, check the signature.
- Imports use the lume-gpt package as `gpt`, and **this checkout (like the lume-gpt one)
  imports openPMD-beamphysics as `beamphysics`** — released versions use `pmd_beamphysics`;
  keep the two checkouts in `reference/` version-matched.
