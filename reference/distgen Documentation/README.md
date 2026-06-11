# distgen

Particle-distribution **generator** for accelerator simulations: builds 6D+time (and spin)
initial particle coordinates from a declarative YAML/JSON input, with pint-quantity units
throughout, and writes them for GPT, ASTRA, SIMION, or openPMD consumers. Source repository:
[distgen on GitHub](https://github.com/ColwynGulliford/distgen).
Available on conda-forge (`conda install -c conda-forge distgen`). It is the standard beam
front end for lume-gpt (`run_gpt_with_distgen`) and the LinacSim-style GPT decks.

This folder is a checkout of the package itself: the code is the reference
([distgen/](distgen/)); the upstream mkdocs site ([docs/](docs/)) is a stub, so this README
is the working index.

---

## Overview

```python
from distgen import Generator
G = Generator('distgen.in.yaml', verbose=1)   # file path, YAML string, or dict
G['n_particle'] = 100_000                      # nested-key access, ':'-separated
G['r_dist:sigma_xy:value'] = 1.5               # quantities keep their units
G.run()                                        # returns/fills G.particles
P = G.particles                                # openPMD-beamphysics ParticleGroup
P.plot('x', 'px')
```

`Generator` ([distgen/generator.py](distgen/generator.py)) parses the input, draws random
numbers (`random_type: hammersley` quasirandom by default, or `random`), samples each
requested coordinate distribution, applies the ordered `transforms`, and assembles a
[`Beam`](distgen/beam.py) → `ParticleGroup`. `G.beam()` returns the intermediate `Beam`
(coordinates `x,y,z` [m], `px,py,pz` [eV/c], `t` [s], `weight` [C], optional spin
`sx,sy,sz`). All scalar inputs are **pint quantities** written as
`{value: 1.0, units: mm}` in YAML.

A CLI exists too: `python -m distgen.command_line -f distgen.in.yaml -v 1` (or the
installed `distgen` entry point), driving [distgen/drivers.py](distgen/drivers.py)
`run_distgen()`.

## Input file anatomy

```yaml
n_particle: 2000
species: electron
random_type: hammersley
total_charge: {value: 10, units: fC}
start:                      # how particles begin
  type: cathode             # 'cathode' | 'time' | 'free'
  MTE: {value: 5, units: meV}
r_dist:                     # one <var>_dist block per sampled coordinate
  type: rg                  # radial gaussian
  sigma_xy: {value: 1.0, units: mm}
  n_sigma_cutoff: {value: 3.0, units: dimensionless}
t_dist:
  type: sg                  # super-gaussian in time
  sigma_t: {value: 1, units: ps}
  alpha: {value: 0.25, units: dimensionless}
transforms:
  order: [t0, t1]
  t0: {type: set_stdxy x:y, sigma_xy: {value: 1.5, units: um}}
  t1: {type: set_avg t, avg_t: {value: 0, units: ps}}
output:
  type: gpt                 # 'gpt' | 'astra' | 'openPMD' | 'simion'
  file: gpt.particles.gdf
```

`start.type` sets the launch convention:
- **`cathode`** — particles begin on the photocathode plane: forward-hemisphere momenta from
  `MTE` (Maxwell-Boltzmann transverse-energy sampling) or an explicit photoemission model,
  emission *time* spread from `t_dist`. GPT output gets cathode-start metadata
  (`status` handling matches `ParticleGroup` cathode conventions).
- **`time`** — a bunch at one instant; sample `z` (or use transforms) instead of `t`.
- **`free`** — fully explicit: every coordinate from its own dist.

Relative `file:` paths inside dist blocks are expanded against the input file's directory
(`parsing.expand_input_filepaths`).

## Distribution catalog (`type:` names, [distgen/dist.py](distgen/dist.py))

1D (`x_dist`, `y_dist`, `z_dist`, `t_dist`, `px_dist`, …):

| `type` (alias) | Class | Notes |
|----------------|-------|-------|
| `uniform` (`u`) | `Uniform` | min/max or avg+sigma |
| `gaussian` (`g`) | `Norm` | optional `n_sigma_cutoff`, asymmetric cuts |
| `super_gaussian` (`sg`) | `SuperGaussian` | `alpha` flat-top↔gaussian morph |
| `tukey` | `Tukey` | tapered flat-top |
| `sech2` | `Sech2` | soliton-like laser pulse |
| `file1d` | `File1d` | tabulated PDF from file |
| `interp` | `Interpolation1d` | interpolated control points |
| `deformable` | `Deformable` | spline-deformable profile |
| `superposition` (`sup`) / `product` (`pro`) | `Superposition`/`Product` | combine sub-dists |
| `maxell_boltzmann` (`mb`) / `maxell_boltzmann_kinetic_energy` (`mbe`) | momentum / KE thermal dists *(sic — upstream spelling)* |
| `crystals` | `TemporalLaserPulseStacking` | birefringent-crystal pulse stacking |

Radial (`r_dist`) and angles:

| `type` (alias) | Class |
|----------------|-------|
| `radial_uniform` (`ru`) | `UniformRad` |
| `radial_gaussian` (`rg`) | `NormRad` |
| `radial_super_gaussian` (`rsg`) | `SuperGaussianRad` |
| `radial_tukey` | `TukeyRad` |
| `radfile` / `radial_interpolation` (`ri`) / `raddeformable` (`dr`) | `RadFile` / `InterpolationRad` / `DeformableRad` |
| `uniform_theta` (`ut`) / `uniform_phi` (`up`) | `UniformTheta` / `UniformPhi` |

(`Linear`/`LinearRad` classes exist but are not dispatched as standalone `type:` names —
they back the `deformable`/`raddeformable` dists.)

2D / ND (`xy_dist`, momentum models):

| `type` | Class | Notes |
|--------|-------|-------|
| `file2d` | `File2d` | gridded 2D density file |
| `image2d` | `Image2d` | image file (PNG/JPG…) → transverse density |
| `uniform_laser_speckle` | `UniformLaserSpeckle` | synthetic speckle ([laser_speckle.py](distgen/laser_speckle.py)) |
| `superposition` / `product` (2D) | `SuperPosition2d` / `Product2d` | |
| `nd_gaussian` | `GaussianNd` | correlated N-dim gaussian |
| `fermi_dirac_3step_barrier_photocathode` (`fd3sb`) | `FermiDirac3StepBarrierMomentumDist` | Dowell-Schmerge-style photoemission momentum model ([fermi_dirac_3step_barrier_photocathode_model.py](distgen/fermi_dirac_3step_barrier_photocathode_model.py)) |

## Transforms ([distgen/transforms.py](distgen/transforms.py))

Applied in `transforms.order` sequence; the `type` string is `"<name> <vars>"` with
`:`-separated variables (e.g. `set_stdxy x:y`, `rotate2d x:y`):

`translate`, `set_avg`, `scale`, `set_std`, `set_stdxy`, `set_avg_and_std`, `rotate2d`,
`shear`, `polynomial`, `cosine`, `matrix2d`, `magnetize` (adds canonical angular momentum —
a θ-dependent kick equivalent to cathode-field immersion), `set_twiss` (set β/α/ε per
plane), `set_min_r_at_theta` / `set_min_r_at_theta_for_round_beam` (hollow-beam shaping).

## Output writers ([distgen/writers.py](distgen/writers.py))

| `output.type` | Function | Notes |
|---------------|----------|-------|
| `gpt` | `write_gpt` | writes ASCII then converts via `$ASCI2GDF_BIN` (`asci2gdf`); GPT columns `x y z GBx GBy GBz t q nmacro` |
| `astra` | `write_astra` | ASTRA conventions incl. reference particle |
| `openPMD` | `write_openPMD` | openPMD-beamphysics HDF5 |
| `simion` | `write_simion` | SIMION `.ion` format |

Programmatically you rarely need these — take `G.particles` (a `ParticleGroup`) and use its
own writers (`P.write(...)`, `P.write_impact(...)`, etc.).

## Supporting modules

| Module | Purpose |
|--------|---------|
| [distgen/beam.py](distgen/beam.py) | `Beam` container (pint arrays, stats, twiss) |
| [distgen/physical_constants.py](distgen/physical_constants.py) | `PHYSICAL_CONSTANTS`, species masses/charges, the shared pint `unit_registry` |
| [distgen/hammersley.py](distgen/hammersley.py) | Hammersley/Halton quasirandom sequences |
| [distgen/parsing.py](distgen/parsing.py) | quantity-dict parsing, file-path expansion |
| [distgen/tools.py](distgen/tools.py) | weighted stats, radial integrals, nested-dict `:` flat keys, image reading |
| [distgen/metrics.py](distgen/metrics.py) | profile-nonuniformity metrics (KL divergence, rms-equivalent) |
| [distgen/archive.py](distgen/archive.py) | HDF5 archive of generator input |
| [distgen/plot.py](distgen/plot.py) | `plot_dist1d/2d`, current profile, radial dist plots |
| [distgen/dist.py](distgen/dist.py) | all distribution classes + `get_dist` dispatch |

## Gotchas

- **Every dimensional scalar must carry units** (`{value: ..., units: ...}`); raw floats are
  only for counts/flags. Wrong-dimension units raise immediately via pint.
- Nested-key access uses **`:` separators** (`G['t_dist:sigma_t:value'] = 2.0`); this is
  also the convention lume-gpt's `settings={'distgen:n_particle': ...}` routing relies on.
- `random_type: hammersley` (default in many examples) is deterministic — identical inputs
  reproduce identical beams; use `random_type: random` + seed control for statistical
  independence.
- For `start: cathode`, GPT/ASTRA outputs encode "born later" via the particle time `t` —
  downstream readers must respect emission-time spread rather than assume a z-slice. This
  matches how this repo's WarpX cathode stage emits over time.
- The thermal-distribution type names are misspelled upstream (`maxell_boltzmann`) — use the
  `mb`/`mbe` aliases to avoid surprises.
- `total_charge: {value: 0, units: pC}` gives equal-weight zero-charge macros (tracking
  without space charge); species charge sign is applied internally (electron weights are
  positive `|q|`).
