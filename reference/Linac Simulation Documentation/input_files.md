# Input File Templates

> **Source:** Adam Bartnik, June 2026 — *"Here's the input files — at least the templates for them with default settings. The magnets / cavities / etc may default to some weird settings, but at least they should be located in the right place."*

These are the template input files the Linac Sim GUI ships and fills in at run time. They define the three-region chain (1D cathode → GPT → BMAD/Tao) described in the [README](README.md) and detailed in [Technical Details](details.md). They are **templates**: the GUI substitutes user/machine settings (and `$…$` placeholders) before each run, so the committed defaults are placeholders — most magnets, cavities, and powers default to `0`. Trust the **element positions and lattice topology**, not the default field strengths.

All six files live in [`input_files/`](input_files/).

> **These are reference artifacts, not a format to match.** The old GUI's input files (cathode 1D / GPT / BMAD / Tao) do **not** need to be mirrored by the WarpX rebuild in this repo. The WarpX stages are configured by Python module-level constants and `config()` overrides (see the root `CLAUDE.md`), not by these files. Use these templates only to read off the original LinacSim physics and geometry — there is no requirement to reproduce their file layout, syntax, or `@GUI` conventions.

| File | Code / consumer | Region | Role |
|------|-----------------|--------|------|
| `cathode_master.in` | Custom 1D cathode code (Java) | Cathode → grid | Cathode/grid geometry and voltage-pulse parameters |
| `gpt_master.in` | GPT | Gun → end of section 1 | Full GPT beamline: gun, prebunchers, lenses/solenoids, section-1 TW linac |
| `bmad_master.in` | BMAD (via Tao) | After section 1 → end of linac | Master lattice: user params, RF/quad conversions, lattice assembly |
| `section_1_4_layout.bmad` | BMAD | Sections 1–4 | Element-by-element layout (drifts, quads, `lcavity`s) |
| `section_5_8_layout.bmad` | BMAD | Sections 5–8 | Element-by-element layout (drifts, quads, `lcavity`s) |
| `tao.init` | Tao | BMAD driver | Tao startup: universe, beam tracking, RF/radiation switches |

## The `@GUI` annotation convention

Lines beginning with `# @GUI` (GPT/cathode) or `! @GUI` (BMAD) are **machine-readable comments** the GUI parses to build its parameter panels — they are inert to the underlying code. The grammar seen across the files:

- `@GUI beamline "<name>"` / `@GUI section "<name>"` — group headings in the left panel.
- `@GUI element <Type> "<label>" <position_m> <length_m> "<code>" [...]` — one beamline element; the trailing code tag (`"Cathode1D"`, `"GPT"`, `"BMAD"`) routes it to the right region. The two numbers are **z-position (m) and length (m)** — this is the "located in the right place" Adam refers to.
- `@GUI parameter "<label>" <var> <default> "<unit>" <editable>` — a tunable bound to a variable used later in the file (e.g. `egun_gun_volt`, `ACC_2_power`, `QH1_current`).
- `@GUI phasedparameter …` — an on-crest phase, paired with its power/relative-phase controls (prebunchers, section-1 linac).
- `@GUI APPENDMARKERS` — directive (in `bmad_master.in`) telling the GUI where to inject viewscreen/output markers.

Because the variables named in `@GUI parameter` are the same ones assigned in the body, editing a value in the GUI just rewrites that assignment before the run.

## Template placeholders (`$…$`)

`gpt_master.in` contains `$…$` tokens the GUI substitutes per run:

- `$INIT_DIST_FILE$` — the cathode code's output distribution, fed to GPT via `settdist(...)` (the cathode → GPT handoff).
- `$FIELD_DIR$` — directory holding the GPT field maps (`CESR_gun.gdf`, `prebuncher_25D.gdf`, `LENS_0A…0E.gdf`, `SOL_0.gdf`, `SLAC-3mLinac-field1/2.gdf`). These are the same maps committed under `fieldmaps/` and rebuilt by each stage's `build_*_field.py`.

`tao.init` likewise references `bmad_master-run.in` — the filled-in copy of `bmad_master.in` the GUI emits (not the raw template).

## Notes per file

### `cathode_master.in`
Pure `@GUI` block — only the cathode/grid panel. Parameters: cathode-grid distance `l` (0.2 mm), pulse voltage `Vpulse` (60 V), off voltage `Voff` (−30 V), max slope `Vp` (30 V/ns), pulse width `twidth` (2 ns), grid transmission `trans` (80 %). These map onto the voltage-pulse model `V(t) = Voff + Vpulse·f(t; Vp, t0)` in [Technical Details](details.md#voltage-pulse-model).

### `gpt_master.in`
The GPT region in full, gated by `if`-switches at the top (`gun_section_on`, `prebunchers_section_on`, `section_1_section_on`, `auto_phase`, `space_charge`, `single_particle`, …), all defaulting to `0`/off. Key physics that mirror the WarpX rebuild:
- **Globals:** `Master_RF = 11.89915e6` Hz; `Linac_RF = Master_RF·240`, `Prebuncher_RF = Master_RF·18`, `CESR_RF = Master_RF·42`.
- **Gun:** `Map2D_E(..., CESR_gun.gdf, ...)` simulated at 1 kV with the wrong sign, then scaled to `egun_gun_volt` (150 kV default) — the negative-field-scale convention noted in `gun/README.md`.
- **Prebunchers:** amplitude from measured loaded Q (`prebuncher1_QL=3000`, `prebuncher2_QL=4300`); prebuncher 2 placed with a reversed `-1,0,0` direction vector (installed backwards).
- **Section-1 linac:** traveling wave as two standing waves 90° apart (`SLAC-3mLinac-field1/2.gdf`, normalized to 0.001 MW).
- **Apertures/`ZSTOP`:** `scatterpipe`/`scatteriris`/`forwardscatter` set beam-pipe boundaries; each subsection sets its own `ZSTOP` (gun 0.5 m, prebunchers 2.1 m, section 1 5.45 m). Section-1 ends at z ≈ 5.4 m, the `@BMADSTART` handoff to BMAD.
- **Initial distribution:** uniform circular spot (radius `egun_cath_diam/2`), thermal spread from `egun_cath_T` (1425 K), time profile from `$INIT_DIST_FILE$`.

### `bmad_master.in` + layout files
`bmad_master.in` is the master: a big `@GUI` block (positions + powers + quad currents for all 8 sections), the user-parameter assignments (all powers/currents default `0`), then `call, file = section_1_4_layout.bmad` and `section_5_8_layout.bmad`, then the conversions "no need to touch below here":
- **RF→energy:** gradient `= sqrt(power)·dEdP/L`, with per-section `dEdP` calibrations (`33e6/sqrt(15e6)` etc.) and `phi0 = rel_phase/(2π)`.
- **Quad current→field:** `b1_gradient = ±conv·current` with placeholder `conv = 0.0125 T/A` for every quad — these are the "weird default" conversions; real A→T calibrations are still unknown (see [Quadrupole Magnets](details.md#quadrupole-magnets)).
- **Assembly:** `sec_1_8: line = (sec_1_4, sec_5_8)`; `use, sec_1_8`.

The two `*_layout.bmad` files define the elements with **measured/estimated geometry** (`l`, elliptical `aperture`s, `lcavity` `n_cell` and `rf_frequency = Linac_RF`). `section_1_4` begins at z = 5.4 m (after section-1 solenoids have decayed, before quad Q1). Both carry the same standing caveats in their headers: linacs still need swapping to traveling wave, lengths are flange-to-flange (not true active length), and pipe radii are guesses.

### `tao.init`
Minimal Tao driver: one universe, `track_type = 'beam'`, `rf_on = T`, radiation damping/fluctuations and CSR all **off**. Loads `bmad_master-run.in` (the GUI-filled lattice) and saves beam at every marker (`beam_saved_at = "marker::*"`).

## Relationship to the WarpX rebuild in this repo

These templates are the **upstream reference** for the from-scratch WarpX pipeline (`cathode/` → `gun/` → `prebuncher/` → `linac_sec1/`). The GPT field-map calls, RF harmonics, gun voltage scaling, prebuncher Q-based scaling, and SLAC two-standing-wave decomposition documented here are exactly what each WarpX stage reproduces from first principles. When a WarpX operating point looks off, these files are the source of truth for the original LinacSim settings — but remember the committed values are GUI defaults (mostly zero/off), not Adam's tuned working point.

### Fidelity vs. the WarpX rebuild

**Important context:** the WarpX stages in this repo were built and their numeric values chosen *before* these original LinacSim input files were available. So where the rebuild diverges from these files, it is generally because the original value was simply **not yet known** — not a deliberate decision to differ. These files are now the source of truth; the rebuild can be reconciled toward them where it matters. The list below reflects the state as of when the files arrived (June 2026); reconcile and update it as stages are brought into line.

**Core physics that already matches** (these were derivable from the field maps and physics, so the rebuild got them right independently): gun voltage (150 kV) and its 1 kV / negative-scale sign convention; the field maps themselves (`CESR_gun.gdf`, `prebuncher_25D.gdf`, `SLAC-3mLinac-field1/2.gdf`); RF harmonics (prebuncher 18× = 214.18 MHz, linac 240× ≈ 2856 MHz); the prebuncher Q-based amplitude scaling `√(1e3·Q·P/(2πf))`; the cos/sin (E/B 90°-apart) RF drive; the SLAC two-standing-wave traveling-wave decomposition with its +π/2 quadrature offset and 0.001 MW (1 kW) normalization; the linac bore (≈9.55 mm) and ~37 MeV-at-15-MW energy target.

| Quantity | Original LinacSim file |
|----------|------------------------|
| Bunch charge | **1 nC** (`gpt_master.in` `total_charge = -1e-9`) |
| Cathode diameter | 16 mm (`egun_cath_diam`) |
| Cathode temperature | 1425 K (`egun_cath_T`) |
| Cathode–grid distance | 0.2 mm (`l`) |
| Prebuncher z-position | 0.534 m (P1), 1.318 m (P2) |
| Prebuncher power / phase | GUI 8 kW @ −70° (P1), 10 kW @ −45° (P2) |
| Linac input power | 11 MW (GUI default) |
| Linac on-crest phase | 183.5° absolute |

**Elements not yet modeled at all**:

- Solenoid lenses **0A–0E** and **Sol 0** — absent from the WarpX prebuncher transport.
- **Prebuncher 2** (installed reversed, `-1,0,0`, Q=4300, z=1.318 m) — the rebuild models a single prebuncher.
- Collimator/iris apertures (9.547 mm) in the injector line.
- **BMAD sections 2–8** (`section_1_4_layout.bmad`, `section_5_8_layout.bmad`, `bmad_master.in`, `tao.init`) — no WarpX counterpart yet. (Note: the BMAD "Section 1" block is drift/quad transport *after* the section-1 linac; the accelerating cavity itself lives only in `gpt_master.in`, which `linac_sec1/` reproduces.)
