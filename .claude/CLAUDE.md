# CLAUDE.md

Guidance for Claude Code when working in this repository.

## Overview

A beam-dynamics simulation of the **Cornell High Energy Synchrotron Source (CHESS)** electron
source, built for a Research Experience for Undergraduates (REU) at the Cornell Center for Bright
Beams (CBB) / CLASSE. It rebuilds the front end of the Cornell Linac (Adam Bartnik's LinacSim:
thermionic source → gun → injector → linac) from first principles. The stages form one
self-consistent chain — each reads the previous stage's openPMD beam. Sections 1–4 of the linac are
**WarpX** (RZ); the **4→5 boundary** (the linac section-4 exit) is the slot for the e⁺/e⁻ converter
target (G4beamline/Geant4); sections 5–8 are **Impact-T**.

```
cathode ─► gun ─► injector ─► linac1 ─► linac2 ─► linac3 ─► linac4 ─►[4→5]─► converter ─► linac5-8
SCL       CESR    2 prebunchers  SLAC 3 m   CEA 2     CEA 3    CU 5    converter  W target   CU/CEA S-band
diode     gun     + 6 solenoids  TW capture                            slot                  (Impact-T)
(2D)      (RZ)    (RZ)           (RZ)       (RZ)      (RZ)     (RZ)               (G4bl)      (sections 5–8)
```

This repository is a **simplifying refactor** of an earlier per-stage-package layout (`../Cornell`):
the old `cathode/`, `gun/`, … `linac_rest/` packages and the `pipeline/` orchestration package were
consolidated into `config/` (one YAML per stage, all tuning hardcoded), `sim/` (one driver per
stage + shared `helpers/` and `plot/`), `logs/` (all generated output), and `docs/` (per-stage
physics notes). The four WarpX linac sections share **one** parametrized driver (`sim/linac1-4.py`,
section chosen by a CLI argument); the four Impact-T sections share `sim/linac5-8.py`. The e⁺/e⁻
converter target between them (`sim/converter.py`) is a G4beamline/Geant4 shower run, not WarpX/Impact-T.

Use the CBB conda environment: `conda activate CBB` (Miniforge at ~/miniforge3; if conda isn't on
PATH, first `source ~/miniforge3/bin/activate`).

## Codebase Standards

Documentation lives in `docs/` and the `README.md`; code comments hold only what a reader of that
line genuinely needs. Keep the code scannable and the prose in one place.

- **Docs live in `docs/<stage>.md` and `README.md`, not code.** Physics rationale, operating
  points, field-map provenance, and tuning narrative belong there (architecture goes here). The
  docs must NOT quote simulation **result** numbers (energies, transmission %, emittances) — those
  change with tuning; state configured operating points only.
- **Module docstrings are short** — what the file does (≤ ~6 lines) + `See docs/<stage>.md`.
- **Keep a comment only if it prevents a silent regression** a reader of that line could not
  otherwise know: sign/unit conventions (negative field scale; `u = γβ`; momentum `= γβ·mₑ·c`),
  WarpX/Impact-T footguns tied to a line (MLMG `dirichlet` outer wall; "last applied field must
  load_E"; θ₀ is *absolute* per Impact-T section; `species` singular vs `electrons` plural;
  field-free-pad sampling), or non-obvious local logic (multi-plane iris scrape; transit-stop
  sizing; "transmission measured before charge re-impose"). Prefer one terse line.
- **No dead code, no history/TODO prose** in comments.
- **Naming:** module constants `UPPER_SNAKE`, functions `lower_snake`, classes `PascalCase`,
  private helpers leading `_`. All tuning is hardcoded in `config/<stage>.yaml` — there is no
  `config()` override layer (a deliberate simplification from the old repo).

## Commands

All commands run from the **repo root** in the `CBB` environment (stage scripts use repo-relative
paths). Keep `OMP_NUM_THREADS=1` (the default) — see *Threads*.

```bash
conda activate CBB
pip install -r requirements.txt                 # pywarpx / impact-t / openpmd-api best via conda

python sim/main.py                              # full chain, per-stage subprocesses + final summary
python sim/cathode.py                           # one stage's simulation
python sim/plot/cathode.py                      # its figures (from existing diagnostics)
python sim/linac1-4.py 2                        # a linac section (argument selects 1/2/3)
python sim/plot/linac1-4.py 2                   # a linac section's figures (from existing diagnostics)
```

- **No `config()`/profile API.** Each WarpX stage reads `config/<stage>.yaml`; the Impact-T stage
  reads `config/linac5-8.yaml`. Retune by editing the YAML. The shipped values are the single
  Balanced operating point.
- **Frozen RF setpoints.** The linac RF crest phases and field scales are **hardcoded** in
  `config/linac2.yaml`, `linac3.yaml`, `linac4.yaml`, `linac5-8.yaml` (`CREST_PHASE_DEG`/`FIELD_SCALE`,
  and the per-section `crest_phase_deg`/`field_scale` table). They were derived once from the beam and
  the drivers simply read+apply them — there is no runtime crest-finding or calibration loop. **If you
  change an upstream knob that shifts the beam, the affected setpoint must be re-derived** (load the
  upstream exit dump + the on-axis field, run the crest math once, paste the result; tooling:
  `sim/autophase.py` for WarpX 1–4, `sim/autophase_impact.py` for Impact-T 5–8). The Impact-T
  per-section crest is *absolute* (`theta0_deg`), so it is only valid for the deck geometry it was
  derived on — keep the real-length zero-K1 inter-section quads unchanged. **linac4 + linac5-8
  setpoints are currently un-derived placeholders** (see the linac docs).
- **Space charge / speed.** The self-field MLMG Poisson solve dominates WarpX runtime. The
  relativistic linac sections 2–4 run with `warpx_do_not_deposit: true` (**SC off**): at γ≳45 the
  self-field is 1/γ²-negligible, so this is byte-identical physics for ~50–80× speedup. SC stays
  **on** where it matters — cathode (the Child–Langmuir limit it exists to show), gun (150 keV
  magnetic pinch), and linac section 1 (150 keV capture). Set `warpx_do_not_deposit: false` for a
  fully self-consistent linac run.
- **Threads:** keep `OMP_NUM_THREADS=1`. These grids are small and the MLMG solve is
  memory-bandwidth bound, so OpenMP threads contend and add overhead with no gain — raising it
  *slows* the small-grid stages (cathode/gun) via oversubscription. `prepare_env()` and
  `sim/main.py` default it to 1.
- **Run one stage off existing upstream output:** each stage reads the previous stage's openPMD
  output from `logs/diags/`, so any unmodified upstream output is reused — run `python
  sim/linac1-4.py 2` directly if `logs/diags/linac1-4/sec1/` already exists.

There is no test suite or linter — validation is physics sanity checks (energy gain, Child–Langmuir
current, bunching, the ~296 MeV exit) printed by each run and inspected in `logs/plots/`.

## Project Architecture

```
config/   one YAML per stage — all tuning hardcoded (operating point + Balanced profile)
sim/      main.py + one driver per stage (cathode, gun, injector, linac1-4, linac5-8)
  helpers/ stage-agnostic plumbing (see below)
  plot/    one plotter per stage + common.py (shared figures)
logs/     diags/<stage>/ (openPMD + injection_summary.json) · plots/<stage>/ · pipeline/log_<date>.log
docs/     per-stage physics notes
fieldmaps/ gdf/ (GPT inputs, committed) · h5/ (built, git-ignored) · rfdata/ (Impact-T TW template)
```

**Shared helpers (`sim/helpers/`)** — the single source of truth the stage drivers import
(consolidated from the old `pipeline/` package):
- `tools.py` — scipy physical constants (`C_LIGHT`/`E_CHARGE`/`M_E`/`MC2_EV`…), emission physics
  (`child_langmuir_current_density`, `thermal_velocity_sigma`), the RF `warpx_*_time_function`
  string builder (`rf_time_functions`), and `prepare_env()` (OMP=1, HDF5 locking, fd-limit raise,
  repo-root chdir — call before importing `warpx`).
- `buildfields.py` — GDF→openPMD thetaMode field builders (`load_cols`/`to_grid`/`pad_r`/
  `write_thetamode_series` + `build_gun_field`/`build_injector_fields`/`build_linac_slac_fields`).
  Reads `fieldmaps/gdf/`, writes `fieldmaps/h5/`. The single SLAC map is shared by linac 1–4.
- `loadparticles.py` — beam-handoff IO (`open_particle_series`/`make_particle_group`/`downsample`/
  `beam_kinematics`/`load_warpx_exit_bunch` (captured-core cut)/`upstream_exit_lab_z` (lab-z chain)/
  the multi-plane iris scrape `pipe_violator_ids`+`survivor_mask`/the Impact-T adapters
  `read_warpx_dump`+`write_openpmd_particles`).
- `metrics.py` — pure-numpy beam moments (`screen_profile`).
- `tqdmwrapper.py` — the Impact-T `fort.18` progress bar (WarpX stages get theirs from lume-warpx's
  `w.run(progress=…)`).

**Stage drivers (`sim/`)** drive lume-warpx (`from warpx import WarpX`) for the WarpX stages and
lume-impact (`import impact`) for `linac5-8`. Each: `prepare_env()` → build field → import the
upstream beam (`WarpX(initial_particles=…)` / `read_warpx_dump`) → apply runtime values (RF time
functions, `dt`, `max_steps`, diag period) or the frozen setpoints → run → write `logs/diags/<stage>`
+ `injection_summary.json`. `sim/linac1-4.py <N>` branches on the section: N=1 is the capture stage
(iris scrape, `scale=√(POWER/RF_NORM)`, on-crest phase); N=2,3,4 read the previous exit's captured
core and apply the frozen `FIELD_SCALE`/`CREST_PHASE_DEG` (section 4's exit is the 4→5 boundary into
the converter).

**Orchestration:** `sim/main.py` runs each stage as a **fresh subprocess** (pywarpx binds one
geometry — 2D/RZ — per interpreter, so the WarpX stages must be isolated; Impact-T runs the same
way for uniformity). The sim and its plotter are separate subprocess calls; sim failures abort,
plot failures warn (a figure bug must not discard completed physics). Subprocess stdout → the run
log `logs/pipeline/log_<date>.log`; the progress bar (stderr) stays on the terminal.

**Inter-stage contract (order-dependent):**

| Stage | Reads | Writes |
|-------|-------|--------|
| `sim/cathode.py` | — (emits at cathode) | `logs/diags/cathode/{fields,particles}` |
| `sim/gun.py` | `logs/diags/cathode/particles` + `fieldmaps/h5/gun_E.h5` | `logs/diags/gun/{fields,particles,handoff}` |
| `sim/injector.py` | `logs/diags/gun/handoff` (else `…/particles`) + `fieldmaps/h5/{preb1_EB,preb2_EB,lens0a…e,sol0}.h5` | `logs/diags/injector/main` |
| `sim/linac1-4.py 1` | `logs/diags/injector/main/particles` (dump near z≈2.03 m, iris-scraped) + `fieldmaps/h5/linac_rf{1,2}.h5` | `logs/diags/linac1-4/sec1/main` |
| `sim/linac1-4.py 2` | `logs/diags/linac1-4/sec1/main/particles` (captured core) + the shared SLAC maps | `logs/diags/linac1-4/sec2/main` |
| `sim/linac1-4.py 3` | `logs/diags/linac1-4/sec2/main/particles` (captured core) | `logs/diags/linac1-4/sec3/main` |
| `sim/linac1-4.py 4` | `logs/diags/linac1-4/sec3/main/particles` (captured core) | `logs/diags/linac1-4/sec4/main` (the 4→5 boundary) |
| `sim/converter.py` | `logs/diags/linac1-4/sec4/main/particles` (e⁻) + g4bl (external) | `logs/diags/converter/main` (e⁺ beam) |
| `sim/linac5-8.py` | `logs/diags/converter/main/particles` (≥ `MIN_KE_MEV` e⁺ core) + `fieldmaps/rfdata/rfdata4–7` | `logs/diags/linac5-8/main` |

**Key gotchas** (detail in `docs/<stage>.md` — read before modifying a stage):
- **Cathode** keeps `warpx_do_not_deposit:false` (SC *is* the Child–Langmuir mechanism); dense-early
  diagnostic union slice.
- **Gun** is RZ electromagnetostatic (the γ² magnetic pinch); the cathode beam is remapped 2D-slab→RZ
  (r-importance resample for the 2πr Jacobian); **timed release only** (the 2 ns grid pulse, via a
  `beforestep` callback — the old "snapshot" path is dropped); the exit handoff is reconstructed by
  id-tracking each particle at its **first appearance in the field-free pad** past the field map
  (sampling an in-field particle inflates ε_n,x ~8×), and the run stops while the beam is in the pad
  (draining the domain aborts MLMG).
- **Injector** is RZ with SC off; solenoids at **native absolute machine-z** (not GUI-argmax-aligned —
  that mis-places the flat-top SOL_0 by +1.08 m); reversed Preb-2 ≡ `PREB2_REV_PHASE=π` in absolute
  phase; B-only solenoid maps must be listed **before** the E-loading RF maps (picmi disables the
  global E init-style for any `load_E:false` field); the 9.547 mm collimator is a **multi-plane id
  scrape** over the 1.922→2.03 m pipe (the beam converges through that tail, so a single 2.03 m cut
  overstates iris transmission ~3×).
- **Linac 1–4** share the SLAC quadrature maps (Re/Im halves of *one* TW section, summed at 90°);
  `RMAX=9.547 mm` is the SLAC bore / iris; keep cells near ≈3:1 aspect or the MLMG self-field solve
  diverges; the captured-core cut (`KE ≥ 0.5·median`) drops the slipping tail; lab-z chaining via
  `injection_summary.json` (`z_handoff_m` for sec1, `z_inject_lab_m` for sec2/3/4). **Section 4 (CU 5)
  is new** — its `config/linac4.yaml` `CREST_PHASE_DEG`/`FIELD_SCALE` are placeholders to be
  re-derived (`sim/autophase.py 4` + a field-scale fit) before it is physical.
- **Linac 5–8** (Impact-T): no field maps exist — the four S-band TW sections reuse the vendored
  `rfdata4–7` shape as a 4-line `solrf` superposition (+0/+30/+90/+0, body `/sin(β₀d)`), all physics
  in the frozen per-section `field_scale`; it accelerates the **converter positron beam** (`q=+e`,
  `Bcharge=+1`); `theta0_deg` is **absolute** so each section's crest is a distinct frozen number
  (and depends on deck geometry — keep the real-length quads); the shipped crests are STALE
  placeholders (re-derive for positrons + the section-5-start deck via `sim/autophase_impact.py`);
  SC off, quads off (K1=0); transmission measured from the **macro count before** re-imposing charge;
  `ParticleGroup.species` is `"electron"` (singular) but openPMD readers key `"electrons"` (plural);
  `ParticleGroup.write()` emits a viewer-incompatible openPMD — the handoff uses
  `loadparticles.write_openpmd_particles`.

**Determinism:** `RNG_SEED=0` everywhere and `OMP_NUM_THREADS=1` make the chain reproducible (the
frozen setpoints depend on it).

**Conventions:** generated output (`logs/`, `fieldmaps/h5/`, `warpx_used_inputs`, `fort.*`) is
git-ignored — regenerate by re-running. Field maps in `fieldmaps/gdf/` and the `rfdata/` template are
committed. To keep a result figure, `git add -f logs/plots/<stage>/<fig>.png`.

## Reference Materials

The `reference/` documentation library (WarpX, IMPACT-T/Z, GPT, BMAD, G4beamline, LinacSim,
lume-impact, lume-gpt, distgen, GPT_tools, openPMD-beamphysics/viewer, easygdf, and `Papers/`) lives
in the original repo at `../Cornell/reference/`. Read the relevant tool docs and papers in full
before writing physics code — under-reading (missing a convention → wrong physics) is far more
costly than over-reading.

- **lume-warpx** (the user's package, import name `warpx`): static YAML → `WarpX(input_file=…,
  path=…)` → `.get/.update` override API, `install_callback`, `initial_particles`→
  `FromInitialParticles`, `AppliedFromFile`→`LoadAppliedField`, `run(progress=…)` + `plot2D`/
  `plot1D`/`plot_fields`. Pinned `lume-warpx==1.0.1`.
- **lume-impact** (import name `impact`): `Impact` class for the `linac5-8` deck (solrf elements,
  `rf_field_scale` ControlGroup, `I.run()`, `I.stat(...)`). Pinned `lume-impact==0.11.0`.

### Key Concepts
- **PIC (Particle-In-Cell):** WarpX tracks particles on a mesh, solves fields on the grid.
- **Space charge:** beam Coulomb self-repulsion, dominant at low energy / high current; 1/γ²-suppressed
  (negligible) for the relativistic linac sections.
- **GDF:** GPT's native binary field-map format (read via `easygdf`).
- **Impact-T `fort.18`:** the longitudinal-position output the progress bar polls.
