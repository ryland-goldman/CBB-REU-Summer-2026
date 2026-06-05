# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

The code in this folder is meant to contain simulations for a Research Experience for Undergraduates (REU) program at the Cornell Center for Bright Beams (CBB) and Cornell Laboratory for Accelerator ScienceS and Education (CLASSE). The project focus is to build a beam simulation for the Cornell High Energy Synchotron Source (CHESS).

Concretely, the repo rebuilds the front end of the Cornell Linac (Adam Bartnik's LinacSim: thermionic source → gun → prebuncher) from first principles in **WarpX** via its Python/PICMI interface (`pywarpx`). The three stages form one self-consistent chain — each reads the previous stage's openPMD beam as input:

```
cathode  ─►  gun  ─►  prebuncher
(SCL diode)  (~148 keV)  (RF velocity bunching)
```

Please use the CBB conda environment by running `conda activate CBB` (Miniforge is installed at ~/miniforge3; if conda isn't on your PATH yet, first run source ~/miniforge3/bin/activate).

## Keeping Documentation in Sync

Whenever you add or significantly change a feature, stage, folder, script, or dependency, **update the docs in the same change** so they never drift from the code:

- **New stage / top-level directory** → add a `README.md` to it, add a row to the component table in the root `README.md`, and update the inter-stage contract table and architecture notes in this file (`CLAUDE.md`).
- **New or changed command / workflow** (run scripts, CLI flags, config toggles, env vars) → update the **Commands** section here and the relevant `README.md`.
- **New inter-stage input/output path** → update the inter-stage contract table in **Project Architecture**.
- **New dependency** → add it (pinned) to `requirements.txt`.
- **New reference doc or paper** → add it to the **Reference Materials** table below and, for papers, to `reference/Papers/README.md` (see *Adding New Papers*).
- **New WarpX gotcha / non-obvious convention** discovered while working on a stage → record it in that stage's `README.md`.
- **New, removed, or renamed result figure** (a `plot_*.py` adds/drops/renames a PNG) → update the corresponding entry in `FIGURES.md` (the visual index of every stage's `results/` figures) in the same change.

When in doubt, treat a doc update as part of "done" — a feature isn't complete until `CLAUDE.md`, the root `README.md`, and the stage `README.md` reflect it.

## Commands

All commands run from the **repo root** in the `CBB` environment. Stage scripts use hard-coded relative paths (e.g. `gun/diags`), so running from anywhere else breaks the inter-stage handoff.

```bash
conda activate CBB
pip install -r requirements.txt                 # pywarpx/openpmd-api are best via mamba

python pipeline/run_pipeline.py                  # full chain, live progress + final-beam summary
```

- **Stage API:** each stage package (`cathode/`, `gun/`, `prebuncher/`) is a top-level facade. `import cathode; cathode.run()` (likewise for `gun` and `prebuncher`) builds the field map (if any), runs the WarpX sim, and generates that stage's plots. Each exposes `config(**kwargs)` (override module-level constants — keys must match the names in each `<stage>/*.py`), `run(plots=True)`, and `plot()` (figures from existing diags). Implementation lives in `pipeline/_runner.py`.
- **Performance knobs:** the runtime-critical parameters are module-level constants in each `<stage>/*.py`, so they are tunable via `config()` from `run_pipeline.py` (which ships an editable "PERFORMANCE KNOBS" block: a Balanced profile active by default plus commented Conservative/Aggressive presets). Per stage: grid (`nx/nz`, `nr/nz`, `NR/NZ`), step count (`MAX_STEPS`/`CFL`/`TRANSIT_MARGIN`/`AVG_SPEED_FRAC`), Poisson solve (`REQUIRED_PRECISION`, `MAX_ITERS`), macroparticles (`PPC` for the cathode, `MAX_PART` downsample cap for the gun and prebuncher), and diagnostic dump count (`N_DIAGS`, `DIAG_PERIOD`). The prebuncher is ≈75% of total runtime. The gun's cells are near-isotropic, so it scales ≈ `nz²` (per-step cost ∝ cells, and `dz = zmax/nz` ⇒ fewer steps as `nz` drops) — coarsening `nz` is a big clean win. The **prebuncher is the opposite** (measured): its long-thin box has anisotropic cells and a convergence-bound MLMG solve, so coarsening `NZ` *slows* the per-step solve faster than it removes cells (1.37× slower per step at `NZ=512`) and under-resolves the ~1 mm bunch — keep `NZ=1024` and speed it via `CFL` (fewer steps) and `MAX_ITERS`/`REQUIRED_PRECISION` (cheaper solve) instead. Lowering knobs trades accuracy for speed; the stage-module defaults reproduce the original (8.8-min) run.
- **Run one stage off existing upstream output:** comment out the unwanted `<stage>.run()` calls in `pipeline/run_pipeline.py`, or just import the one you want (`import prebuncher; prebuncher.run()`) — each stage reads the previous stage's openPMD output from disk, so any unmodified upstream output is reused.
- **Prebuncher power/phase scan:** call `prebuncher.run()` once per operating point in a Python loop, e.g. `for p in (160, 300, 500, 800): prebuncher.config(POWER_W=p, OUTDIR=f"prebuncher/diags/P{p}_zc"); prebuncher.run(plots=False)`. Then a single `prebuncher.plot()` aggregates every `diags/P*` directory (see `prebuncher/README.md`).
- **Plots:** `<stage>.plot()` reads its `diags/` and writes PNGs to `results/`. `run()` calls `plot()` by default; pass `plots=False` to skip.
- **Threads:** `OMP_THREADS` (default 6) — the MLMG Poisson solve is memory-bandwidth bound, so using all cores is *slower*. Override via the `OMP_THREADS` env var (set before any pywarpx import; `config()` cannot set it).

There is no test suite, linter, or build step — validation is physics sanity checks (energy gain, Child–Langmuir current, bunching) printed by each run and inspected in the `results/` plots.

## Project Architecture

Each stage lives in its own `<stage>/` directory and follows the same script layout:

- `build_*_field.py` — converts a GPT `.gdf` field map from `fieldmaps/` into an openPMD `.h5` field mesh (via `easygdf` to read + `openPMD-api` to write) that WarpX loads as an external field. (The cathode has no field map; its field is self-consistent.)
- `*_sim.py` — the WarpX/PICMI run. Reads the upstream beam with `openPMD-viewer`, injects it, tracks through the stage, writes openPMD particle diagnostics to its own `diags/`.
- `plot_*.py` — reads `diags/`, writes figures to `results/`.
- `README.md` — the stage's physics, field map, operating point, and outputs.

**Inter-stage contract (the chain is order-dependent):**

| Stage | Reads | Writes |
|-------|-------|--------|
| `cathode/cathode_diode.py` | — (emits at cathode) | `cathode/diags/particles` |
| `gun/gun_sim.py` | `cathode/diags/particles` + `gun/gun_field/gun_E.h5` | `gun/diags` |
| `prebuncher/prebuncher_sim.py` | `gun/diags/particles` + `prebuncher/prebuncher_field/prebuncher_EB.h5` | `prebuncher/diags/<P..._...>` |

`pipeline/run_pipeline.py` orchestrates the whole chain by calling `cathode.run()`, `gun.run()`, `prebuncher.run()` in order. The shared runner in `pipeline/_runner.py` builds field maps and generates plots **in-process**, but spawns a **fresh Python subprocess** (`pipeline/_launch_sim.py`) for each sim — pywarpx binds globally to one geometry (2D/RZ/3D) at first `.so` load and caches diagnostic state by name, so chaining cathode (2D) → gun (RZ) → prebuncher (RZ) in one interpreter would trip `AssertionError: Diagnostic attributes not consistent`. Inside each sim, `run_step(...)` installs a `pywarpx.callbacks.installcallback("afterstep", …)` hook to drive a tqdm progress bar and redirects WarpX's noisy per-step stdout to the pipeline log file, so the bar updates on a clean terminal line. A structured DEBUG log lands in `pipeline/logs/pipeline_<timestamp>.log`. The cathode is 2D x–z; the gun and prebuncher are RZ (cylindrically symmetric).

Field maps (`CESR_gun.gdf`, `prebuncher_25D.gdf`) live in `fieldmaps/`; field-map paths are set near the top of each `build_*_field.py`. WarpX-specific gotchas accumulated per stage (negative field-scale conventions, thetaMode openPMD axis order, charge renormalization, where to stop the run to avoid MLMG aborts) are documented in each stage's `README.md` — read it before modifying a stage.

**Conventions:**

- Simulation outputs are git-ignored (`diags/`, `results/`, `*.h5`, `*.gdf`, logs); regenerate by re-running. Field maps in `fieldmaps/` are committed.
- Commit convention (matches existing history): for a stage, commit its `*.py` scripts + `README.md`, and `git add -f <stage>/results/*.png` to include the result figures (since `results/` is git-ignored). Do **not** commit `diags/`, `.h5`, or logs.

## Reference Materials

**Read aggressively and up front — the context window is a resource to be used, not conserved.** Before writing or modifying any code, load the documentation and papers into context generously rather than minimally. A typical task should begin by reading *several* relevant files in full, not skimming one. Specifically, at the start of a task:

- Read the `README.md` of **every** simulation stage involved in the task, plus the stages immediately upstream and downstream of it (the chain is order-dependent, so neighboring stages' conventions matter).
- Read the **full** reference doc(s) for each tool the task touches — not just the section you think is relevant. The relevant detail (a field-scale sign, an axis-order convention, a solver flag) is often elsewhere in the doc.
- Read the relevant `reference/Papers/` entries in full when the task involves the underlying physics (emission, space charge, RF bunching, beam optics).
- When unsure whether a doc is relevant, **read it anyway.** Under-reading (missing a convention and producing wrong physics) is far more costly here than over-reading. Err toward filling the context window with primary sources before you start coding.

This applies on every model, and is mandatory on Opus: load all relevant documentation and papers into context immediately after the task is specified, before planning or editing.

The tables below index what's available; `reference/Papers/README.md` indexes the papers.

### Simulation Codes

| Tool | Location | Purpose |
|------|----------|---------|
| **IMPACT-T** | `reference/Impact-T Documentation/README.md` | 3D relativistic particle tracking with space charge, wakefields, and CSR. Parallel implementation, used in photoinjector design. |
| **IMPACT-Z** | `reference/Impact-Z Documentation/README.md` | 3D parallel PIC code for intense beams through drifts, quadrupoles, solenoids, bending magnets, multipoles, and RF cavities. |
| **GPT** | `reference/GPT Documentation/README.md` | General Particle Tracer — 3D charged particle dynamics including space charge. Uses GDF file format for I/O. |
| **WarpX** | `reference/WarpX Documentation/README.md` | Massively parallel PIC code (EM and electrostatic). Supports GPU backends (CUDA/HIP/SYCL), adaptive mesh refinement, Python interface via `pywarpx`/PICMI. |
| **G4beamline** | `reference/G4beamline Documentation/README.md` | Geant4-based beamline simulation — command-driven input file, full physics lists, virtual detectors, NTuples, and 3D visualization. |
| **BMAD** | `reference/BMAD Documentation/README.md` | Fortran90 subroutine library for reading MAD-format lattice files, computing Twiss parameters, and tracking particles. Developed at Cornell (CESR/CLASSE). Supports Taylor maps, Runge-Kutta, symplectic integrators, and PTC interface. MAD-X User Manual at `reference/BMAD Documentation/MAD-X User Manual/README.md`. |
| **Linac Sim GUI** | `reference/Linac Simulation Documentation/README.md` | Adam Bartnik's CESR Linac simulation GUI (Java). Chains a custom 1D cathode code → GPT (space charge, cylindrical symmetry) → BMAD (high-energy, 3D). Includes fieldmaps for the thermionic gun, prebunchers, solenoid lenses, and SLAC-design linac cavities. |
| **LUME-Impact** | `reference/lume-impact Documentation/README.md` | Python interface for IMPACT-T and IMPACT-Z. Provides `Impact` and `ImpactZ` classes for input configuration, execution, output parsing, and plotting. Integrates with openPMD-beamphysics and BMAD. |
| **openPMD-beamphysics** | `reference/openPMD-beamphysics Documentation/README.md` | Python tools for particle/field data in the openPMD beamphysics standard. `ParticleGroup` (particle data) and `FieldMesh` (field maps, e.g. `from_onaxis` for the DC gun field). |
| **openPMD-viewer** | `reference/openPMD-viewer Documentation/README.md` | Python API + Jupyter GUI for reading/visualizing openPMD file series. `OpenPMDTimeSeries.get_field`/`get_particle`. Used to read WarpX diagnostics in `cathode/plot_cathode.py`. |
| **easygdf** | `reference/easygdf Documentation/README.md` | Pure-Python reader/writer for GPT's GDF binary format. `load`/`save` for raw blocks; `load_screens_touts`/`save_screens_touts` for GPT output; `load`/`save_initial_distribution` for GPT input distributions. |

### Key Concepts

- **PIC (Particle-In-Cell)**: The computational method used by IMPACT-Z and WarpX — particles tracked on a mesh, fields solved on grid.
- **Space charge**: Coulomb self-repulsion of the beam, dominant at low energy and high current; all four codes model it.
- **CSR (Coherent Synchrotron Radiation)**: Wakefield from relativistic bunches in bending magnets; modeled in IMPACT-T.
- **GDF**: GPT's native binary data format; convert to/from ASCII with `GDF2A`/`ASCI2GDF`.
- **IMPACT-Z output files**: Named `fort.18`, `fort.24`–`fort.30`, `fort.32` — see `reference/Impact-Z Documentation/output_files/` for column definitions.

### Adding New Papers

When saving a new paper to `reference/Papers/`, add a summary entry to `reference/Papers/README.md` following the existing format (title, file, author, abstract summary).
