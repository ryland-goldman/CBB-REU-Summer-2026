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

All commands run from the **repo root** in the `CBB` environment. Stage scripts use hard-coded relative paths (e.g. `warpx_gun/diags`), so running from anywhere else breaks the inter-stage handoff.

```bash
conda activate CBB
pip install -r requirements.txt                 # pywarpx/openpmd-api are best via mamba

python pipeline/run_pipeline.py                  # full chain, live progress + final-beam summary
```

- **Run one stage off existing upstream output:** toggle `RUN_CATHODE / RUN_GUN / RUN_PREBUNCHER / MAKE_PLOTS` in the `CONFIG` block at the top of `pipeline/run_pipeline.py` (e.g. set cathode+gun `False` to re-run only the prebuncher against the saved gun beam). `PREBUNCHER_POWER_W` / `PREBUNCHER_PHASE` set the prebuncher operating point.
- **Run a stage directly:** `python warpx_gun/gun_sim.py`, etc. The prebuncher takes CLI args: `python warpx_prebuncher/prebuncher_sim.py --power 800 --phase zc --outdir warpx_prebuncher/diags/P800_zc` (`--phase` is `zc` = zero-crossing bunching or `crest` = max energy gain; `--power 0` = drift-only baseline).
- **Prebuncher power/phase scan:** `python warpx_prebuncher/run_scan.py`.
- **Plots:** each stage has a `plot_*.py` that reads its `diags/` and writes PNGs to `results/`.
- **Threads:** `OMP_THREADS` (default 6) — the MLMG Poisson solve is memory-bandwidth bound, so using all cores is *slower*. Override via env var or the `CONFIG` block.

There is no test suite, linter, or build step — validation is physics sanity checks (energy gain, Child–Langmuir current, bunching) printed by each run and inspected in the `results/` plots.

## Project Architecture

Each stage lives in its own `warpx_<stage>/` directory and follows the same script layout:

- `build_*_field.py` — converts a GPT `.gdf` field map from `fieldmaps/` into an openPMD `.h5` field mesh (via `easygdf` + `openPMD-beamphysics`) that WarpX loads as an external field. (The cathode has no field map; its field is self-consistent.)
- `*_sim.py` — the WarpX/PICMI run. Reads the upstream beam with `openPMD-viewer`, injects it, tracks through the stage, writes openPMD particle diagnostics to its own `diags/`.
- `plot_*.py` — reads `diags/`, writes figures to `results/`.
- `README.md` — the stage's physics, field map, operating point, and outputs.

**Inter-stage contract (the chain is order-dependent):**

| Stage | Reads | Writes |
|-------|-------|--------|
| `warpx_cathode/cathode_diode.py` | — (emits at cathode) | `warpx_cathode/diags/particles` |
| `warpx_gun/gun_sim.py` | `warpx_cathode/diags/particles` + `warpx_gun/gun_field/gun_E.h5` | `warpx_gun/diags` |
| `warpx_prebuncher/prebuncher_sim.py` | `warpx_gun/diags/particles` + `warpx_prebuncher/prebuncher_field/prebuncher_EB.h5` | `warpx_prebuncher/diags/<P..._...>` |

`pipeline/run_pipeline.py` orchestrates the whole chain as a sequence of subprocesses (cd's to repo root, runs each `build`/`sim`/`plot` in order), surfacing key physics lines and a live progress bar to the terminal while writing a full DEBUG log to `pipeline/logs/pipeline_<timestamp>.log`. The cathode is 2D x–z; the gun and prebuncher are RZ (cylindrically symmetric).

Field maps (`CESR_gun.gdf`, `prebuncher_25D.gdf`) live in `fieldmaps/`; field-map paths are set near the top of each `build_*_field.py`. WarpX-specific gotchas accumulated per stage (negative field-scale conventions, thetaMode openPMD axis order, charge renormalization, where to stop the run to avoid MLMG aborts) are documented in each stage's `README.md` — read it before modifying a stage.

**Conventions:**

- Simulation outputs are git-ignored (`diags/`, `results/`, `*.h5`, `*.gdf`, logs); regenerate by re-running. Field maps in `fieldmaps/` are committed.
- Commit convention (matches existing history): for a stage, commit its `*.py` scripts + `README.md`, and `git add -f warpx_<stage>/results/*.png` to include the result figures (since `results/` is git-ignored). Do **not** commit `diags/`, `.h5`, or logs.

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
| **openPMD-viewer** | `reference/openPMD-viewer Documentation/README.md` | Python API + Jupyter GUI for reading/visualizing openPMD file series. `OpenPMDTimeSeries.get_field`/`get_particle`. Used to read WarpX diagnostics in `warpx_cathode/plot_cathode.py`. |
| **easygdf** | `reference/easygdf Documentation/README.md` | Pure-Python reader/writer for GPT's GDF binary format. `load`/`save` for raw blocks; `load_screens_touts`/`save_screens_touts` for GPT output; `load`/`save_initial_distribution` for GPT input distributions. |

### Key Concepts

- **PIC (Particle-In-Cell)**: The computational method used by IMPACT-Z and WarpX — particles tracked on a mesh, fields solved on grid.
- **Space charge**: Coulomb self-repulsion of the beam, dominant at low energy and high current; all four codes model it.
- **CSR (Coherent Synchrotron Radiation)**: Wakefield from relativistic bunches in bending magnets; modeled in IMPACT-T.
- **GDF**: GPT's native binary data format; convert to/from ASCII with `GDF2A`/`ASCI2GDF`.
- **IMPACT-Z output files**: Named `fort.18`, `fort.24`–`fort.30`, `fort.32` — see `reference/Impact-Z Documentation/output_files/` for column definitions.

### Adding New Papers

When saving a new paper to `reference/Papers/`, add a summary entry to `reference/Papers/README.md` following the existing format (title, file, author, abstract summary).
