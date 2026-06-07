# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

The code in this folder is meant to contain simulations for a Research Experience for Undergraduates (REU) program at the Cornell Center for Bright Beams (CBB) and Cornell Laboratory for Accelerator ScienceS and Education (CLASSE). The project focus is to build a beam simulation for the Cornell High Energy Synchotron Source (CHESS).

Concretely, the repo rebuilds the front end of the Cornell Linac (Adam Bartnik's LinacSim: thermionic source → gun → injector → linac) from first principles in **WarpX** via its Python/PICMI interface (`pywarpx`). The four stages form one self-consistent chain — each reads the previous stage's openPMD beam as input:

```
cathode  ─►  gun  ─►  injector  ─►  linac_sec1
(SCL diode)  (~146 keV)  (2 prebunchers + 3 solenoids,  (~16 MeV captured, SLAC TW section)
              velocity bunching + focusing, ~2 m)
```

The `injector` stage is the full LinacSim injector subsection in one self-consistent RZ
space-charge run — Lens 0A → Prebuncher 1 → Prebuncher 2 (reversed) → Sol 0 / Lens 0E —
handing a focused, velocity-bunched beam through the 9.547 mm collimator to `linac_sec1`
at the true linac entrance z ≈ 2.03 m. (It replaced the earlier single-cavity `prebuncher/`
stage.)

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

- **Stage API:** each stage package (`cathode/`, `gun/`, `injector/`, `linac_sec1/`) is a top-level facade. `import cathode; cathode.run()` (likewise for `gun`, `injector`, `linac_sec1`) builds the field map (if any), runs the WarpX sim, and generates that stage's plots. Each exposes `config(**kwargs)` (override module-level constants — keys must match the names in each `<stage>/*.py`), `run(plots=True)`, and `plot()` (figures from existing diags). Implementation lives in `pipeline/_runner.py`.
- **Performance knobs:** the runtime-critical parameters are module-level constants in each `<stage>/*.py`, so they are tunable via `config()` from `run_pipeline.py` (which ships an editable "PERFORMANCE KNOBS" block: a Balanced profile active by default plus commented Conservative/Aggressive presets). Per stage: grid (`nx/nz`, `nr/nz`, `NR/NZ`), step count (`MAX_STEPS`/`CFL`/`TRANSIT_MARGIN`/`AVG_SPEED_FRAC`), Poisson solve (`REQUIRED_PRECISION`, `MAX_ITERS`), macroparticles (`PPC` for the cathode, `MAX_PART` downsample cap for the gun and injector), and diagnostic dump count (`N_DIAGS`, `DIAG_PERIOD`). The **injector dominates total runtime** (its self-field MLMG solve over the ~2 m long-thin box). The gun's cells are near-isotropic, so it scales ≈ `nz²` (per-step cost ∝ cells, and `dz = zmax/nz` ⇒ fewer steps as `nz` drops) — coarsening `nz` is a big clean win. The **injector is the opposite** (measured): its long-thin box (2.10 m × 36 mm) has anisotropic cells and a convergence-bound MLMG solve, so coarsening `NZ` *slows* the per-step solve faster than it removes cells and under-resolves the ~1 mm bunch — keep `NZ=1664` (dz≈1.26 mm ⇒ 2.80:1 aspect at `NR=80`) and speed it via `CFL` (fewer steps) and `MAX_ITERS`/`REQUIRED_PRECISION` (cheaper solve) instead. The injector run is **convergence-bound, so its cost over the longer 2.10 m box rose >2×** vs the old 1.30 m prebuncher (~60 s vs ~24 s); the linac got cheaper (dropped the solenoid + heavy radial scrape), partly offsetting. Lowering knobs trades accuracy for speed. The `linac_sec1` stage is a single ~3.5 m RZ run (≈45 s at the default `NZ=1664`, `NR=16` over the 9.547 mm bore → ~3.5:1 cells); like the injector it is a long-thin box, so its cells must stay near the ≈3:1 rule or the MLMG self-field solve diverges (`MLMG failed`).
- **Run one stage off existing upstream output:** comment out the unwanted `<stage>.run()` calls in `pipeline/run_pipeline.py`, or just import the one you want (`import injector; injector.run()`) — each stage reads the previous stage's openPMD output from disk, so any unmodified upstream output is reused.
- **Injector power/phase scan:** call `injector.run()` once per operating point in a Python loop with an explicit `OUTDIR` per case, e.g. `for p in (160, 300, 500, 800): injector.config(PREB1_KW=p, OUTDIR=f"injector/diags/P{p}_zc"); injector.run(plots=False)`. A single `injector.plot()` then aggregates every `diags/P*` directory alongside the default `diags/main` (see `injector/README.md`). The faithful default is the 8 kW / 10 kW two-cavity point — the scan is exploration only. **Caveat:** a hard Preb-1 power scan desyncs the Preb-2 phase reference (Preb-2 is timed from the injection β + Preb-1's faithful kick); a hardened-Preb-1 study needs a two-pass run (see `injector/README.md`).
- **Plots:** `<stage>.plot()` reads its `diags/` and writes PNGs to `results/`. `run()` calls `plot()` by default; pass `plots=False` to skip.
- **Threads:** `OMP_THREADS` (**default 1 — keep this pipeline single-threaded**). These stages run fastest on a single thread: the grids are small and the MLMG Poisson solve is memory-bandwidth bound, so OpenMP threads contend for the same memory bus and add fork/join + barrier overhead without speeding the solve (measured on an M4 Pro: full Balanced chain ~1.1 min at `OMP_THREADS=1`; `OMP_THREADS=14` ran at only ~450% CPU with no gain, and MPI was far worse — see [[project-warpx-mpi-build]]). Only raise `OMP_THREADS` for the much larger *original-config* grids, where per-thread work outgrows the overhead. Override via the `OMP_THREADS` env var (set before any pywarpx import; `config()` cannot set it).

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
| `injector/injector_sim.py` | `gun/diags/particles` + `injector/injector_field/{preb1_EB,preb2_EB,lens0a,sol0,lens0e}.h5` | `injector/diags/main` |
| `linac_sec1/linac_sec1_sim.py` | `injector/diags/main/particles` (snapshot nearest the z≈2.03 m handoff) + `linac_sec1/linac_sec1_field/linac_{rf1,rf2}.h5` | `linac_sec1/diags/main` |

`pipeline/run_pipeline.py` orchestrates the whole chain by calling `cathode.run()`, `gun.run()`, `injector.run()`, `linac_sec1.run()` in order, then `pipeline.plot_chain.main()` for the cross-stage figures (`results/`). The shared runner in `pipeline/_runner.py` builds field maps and generates plots **in-process**, but spawns a **fresh Python subprocess** (`pipeline/_launch_sim.py`) for each sim — pywarpx binds globally to one geometry (2D/RZ/3D) at first `.so` load and caches diagnostic state by name, so chaining cathode (2D) → gun (RZ) → injector (RZ) in one interpreter would trip `AssertionError: Diagnostic attributes not consistent`. Inside each sim, `run_step(...)` installs a `pywarpx.callbacks.installcallback("afterstep", …)` hook to drive a tqdm progress bar and redirects WarpX's noisy per-step stdout to the pipeline log file, so the bar updates on a clean terminal line. A structured DEBUG log lands in `pipeline/logs/pipeline_<timestamp>.log`. The cathode is 2D x–z; the gun, injector, and linac_sec1 are RZ (cylindrically symmetric).

- The default chain writes `injector/diags/main`; the power/phase scan facility survives as an optional per-case `OUTDIR` override (`injector/diags/<case>`). `injector.resolve_outdir()` returns `injector/diags/main` by default.
- `linac_sec1` no longer builds `linac_sol.h5` or carries `I_SOL`/`SOL_Z`/`SOL_MAP` — transverse focusing moved upstream into the injector (the three real lenses at their true lab z). Its only applied fields are the two SLAC quadrature RF maps; its `RMAX=9.547 mm` is the SLAC bore / injector→linac collimator iris, and the `r ≤ RMAX` injection cut IS the physical 9.547 mm collimation.

Field maps (`CESR_gun.gdf`, `prebuncher_25D.gdf`, `SLAC-3mLinac-field1/field2.gdf`, `SOL_0.gdf`, `LENS_0A…0E.gdf`) live in `fieldmaps/`; field-map paths are set near the top of each `build_*_field.py`. (The two `SLAC-3mLinac` files are the quadrature Re/Im halves of **one** traveling-wave section, summed at 90° — not two sections. The injector builds **both** prebunchers from the *one* `prebuncher_25D.gdf`: `preb1_EB.h5` at z=0.534 m and `preb2_EB.h5` at z=1.318 m — same forward field, different `grid_global_offset`; the reversed Preb-2 install is encoded as a run-time phase, not a negated map.) WarpX-specific gotchas accumulated per stage (negative field-scale conventions, thetaMode openPMD axis order, charge renormalization, where to stop the run to avoid MLMG aborts, the ≈3:1 cell-aspect requirement for the long-thin boxes, that picmi `LoadAppliedField` disables the global E init-style for any `load_E=False` field so B-only maps must be added *before* the accelerating E maps, the reversed-cavity `-1,0,0` install ≡ +π in **absolute** phase but absorbed by crest-referencing the loaded field so the applied `PREB2_REV_PHASE=0`, and that this pywarpx RZ build's particle-position SoA accessors raise "Component x does not exist" so an in-run radial scrape isn't available — the 9.547 mm collimator is a post-hoc radial cut at the handoff instead) are documented in each stage's `README.md` — read it before modifying a stage.

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
