# End-to-end pipeline

`run_pipeline.py` runs the full Cornell Linac beam simulation in WarpX, in order, by importing
each stage's top-level facade module and calling `.run()`:

```
cathode.run()      # SCL emission (2D Child–Langmuir diode) + plots
gun.run()          # build gun field map + RZ acceleration (~146 keV) + plots
injector.run()     # build injector fields (2 cavities + 3 solenoids) + RZ run + plots
linac_sec1.run()   # build linac field maps + RZ SLAC TW section (~16 MeV captured) + plots
pipeline.plot_chain.main()   # cross-stage figures into the repo-root results/
```

Each `run()` reads the previous stage's openPMD output, runs the WarpX sim, and (by default)
generates that stage's figures. After the four stages, `run_pipeline.main()` calls
`pipeline.plot_chain.main()` for the cross-stage figures (also exposed as `pipeline.plot_chain()`).

Each downstream stage reads the previous stage's openPMD output, so the order is fixed. The
**injector → linac handoff plane is z ≈ 2.03 m** (Z_acc_1): the linac selects the injector dump
whose ⟨z⟩ is nearest 2.03 m and applies the 9.547 mm collimator cut at injection.

Within a single `.run()`, the field-map build and the plots happen in-process, but the WarpX
sim itself is spawned in a **fresh Python subprocess** (`pipeline._launch_sim`). This sidesteps
pywarpx's per-process geometry binding — pywarpx binds globally to one geometry (2D/RZ/3D) at
first `.so` load and caches diagnostic state per name, so chaining cathode (2D) → gun (RZ) →
injector (RZ) → linac_sec1 (RZ) in one interpreter would trip `AssertionError: Diagnostic
attributes not consistent`. The tqdm progress bar is driven by a WarpX `afterstep` callback, with WarpX's
per-step stdout redirected to the log file so the bar stays on a clean terminal line.

## Run

```
conda activate CBB
python pipeline/run_pipeline.py
```

## Configuration

`run_pipeline.py` is four `.run()` calls plus a `config()` header. To change behaviour, either:

- **Comment out** the stages you don't want to re-run (each downstream stage still reads the
  previous stage's saved openPMD output from disk).
- **Override physics inputs** at the top of `run_pipeline.py` using each stage's `config()`:

  ```python
  cathode.config(V_anode=60.0, gap_d=200e-6)
  gun.config(GUN_VOLTAGE=150e3, BUNCH_CHARGE=1.0e-9)
  injector.config(PREB1_KW=8, PREB2_KW=10, PHASE="crest")   # default OUTDIR -> diags/main
  linac_sec1.config(POWER_MW=11.0)                            # no I_SOL — focusing is upstream
  ```

  Keys must match the module-level constants in each `<stage>/*.py` — `config()` writes them
  through with `setattr`.
- **Performance knobs (accuracy ↔ speed):** `run_pipeline.py` has a **PERFORMANCE KNOBS**
  block — a *Balanced* profile active by default, plus commented *Conservative* and *Aggressive*
  presets. The **injector dominates total runtime** (its self-field MLMG solve over the ~2 m
  long-thin box; ~60 s, >2× the old 1.30 m prebuncher because the convergence-bound solve over
  the longer box rises super-linearly). The two RZ-cavity stages differ (measured): the **gun**
  has near-isotropic cells and scales ≈ `nz²` (smaller `nz` ⇒ fewer cells *and*, via `dz=zmax/nz`,
  fewer steps — halving `nz` ≈ 4× faster). The **injector** is the opposite — its long-thin box
  has anisotropic cells and a convergence-bound MLMG solve, so coarsening `NZ` *slows* the
  per-step solve and under-resolves the ~1 mm bunch; keep `NZ=1664` (dz≈1.26 mm ⇒ 2.80:1 at
  `NR=80`) and speed it via `CFL` (fewer steps) and `MAX_ITERS`/`REQUIRED_PRECISION` (cheaper
  solve). The **linac_sec1** stage is a single ~3.5 m on-crest RZ run (≈45 s at the default
  `NZ=1664`, `NR=16` over the 9.547 mm bore); like the injector it is a long-thin box, so its
  cells must stay near the ≈3:1 rule or the MLMG self-field solve diverges (`MLMG failed`).
  Per-stage knobs (all `config()`-overridable module constants): grid (`nx/nz`, `nr/nz`, `NR/NZ`),
  run length (`MAX_STEPS`, `CFL`, `TRANSIT_MARGIN`, `AVG_SPEED_FRAC`), Poisson solve
  (`REQUIRED_PRECISION`, `MAX_ITERS`), macroparticles (`PPC` for the cathode, `MAX_PART` downsample
  for gun/injector), and diagnostic dumps (`N_DIAGS`, `DIAG_PERIOD`).
- **OMP threads:** set the `OMP_THREADS` environment variable (default 1 — keep this pipeline
  single-threaded). The grids are small and the MLMG solve is memory-bandwidth-bound, so threads
  add fork/join + barrier overhead with no gain; only raise `OMP_THREADS` for the much larger
  original-config grids. The shim sets `OMP_NUM_THREADS` before any pywarpx import.

## Output

**Terminal:** a live tqdm progress bar per simulation step, the key prints from each component
(beam energy, bunch charge, field magnitudes, …), per-stage timing, two final-beam summaries
(⟨z⟩, σ_z, ⟨KE⟩, σ_KE, charge) — one read from the injector exit (in keV) and one from the
linac_sec1 exit (in MeV, with the captured-charge fraction vs the true injected charge) — and the
cross-stage scorecard table (per-stage entry/exit moments) from `plot_chain`.

**Log file** (`pipeline/logs/pipeline_<timestamp>.log`, path printed at start/end): the same
information with the progress bar replaced by periodic `step N/total (%) — elapsed / rate / ETA`
lines, per-stage durations and `ok=true/false`, and DEBUG-level WarpX output (the full text the
progress bar replaces on the terminal; MLMG per-iteration convergence spam is dropped). ANSI
colours are written only to an interactive terminal, so the log stays plain text.

Per-stage figures land in `{cathode,gun,injector,linac_sec1}/results/`. The **cross-stage**
figures land in the repo-root `results/` (from `pipeline.plot_chain`): `chain_evolution.png`
(3×2 panels of ⟨KE⟩, ε_n,x, σ_x/σ_r, σ_z, charge fraction, I_peak vs lab ⟨z⟩),
`emittance_budget.png`, `transmission_waterfall.png`, and `chain_scorecard.png`. All are
git-ignored; commit with `git add -f results/*.png`. `plot_chain` reads each stage's existing
openPMD series and works even if only some stages have run.

## Cross-stage figures (`pipeline/plot_chain.py`)

`pipeline.plot_chain()` (a thin wrapper for `pipeline.plot_chain.main()`) builds ONE per-dump
moment table per stage and renders four figures — all views of that table. It is in-process
(no pywarpx) and reads `cathode/diags/particles`, `gun/diags/particles`,
`injector/diags/main/particles`, `linac_sec1/diags/main/particles`. Notes: the cathode is 2D
(x–z; only x/ux requested) so its ε_n,x is the slab x-emittance — the cathode→gun ε_n step is a
2D→RZ **definitional** discontinuity, annotated as such. Longitudinal ε_n,z is the z–(γβ_z)
emittance in mm (NOT mm·mrad). Capture is reported vs the TRUE injected charge
(`linac_sec1/diags/main/injection_summary.json`), and the σ_r / capture panels carry the
γ²≈1.7× ES-self-field conservative-lower-bound caveat.
