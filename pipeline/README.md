# End-to-end pipeline

`run_pipeline.py` runs the full Cornell Linac beam simulation in WarpX, in order, by importing
each stage's top-level facade module and calling `.run()`:

```
cathode.run()      # SCL emission (2D Child–Langmuir diode) + plots
gun.run()          # build gun field map + RZ acceleration (~148 keV) + plots
prebuncher.run()   # build prebuncher field map + RZ RF cavity + plots
linac_sec1.run()   # build linac field maps + RZ SLAC TW section (~37 MeV) + plots
```

Each `run()` reads the previous stage's openPMD output, runs the WarpX sim, and (by default)
generates that stage's figures.

Each downstream stage reads the previous stage's openPMD output, so the order is fixed.

Within a single `.run()`, the field-map build and the plots happen in-process, but the WarpX
sim itself is spawned in a **fresh Python subprocess** (`pipeline._launch_sim`). This sidesteps
pywarpx's per-process geometry binding — pywarpx binds globally to one geometry (2D/RZ/3D) at
first `.so` load and caches diagnostic state per name, so chaining cathode (2D) → gun (RZ) →
prebuncher (RZ) → linac_sec1 (RZ) in one interpreter would trip `AssertionError: Diagnostic
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
  cathode.config(V_anode=50.0, gap_d=100e-6)
  gun.config(GUN_VOLTAGE=150e3, BUNCH_CHARGE=0.1e-9)
  prebuncher.config(POWER_W=800, PHASE="zc",
                    OUTDIR="prebuncher/diags/P800_zc")
  linac_sec1.config(POWER_MW=15.0)
  ```

  Keys must match the module-level constants in each `<stage>/*.py` — `config()` writes them
  through with `setattr`.
- **Performance knobs (accuracy ↔ speed):** `run_pipeline.py` has a **PERFORMANCE KNOBS**
  block — a *Balanced* profile (~1.7× faster, ~5 min) active by default, plus commented
  *Conservative* (~1.3×) and *Aggressive* (~2.2×) presets. Measured runtime split is cathode 7%,
  gun 17%, **prebuncher 75%** (its self-field MLMG solve dominates). The two RZ stages differ
  (measured): the **gun** has near-isotropic cells and scales ≈ `nz²` (smaller `nz` ⇒ fewer cells
  *and*, via `dz=zmax/nz`, fewer steps — halving `nz` ≈ 4× faster). The **prebuncher** is the
  opposite — its long-thin box has anisotropic cells and a convergence-bound MLMG solve, so
  coarsening `NZ` *slows* the per-step solve and under-resolves the ~1 mm bunch; keep `NZ=1024`
  and speed it via `CFL` (fewer steps) and `MAX_ITERS`/`REQUIRED_PRECISION` (cheaper solve). The
  prebuncher's space-charge solve over the full 1.3 m transit is the irreducible runtime floor.
  The **linac_sec1** stage is a single ~3.5 m on-crest RZ run (≈40 s at the default `NZ=1664`);
  like the prebuncher it is a long-thin box, so its cells must stay ≈3:1 (`NR=16`) or the MLMG
  self-field solve diverges (`MLMG failed`).
  Per-stage knobs (all `config()`-overridable module constants): grid (`nx/nz`, `nr/nz`, `NR/NZ`),
  run length (`MAX_STEPS`, `CFL`, `TRANSIT_MARGIN`, `AVG_SPEED_FRAC`), Poisson solve
  (`REQUIRED_PRECISION`, `MAX_ITERS`), macroparticles (`PPC` for the cathode, `MAX_PART` downsample
  for gun/prebuncher), and diagnostic dumps (`N_DIAGS`, `DIAG_PERIOD`). Comment the Balanced block
  to restore the exact original (8.8-min) physics.
- **OMP threads:** set the `OMP_THREADS` environment variable (default 6; the MLMG solve is
  bandwidth-bound, so all cores is slower). The shim sets `OMP_NUM_THREADS` before any pywarpx
  import.

## Output

**Terminal:** a live tqdm progress bar per simulation step, the key prints from each component
(beam energy, bunch charge, field magnitudes, …), per-stage timing, and two final-beam summaries
(⟨z⟩, σ_z, ⟨KE⟩, σ_KE, charge) — one read from the prebuncher exit (in keV) and one from the
linac_sec1 exit (in MeV, with the captured-charge fraction).

**Log file** (`pipeline/logs/pipeline_<timestamp>.log`, path printed at start/end): the same
information with the progress bar replaced by periodic `step N/total (%) — elapsed / rate / ETA`
lines, per-stage durations and `ok=true/false`, and DEBUG-level WarpX output (the full text the
progress bar replaces on the terminal; MLMG per-iteration convergence spam is dropped). ANSI
colours are written only to an interactive terminal, so the log stays plain text.

Figures land in `{cathode,gun,prebuncher,linac_sec1}/results/`.
