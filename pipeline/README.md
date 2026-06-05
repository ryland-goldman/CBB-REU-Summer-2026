# End-to-end pipeline

`run_pipeline.py` runs the full Cornell Linac beam simulation in WarpX, in order, by importing
each stage's top-level facade module and calling `.run()`:

```
cathode.run()      # SCL emission (2D Child–Langmuir diode) + plots
gun.run()          # build gun field map + RZ acceleration (~148 keV) + plots
prebuncher.run()   # build prebuncher field map + RZ RF cavity + plots
```

Each `run()` reads the previous stage's openPMD output, runs the WarpX sim, and (by default)
generates that stage's figures.

Each downstream stage reads the previous stage's openPMD output, so the order is fixed.

Within a single `.run()`, the field-map build and the plots happen in-process, but the WarpX
sim itself is spawned in a **fresh Python subprocess** (`pipeline._launch_sim`). This sidesteps
pywarpx's per-process geometry binding — pywarpx binds globally to one geometry (2D/RZ/3D) at
first `.so` load and caches diagnostic state per name, so chaining cathode (2D) → gun (RZ) →
prebuncher (RZ) in one interpreter would trip `AssertionError: Diagnostic attributes not
consistent`. The tqdm progress bar is driven by a WarpX `afterstep` callback, with WarpX's
per-step stdout redirected to the log file so the bar stays on a clean terminal line.

## Run

```
conda activate CBB
python pipeline/run_pipeline.py
```

## Configuration

`run_pipeline.py` is just three `.run()` calls. To change behaviour, either:

- **Comment out** the stages you don't want to re-run (each downstream stage still reads the
  previous stage's saved openPMD output from disk).
- **Override physics inputs** at the top of `run_pipeline.py` using each stage's `config()`:

  ```python
  cathode.config(V_anode=50.0, gap_d=100e-6)
  gun.config(GUN_VOLTAGE=150e3, BUNCH_CHARGE=0.1e-9)
  prebuncher.config(POWER_W=800, PHASE="zc",
                    OUTDIR="prebuncher/diags/P800_zc")
  ```

  Keys must match the module-level constants in each `<stage>/*.py` — `config()` writes them
  through with `setattr`.
- **OMP threads:** set the `OMP_THREADS` environment variable (default 6; the MLMG solve is
  bandwidth-bound, so all cores is slower). The shim sets `OMP_NUM_THREADS` before any pywarpx
  import.

## Output

**Terminal:** a live tqdm progress bar per simulation step, the key prints from each component
(beam energy, bunch charge, field magnitudes, …), per-stage timing, and a final-beam summary
(⟨z⟩, σ_z, ⟨KE⟩, σ_KE, charge) read from the prebuncher exit.

**Log file** (`pipeline/logs/pipeline_<timestamp>.log`, path printed at start/end): the same
information with the progress bar replaced by periodic `step N/total (%) — elapsed / rate / ETA`
lines, per-stage durations and `ok=true/false`, and DEBUG-level WarpX output (the full text the
progress bar replaces on the terminal; MLMG per-iteration convergence spam is dropped). ANSI
colours are written only to an interactive terminal, so the log stays plain text.

Figures land in `{cathode,gun,prebuncher}/results/`.
