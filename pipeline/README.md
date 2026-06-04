# End-to-end pipeline

`run_pipeline.py` runs the full Cornell Linac beam simulation in WarpX, in order, as a chain of
subprocesses with a live progress bar and the key physics output for each component:

```
1. Cathode diode        warpx_cathode/cathode_diode.py        (space-charge-limited emission)
2. Gun field map        warpx_gun/build_gun_field.py          (CESR_gun.gdf -> openPMD)
3. Gun acceleration     warpx_gun/gun_sim.py                  (RZ, -> ~148 keV)
4. Prebuncher field map warpx_prebuncher/build_prebuncher_field.py
5. Prebuncher bunching  warpx_prebuncher/prebuncher_sim.py    (RZ RF cavity)
6. Plots                each stage's plot_*.py
```

Each downstream stage reads the previous stage's openPMD output, so the order is fixed.

## Run

```
conda activate CBB
python pipeline/run_pipeline.py
```

## Configuration

Edit the `CONFIG` block at the top of `run_pipeline.py`:

- `RUN_CATHODE` / `RUN_GUN` / `RUN_PREBUNCHER` / `MAKE_PLOTS` — toggle stage groups (e.g. set
  `RUN_CATHODE = RUN_GUN = False` to re-run only the prebuncher off existing gun output).
- `PREBUNCHER_POWER_W`, `PREBUNCHER_PHASE` — the prebuncher operating point (passed to
  `prebuncher_sim.py`; `0 W` = drift-only baseline).
- `OMP_THREADS` — OpenMP threads per WarpX run (default 6; the MLMG solve is bandwidth-bound, so
  all cores is slower). Also overridable via the `OMP_THREADS` environment variable.

## Output

**Terminal:** a live progress bar per simulation stage, the key prints from each component
(beam energy, bunch charge, field magnitudes, …), per-stage timing, and a final-beam summary
(⟨z⟩, σ_z, ⟨KE⟩, σ_KE, charge) read from the prebuncher exit.

**Log file** (`pipeline/logs/pipeline_<timestamp>.log`, path printed at start/end): the same
information with the progress bar replaced by periodic `step N/total (%) — elapsed / rate / ETA`
lines, plus per-stage return codes and durations and extra `DEBUG` detail — the full command and
environment, WarpX init, per-step timing, and the run's profiler table. (MLMG per-iteration
convergence spam is dropped.) ANSI colours are written only to an interactive terminal, so the
log stays plain text.

Figures land in `warpx_{cathode,gun,prebuncher}/results/`.
