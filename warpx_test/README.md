# WarpX Positron Beam Simulations

Introductory demos of [WarpX](https://warpx.readthedocs.io) for positron beam physics, using the Python/PICMI interface (`pywarpx`).

Run any script with:
```bash
conda activate CBB
python warpx_test/<script>.py
```

---

## Scripts

### `01_single_positron.py` — Single relativistic positron
Traces a single 100 MeV positron through free space for 50 time steps.

- Uses `ParticleListDistribution` to place one particle at the origin with momentum along z
- Confirms WarpX is working: checks γ ≈ 197, uz ≈ 5.9 × 10¹⁰ m/s
- Writes full particle state (position + momentum) at every step to `diags_single/`

### `02_positron_bunch.py` — Gaussian bunch, negligible space charge
Propagates a 5 MeV, **1 pC** Gaussian positron bunch for 500 steps (~12 cm).

- With only 1 pC the Coulomb self-field is ~1000× smaller than the 1 nC case
- σ_x stays constant → baseline reference for space-charge comparison
- Writes beam statistics (RMS size, emittance, etc.) to `diags_bunch/beam_stats.txt`

### `03_positron_space_charge.py` — Gaussian bunch, strong space charge
Identical setup to Script 02 except charge = **1 nC**.

- At 5 MeV (γ ≈ 11) the relativistic suppression 1/γ² ≈ 1/120 is small enough
  that the Coulomb repulsion drives **~7% transverse growth** in 500 steps
- Uses the FFT-based (Integrated Green Function) relativistic Poisson solver
- Writes beam statistics to `diags_space_charge/beam_stats_sc.txt`

### `plot_results.py` — Comparison figure
Reads the BeamRelevant CSV outputs from Scripts 02 and 03 and produces
`results/beam_comparison.png`: σ_x and transverse emittance vs. simulation step.

```bash
python warpx_test/plot_results.py
```

---

## Key physics

| Parameter | Value |
|-----------|-------|
| Beam energy | 5 MeV (γ ≈ 10.78, β ≈ 0.994) |
| σ_x = σ_y (initial) | 0.5 mm |
| σ_z (initial) | 1 mm |
| Low-charge reference | 1 pC |
| High-charge (space charge) | 1 nC |
| Grid | 32 × 32 × 32 cells, ±4 mm, dx = 0.25 mm |
| Time step | ≈ 0.83 ps |
| Steps | 500 |
| Beam travel | ≈ 12 cm |

**Why 5 MeV?** At GeV energies, space charge is suppressed by 1/γ² ≈ 1/4 × 10⁶ —
completely negligible. At 5 MeV the suppression is only ~120×, so a 1 nC bunch produces
visible beam blow-up within a short simulation run. This energy range is also relevant to
Cornell's ERL photoinjector, where space charge is a dominant concern.

**Why electrostatic solver?** The standard Yee FDTD (electromagnetic) solver suffers from
the **Numerical Cherenkov Instability** (NCI) for relativistic beams (γ >> 1), causing
unphysical exponential energy growth. The FFT-based Poisson solver avoids this entirely
and is the correct tool for space-charge studies.

---

## Output structure

```
warpx_test/
├── diags_single/          # per-step particle plotfiles from Script 01
├── diags_bunch/
│   ├── beam_stats.txt     # BeamRelevant CSV: RMS size, emittance, etc. (Script 02)
│   └── part*/             # particle plotfiles every 50 steps
├── diags_space_charge/
│   ├── beam_stats_sc.txt  # same diagnostic for the 1 nC run (Script 03)
│   └── part*/
└── results/
    └── beam_comparison.png
```

Plotfiles can be read with [`yt`](https://yt-project.org) or
[`openpmd-viewer`](https://github.com/openPMD/openPMD-viewer):

```python
import yt
ds = yt.load("warpx_test/diags_single/part000025")
```

BeamRelevant CSVs are plain whitespace-separated text (columns documented in `plot_results.py`).

---

## References

- WarpX docs: `reference/WarpX Documentation/README.md`
- Gaussian beam example: `reference/WarpX Documentation/usage/examples/gaussian_beam/README.md`
- Cornell ERL injector paper: `reference/Papers/PhysRevSTAB.16.073401.pdf`
