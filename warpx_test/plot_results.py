"""
Plot beam size evolution: free space vs. space charge.

Reads BeamRelevant reduced diagnostic files produced by:
  - 02_positron_bunch.py        → warpx_test/diags_bunch/beam_stats.txt
  - 03_positron_space_charge.py → warpx_test/diags_space_charge/beam_stats_sc.txt

Produces warpx_test/results/beam_comparison.png with:
  - Left:  Transverse RMS beam size (x and y) vs. step
  - Right: Normalized transverse emittance vs. step

WarpX BeamRelevant column layout (0-indexed):
  0  step
  1  time (s)
  2  x_mean (m)     3  y_mean (m)    4  z_mean (m)
  5  px_mean        6  py_mean       7  pz_mean        8  gamma_mean
  9  x_rms (m)     10  y_rms (m)   11  z_rms (m)
  12 px_rms        13 py_rms       14 pz_rms          15 gamma_rms
  16 emittance_x   17 emittance_y  18 emittance_z
  19 alpha_x       20 alpha_y      21 beta_x          22 beta_y
  23 charge (C)

Run with:
    conda run -n CBB python warpx_test/plot_results.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")           # non-interactive: works without a display
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ── Column indices ─────────────────────────────────────────────────────────────
COL_STEP   = 0
COL_TIME   = 1
COL_X_RMS  = 9
COL_Y_RMS  = 10
COL_Z_RMS  = 11
COL_EMT_X  = 16
COL_EMT_Y  = 17

# ── Load data ─────────────────────────────────────────────────────────────────
free_path = "warpx_test/diags_bunch/beam_stats.txt"
sc_path   = "warpx_test/diags_space_charge/beam_stats_sc.txt"

for p in [free_path, sc_path]:
    if not os.path.exists(p):
        raise FileNotFoundError(
            f"Missing: {p}\nRun 02_positron_bunch.py and "
            "03_positron_space_charge.py first."
        )

free = np.genfromtxt(free_path, comments="#")
sc   = np.genfromtxt(sc_path,   comments="#")

# Guard against single-row arrays
if free.ndim == 1:
    free = free[np.newaxis, :]
if sc.ndim == 1:
    sc = sc[np.newaxis, :]

steps_free = free[:, COL_STEP]
steps_sc   = sc[:,   COL_STEP]

x_rms_free_mm = free[:, COL_X_RMS] * 1e3
y_rms_free_mm = free[:, COL_Y_RMS] * 1e3
x_rms_sc_mm   = sc[:,   COL_X_RMS] * 1e3
y_rms_sc_mm   = sc[:,   COL_Y_RMS] * 1e3

emt_x_free = free[:, COL_EMT_X]
emt_y_free = free[:, COL_EMT_Y]
emt_x_sc   = sc[:,   COL_EMT_X]
emt_y_sc   = sc[:,   COL_EMT_Y]

print(f"Free-space  — initial sigma_x: {x_rms_free_mm[0]:.3f} mm  "
      f"final: {x_rms_free_mm[-1]:.3f} mm")
print(f"Space-charge — initial sigma_x: {x_rms_sc_mm[0]:.3f} mm  "
      f"final: {x_rms_sc_mm[-1]:.3f} mm")

# ── Plot ──────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(12, 5))
gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.35)

# --- Panel 1: RMS beam size --------------------------------------------------
ax1 = fig.add_subplot(gs[0])

ax1.plot(steps_free, x_rms_free_mm, "b-",  lw=2, label=r"Free space $\sigma_x$")
ax1.plot(steps_free, y_rms_free_mm, "b--", lw=2, label=r"Free space $\sigma_y$", alpha=0.7)
ax1.plot(steps_sc,   x_rms_sc_mm,   "r-",  lw=2, label=r"Space charge $\sigma_x$")
ax1.plot(steps_sc,   y_rms_sc_mm,   "r--", lw=2, label=r"Space charge $\sigma_y$", alpha=0.7)

ax1.set_xlabel("Simulation step")
ax1.set_ylabel(r"RMS beam size [mm]")
ax1.set_title("Transverse RMS Beam Size")
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.3)

# --- Panel 2: Normalized emittance -------------------------------------------
ax2 = fig.add_subplot(gs[1])

ax2.plot(steps_free, emt_x_free * 1e6, "b-",  lw=2, label=r"Free space $\varepsilon_x$")
ax2.plot(steps_free, emt_y_free * 1e6, "b--", lw=2, label=r"Free space $\varepsilon_y$", alpha=0.7)
ax2.plot(steps_sc,   emt_x_sc   * 1e6, "r-",  lw=2, label=r"Space charge $\varepsilon_x$")
ax2.plot(steps_sc,   emt_y_sc   * 1e6, "r--", lw=2, label=r"Space charge $\varepsilon_y$", alpha=0.7)

ax2.set_xlabel("Simulation step")
ax2.set_ylabel(r"Emittance [$\mu$m]")
ax2.set_title("Normalized Transverse Emittance")
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)

fig.suptitle(
    "5 MeV Positron Bunch ($\\gamma\\approx11$): Negligible vs. Strong Space Charge\n"
    r"$\sigma_{x0}=0.5\,$mm, $\sigma_{z0}=1\,$mm, zero initial divergence"
    "\n1 pC (reference) vs 1 nC (space-charge blow-up)",
    fontsize=10,
)

os.makedirs("warpx_test/results", exist_ok=True)
out_path = "warpx_test/results/beam_comparison.png"
fig.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"\nFigure saved → {out_path}")
