"""
Plot the single-positron trajectory from Script 01.

Reads the per-step AMReX plotfiles in diags_single/ with yt, extracts the
positron's z position and z momentum at each step, and produces
results/single_positron.png:

  (left)  z vs time   -- should be a straight line at v = beta*c
  (right) gamma vs time -- should stay flat at ~196.8 (free space, no fields)

Run with:
    conda run -n CBB python warpx_test/plot_single.py
"""

import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import yt

yt.set_log_level("error")  # quiet the per-file INFO spam

HERE = os.path.dirname(os.path.abspath(__file__))
DIAG_DIR = os.path.join(HERE, "diags_single")
OUT_DIR = os.path.join(HERE, "results")

# Physical constants (SI)
c = 2.99792458e8        # m/s
m_e = 9.1093837015e-31  # kg

# Match only the canonical plotfiles part000000 .. partNNNNNN
# (ignores stray copies like "part000010 3").
plotfiles = sorted(
    p for p in glob.glob(os.path.join(DIAG_DIR, "part*"))
    if re.fullmatch(r"part\d+", os.path.basename(p))
)
if not plotfiles:
    raise SystemExit(f"No plotfiles found in {DIAG_DIR}. Run 01_single_positron.py first.")

t, z, gamma = [], [], []
for pf in plotfiles:
    ds = yt.load(pf)
    ad = ds.all_data()
    # Single particle -> take element [0]
    pz = float(ad[("positrons", "particle_momentum_z")][0])  # kg*m/s
    zz = float(ad[("positrons", "particle_position_z")][0])  # m
    # gamma from gamma*beta*c = p/m  ->  gamma = sqrt(1 + (p/(m c))^2)
    g = np.sqrt(1.0 + (pz / (m_e * c)) ** 2)
    t.append(float(ds.current_time))
    z.append(zz)
    gamma.append(g)

t = np.array(t) * 1e9   # ns
z = np.array(z)
gamma = np.array(gamma)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

ax1.plot(t, z * 100, "o-", ms=3, color="C0")
ax1.set_xlabel("time [ns]")
ax1.set_ylabel("z position [cm]")
ax1.set_title("Positron trajectory (free space)")
ax1.grid(alpha=0.3)

ax2.plot(t, gamma, "o-", ms=3, color="C1")
ax2.axhline(196.8, ls="--", color="gray", lw=1, label="expected γ ≈ 196.8")
ax2.set_xlabel("time [ns]")
ax2.set_ylabel("Lorentz factor γ")
ax2.set_title("Energy (should be constant)")
ax2.legend()
ax2.grid(alpha=0.3)

fig.tight_layout()
os.makedirs(OUT_DIR, exist_ok=True)
out = os.path.join(OUT_DIR, "single_positron.png")
fig.savefig(out, dpi=150)
print(f"Saved {out}")

# Quick numeric sanity check
v_fit = np.polyfit(t * 1e-9, z, 1)[0]  # m/s
print(f"Fitted velocity: {v_fit:.4e} m/s  ({v_fit / c:.5f} c)")
print(f"Mean gamma: {gamma.mean():.2f}  (spread {gamma.max() - gamma.min():.2e})")
