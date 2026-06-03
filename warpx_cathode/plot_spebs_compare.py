"""
Compare the self-consistent WarpX cathode field to the SPEBS analytic vacuum gun.

Overlays the on-axis longitudinal field E_z(z) from two RZ WarpX runs (the SPEBS
~1 A operating point and the Child–Langmuir-limited case) against:
  - the vacuum parallel-plate value  E_z = -V/d,
  - the planar Child–Langmuir law     E_z = -(4V/3d)(z/d)^{1/3},
  - the SPEBS reported peak gun field  -2.48 MV/m (from the "Gun Fieldmap" slide).

Writes warpx_cathode/results/spebs_gun_compare.png.

Run with:
    conda run -n CBB python warpx_cathode/plot_spebs_compare.py
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from openpmd_viewer import OpenPMDTimeSeries

# ── SPEBS gun parameters (match spebs_gun.py) ───────────────────────────────────
V_anode = 50.0e3
gap_d   = 23.1e-3
SPEBS_PEAK = -2.48e6        # MV/m, reported analytic-vacuum gun peak at anode

RESULTS = "warpx_cathode/results"
os.makedirs(RESULTS, exist_ok=True)


def on_axis_Ez(tag):
    ts = OpenPMDTimeSeries(f"warpx_cathode/diags_{tag}/fields")
    ez, meta = ts.get_field("E", "z", iteration=ts.iterations[-1])
    ir0 = np.argmin(np.abs(meta.r))         # on-axis column
    return meta.z, ez[:, ir0]


z, ez_spebs = on_axis_Ez("spebs")
_, ez_cl    = on_axis_Ez("cl")

# Analytic references on the same z grid
ez_vac = -(V_anode / gap_d) * np.ones_like(z)
ez_clt = -(4.0 / 3.0) * (V_anode / gap_d) * (z / gap_d) ** (1.0 / 3.0)

# ── Figure ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8.5, 5.5), constrained_layout=True)

ax.plot(z * 1e3, ez_spebs * 1e-6, "-", color="C0", lw=2.2,
        label="WarpX self-consistent, SPEBS 1 A (0.41× $J_{CL}$)")
ax.plot(z * 1e3, ez_cl * 1e-6, "-", color="C3", lw=2.2,
        label="WarpX self-consistent, Child–Langmuir limit")
ax.plot(z * 1e3, ez_vac * 1e-6, "--", color="gray", lw=1.8,
        label=r"vacuum parallel plate  $-V/d$")
ax.plot(z * 1e3, ez_clt * 1e-6, ":", color="k", lw=1.8,
        label=r"planar Child–Langmuir  $-\frac{4V}{3d}(z/d)^{1/3}$")
ax.scatter([gap_d * 1e3], [SPEBS_PEAK * 1e-6], color="purple", zorder=5, s=60,
           label="SPEBS analytic-vacuum gun peak (−2.48 MV/m)")

ax.axvline(0, color="0.7", lw=0.8)
ax.text(0.3, ax.get_ylim()[0] * 0.5, "cathode", rotation=90,
        va="center", color="0.4", fontsize=9)
ax.set_xlabel("z  [mm]   (cathode at 0 → anode at 23.1 mm)")
ax.set_ylabel("on-axis $E_z$  [MV/m]")
ax.set_title("SPEBS gun: self-consistent space-charge field vs. vacuum model\n"
             "(50 kV, 23.1 mm gap, 4 mm cathode, RZ)")
ax.legend(loc="lower center", fontsize=8.5)
fig.savefig(f"{RESULTS}/spebs_gun_compare.png", dpi=140)
print(f"wrote {RESULTS}/spebs_gun_compare.png")

# ── Quantitative cathode-field depression ───────────────────────────────────────
Ez_vac_cathode = V_anode / gap_d
print(f"\nVacuum parallel-plate field            : {Ez_vac_cathode/1e6:6.3f} MV/m")
print(f"SPEBS reported vacuum gun peak (anode)  : {abs(SPEBS_PEAK)/1e6:6.3f} MV/m")
for tag, ez in (("SPEBS 1 A (0.41x CL)", ez_spebs), ("CL limit", ez_cl)):
    Ec = abs(ez[0])         # field at the cathode surface (first cell)
    print(f"WarpX cathode-surface |Ez|, {tag:20s}: {Ec/1e6:6.3f} MV/m "
          f"({100*Ec/Ez_vac_cathode:4.0f}% of vacuum)")
