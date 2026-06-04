"""
Figures for the finite-cathode space-charge-limited diode (cathode_diode.py).

Reads the openPMD output under warpx_cathode/diags/ and writes four figures to
warpx_cathode/results/:

  1. cathode_2d.png       — 2D maps of charge density, potential, |E|: shows the
                            beam emitted from the finite cathode and the field
                            enhancement at the cathode edges.
  2. child_langmuir.png   — on-axis phi(z), Ez(z) vs. the Child–Langmuir laws.
  3. current_saturation.png — transmitted current vs. time, saturating at J_CL
                            even though we inject 2× J_CL.
  4. rho_z_time.png       — on-axis charge density rho(z, t): build-up of the
                            space-charge cloud filling the gap during turn-on.

Run with:
    conda run -n CBB python warpx_cathode/plot_cathode.py
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from openpmd_viewer import OpenPMDTimeSeries
from scipy.constants import e as q_e, epsilon_0, m_e

# ── Diode parameters (must match cathode_diode.py) ──────────────────────────────
V_anode   = 50.0
gap_d     = 100.0e-6
R_cathode = 6.0e-3
over_inject = 2.0

J_CL = (4.0 / 9.0) * epsilon_0 * np.sqrt(2.0 * q_e / m_e) * V_anode**1.5 / gap_d**2

RESULTS = "warpx_cathode/results"
os.makedirs(RESULTS, exist_ok=True)

ts  = OpenPMDTimeSeries("warpx_cathode/diags/fields")
it  = ts.iterations[-1]            # final (steady-state) snapshot

# Fields are stored with shape (nz, nx); meta.z, meta.x give the axes.
phi, meta = ts.get_field("phi", iteration=it)
ex,  _    = ts.get_field("E", "x", iteration=it)
ez,  _    = ts.get_field("E", "z", iteration=it)
rho, _    = ts.get_field("rho", iteration=it)
z, x = meta.z, meta.x
ix0 = np.argmin(np.abs(x))          # column nearest the axis x = 0
Emag = np.sqrt(ex**2 + ez**2)

extent = [x.min() * 1e3, x.max() * 1e3, z.min() * 1e3, z.max() * 1e3]  # mm

# ════════════════════════════════════════════════════════════════════════════════
# Figure 1 — 2D maps
# ════════════════════════════════════════════════════════════════════════════════
fig, axs = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True)

# |ρ| on a square-root (power-law) scale: electrons make ρ < 0, so plot the
# magnitude.  γ=½ compresses the bright cathode layer to expose the bulk gradient
# without log's amplification of the near-zero noise floor.
absrho = np.abs(rho) * 1e6
rho_norm = PowerNorm(gamma=0.5, vmin=0.0, vmax=absrho.max())

panels = [
    (absrho, "|Charge density|  |ρ|  [µC/m³]  (√)", "viridis", rho_norm),
    (phi,       "Potential  φ  [V]",          "plasma",  None),
    (Emag * 1e-3, "Field magnitude  |E|  [kV/m]", "inferno", None),
]
for ax, (data, title, cmap, norm) in zip(axs, panels):
    im = ax.imshow(data, extent=extent, origin="lower", aspect="auto",
                   cmap=cmap, norm=norm)
    fig.colorbar(im, ax=ax, shrink=0.9)
    # mark the emitting cathode patch (z = 0, |x| < R)
    ax.plot([-R_cathode * 1e3, R_cathode * 1e3], [0, 0], "w-", lw=3,
            solid_capstyle="butt")
    ax.set_title(title)
    ax.set_xlabel("x  [mm]")
axs[0].set_ylabel("z  [mm]   (cathode → anode)")
fig.suptitle("Finite thermionic cathode in WarpX — emission from |x| < "
             f"{R_cathode*1e3:.0f} mm (white bar); note edge field enhancement",
             fontsize=12)
fig.savefig(f"{RESULTS}/cathode_2d.png", dpi=140)
print(f"wrote {RESULTS}/cathode_2d.png")

# ════════════════════════════════════════════════════════════════════════════════
# Figure 2 — on-axis profiles vs. Child–Langmuir theory
# ════════════════════════════════════════════════════════════════════════════════
phi_axis = phi[:, ix0]
ez_axis  = ez[:, ix0]

# Child–Langmuir (1D, planar) reference for the same V and gap:
zt = z
phi_cl = V_anode * (zt / gap_d) ** (4.0 / 3.0)
ez_cl  = -(4.0 / 3.0) * (V_anode / gap_d) * (zt / gap_d) ** (1.0 / 3.0)
# vacuum (no space charge) reference: linear potential, uniform field
phi_vac = V_anode * (zt / gap_d)
ez_vac  = -(V_anode / gap_d) * np.ones_like(zt)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)

a1.plot(zt * 1e3, phi_axis, "o-", color="C0", ms=3, label="WarpX (on axis)")
a1.plot(zt * 1e3, phi_cl, "k--", label=r"Child–Langmuir  $V(z/d)^{4/3}$")
a1.plot(zt * 1e3, phi_vac, ":", color="gray", label="vacuum (no space charge)")
a1.set_xlabel("z  [mm]"); a1.set_ylabel("φ  [V]")
a1.set_title("On-axis potential"); a1.legend()

a2.plot(zt * 1e3, ez_axis * 1e-3, "o-", color="C3", ms=3, label="WarpX (on axis)")
a2.plot(zt * 1e3, ez_cl * 1e-3, "k--",
        label=r"Child–Langmuir  $-\frac{4V}{3d}(z/d)^{1/3}$")
a2.plot(zt * 1e3, ez_vac * 1e-3, ":", color="gray", label="vacuum")
a2.set_xlabel("z  [mm]"); a2.set_ylabel("$E_z$  [kV/m]")
a2.set_title("On-axis longitudinal field"); a2.legend()

fig.suptitle("Space-charge depression of the field at the cathode follows "
             "the Child–Langmuir law", fontsize=12)
fig.savefig(f"{RESULTS}/child_langmuir.png", dpi=140)
print(f"wrote {RESULTS}/child_langmuir.png")

# ════════════════════════════════════════════════════════════════════════════════
# Figure 3 — transmitted current saturates at J_CL despite 2× over-injection
# ════════════════════════════════════════════════════════════════════════════════
# Measure the transmitted current just upstream of the anode.  To be robust to
# the transverse spreading of the finite beam, integrate the total current
# ∫ jz dx across the full domain (charge conservation) and reference it to the
# cathode width 2R — this is the current density the cathode actually delivers.
# jz is negative (electrons moving +z carry negative charge); report its magnitude.
iz_anode = -2                       # row just inside the anode
dx       = x[1] - x[0]

times, J_trans, rho_zt = [], [], []
for i, itr in enumerate(ts.iterations):
    jz, _ = ts.get_field("j", "z", iteration=itr)
    I_line = np.abs(jz[iz_anode, :].sum() * dx)     # total current [A/m depth]
    J_trans.append(I_line / (2.0 * R_cathode))      # referenced to cathode width
    rho_it, _ = ts.get_field("rho", iteration=itr)
    rho_zt.append(rho_it[:, ix0])                   # on-axis charge density column
    times.append(ts.t[i])
times = np.array(times) * 1e9       # ns
J_trans = np.array(J_trans)
rho_zt = np.array(rho_zt).T         # shape (nz, n_times): rho(z, t) on axis

fig, ax = plt.subplots(figsize=(7, 4.8), constrained_layout=True)
ax.plot(times, J_trans, "o-", color="C2", label="WarpX transmitted current")
ax.axhline(J_CL, color="k", ls="--", label=r"Child–Langmuir limit $J_{CL}$")
ax.axhline(over_inject * J_CL, color="r", ls=":",
           label=f"injected current ({over_inject:.0f}× $J_{{CL}}$)")
ax.set_xlabel("time  [ns]")
ax.set_ylabel(r"current density at anode  $|J_z|$  [A/m²]")
ax.set_title("Emission self-limits to the Child–Langmuir value")
# Linear y-axis anchored at the origin so the turn-on ramp and the plateau
# relative to J_CL are both visible to scale.
ax.set_xlim(0, 0.15)
ax.set_ylim(0, J_CL * 1.4)
ax.legend()
fig.savefig(f"{RESULTS}/current_saturation.png", dpi=140)
print(f"wrote {RESULTS}/current_saturation.png")

# ════════════════════════════════════════════════════════════════════════════════
# Figure 4 — on-axis charge density ρ(z, t): build-up of the space-charge cloud
# ════════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
abs_rho_zt = np.abs(rho_zt) * 1e6                   # |ρ| in µC/m³
# Time sampling is non-uniform (dense through the gap-fill transient, sparse in
# steady state), so use pcolormesh with the true time coordinates — imshow would
# force equal-width columns and distort the time axis.
im = ax.pcolormesh(
    times, z * 1e3, abs_rho_zt,
    shading="nearest", cmap="viridis",
    norm=PowerNorm(gamma=0.5, vmin=0.0, vmax=abs_rho_zt.max()),
)
fig.colorbar(im, ax=ax, label="|charge density|  |ρ|  [µC/m³]  (√)")
ax.set_xlim(0, 0.15)            # densely-sampled transient; cloud is steady after
ax.set_xlabel("time  [ns]")
ax.set_ylabel("z  [mm]   (cathode → anode)")
ax.set_title("On-axis charge density vs. time — space-charge cloud fills the gap")
fig.savefig(f"{RESULTS}/rho_z_time.png", dpi=140)
print(f"wrote {RESULTS}/rho_z_time.png")

# ── Quantitative sanity check ───────────────────────────────────────────────────
J_final = J_trans[-1]
print(f"\nJ_CL (theory)          = {J_CL:8.1f} A/m²")
print(f"injected               = {over_inject*J_CL:8.1f} A/m²")
print(f"transmitted (steady)   = {J_final:8.1f} A/m²  "
      f"({100*J_final/J_CL:.0f}% of J_CL)")
