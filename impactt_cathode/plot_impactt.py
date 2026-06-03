"""
Render the SPEBS-gun IMPACT-T figures from the stats_*.json written by
spebs_gun_impactt.py. Run the driver for both settings first:

    python impactt_cathode/spebs_gun_impactt.py img
    python impactt_cathode/spebs_gun_impactt.py noimg
    python impactt_cathode/plot_impactt.py

Produces three PNGs in impactt_cathode/results/:
  gun_dynamics.png    energy gain, beam size, and emittance vs z (image-charge run)
  final_phasespace.png  x-px and z-pz of the final bunch
  image_charge.png    image-charge ON vs OFF: emittance and exit energy
"""

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")


def load(tag):
    path = os.path.join(RESULTS, f"stats_{tag}.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


img = load("img")
noimg = load("noimg")
if img is None:
    raise SystemExit("Missing results/stats_img.json — run spebs_gun_impactt.py img first.")

z = np.array(img["z"]) * 1e3          # mm
E0 = img["E0"]
gap = img["gap"]

# ── Figure 1: gun dynamics (energy, size, emittance) ─────────────────────────────
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))

ax[0].plot(z, np.array(img["KE"]) / 1e3, lw=2, label="IMPACT-T")
zc = np.linspace(0, gap * 1e3, 50)
ax[0].plot(zc, E0 * (zc / 1e3) / 1e3, "k--", label=r"uniform field $eE_0 z$")
ax[0].axhline(img["V"] / 1e3, color="gray", ls=":", label=f"anode {img['V']/1e3:.0f} kV")
ax[0].set_xlabel("z [mm]"); ax[0].set_ylabel("mean kinetic energy [keV]")
ax[0].set_title("Energy gain across the gap"); ax[0].legend(); ax[0].grid(alpha=0.3)

ax[1].plot(z, np.array(img["sigma_x"]) * 1e3, lw=2, label=r"$\sigma_x$")
ax[1].plot(z, np.array(img["sigma_y"]) * 1e3, lw=2, ls="--", label=r"$\sigma_y$")
ax[1].set_xlabel("z [mm]"); ax[1].set_ylabel("RMS beam size [mm]")
ax[1].set_title("Transverse beam size"); ax[1].legend(); ax[1].grid(alpha=0.3)

ax[2].plot(z, np.array(img["norm_emit_x"]) * 1e6, lw=2, label=r"$\varepsilon_{n,x}$")
ax[2].plot(z, np.array(img["norm_emit_y"]) * 1e6, lw=2, ls="--", label=r"$\varepsilon_{n,y}$")
ax[2].set_xlabel("z [mm]"); ax[2].set_ylabel(r"norm. emittance [$\mu$m]")
ax[2].set_title("Emittance growth (1500 K cathode)"); ax[2].legend(); ax[2].grid(alpha=0.3)

fig.suptitle(f"SPEBS gun in IMPACT-T — {img['V']/1e3:.0f} kV, {gap*1e3:.1f} mm gap, "
             f"image charge ON, transmission {img['transmission']*100:.0f}%", y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(RESULTS, "gun_dynamics.png"), dpi=130, bbox_inches="tight")
print("wrote results/gun_dynamics.png")

# ── Figure 2: final phase space ──────────────────────────────────────────────────
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].scatter(np.array(img["fp_x"]) * 1e3, np.array(img["fp_px"]) / 1e3, s=3, alpha=0.4)
ax[0].set_xlabel("x [mm]"); ax[0].set_ylabel(r"$p_x$ [keV/c]")
ax[0].set_title("Transverse phase space (final)"); ax[0].grid(alpha=0.3)

ax[1].scatter(np.array(img["fp_z"]) * 1e3, np.array(img["fp_pz"]) / 1e3, s=3, alpha=0.4)
ax[1].set_xlabel("z [mm]"); ax[1].set_ylabel(r"$p_z$ [keV/c]")
ax[1].set_title("Longitudinal phase space (final)"); ax[1].grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(RESULTS, "final_phasespace.png"), dpi=130, bbox_inches="tight")
print("wrote results/final_phasespace.png")

# ── Figure 3: image-charge comparison ────────────────────────────────────────────
if noimg is not None:
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for d, c, lab in [(img, "C0", "image charge ON"), (noimg, "C1", "image charge OFF")]:
        zz = np.array(d["z"]) * 1e3
        ax[0].plot(zz, np.array(d["norm_emit_x"]) * 1e6, color=c, lw=2, label=lab)
        ax[1].plot(zz, np.array(d["KE"]) / 1e3, color=c, lw=2, label=lab)
    ax[0].set_xlabel("z [mm]"); ax[0].set_ylabel(r"norm. emittance $\varepsilon_{n,x}$ [$\mu$m]")
    ax[0].set_title("Cathode image-charge effect on emittance")
    ax[0].legend(); ax[0].grid(alpha=0.3)
    ax[0].set_xlim(0, 3)  # zoom near the cathode where image charge matters
    ax[1].set_xlabel("z [mm]"); ax[1].set_ylabel("mean kinetic energy [keV]")
    ax[1].set_title("Energy gain: image charge ON vs OFF")
    ax[1].legend(); ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "image_charge.png"), dpi=130, bbox_inches="tight")
    print("wrote results/image_charge.png")

    print(f"\nimage ON : exit KE={img['final_KE_eV']/1e3:.2f} keV, "
          f"eps_n,x={img['final_norm_emit_x']*1e6:.3f} um, "
          f"transmission={img['transmission']*100:.0f}%")
    print(f"image OFF: exit KE={noimg['final_KE_eV']/1e3:.2f} keV, "
          f"eps_n,x={noimg['final_norm_emit_x']*1e6:.3f} um, "
          f"transmission={noimg['transmission']*100:.0f}%")
else:
    print("(run 'spebs_gun_impactt.py noimg' to also produce image_charge.png)")
