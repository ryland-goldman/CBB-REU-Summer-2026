"""
Figures for the e+/e- converter stage (sim/converter.py). Reads
logs/diags/converter/main/{particles, injection_summary.json}, writes PNGs to logs/plots/converter/.

main() runs ONLY plotting (sim/converter.py must have been run first). Run as
`python sim/plot/converter.py`. See docs/converter.md.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sim.helpers.tools import MC2_EV
from sim.helpers.loadparticles import read_warpx_dump
from sim.plot import common as px

DIAG_DIR = "logs/diags/converter/main"
RESULTS = "logs/plots/converter"


def _load_summary(diag):
    path = os.path.join(diag, "injection_summary.json")
    if os.path.isfile(path):
        with open(path) as fh:
            return json.load(fh)
    return {}


def main():
    diag = DIAG_DIR
    os.makedirs(RESULTS, exist_ok=True)
    summ = _load_summary(diag)

    parts = os.path.join(diag, "particles")
    if not os.path.isdir(parts):
        print(f"plot converter: no particle dump in {parts} -- skipping.", flush=True)
        return

    pg = read_warpx_dump(parts, species="positrons")
    if pg.n_particle == 0:
        print(f"plot converter: empty positron dump in {parts} -- skipping.", flush=True)
        return

    fig = px.energy_spectrum(pg, use_ke=True, e_unit="MeV")
    fig.savefig(os.path.join(RESULTS, "positron_spectrum.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.6), constrained_layout=True)
    labels = ["positron", "electron", "gamma"]
    vals = [summ.get("yield_positron", np.nan), summ.get("yield_electron", np.nan),
            summ.get("yield_gamma", np.nan)]
    bars = ax.bar(labels, vals, color=["C3", "C0", "C2"], alpha=0.85)
    for b, v in zip(bars, vals):
        if np.isfinite(v):
            ax.annotate(f"{v:.3g}", (b.get_x() + b.get_width() / 2, v),
                        ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("per incident e$^-$")
    tlen = summ.get("target_length_mm", "?")
    tmat = summ.get("target_material", "?")
    ax.set_title(f"Converter secondary yield  ({tlen} mm {tmat} target)")
    ax.grid(alpha=0.3, axis="y")
    fig.savefig(os.path.join(RESULTS, "yield_bars.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)

    # pg.px is momentum in eV/c (= ux*MC2_EV); ux = px/MC2_EV.
    x = np.asarray(pg.x, float)
    y = np.asarray(pg.y, float)
    ux = np.asarray(pg.px, float) / MC2_EV
    uy = np.asarray(pg.py, float) / MC2_EV
    w = np.asarray(pg.weight, float)
    div_rms = summ.get("div_pos_rms_mrad")
    title = "Converter positron transverse phase space  (r, p_r)"
    if div_rms is not None:
        title += f"\nRMS divergence {div_rms:.0f} mrad (~{div_rms / 1e3:.2g} rad) -- capture challenge"
    fig = px.transverse_rpr(x, y, ux, uy, w, title=title, p_unit="MeV")
    fig.savefig(os.path.join(RESULTS, "positron_divergence.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)

    # x-axis clipped at the 99th percentile -- the wide-angle tail otherwise dominates the plot range.
    uz = np.asarray(pg.pz, float) / MC2_EV
    rdiv_mrad = 1e3 * np.hypot(ux, uy) / uz
    hi = float(np.percentile(rdiv_mrad, 99))
    med = float(np.median(rdiv_mrad))
    fig, ax = plt.subplots(figsize=(7.0, 4.8), constrained_layout=True)
    ax.hist(np.clip(rdiv_mrad, 0.0, hi), bins=80, weights=w * 1e15, color="C3", alpha=0.85)
    ax.axvline(med, ls="--", color="k", lw=1.2,
               label=f"median {med:.0f} mrad\nRMS {div_rms:.0f} mrad" if div_rms is not None
               else f"median {med:.0f} mrad")
    ax.set_xlabel(r"radial divergence  $p_r/p_z$  [mrad]")
    ax.set_ylabel("charge / bin  [fC]")
    ax.set_title("Converter positron radial divergence  (99th-pct x-clip)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, axis="y")
    fig.savefig(os.path.join(RESULTS, "positron_radial_divergence.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)

    ke_mev = pg["kinetic_energy"] / 1e6
    z_mm = (np.asarray(pg.z, float) - np.average(pg.z, weights=w)) * 1e3
    fig, ax = plt.subplots(figsize=(7.0, 4.8), constrained_layout=True)
    hb = ax.hexbin(z_mm, ke_mev, gridsize=70, cmap="viridis", mincnt=1)
    fig.colorbar(hb, ax=ax, label="macroparticles / bin")
    ax.set_xlabel("z - <z>  [mm]")
    ax.set_ylabel("kinetic energy  [MeV]")
    ax.set_title("Converter positron longitudinal phase space  (z, KE)")
    fig.savefig(os.path.join(RESULTS, "z_ke.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)

    print(f"plot converter: wrote 5 figures to {RESULTS}/ "
          f"({pg.n_particle} positrons, <KE> {np.average(ke_mev, weights=w):.1f} MeV).",
          flush=True)


if __name__ == "__main__":
    main()
