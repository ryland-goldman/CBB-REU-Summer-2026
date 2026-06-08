"""
Figures for the Cornell Linac Sections 2–8 stage (linac_rest_sim.py, Impact-T).

Reads the stage's openPMD beam diagnostics (linac_rest/diags/main/particles — the
along-z slices the sim writes via pipeline.impact_io) and the per-section calibration
table recorded in linac_rest/diags/main/injection_summary.json. Writes to
linac_rest/results/:

  1. energy_gain.png      — ⟨KE⟩ (± σ_KE band) vs ⟨z⟩ across the 7 sections, with the
                            per-section ΔE_target annotated; the on-crest cumulative rise
                            to ≈308 MeV at 11 MW (307.97 survivors).
  2. energy_spread.png    — σ_KE and the relative spread σ_KE/⟨KE⟩ vs ⟨z⟩ (the relative
                            spread adiabatically shrinks as ⟨KE⟩ grows).
  3. emittance.png        — normalized emittance εn,x / εn,y vs ⟨z⟩ (quads OFF ⇒ near-
                            conserved; flags numerical growth). The headline beam.
  4. section_gains.png     — per-section achieved vs target ΔE bar chart (from the
                            calibration table) — the §5 gate-1 visual.
  5. fodo_optics.png       — EXPLORATORY (labeled "placeholder optics — not predictive"):
                            σ_x vs ⟨z⟩; only meaningful with QUADS_ON. Always written so
                            the figure set is stable; titled to flag the quad state.

Run with:
    conda run -n CBB python -c "import linac_rest; linac_rest.plot()"
"""

import os
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from openpmd_viewer import OpenPMDTimeSeries

from . import DEFAULT_OUTDIR
from . import build_linac_rest_lattice as L

MC2 = 0.51099895069              # electron rest energy [MeV]
OUTDIR = None                    # config(OUTDIR=...) sets this; None → DEFAULT_OUTDIR
RESULTS = "linac_rest/results"


def _wstat(a, w):
    """Weighted mean and standard deviation."""
    m = np.average(a, weights=w)
    return m, np.sqrt(np.average((a - m) ** 2, weights=w))


def _norm_emit(x, ux, w):
    """Normalized RMS emittance ε_n = √(⟨x²⟩⟨ux²⟩ − ⟨x·ux⟩²) (ux = γβ_x)."""
    xm = np.average(x, weights=w)
    um = np.average(ux, weights=w)
    xx = np.average((x - xm) ** 2, weights=w)
    uu = np.average((ux - um) ** 2, weights=w)
    xu = np.average((x - xm) * (ux - um), weights=w)
    return float(np.sqrt(max(0.0, xx * uu - xu * xu)))


def _read_slices(diag):
    """Return per-dump beam moments sorted by ⟨z⟩: list of dicts (z, ke, dke, enx, eny, q)."""
    ts = OpenPMDTimeSeries(os.path.join(diag, "particles"))
    rows = []
    for it in ts.iterations:
        x, y, z, ux, uy, uz, w = ts.get_particle(
            ["x", "y", "z", "ux", "uy", "uz", "w"], species="electrons", iteration=it)
        if len(z) < 50:
            continue
        g = np.sqrt(1.0 + ux**2 + uy**2 + uz**2)
        ke = (g - 1.0) * MC2
        km, dk = _wstat(ke, w)
        zm = np.average(z, weights=w)
        rows.append(dict(
            z=zm, ke=km, dke=dk,
            enx=_norm_emit(x, ux, w), eny=_norm_emit(y, uy, w),
            sx=_wstat(x, w)[1], q=float(w.sum()),
        ))
    rows.sort(key=lambda r: r["z"])
    return rows


def _load_summary(diag):
    path = os.path.join(diag, "injection_summary.json")
    if os.path.isfile(path):
        with open(path) as fh:
            return json.load(fh)
    return {}


def _arr(rows, key):
    return np.array([r[key] for r in rows])


def main():
    diag = OUTDIR or DEFAULT_OUTDIR
    os.makedirs(RESULTS, exist_ok=True)
    summ = _load_summary(diag)
    calib = summ.get("calibration", [])
    # Prefer Impact-T's continuous I.stat() vs-z table (written to the summary) for the
    # evolution panels; fall back to the sparse openPMD particle slices if it's absent.
    svz = summ.get("stat_vs_z", {})
    if svz.get("z_m"):
        z = np.array(svz["z_m"])
        ke = np.array(svz["ke_mev"])
        dke = np.array(svz["sigma_ke_mev"])
        enx = np.array(svz["norm_emit_x"])
        eny = np.array(svz["norm_emit_y"])
        sx = np.array(svz["sigma_x_m"])
    else:
        rows = _read_slices(diag)
        if not rows:
            print(f"plot_linac_rest: no stat_vs_z and no usable dumps in {diag} — skipping.",
                  flush=True)
            return
        z, ke, dke = _arr(rows, "z"), _arr(rows, "ke"), _arr(rows, "dke")
        enx, eny, sx = _arr(rows, "enx"), _arr(rows, "eny"), _arr(rows, "sx")

    # 1) energy gain vs z
    fig, ax = plt.subplots(figsize=(9.2, 4.8), constrained_layout=True)
    ax.plot(z, ke, "-o", ms=3, color="C5", label="⟨KE⟩")
    ax.fill_between(z, ke - dke, ke + dke, color="C5", alpha=0.18, label="± σ_KE")
    ax.set_xlabel("⟨z⟩ (local Impact-T frame) [m]")
    ax.set_ylabel("kinetic energy [MeV]")
    exp = summ.get("expected_ke_out_mev")
    ttl = "linac_rest: cumulative energy gain (sections 2–8, on-crest"
    ttl += f", {summ.get('power_mw', L.POWER_MW):g} MW)"
    if exp:
        ax.axhline(exp, color="0.5", ls=":", lw=1,
                   label=f"expected exit {exp:.0f} MeV")
    ax.set_title(ttl)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.savefig(os.path.join(RESULTS, "energy_gain.png"), dpi=130)
    plt.close(fig)

    # 2) energy spread vs z
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.6, 6.2), constrained_layout=True,
                                 sharex=True)
    a1.plot(z, dke, "-o", ms=3, color="C5")
    a1.set_ylabel("σ_KE [MeV]")
    a1.set_title("linac_rest: energy spread (absolute grows ~3.9×, relative shrinks)")
    a1.grid(alpha=0.3)
    rel = np.where(ke > 0, dke / ke * 100.0, np.nan)
    a2.plot(z, rel, "-o", ms=3, color="C6")
    a2.set_ylabel("σ_KE/⟨KE⟩ [%]")
    a2.set_xlabel("⟨z⟩ [m]")
    a2.grid(alpha=0.3)
    fig.savefig(os.path.join(RESULTS, "energy_spread.png"), dpi=130)
    plt.close(fig)

    # 3) normalized emittance vs z (headline, quads OFF ⇒ ~conserved)
    fig, ax = plt.subplots(figsize=(9.2, 4.8), constrained_layout=True)
    ax.plot(z, enx * 1e6, "-", color="C5", label="ε_n,x")
    ax.plot(z, eny * 1e6, "-", color="C6", label="ε_n,y")
    ax.set_xlabel("⟨z⟩ [m]")
    ax.set_ylabel("normalized emittance [mm·mrad]")
    ax.set_title("linac_rest: normalized emittance (quads OFF ⇒ RF + drift only)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.savefig(os.path.join(RESULTS, "emittance.png"), dpi=130)
    plt.close(fig)

    # 4) per-section achieved vs target ΔE (calibration table)
    fig, ax = plt.subplots(figsize=(9.2, 4.8), constrained_layout=True)
    if calib:
        names = [c.get("name", f"sec{c['index']+2}") for c in calib]
        tgt = [c["target_de_mev"] for c in calib]
        ach = [c["achieved_de_mev"] for c in calib]
        xi = np.arange(len(calib))
        ax.bar(xi - 0.2, tgt, 0.4, label="target ΔE", color="0.6")
        ax.bar(xi + 0.2, ach, 0.4, label="achieved ΔE", color="C5")
        ax.set_xticks(xi)
        ax.set_xticklabels(names, rotation=30, ha="right", fontsize=8)
        ax.set_ylabel("ΔE per section [MeV]")
        ax.set_title("linac_rest: per-section gain — target vs achieved (±3% gate)")
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "no calibration table in injection_summary.json",
                ha="center", va="center", transform=ax.transAxes)
    ax.grid(alpha=0.3, axis="y")
    fig.savefig(os.path.join(RESULTS, "section_gains.png"), dpi=130)
    plt.close(fig)

    # 5) exploratory FODO optics (σ_x vs z) — labeled placeholder
    fig, ax = plt.subplots(figsize=(9.2, 4.8), constrained_layout=True)
    ax.plot(z, sx * 1e3, "-", color="C7")
    ax.set_xlabel("⟨z⟩ [m]")
    ax.set_ylabel("σ_x [mm]")
    quads_on = summ.get("quads_on", False)
    ax.set_title("linac_rest: transverse envelope σ_x  "
                 + ("(QUADS ON — exploratory)" if quads_on
                    else "(quads OFF — placeholder optics, NOT predictive)"))
    ax.grid(alpha=0.3)
    fig.savefig(os.path.join(RESULTS, "fodo_optics.png"), dpi=130)
    plt.close(fig)

    print(f"plot_linac_rest: wrote 5 figures to {RESULTS}/ "
          f"({len(z)} vs-z points, exit ⟨KE⟩ {ke[-1]:.1f} MeV).", flush=True)


if __name__ == "__main__":
    main()
