"""
Figures for the Cornell Linac sections 4-8 stage (Impact-T, sim/linac4-8.py).

Reads logs/diags/linac4-8/main/{particles, injection_summary.json} and writes five PNGs to
logs/plots/linac4-8/: energy_gain, energy_spread, emittance, section_gains (per-section achieved
vs frozen-target ΔE bars), fodo_optics (quads-OFF sigma_x/y vs z). The vs-z curves come from the
summary's stat_vs_z table (Impact-T I.stat); a sparse particle-slice fallback covers legacy dumps.

main() runs ONLY plotting (sim/linac4-8.py must have been run first). Run as
`python sim/plot/linac4-8.py` (hyphenated name is not importable).
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

MC2 = MC2_EV / 1e6                  # electron rest energy [MeV]
DIAG_DIR = "logs/diags/linac4-8/main"
RESULTS = "logs/plots/linac4-8"


def _wstat(a, w):
    """Weighted mean and standard deviation."""
    m = np.average(a, weights=w)
    return m, np.sqrt(np.average((a - m) ** 2, weights=w))


def _norm_emit(x, ux, w):
    """Normalized RMS emittance eps_n = sqrt(<x^2><ux^2> - <x*ux>^2) (ux = gamma*beta_x)."""
    xm = np.average(x, weights=w)
    um = np.average(ux, weights=w)
    xx = np.average((x - xm) ** 2, weights=w)
    uu = np.average((ux - um) ** 2, weights=w)
    xu = np.average((x - xm) * (ux - um), weights=w)
    return float(np.sqrt(max(0.0, xx * uu - xu * xu)))


def _read_slices(diag):
    """Per-dump beam moments sorted by <z>: list of dicts (z, ke, dke, enx, eny, sx, sy, q).

    Fallback for summaries without a stat_vs_z table (reads the openPMD particle slices).
    """
    from openpmd_viewer import OpenPMDTimeSeries
    ts = OpenPMDTimeSeries(os.path.join(diag, "particles"))
    rows = []
    for it in ts.iterations:
        x, y, z, ux, uy, uz, w = ts.get_particle(
            ["x", "y", "z", "ux", "uy", "uz", "w"], species="electrons", iteration=it)
        if len(z) < 50:
            continue
        g = np.sqrt(1.0 + ux ** 2 + uy ** 2 + uz ** 2)
        ke = (g - 1.0) * MC2
        km, dk = _wstat(ke, w)
        rows.append(dict(
            z=np.average(z, weights=w), ke=km, dke=dk,
            enx=_norm_emit(x, ux, w), eny=_norm_emit(y, uy, w),
            sx=_wstat(x, w)[1], sy=_wstat(y, w)[1], q=float(w.sum())))
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


def _vs_z(diag, summ):
    """(z, ke, dke, enx, eny, sx, sy) [m / MeV / mm-as-m units kept SI] from the summary's
    stat_vs_z (preferred) or the particle slices. Returns None if neither is usable.
    """
    svz = summ.get("stat_vs_z", {})
    if svz.get("z_m"):
        z = np.array(svz["z_m"])
        return (z, np.array(svz["ke_mev"]), np.array(svz["sigma_ke_mev"]),
                np.array(svz["norm_emit_x"]), np.array(svz["norm_emit_y"]),
                np.array(svz["sigma_x_m"]),
                np.array(svz.get("sigma_y_m", svz["sigma_x_m"])))
    rows = _read_slices(diag)
    if not rows:
        return None
    return (_arr(rows, "z"), _arr(rows, "ke"), _arr(rows, "dke"),
            _arr(rows, "enx"), _arr(rows, "eny"), _arr(rows, "sx"), _arr(rows, "sy"))


def _achieved_de(calib, z, ke):
    """Per-section achieved ΔE [MeV] from the vs-z KE curve: KE at each section's exit z minus
    the entry KE. Section exit z's are reconstructed by cumulating the calibration order onto the
    z-grid extent (the frozen calib carries no z, so split the line evenly across the sections --
    a coarse but monotone read of the per-section gain for the target-vs-achieved bar).
    """
    n = len(calib)
    if n == 0 or len(z) < 2:
        return []
    z0, z1 = float(z[0]), float(z[-1])
    edges = np.linspace(z0, z1, n + 1)
    ke_at = lambda zz: float(np.interp(zz, z, ke))
    return [ke_at(edges[i + 1]) - ke_at(edges[i]) for i in range(n)]


def main():
    diag = DIAG_DIR
    os.makedirs(RESULTS, exist_ok=True)
    summ = _load_summary(diag)
    calib = summ.get("calibration", [])

    vs = _vs_z(diag, summ)
    if vs is None:
        print(f"plot linac4-8: no stat_vs_z and no usable dumps in {diag} -- skipping.",
              flush=True)
        return
    z, ke, dke, enx, eny, sx, sy = vs
    power_mw = summ.get("power_mw", 11.0)

    # 1) energy gain vs z
    fig, ax = plt.subplots(figsize=(9.2, 4.8), constrained_layout=True)
    ax.plot(z, ke, "-o", ms=3, color="C5", label="<KE>")
    ax.fill_between(z, ke - dke, ke + dke, color="C5", alpha=0.18, label="+/- sigma_KE")
    ax.set_xlabel("<z> (local Impact-T frame) [m]")
    ax.set_ylabel("kinetic energy [MeV]")
    exp = summ.get("ke_out_mev")
    if exp:
        ax.axhline(exp, color="0.5", ls=":", lw=1, label=f"exit {exp:.0f} MeV")
    ax.set_title(f"linac4-8: cumulative energy gain (sections 4-8, on-crest, {power_mw:g} MW)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.savefig(os.path.join(RESULTS, "energy_gain.png"), dpi=130)
    plt.close(fig)

    # 2) energy spread vs z (absolute grows on the cosine crest curvature; relative shrinks)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.6, 6.2), constrained_layout=True, sharex=True)
    a1.plot(z, dke, "-o", ms=3, color="C5")
    a1.set_ylabel("sigma_KE [MeV]")
    a1.set_title("linac4-8: energy spread (absolute grows, relative shrinks)")
    a1.grid(alpha=0.3)
    rel = np.where(ke > 0, dke / ke * 100.0, np.nan)
    a2.plot(z, rel, "-o", ms=3, color="C6")
    a2.set_ylabel("sigma_KE/<KE> [%]")
    a2.set_xlabel("<z> [m]")
    a2.grid(alpha=0.3)
    fig.savefig(os.path.join(RESULTS, "energy_spread.png"), dpi=130)
    plt.close(fig)

    # 3) normalized emittance vs z (quads OFF: the ~2.4x rise is a fort.10N diagnostic artifact)
    fig, ax = plt.subplots(figsize=(9.2, 4.8), constrained_layout=True)
    ax.plot(z, enx * 1e6, "-", color="C5", label="eps_n,x")
    ax.plot(z, eny * 1e6, "-", color="C6", label="eps_n,y")
    ax.set_xlabel("<z> [m]")
    ax.set_ylabel("normalized emittance [mm.mrad]")
    ax.set_title("linac4-8: normalized emittance\n"
                 "quads OFF -- eps_n rises ~2.4x: a fort.10N diagnostic artifact, not physical",
                 fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.savefig(os.path.join(RESULTS, "emittance.png"), dpi=130)
    plt.close(fig)

    # 4) per-section achieved vs frozen-target ΔE
    fig, ax = plt.subplots(figsize=(9.2, 4.8), constrained_layout=True)
    if calib:
        names = [c.get("name", f"sec{c['index'] + 4}") for c in calib]
        tgt = [c.get("target_de_mev", np.nan) for c in calib]
        ach = _achieved_de(calib, z, ke)
        xi = np.arange(len(calib))
        ax.bar(xi - 0.2, tgt, 0.4, label="frozen target dE", color="0.6")
        if ach:
            ax.bar(xi + 0.2, ach, 0.4, label="achieved dE (from vs-z KE)", color="C5")
        ax.set_xticks(xi)
        ax.set_xticklabels(names, rotation=30, ha="right", fontsize=8)
        ax.set_ylabel("dE per section [MeV]")
        ax.set_title("linac4-8: per-section gain -- frozen target vs achieved")
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "no calibration table in injection_summary.json",
                ha="center", va="center", transform=ax.transAxes)
    ax.grid(alpha=0.3, axis="y")
    fig.savefig(os.path.join(RESULTS, "section_gains.png"), dpi=130)
    plt.close(fig)

    # 5) FODO optics (quads OFF): sigma_x / sigma_y vs z -- placeholder optics, NOT predictive
    fig, ax = plt.subplots(figsize=(9.2, 4.8), constrained_layout=True)
    ax.plot(z, sx * 1e3, "-", color="C7", label="sigma_x")
    ax.plot(z, sy * 1e3, "-", color="C8", label="sigma_y")
    ax.set_xlabel("<z> [m]")
    ax.set_ylabel("transverse RMS size [mm]")
    ax.set_title("linac4-8: transverse envelope sigma_x / sigma_y\n"
                 "quads OFF -- no focusing, placeholder optics, NOT predictive", fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.savefig(os.path.join(RESULTS, "fodo_optics.png"), dpi=130)
    plt.close(fig)

    print(f"plot linac4-8: wrote 5 figures to {RESULTS}/ "
          f"({len(z)} vs-z points, exit <KE> {ke[-1]:.1f} MeV).", flush=True)


if __name__ == "__main__":
    main()
