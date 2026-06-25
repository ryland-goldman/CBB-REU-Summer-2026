"""
Figures for the Cornell Linac sections 5-8 stage (Impact-T, sim/linac5-8.py).

Reads logs/diags/linac5-8/main/{particles, injection_summary.json} and writes PNGs to
logs/plots/linac5-8/: evolution_vs_z (mean KE / eps_n,x / sigma_x / surviving charge), energy_spread,
section_gains (per-section achieved vs frozen-target ΔE bars), and from the exit particle slice
energy_spectrum, phase_space_z_KE (longitudinal) and transverse_r_pr. The vs-z curves come from the
summary's stat_vs_z table (Impact-T I.stat); a sparse particle-slice fallback covers legacy dumps.

main() runs ONLY plotting (sim/linac5-8.py must have been run first). Run as
`python sim/plot/linac5-8.py` (hyphenated name is not importable).
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
from sim.helpers.loadparticles import make_particle_group
from sim.plot import common as px

MC2 = MC2_EV / 1e6                  # electron rest energy [MeV]
DIAG_DIR = "logs/diags/linac5-8/main"
RESULTS = "logs/plots/linac5-8"


def _species(ts):
    """The openPMD particle-group key to read (positron handoff => "positrons")."""
    return ts.avail_species[0] if getattr(ts, "avail_species", None) else "electrons"


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
    try:
        ts = OpenPMDTimeSeries(os.path.join(diag, "particles"))
    except Exception:
        return []
    sp = _species(ts)
    rows = []
    for it in ts.iterations:
        x, y, z, ux, uy, uz, w = ts.get_particle(
            ["x", "y", "z", "ux", "uy", "uz", "w"], species=sp, iteration=it)
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
    """(z, ke, dke, enx, eny, sx, sy, charge_pc) [m / MeV / mm-as-m units kept SI; charge pC] from
    the summary's stat_vs_z (preferred) or the particle slices. `charge_pc` is None when the summary
    predates the surviving-charge column. Returns None if neither source is usable.
    """
    svz = summ.get("stat_vs_z", {})
    if svz.get("z_m"):
        z = np.array(svz["z_m"])
        charge = np.array(svz["charge_pc"]) if svz.get("charge_pc") else None
        return (z, np.array(svz["ke_mev"]), np.array(svz["sigma_ke_mev"]),
                np.array(svz["norm_emit_x"]), np.array(svz["norm_emit_y"]),
                np.array(svz["sigma_x_m"]),
                np.array(svz.get("sigma_y_m", svz["sigma_x_m"])), charge)
    rows = _read_slices(diag)
    if not rows:
        return None
    return (_arr(rows, "z"), _arr(rows, "ke"), _arr(rows, "dke"),
            _arr(rows, "enx"), _arr(rows, "eny"), _arr(rows, "sx"), _arr(rows, "sy"),
            _arr(rows, "q") * 1e12)


def _achieved_de(calib, z, ke):
    """Per-section achieved ΔE [MeV] from the vs-z KE curve: KE at each section's exit z minus KE
    at its entry z. The section z-edges come from the calibration table (`z_entry_m`/`z_exit_m`,
    the real deck geometry written by sim/linac5-8.py); a legacy summary without them falls back to
    an even split of the z-grid (coarse -- it mis-attributes gain across the inter-section drifts).
    """
    n = len(calib)
    if n == 0 or len(z) < 2:
        return []
    ke_at = lambda zz: float(np.interp(zz, z, ke))
    if all("z_entry_m" in c and "z_exit_m" in c for c in calib):
        return [ke_at(c["z_exit_m"]) - ke_at(c["z_entry_m"]) for c in calib]
    edges = np.linspace(float(z[0]), float(z[-1]), n + 1)
    return [ke_at(edges[i + 1]) - ke_at(edges[i]) for i in range(n)]


def _exit_slice(diag):
    """Particle arrays (x, y, z, ux, uy, uz, w) of the exit dump (largest <z>), or None."""
    from openpmd_viewer import OpenPMDTimeSeries
    parts = os.path.join(diag, "particles")
    if not os.path.isdir(parts):
        return None
    try:
        ts = OpenPMDTimeSeries(parts)
    except Exception:
        return None
    if not list(ts.iterations):
        return None
    sp = _species(ts)
    best = None
    for it in ts.iterations:
        x, y, z, ux, uy, uz, w = ts.get_particle(
            ["x", "y", "z", "ux", "uy", "uz", "w"], species=sp, iteration=it)
        if len(z) < 50:
            continue
        zc = float(np.average(z, weights=w))
        if best is None or zc > best[0]:
            best = (zc, (x, y, z, ux, uy, uz, w))
    return best[1] if best else None


def _save_exit_figures(x, y, z, ux, uy, uz, w):
    """Energy spectrum, longitudinal (z, KE) and transverse (r, p_r) phase space at the exit."""
    pg = make_particle_group(x, y, z, ux, uy, uz, w)
    fig = px.energy_spectrum(pg, use_ke=True, e_unit="MeV")
    fig.savefig(os.path.join(RESULTS, "energy_spectrum.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)

    ke = px.ke_kev_from_u(ux, uy, uz) / 1e3                  # keV -> MeV
    zc_mm = (z - np.average(z, weights=w)) * 1e3
    fig, ax = plt.subplots(figsize=(7.0, 4.8), constrained_layout=True)
    hb = ax.hexbin(zc_mm, ke, gridsize=70, cmap="viridis", mincnt=1)
    fig.colorbar(hb, ax=ax, label="macroparticles / bin")
    ax.set_xlabel("z - <z>  [mm]"); ax.set_ylabel("kinetic energy  [MeV]")
    ax.set_title("linac5-8 exit longitudinal phase space  (z, KE)")
    fig.savefig(os.path.join(RESULTS, "phase_space_z_KE.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)

    fig = px.transverse_rpr(x, y, ux, uy, w,
                            title="linac5-8 exit transverse phase space  (r, p_r)",
                            p_unit="MeV")
    fig.savefig(os.path.join(RESULTS, "transverse_r_pr.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)


def main():
    diag = DIAG_DIR
    os.makedirs(RESULTS, exist_ok=True)
    summ = _load_summary(diag)
    calib = summ.get("calibration", [])

    vs = _vs_z(diag, summ)
    if vs is None:
        print(f"plot linac5-8: no stat_vs_z and no usable dumps in {diag} -- skipping.",
              flush=True)
        return
    z, ke, dke, enx, eny, sx, sy, charge = vs
    power_mw = summ.get("power_mw", 11.0)

    # 1) beam evolution vs z: mean KE / eps_n,x / sigma_x (+ surviving charge when available)
    fig = px.evolution_vs_z(
        z, ke, enx * 1e6, sx * 1e3, charge_pc=charge, ke_unit="MeV",
        title=f"linac5-8 beam evolution (sections 5-8, on-crest, {power_mw:g} MW)",
        notes={"emit": "quads OFF: eps_n rises ~2.4x -- a fort.10N diagnostic artifact, not physical",
               "sigma": "quads OFF: no focusing, placeholder optics, NOT predictive",
               "charge": "surviving core charge (macro count x q/macro); quads OFF -> aperture loss"})
    fig.savefig(os.path.join(RESULTS, "evolution_vs_z.png"), dpi=130)
    plt.close(fig)

    # 2) energy spread vs z (absolute grows on the cosine crest curvature; relative shrinks)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.6, 6.2), constrained_layout=True, sharex=True)
    a1.plot(z, dke, "-o", ms=3, color="C5")
    a1.set_ylabel("sigma_KE [MeV]")
    a1.set_title("linac5-8: energy spread (absolute grows, relative shrinks)")
    a1.grid(alpha=0.3)
    rel = np.where(ke > 0, dke / ke * 100.0, np.nan)
    a2.plot(z, rel, "-o", ms=3, color="C6")
    a2.set_ylabel("sigma_KE/<KE> [%]")
    a2.set_xlabel("<z> [m]")
    a2.grid(alpha=0.3)
    fig.savefig(os.path.join(RESULTS, "energy_spread.png"), dpi=130)
    plt.close(fig)

    # 3) per-section achieved vs frozen-target ΔE
    fig, ax = plt.subplots(figsize=(9.2, 4.8), constrained_layout=True)
    if calib:
        names = [c.get("name", f"sec{c['index'] + 5}") for c in calib]
        tgt = [c.get("target_de_mev", np.nan) for c in calib]
        ach = _achieved_de(calib, z, ke)
        xi = np.arange(len(calib))
        ax.bar(xi - 0.2, tgt, 0.4, label="frozen target dE", color="0.6")
        if ach:
            ax.bar(xi + 0.2, ach, 0.4, label="achieved dE (from vs-z KE)", color="C5")
        ax.set_xticks(xi)
        ax.set_xticklabels(names, rotation=30, ha="right", fontsize=8)
        ax.set_ylabel("dE per section [MeV]")
        ax.set_title("linac5-8: per-section gain -- frozen target vs achieved")
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "no calibration table in injection_summary.json",
                ha="center", va="center", transform=ax.transAxes)
    ax.grid(alpha=0.3, axis="y")
    fig.savefig(os.path.join(RESULTS, "section_gains.png"), dpi=130)
    plt.close(fig)

    # 4-6) exit-slice figures: energy spectrum, longitudinal (z, KE) and transverse (r, p_r)
    n_slice = _exit_slice(diag)
    if n_slice is not None:
        _save_exit_figures(*n_slice)
    else:                                          # no exit beam (e.g. near-total loss): drop stale
        for f in ("energy_spectrum", "phase_space_z_KE", "transverse_r_pr"):
            p = os.path.join(RESULTS, f + ".png")
            if os.path.exists(p):
                os.remove(p)

    last_ke = f"{ke[-1]:.1f} MeV" if len(ke) else "n/a"
    print(f"plot linac5-8: wrote figures to {RESULTS}/ "
          f"({len(z)} vs-z points, exit <KE> {last_ke}).", flush=True)


if __name__ == "__main__":
    main()
