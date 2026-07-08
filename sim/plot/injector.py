"""Figures for the WarpX RZ CESR injector (sim/injector.py) over logs/diags/injector/main/.
Writes PNGs to logs/plots/injector/. See docs/injector.md for the physics each figure shows.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import openpmd_api as io
from openpmd_viewer import OpenPMDTimeSeries

from sim.plot import common as px
from sim.helpers.tools import C_LIGHT, E_CHARGE, MC2_KEV, prepare_env
from sim.helpers.buildfields import Z_GAP_CENTER_1, Z_GAP_CENTER_2, V1J_KEV, INJ_Z_HANDOFF as Z_HANDOFF

CONFIG = "config/injector.yaml"
RESULTS = "logs/plots/injector"
DIAG = "logs/diags/injector/main/particles"
PREB1_FIELD = "fieldmaps/h5/preb1_EB.h5"
PREB2_FIELD = "fieldmaps/h5/preb2_EB.h5"


def _last_populated(diag, species="electrons"):
    ts = OpenPMDTimeSeries(diag)
    for it in reversed([int(i) for i in ts.iterations]):
        try:
            x, = ts.get_particle(["x"], species=species, iteration=it)
        except Exception:
            continue
        if len(x):
            return it
    return int(ts.iterations[-1])


def _save(fig, name):
    fig.savefig(f"{RESULTS}/injector_{name}.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {RESULTS}/injector_{name}.png")


def _wmean_std(v, w):
    mean = np.average(v, weights=w)
    return mean, np.sqrt(np.average((v - mean) ** 2, weights=w))


def _peak_current(z, w, v_beam, n_bins=400):
    """Peak longitudinal current I = max(lambda)*v_beam from the line-charge density lambda(z)."""
    if z.max() <= z.min():
        return 0.0
    charge, edges = np.histogram(z, bins=n_bins, weights=w * E_CHARGE)
    return float(charge.max() / (edges[1] - edges[0]) * v_beam)


def rf_scale(power_kw, q_l, f_rf):
    """Field scale sqrt(stored energy / 1 J), stored energy = 1e3*Q*P/(2pi f_RF)."""
    return float(np.sqrt(1e3 * q_l * power_kw / (2.0 * np.pi * f_rf))) if power_kw > 0 else 0.0


def analyse(diag):
    """Sorted by <z>: station_picks/line_figure rely on ascending order."""
    ts = OpenPMDTimeSeries(diag)
    rec = {k: [] for k in ("zmean", "sigz", "ke", "dke", "ipk", "it")}
    snaps = {}
    v_beam = None
    for it in ts.iterations:
        z, ux, uy, uz, w = ts.get_particle(
            ["z", "ux", "uy", "uz", "w"], species="electrons", iteration=it)
        if len(z) < 50:                                     # skip near-empty boundary dumps
            continue
        gamma = px.gamma_from_u(ux, uy, uz)
        ke = (gamma - 1.0) * MC2_KEV
        if v_beam is None:
            v_beam = float(np.average(uz / gamma, weights=w) * C_LIGHT)
        zmean, sigz = _wmean_std(z, w)
        kmean, dke = _wmean_std(ke, w)
        rec["zmean"].append(zmean); rec["sigz"].append(sigz)
        rec["ke"].append(kmean); rec["dke"].append(dke)
        rec["ipk"].append(_peak_current(z, w, v_beam)); rec["it"].append(it)
        snaps[it] = (z, ke, w)
    if not rec["zmean"]:
        return None, None, None
    order = np.argsort(rec["zmean"])
    for k in ("zmean", "sigz", "ke", "dke", "ipk"):
        rec[k] = np.asarray(rec[k])[order]
    rec["it"] = [rec["it"][i] for i in order]
    return rec, snaps, v_beam


def station_picks(rec):
    its, zmean = rec["it"], rec["zmean"]
    nearest = lambda z: its[int(np.argmin(np.abs(zmean - z)))]
    post = zmean > Z_GAP_CENTER_2
    z_focus = zmean[np.where(post)[0][np.argmin(rec["sigz"][post])]] if post.any() else zmean[-1]
    picks = [its[0], nearest(Z_GAP_CENTER_2 + 0.06), nearest(z_focus), nearest(Z_HANDOFF)]
    titles = ["injection", "after Preb 2", "best focus (min sigma_z)", "handoff (2.03 m)"]
    return picks, titles


def line_figure(rec):
    z_mm = rec["zmean"] * 1e3
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.4), constrained_layout=True)
    a1.plot(z_mm, rec["sigz"] * 1e3, "o-", ms=3, color="C0")
    for z_gap, c in ((Z_GAP_CENTER_1, "C3"), (Z_GAP_CENTER_2, "C5")):
        a1.axvline(z_gap * 1e3, color=c, ls=":")
    a1.set_xlabel("<z>  [mm]"); a1.set_ylabel(r"$\sigma_z$  [mm]")
    a1.set_title("Bunch length along the line (velocity bunching -> waist)")

    a2b = a2.twinx()
    a2.plot(z_mm, rec["ipk"], "o-", ms=3, color="C2")
    a2b.plot(z_mm, rec["ke"], "s--", ms=3, color="C4")
    for z_gap in (Z_GAP_CENTER_1, Z_GAP_CENTER_2):
        a2.axvline(z_gap * 1e3, color="0.6", ls=":")
    a2.set_xlabel("<z>  [mm]")
    a2.set_ylabel("peak current  [A]", color="C2"); a2.tick_params(axis="y", labelcolor="C2")
    a2b.set_ylabel("mean KE  [keV]", color="C4"); a2b.tick_params(axis="y", labelcolor="C4")
    a2.set_title("Peak current and mean energy")
    _save(fig, "line")


def bunch_profile_figure(rec, snaps, picks, titles):
    fig, axs = plt.subplots(1, len(picks), figsize=(3.4 * len(picks), 4.0), constrained_layout=True)
    for ax, it, title in zip(np.atleast_1d(axs), picks, titles):
        z, _, w = snaps[it]
        zc = z - np.average(z, weights=w)
        sigz = np.sqrt(np.average(zc ** 2, weights=w))
        span = max(4.0 * sigz, 5e-4)
        charge, edges = np.histogram(zc, bins=np.linspace(-span, span, 121), weights=w * E_CHARGE)
        lam = charge / (edges[1] - edges[0]) * 1e9         # C/m -> nC/m
        centres = 0.5 * (edges[:-1] + edges[1:])
        ax.fill_between(centres * 1e3, lam, color="C0", alpha=0.25)
        ax.plot(centres * 1e3, lam, color="C0", lw=1.4)
        zmean = rec["zmean"][rec["it"].index(it)]
        ax.set_title(f"{title}\n(<z>={zmean*1e3:.0f} mm)")
        ax.set_xlabel("z - <z>  [mm]"); ax.set_ylabel("lambda  [nC/m]")
        ax.annotate(f"peak {lam.max():.2f} nC/m\nsigma_z = {sigz*1e3:.2f} mm",
                    xy=(0.96, 0.95), xycoords="axes fraction", ha="right", va="top",
                    fontsize=8, bbox=dict(boxstyle="round", fc="white", alpha=0.8))
    fig.suptitle("Longitudinal line-charge density lambda(z)", fontsize=12)
    _save(fig, "bunch_profile")


def phasespace_figure(rec, snaps, picks, titles):
    def wpercentile(v, w, q):
        order = np.argsort(v)
        cdf = np.cumsum(w[order]); cdf /= cdf[-1]
        return np.interp(q, cdf, v[order])

    fig, axs = plt.subplots(1, len(picks), figsize=(3.6 * len(picks), 4.0), constrained_layout=True)
    for ax, it, title in zip(np.atleast_1d(axs), picks, titles):
        z, ke, w = snaps[it]
        zc = (z - np.average(z, weights=w)) * 1e3
        kc = ke - np.average(ke, weights=w)
        h = ax.hist2d(zc, kc, bins=120, weights=w * E_CHARGE * 1e9,
                      cmap="inferno", cmin=np.finfo(float).tiny)
        fig.colorbar(h[3], ax=ax, label="charge  [nC/bin]", fraction=0.046, pad=0.02)
        for set_lim, vals in ((ax.set_xlim, zc), (ax.set_ylim, kc)):
            lo, hi = wpercentile(vals, w, 1e-3), wpercentile(vals, w, 1 - 1e-3)
            pad = 0.05 * (hi - lo) or 1.0
            set_lim(lo - pad, hi + pad)
        zmean = rec["zmean"][rec["it"].index(it)]
        ax.set_title(f"{title}  (<z>={zmean*1e3:.0f} mm)")
        ax.set_xlabel("z - <z>  [mm]"); ax.set_ylabel("KE - <KE>  [keV]")
    fig.suptitle("Longitudinal phase space along the injector", fontsize=12)
    _save(fig, "phasespace")


def _on_axis_ez(field_path):
    """Field map is stored normalized to 1 J; caller must scale by rf_scale()."""
    s = io.Series(field_path, io.Access.read_only)
    mesh = s.iterations[0].meshes["E"]
    ez = mesh["z"].load_chunk()
    s.flush()
    dz = mesh.grid_spacing[1]
    z0 = mesh.grid_global_offset[1] if mesh.grid_global_offset else 0.0
    return z0 + np.arange(ez[0].shape[1]) * dz, ez[0][0]


def cavity_figure(w):
    p1, p2 = w.get("params/PREB1_KW"), w.get("params/PREB2_KW")
    f_rf, q1, q2 = w.get("params/F_RF"), w.get("params/Q_L_1"), w.get("params/Q_L_2")
    scale1, scale2 = rf_scale(p1, q1, f_rf), rf_scale(p2, q2, f_rf)
    z1, ez1 = _on_axis_ez(PREB1_FIELD)
    z2, ez2 = _on_axis_ez(PREB2_FIELD)

    fig, ax = plt.subplots(figsize=(8.6, 4.4), constrained_layout=True)
    ax.plot(z1 * 1e3, ez1 * scale1 / 1e6, color="C3",
            label=f"Preb 1 ({p1:g} kW, $V_g$~={scale1*V1J_KEV:.0f} kV)")
    if scale2 > 0:
        ax.plot(z2 * 1e3, ez2 * scale2 / 1e6, color="C4",
                label=f"Preb 2 reversed ({p2:g} kW, $V_g$~={scale2*V1J_KEV:.0f} kV)")
    ax.axhline(0, color="k", lw=0.6)
    ax.axvline(Z_HANDOFF * 1e3, color="C2", ls="--", label="handoff (2.03 m)")
    ax.set_xlabel("lab z  [mm]"); ax.set_ylabel(r"on-axis $E_z \times$ scale  [MV/m]")
    ax.set_title("Prebuncher RF field lobes")
    ax.legend(fontsize=9)
    _save(fig, "cavity")


def main():
    prepare_env()
    from warpx import WarpX
    os.makedirs(RESULTS, exist_ok=True)

    w = WarpX(input_file=CONFIG, path="logs/diags/injector")
    w.load_output(diag_dir=DIAG)
    it = _last_populated(DIAG)
    pg = w._particle_group(iteration=it)

    for name, fig in [
        ("phase_space_z_KE", w.plot2D("z", "kinetic_energy", iteration=it)),
        ("energy_spectrum",  px.energy_spectrum(pg)),
    ]:
        _save(fig, name)

    ts = OpenPMDTimeSeries(DIAG)
    x, y, ux, uy, wgt = ts.get_particle(["x", "y", "ux", "uy", "w"],
                                        species="electrons", iteration=it)
    _save(px.transverse_rpr(x, y, ux, uy, wgt,
                            title="Injector handoff transverse phase space  (r, p_r)"),
          "transverse_r_pr")

    z_m, ke, emit, sigma, q_pc = px.evolution_screens(px.pool_trajectories(ts, ts.iterations))
    _save(px.evolution_vs_z(z_m, ke, emit, sigma, charge_pc=q_pc,
                            title="Beam evolution along the injector  (fixed-z virtual screens)"),
          "evolution_vs_z")

    cavity_figure(w)
    rec, snaps, _ = analyse(DIAG)
    if rec is not None:
        picks, titles = station_picks(rec)
        line_figure(rec)
        bunch_profile_figure(rec, snaps, picks, titles)
        phasespace_figure(rec, snaps, picks, titles)


if __name__ == "__main__":
    main()
