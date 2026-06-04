"""
Figures and summary for the WarpX prebuncher scan (repeated prebuncher_sim.py
runs, one --outdir per power/phase).

The gun-exit bunch is short (σ_z ≈ 1 mm ≈ 0.1 % of the 214 MHz RF wavelength) and
carries an intrinsic +1.4 keV/mm (debunching) energy chirp, and at 0.1 nC it is
space-charge dense. In a free drift it therefore *expands*. The prebuncher acts as
a ballistic buncher: at the zero-crossing it flips the chirp negative and
compresses the bunch — but space charge limits how far. The clearest way to show
this is to compare each powered run against a **drift-only baseline** (P = 0):

  * σ_z(z)        — bunch length along the line, cavity runs vs. the drift baseline
  * ratio(z)      — σ_z,drift(z) / σ_z,cavity(z)  (>1 ⇒ cavity is bunching)
  * I_peak(z)     — peak current (peak line density × beam velocity)
  * z–KE space    — the chirp flipping/rotating through the cavity

(The bunching factor at the RF fundamental is ≈1 and flat here, so it is not used.)

Run with:
    conda run -n CBB python warpx_prebuncher/plot_prebuncher.py
"""

import os
import glob
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from openpmd_viewer import OpenPMDTimeSeries

c = 299792458.0
MC2 = 0.51099895e3                # electron rest energy [keV]
DIAG_ROOT = "warpx_prebuncher/diags"
RESULTS = "warpx_prebuncher/results"
Z_GAP_CENTER = 0.20              # [m] cavity gap (for marking plots)
os.makedirs(RESULTS, exist_ok=True)


def wstats(v, w):
    m = np.average(v, weights=w)
    return m, np.sqrt(np.average((v - m) ** 2, weights=w))


def peak_current(z, w, v_beam, nbins=400):
    zlo, zhi = z.min(), z.max()
    if zhi <= zlo:
        return 0.0
    edges = np.linspace(zlo, zhi, nbins + 1)
    dz = edges[1] - edges[0]
    q_e = 1.602176634e-19
    lam, _ = np.histogram(z, bins=edges, weights=w * q_e)
    return float(lam.max() / dz * v_beam)


def analyse_case(path):
    """Per-snapshot metric arrays for one case directory (sorted by ⟨z⟩)."""
    ts = OpenPMDTimeSeries(os.path.join(path, "particles"))
    rec = dict(zmean=[], sigz=[], ke=[], dke=[], ipk=[], it=[])
    snaps = {}
    v_beam = None
    for i, it in enumerate(ts.iterations):
        z, ux, uy, uz, w = ts.get_particle(
            ["z", "ux", "uy", "uz", "w"], species="electrons", iteration=it)
        if len(z) < 50:                       # skip near-empty (boundary) dumps
            continue
        gam = np.sqrt(1.0 + ux**2 + uy**2 + uz**2)
        ke = (gam - 1.0) * MC2
        if v_beam is None:
            v_beam = float(np.average(uz / gam, weights=w) * c)
        zm, sz = wstats(z, w)
        km, dk = wstats(ke, w)
        rec["zmean"].append(zm); rec["sigz"].append(sz)
        rec["ke"].append(km); rec["dke"].append(dk)
        rec["ipk"].append(peak_current(z, w, v_beam)); rec["it"].append(it)
        snaps[it] = (z, ke, w)
    order = np.argsort(rec["zmean"])
    for k in ("zmean", "sigz", "ke", "dke", "ipk"):
        rec[k] = np.asarray(rec[k])[order]
    rec["it"] = [rec["it"][i] for i in order]
    rec["v_beam"] = v_beam or 0.632 * c
    rec["snaps"] = snaps
    return rec


def case_label(name):
    if name == "P0_drift":
        return (0, "drift")
    m = re.match(r"P(\d+)_(zc|crest)$", name)
    return (int(m.group(1)), m.group(2)) if m else None


def per_case_figure(name, rec, base):
    zmm = rec["zmean"] * 1e3
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.4), constrained_layout=True)
    a1.plot(zmm, rec["sigz"] * 1e3, "o-", ms=3, color="C0", label="cavity")
    if base is not None:
        sd = np.interp(rec["zmean"], base["zmean"], base["sigz"])
        a1.plot(zmm, sd * 1e3, "k--", lw=1, label="drift baseline")
        ratio = sd / rec["sigz"]
        ib = int(np.argmax(ratio))
        a1.plot(zmm[ib], rec["sigz"][ib] * 1e3, "*", ms=14, color="C1",
                label=f"max bunching {ratio[ib]:.2f}× @ {zmm[ib]:.0f} mm")
    a1.axvline(Z_GAP_CENTER * 1e3, color="C3", ls=":", label="gap")
    a1.set_xlabel("⟨z⟩  [mm]"); a1.set_ylabel("σ_z  [mm]")
    a1.set_title(f"{name}: bunch length vs. drift"); a1.legend(fontsize=8)

    a2b = a2.twinx()
    a2.plot(zmm, rec["ipk"], "o-", ms=3, color="C2")
    a2b.plot(zmm, rec["ke"], "s--", ms=3, color="C4")
    a2.axvline(Z_GAP_CENTER * 1e3, color="C3", ls=":")
    a2.set_xlabel("⟨z⟩  [mm]")
    a2.set_ylabel("peak current  [A]", color="C2")
    a2b.set_ylabel("mean KE  [keV]", color="C4")
    a2.set_title(f"{name}: peak current & mean energy")
    fig.savefig(f"{RESULTS}/{name}_line.png", dpi=140); plt.close(fig)

    # z–KE phase space at injection / gap exit / best-bunching point.
    snaps, its = rec["snaps"], rec["it"]
    def nearest_it(ztarget):
        return its[int(np.argmin(np.abs(rec["zmean"] - ztarget)))]
    zbest = zmm[int(np.argmax(np.interp(rec["zmean"], base["zmean"], base["sigz"])
                              / rec["sigz"]))] / 1e3 if base is not None else rec["zmean"][-1]
    picks = [its[0], nearest_it(Z_GAP_CENTER + 0.06), nearest_it(zbest)]
    titles = ["injection", "cavity exit", "max bunching"]
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.0), constrained_layout=True)
    for ax, it, ti in zip(axs, picks, titles):
        z, ke, w = snaps[it]
        ax.scatter((z - z.mean()) * 1e3, ke - ke.mean(), s=2, alpha=0.15, color="C0")
        zi = rec["it"].index(it)
        ax.set_xlabel("z − ⟨z⟩  [mm]"); ax.set_ylabel("KE − ⟨KE⟩  [keV]")
        ax.set_title(f"{ti}  (⟨z⟩={rec['zmean'][zi]*1e3:.0f} mm)")
    fig.suptitle(f"{name}: longitudinal phase space", fontsize=12)
    fig.savefig(f"{RESULTS}/{name}_phasespace.png", dpi=140); plt.close(fig)

    out = dict(szmin=float(rec["sigz"].min()),
               sz0=float(rec["sigz"][0]),
               ke_end=float(rec["ke"][-1]), ipk_max=float(rec["ipk"].max()))
    if base is not None:
        sd = np.interp(rec["zmean"], base["zmean"], base["sigz"])
        ratio = sd / rec["sigz"]; ib = int(np.argmax(ratio))
        out.update(ratio=float(ratio[ib]), zbest=float(rec["zmean"][ib]),
                   sz_best=float(rec["sigz"][ib]))
    return out


def main():
    import argparse
    ap = argparse.ArgumentParser(
        description="Plot prebuncher case(s). With no arguments, plots every "
                    f"{DIAG_ROOT}/P* directory; pass case names to restrict.")
    ap.add_argument("cases", nargs="*",
                    help="case dir names to plot, e.g. P800_zc (default: all)")
    a = ap.parse_args()

    if a.cases:
        dirs = [os.path.join(DIAG_ROOT, c) for c in a.cases]
    else:
        dirs = sorted(glob.glob(f"{DIAG_ROOT}/P*"))
    cases = [(d, case_label(os.path.basename(d))) for d in dirs if os.path.isdir(d)]
    cases = [(d, lab) for d, lab in cases if lab]
    if not cases:
        print(f"No case directories found under {DIAG_ROOT}/ "
              f"(run prebuncher_sim.py first).")
        return

    base = None
    for d, lab in cases:
        if lab[1] == "drift":
            print("analysing drift baseline …", flush=True)
            base = analyse_case(d)

    summary = []
    for d, (power, phase) in sorted(cases, key=lambda x: (x[1][1], x[1][0])):
        if phase == "drift":
            continue
        name = os.path.basename(d)
        print(f"analysing {name} …", flush=True)
        rec = analyse_case(d)
        s = per_case_figure(name, rec, base)
        s.update(name=name, power=power, phase=phase, rec=rec)
        summary.append(s)

    # ── Headline comparison (only meaningful with several cases / a baseline) ──
    if len(summary) > 1 or base is not None:
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
        if base is not None:
            a1.plot(base["zmean"] * 1e3, base["sigz"] * 1e3, "k--", lw=2,
                    label="drift (P=0)")
        for s in sorted([s for s in summary if s["phase"] == "zc"],
                        key=lambda s: s["power"]):
            r = s["rec"]
            a1.plot(r["zmean"] * 1e3, r["sigz"] * 1e3, "-", label=f"{s['power']:.0f} W")
        a1.axvline(Z_GAP_CENTER * 1e3, color="C3", ls=":")
        a1.set_xlabel("⟨z⟩  [mm]"); a1.set_ylabel("σ_z  [mm]")
        a1.set_title("Bunch length: drift vs. zero-crossing cavity"); a1.legend(fontsize=8)

        for phase, col in (("zc", "C0"), ("crest", "C3")):
            pts = sorted([s for s in summary if s["phase"] == phase and "ratio" in s],
                         key=lambda s: s["power"])
            if pts:
                a2.plot([s["power"] for s in pts], [s["ratio"] for s in pts], "o-",
                        color=col, label="zero-crossing" if phase == "zc" else "on-crest")
        a2.axhline(1.0, color="k", lw=0.6)
        a2.set_xlabel("RF power P  [W]")
        a2.set_ylabel("max bunching  σ_drift / σ_cavity")
        a2.set_title("Bunching vs. power"); a2.legend(fontsize=8)
        fig.savefig(f"{RESULTS}/compare_power_phase.png", dpi=140); plt.close(fig)
        print(f"wrote {RESULTS}/compare_power_phase.png")
    else:
        print("single case, no drift baseline -> skipping cross-case comparison figure")

    # ── Summary table ─────────────────────────────────────────────────────────
    print("\n" + "=" * 96)
    print(f"{'case':>10} {'P[W]':>5} {'phase':>6} {'σz0[mm]':>8} {'σzmin[mm]':>10} "
          f"{'bunch×':>7} {'zbest[mm]':>9} {'Ipk[A]':>7} {'KEend[keV]':>11}")
    print("-" * 96)
    if base is not None:
        print(f"{'P0_drift':>10} {0:5d} {'drift':>6} {base['sigz'][0]*1e3:8.3f} "
              f"{base['sigz'].min()*1e3:10.3f} {'—':>7} {'—':>9} "
              f"{base['ipk'].max():7.2f} {base['ke'][-1]:11.1f}")
    for s in sorted(summary, key=lambda s: (s["phase"], s["power"])):
        ratio = f"{s['ratio']:7.2f}" if "ratio" in s else f"{'—':>7}"
        zbest = f"{s['zbest']*1e3:9.0f}" if "zbest" in s else f"{'—':>9}"
        print(f"{s['name']:>10} {s['power']:5.0f} {s['phase']:>6} {s['sz0']*1e3:8.3f} "
              f"{s['szmin']*1e3:10.3f} {ratio} {zbest} {s['ipk_max']:7.2f} "
              f"{s['ke_end']:11.1f}")
    print("=" * 96)


if __name__ == "__main__":
    main()
