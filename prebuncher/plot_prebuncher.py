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

Per-case figures use CONFIG-INDEPENDENT filenames (the power/phase lives in the
diags/<case> input dir and the figure titles, not the filename) so changing the
operating point overwrites the same files instead of leaving orphans:

  * prebuncher_line.png        — σ_z(z) (vs. drift) and peak current / mean energy
  * prebuncher_phasespace.png  — z–KE at injection / cavity exit / best focus
  * prebuncher_cavity.png      — the RF DRIVE: on-axis Ez(z) of the scaled 1-J map in
                             the lab frame, plus the cos/sin RF waveform vs. time
                             bracketing the gap arrival, showing whether the bunch
                             centre lands on the field zero-crossing (zc) or crest.
  * prebuncher_bunch_profile.png — the real longitudinal line-charge density λ(z) at
                             the same three snapshots as the phase-space figure,
                             exposing the compression and any space-charge
                             spike/filamentation the scalar σ_z curve cannot show.

With several diags/P* cases present these per-case figures are overwritten (last
case wins) — use compare_power_phase.png for the cross-case scan summary.

Run with:
    conda run -n CBB python prebuncher/plot_prebuncher.py
"""

import os
import glob
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import openpmd_api as io
from openpmd_viewer import OpenPMDTimeSeries

c = 299792458.0
MC2 = 0.51099895e3                # electron rest energy [keV]
Q_E = 1.602176634e-19            # elementary charge [C]
DIAG_ROOT = "prebuncher/diags"
RESULTS = "prebuncher/results"
PREBUNCH_FIELD = "prebuncher/prebuncher_field/prebuncher_EB.h5"
Z_GAP_CENTER = 0.20              # [m] cavity gap (for marking plots)

# ── RF-drive constants (must match prebuncher_sim.py) ─────────────────────────
# The cavity map is stored 1-J-normalised; the run multiplies it by `scale` and
# modulates E ∝ cos(ω t + φ), B ∝ sin(ω t + φ). To draw the field the beam
# actually sees we re-derive scale and φ here from the same formulae the sim uses.
F_RF = 499.7645e6 / 42 * 18      # 18 × master RF = 214.18 MHz (details.md)
OMEGA = 2.0 * np.pi * F_RF       # RF angular frequency [rad/s]
Q_L = 3000                       # loaded Q of prebuncher 1 (details.md)
Z_INJECT = 0.005                 # [m] lab z where the bunch head is launched
V1J_KEV = 430.2                  # 1-J effective gap voltage [keV] (for V_gap label)
os.makedirs(RESULTS, exist_ok=True)


def rf_scale(power):
    """Field scale = sqrt(stored_energy / 1 J), stored_energy = 1e3·Q·P/(2π f_RF).

    Same expression as prebuncher_sim.py: at P=0 (drift baseline) there is no
    cavity field, so the scale is 0.
    """
    if power <= 0:
        return 0.0
    return float(np.sqrt(1e3 * Q_L * power / (2.0 * np.pi * F_RF)))


def rf_phase(phase, t_gap):
    """RF phase φ that puts the bunch-centre gap arrival at zero-crossing/crest.

    Mirrors prebuncher_sim.py: zc → φ = -ω t_gap + π/2 (max +slope, zero net
    kick → velocity bunching); crest → φ = -ω t_gap + π (-cos = +1 → max gain).
    """
    if phase == "crest":
        return -OMEGA * t_gap + np.pi
    return -OMEGA * t_gap + np.pi / 2.0


def load_cavity_axis():
    """On-axis Ez(z) of the raw 1-J cavity map, in LAB z.

    Read the way gun/plot_gun.py reads gun_E.h5: io.Series, mesh 'E',
    component 'z', grid_spacing, and grid_global_offset (set by
    build_prebuncher_field.py so the gap-centred map lands at Z_GAP_CENTER).
    Returns (z_lab [m], ez_axis [V/m]) for the r = 0 row of the (1, nr, nz) mesh.
    """
    s = io.Series(PREBUNCH_FIELD, io.Access.read_only)
    E = s.iterations[0].meshes["E"]
    ez = E["z"].load_chunk()
    s.flush()
    ez = ez[0]                                   # (nr, nz), mode 0
    dz = E.grid_spacing[1]
    # axis_labels = ["r", "z"] -> grid_global_offset[1] is the z origin in the lab.
    z0 = E.grid_global_offset[1] if E.grid_global_offset else (Z_GAP_CENTER - 0.1524)
    z_lab = z0 + np.arange(ez.shape[1]) * dz
    return z_lab, ez[0]                          # r = 0 row


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


def snapshot_picks(rec, base):
    """The three snapshots used by the phase-space and bunch-profile figures:
    injection / cavity exit / best bunching. Factored out so both figures sample
    the SAME iterations.

    With a drift baseline the third pick is the σ_drift/σ_cavity maximum ("max
    bunching"). Without one (the only data on disk here), use the *post-cavity*
    σ_z minimum — the true ballistic focus — not the last dump, which is the
    re-expanded exit beam and no focus at all. (A global σ_z minimum would land
    back at injection, since the bunch is shortest *before* the cavity kick
    expands it; the meaningful focus is the σ_z dip downstream of the gap.)
    """
    its, zmean = rec["it"], rec["zmean"]
    def nearest_it(ztarget):
        return its[int(np.argmin(np.abs(zmean - ztarget)))]
    if base is not None:
        ib = int(np.argmax(np.interp(zmean, base["zmean"], base["sigz"])
                           / rec["sigz"]))
        zbest, third = zmean[ib], "max bunching"
    else:
        post = zmean > Z_GAP_CENTER                # downstream of the cavity gap
        if post.any():
            ib = np.where(post)[0][int(np.argmin(rec["sigz"][post]))]
            zbest, third = zmean[ib], "best focus (min σ_z)"
        else:
            zbest, third = zmean[-1], "exit"
    picks = [its[0], nearest_it(Z_GAP_CENTER + 0.06), nearest_it(zbest)]
    return picks, ["injection", "cavity exit", third]


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE A — prebuncher_cavity.png : the RF DRIVE the bunch actually sees.
#
# The σ_z / phase-space figures show the beam's *response* but never the cavity
# field itself — there is no prebuncher analogue of the gun's gun_field.png. This
# makes the drive visible for the case's (power, phase):
#   left  — the on-axis Ez(z) spatial profile of the 1-J map × scale, placed at
#           the lab gap (Z_GAP_CENTER) via the map's grid_global_offset, in MV/m.
#   right — the temporal RF waveform E ∝ cos(ω t+φ), B ∝ sin(ω t+φ) over ~2 RF
#           periods around the bunch-centre gap-arrival time t_gap, each normalised
#           to ±1, with the bunch σ_t width shaded. Whether the bunch centre lands
#           on the field ZERO-CROSSING (zc) or the CREST (crest) — i.e. velocity
#           bunching vs. pure acceleration — is exactly what this panel makes visual.
# ══════════════════════════════════════════════════════════════════════════════
def cavity_figure(name, rec, power, phase):
    v_beam = rec["v_beam"]
    scale = rf_scale(power)
    t_gap = (Z_GAP_CENTER - Z_INJECT) / v_beam     # bunch-centre gap arrival [s]
    phi = rf_phase(phase, t_gap)

    z_lab, ez_axis = load_cavity_axis()            # raw 1-J on-axis Ez [V/m]
    ez_scaled = ez_axis * scale                    # field amplitude this case sees

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.4), constrained_layout=True)

    # ── Left: spatial Ez(z) of the scaled cavity map, in the lab frame ──────────
    a1.plot(z_lab * 1e3, ez_scaled / 1e6, color="C3")
    a1.axhline(0, color="k", lw=0.6)
    a1.axvline(Z_GAP_CENTER * 1e3, color="C0", ls=":", label="gap centre")
    a1.axvline(Z_INJECT * 1e3, color="C2", ls="--", label="injection")
    v_gap = scale * V1J_KEV                         # physical gap voltage [kV]
    a1.set_xlabel("lab z  [mm]")
    a1.set_ylabel(r"on-axis $E_z \times$ scale  [MV/m]")
    a1.set_title(f"{name}: cavity field (scale={scale:.3f}, "
                 f"$V_{{gap}}$≈{v_gap:.0f} kV)")
    a1.legend(fontsize=8)

    # ── Right: temporal RF waveform around t_gap, normalised to ±1 ──────────────
    # E ∝ cos(ω t+φ), B ∝ sin(ω t+φ); the on-axis Ez is single-signed, so the
    # *energy kick* of an electron ∝ -cos(ω t+φ). zc phases the bunch onto the
    # zero of that kick (max slope → head decelerated / tail accelerated); crest
    # phases it onto the kick maximum (pure acceleration).
    T = 1.0 / F_RF
    tt = np.linspace(t_gap - 1.0 * T, t_gap + 1.0 * T, 800)
    a2.plot((tt - t_gap) * 1e9, np.cos(OMEGA * tt + phi), color="C3",
            label=r"$E \propto \cos(\omega t+\varphi)$")
    a2.plot((tt - t_gap) * 1e9, np.sin(OMEGA * tt + phi), color="C0",
            label=r"$B \propto \sin(\omega t+\varphi)$")
    a2.axhline(0, color="k", lw=0.6)
    a2.axvline(0.0, color="C1", lw=1.2, label="bunch centre @ gap")
    # Shade the bunch temporal width σ_t = σ_z / v_beam at the snapshot nearest
    # the gap, so the reader sees how much RF phase the whole bunch samples.
    igap = int(np.argmin(np.abs(rec["zmean"] - Z_GAP_CENTER)))
    sigma_t = rec["sigz"][igap] / v_beam
    a2.axvspan(-sigma_t * 1e9, sigma_t * 1e9, color="C1", alpha=0.15,
               label=r"$\pm\sigma_t$ bunch")
    # The E field the bunch centre sits on tells zc (≈0) from crest (≈±1 extremum).
    e_at_gap = np.cos(OMEGA * t_gap + phi)
    landing = "ZERO-CROSSING" if phase == "zc" else "CREST"
    a2.annotate(f"bunch centre on field {landing}\n"
                f"(E≈{e_at_gap:+.2f} of peak)",
                xy=(0.0, e_at_gap), xytext=(0.05, 0.06),
                textcoords="axes fraction", fontsize=8,
                arrowprops=dict(arrowstyle="->", color="C1", lw=1.0))
    a2.set_xlabel(r"$t - t_{gap}$  [ns]")
    a2.set_ylabel("normalised field")
    a2.set_title(f"{name}: RF waveform at the gap (φ={phi:.2f} rad)")
    a2.legend(fontsize=8, loc="upper right")

    path = f"{RESULTS}/prebuncher_cavity.png"
    fig.savefig(path, dpi=140); plt.close(fig)
    print(f"wrote {path}")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE B — prebuncher_bunch_profile.png : the REAL longitudinal bunch shape λ(z).
#
# The σ_z(z) curve is a single scalar per snapshot; it cannot show *how* the bunch
# is shaped — the compression, and the space-charge spike/filamentation that grows
# near the focus. This histograms the actual line-charge density λ(z−⟨z⟩) =
# (Σ w q_e in bin)/Δz at the SAME three snapshots the phase-space figure uses
# (injection / cavity exit / max bunching), so the bunch's true profile and its
# peak λ / σ_z are visible. A drift baseline, if present, is overlaid (guarded).
# ══════════════════════════════════════════════════════════════════════════════
def bunch_profile_figure(name, rec, base):
    snaps = rec["snaps"]
    picks, titles = snapshot_picks(rec, base)
    base_snaps = base["snaps"] if base is not None else None
    base_its = base["it"] if base is not None else None

    fig, axs = plt.subplots(1, 3, figsize=(13, 4.0), constrained_layout=True)
    for ax, it, ti in zip(axs, picks, titles):
        z, ke, w = snaps[it]
        zc = z - np.average(z, weights=w)          # centre on ⟨z⟩
        sz = np.sqrt(np.average(zc**2, weights=w))
        # λ(z) [nC/m] = (charge per bin)/Δz; bins span a few σ each side so the
        # tails/spike are resolved without over-binning the sparse wings.
        span = max(4.0 * sz, 5e-4)
        edges = np.linspace(-span, span, 121)
        dzbin = edges[1] - edges[0]
        q_bin, _ = np.histogram(zc, bins=edges, weights=w * Q_E)
        lam = q_bin / dzbin * 1e9                   # C/m -> nC/m
        ctr = 0.5 * (edges[:-1] + edges[1:])
        ax.plot(ctr * 1e3, lam, color="C0", lw=1.4, label="cavity")
        ax.fill_between(ctr * 1e3, lam, color="C0", alpha=0.20)

        # Overlay the drift baseline profile at the matched ⟨z⟩ snapshot (if any).
        if base_snaps is not None:
            jb = int(np.argmin(np.abs(np.asarray(base["zmean"])
                                      - rec["zmean"][rec["it"].index(it)])))
            zb, _, wb = base_snaps[base_its[jb]]
            zbc = zb - np.average(zb, weights=wb)
            qb, _ = np.histogram(zbc, bins=edges, weights=wb * Q_E)
            ax.plot(ctr * 1e3, qb / dzbin * 1e9, "k--", lw=1.0,
                    label="drift baseline")

        zi = rec["it"].index(it)
        ax.set_title(f"{ti}  (⟨z⟩={rec['zmean'][zi]*1e3:.0f} mm)")
        ax.set_xlabel("z − ⟨z⟩  [mm]")
        ax.set_ylabel("λ  [nC/m]")
        ax.annotate(f"peak λ = {lam.max():.2f} nC/m\nσ_z = {sz*1e3:.2f} mm",
                    xy=(0.97, 0.95), xycoords="axes fraction",
                    ha="right", va="top", fontsize=8,
                    bbox=dict(boxstyle="round", fc="white", alpha=0.8))
        if base_snaps is not None:
            ax.legend(fontsize=8, loc="upper left")
    fig.suptitle(f"{name}: longitudinal line-charge density λ(z)", fontsize=12)
    path = f"{RESULTS}/prebuncher_bunch_profile.png"
    fig.savefig(path, dpi=140); plt.close(fig)
    print(f"wrote {path}")


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
    fig.savefig(f"{RESULTS}/prebuncher_line.png", dpi=140); plt.close(fig)

    # z–KE phase space at injection / gap exit / best-bunching point.
    snaps, its = rec["snaps"], rec["it"]
    picks, titles = snapshot_picks(rec, base)
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.0), constrained_layout=True)
    for ax, it, ti in zip(axs, picks, titles):
        z, ke, w = snaps[it]
        ax.scatter((z - z.mean()) * 1e3, ke - ke.mean(), s=2, alpha=0.15, color="C0")
        zi = rec["it"].index(it)
        ax.set_xlabel("z − ⟨z⟩  [mm]"); ax.set_ylabel("KE − ⟨KE⟩  [keV]")
        ax.set_title(f"{ti}  (⟨z⟩={rec['zmean'][zi]*1e3:.0f} mm)")
    fig.suptitle(f"{name}: longitudinal phase space", fontsize=12)
    fig.savefig(f"{RESULTS}/prebuncher_phasespace.png", dpi=140); plt.close(fig)

    out = dict(szmin=float(rec["sigz"].min()),
               sz0=float(rec["sigz"][0]),
               ke_end=float(rec["ke"][-1]), ipk_max=float(rec["ipk"].max()))
    if base is not None:
        sd = np.interp(rec["zmean"], base["zmean"], base["sigz"])
        ratio = sd / rec["sigz"]; ib = int(np.argmax(ratio))
        out.update(ratio=float(ratio[ib]), zbest=float(rec["zmean"][ib]),
                   sz_best=float(rec["sigz"][ib]))
    return out


def main(cases=None):
    """Plot prebuncher case(s).

    With `cases=None`, plots every `{DIAG_ROOT}/P*` directory; pass a list of
    case dir names (e.g. ["P800_zc"]) to restrict.
    """
    if cases:
        dirs = [os.path.join(DIAG_ROOT, c) for c in cases]
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

    n_powered = sum(1 for _, lab in cases if lab[1] != "drift")
    if n_powered > 1:
        print(f"note: {n_powered} powered cases present — the prebuncher_*.png "
              f"per-case figures use fixed filenames and are OVERWRITTEN (last "
              f"case wins); see compare_power_phase.png for the cross-case scan.",
              flush=True)

    summary = []
    for d, (power, phase) in sorted(cases, key=lambda x: (x[1][1], x[1][0])):
        if phase == "drift":
            continue
        name = os.path.basename(d)
        print(f"analysing {name} …", flush=True)
        rec = analyse_case(d)
        s = per_case_figure(name, rec, base)
        cavity_figure(name, rec, power, phase)        # FIGURE A: the RF drive
        bunch_profile_figure(name, rec, base)         # FIGURE B: real λ(z) shape
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
