"""
Figures and summary for the WarpX injector stage.

Default run writes one case to `injector/diags/main`; an optional power/phase scan
writes `injector/diags/P<power>_<phase>` dirs (and a `P0_drift` baseline). This
plotter reads `diags/main` plus any `diags/P*` cases present.

The gun-exit bunch is short (σ_z ≈ a few mm vs the 214 MHz RF wavelength) and
carries an intrinsic +1.4 keV/mm (debunching) energy chirp, and at ~1 nC it is
space-charge dense. In a free drift it therefore *expands*. The prebunchers act as
ballistic bunchers: at the zero-crossing they flip the chirp negative and compress
the bunch — but space charge limits how far. The clearest way to show this is to
compare each powered run against a **drift-only baseline** (P = 0):

  * σ_z(z)        — bunch length along the line, cavity runs vs. the drift baseline
  * ratio(z)      — σ_z,drift(z) / σ_z,cavity(z)  (>1 ⇒ cavity is bunching)
  * I_peak(z)     — peak current (peak line density × beam velocity)
  * z–KE space    — the chirp flipping/rotating through the cavities

Per-case figures use CONFIG-INDEPENDENT filenames (the operating point lives in the
diags/<case> input dir and the figure titles, not the filename) so a new run
overwrites the same files instead of leaving orphans:

  * injector_line.png          — σ_z(z) (vs. drift) and peak current / mean energy
  * injector_phasespace.png    — z–KE at injection / cavity exit / best focus
  * injector_cavity.png        — the RF DRIVE: on-axis Ez(z) of the scaled 1-J map(s)
                                 in the lab frame, plus the cos/sin RF waveform vs time
  * injector_bunch_profile.png — the real longitudinal line-charge density λ(z)
  * compare_power_phase.png    — σ_z(z) for the baseline vs all powers (scan only)

Run with:
    conda run -n CBB python injector/plot_injector.py
"""

import os
import glob
import re
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import openpmd_api as io
from openpmd_viewer import OpenPMDTimeSeries


def _retry_io(fn, *args, tries=6, base=0.25, **kwargs):
    """Call an openPMD read, retrying a transient HDF5 "Inaccessible" open error.

    NOTE: the production "OPEN_FILE failed ... Inaccessible" failure on this
    stage was fd exhaustion (openpmd-viewer leaks an fd per get_particle vs
    macOS's 256-fd default), now fixed by raising RLIMIT_NOFILE — see
    _runner._raise_fd_limit. This retry does NOT help that case (the fds stay
    spent); it is a backstop only for a genuinely transient open, e.g. the
    in-sim handoff report opening a Series while WarpX's own diagnostic Series is
    still releasing a just-flushed file. Re-raises after the last try so a
    genuinely missing file (or unfixed fd exhaustion) still surfaces.
    """
    for i in range(tries):
        try:
            return fn(*args, **kwargs)
        except io.Error:
            if i == tries - 1:
                raise
            time.sleep(base * 2 ** i)

# RF drive constants (F_RF, Q_L, V1J_KEV, gap centres) come from the single source of
# truth in the build module so this figure's re-derived scale/phase cannot drift from
# the actual run.
from .build_injector_field import V1J_KEV, F_RF, Q_L_1, Z_GAP_CENTER_1, Z_GAP_CENTER_2

c = 299792458.0
MC2 = 0.51099895e3                # electron rest energy [keV]
Q_E = 1.602176634e-19            # elementary charge [C]
DIAG_ROOT = "injector/diags"
RESULTS = "injector/results"
PREB1_FIELD = "injector/injector_field/preb1_EB.h5"
# Cavity-1 gap lab-z (for marking plots) is imported from build_injector_field above.
Z_GAP_CENTER = Z_GAP_CENTER_1

# ── RF-drive constants ────────────────────────────────────────────────────────
# The cavity map is stored 1-J-normalised; the run multiplies it by `scale` and
# modulates E ∝ cos(ω t + φ), B ∝ sin(ω t + φ). To draw the field the beam actually
# sees we re-derive scale and φ here from the same formulae the sim uses.
OMEGA = 2.0 * np.pi * F_RF       # RF angular frequency [rad/s]
Z_INJECT = 0.005                 # [m] lab z where the bunch tail (smallest z) is launched
os.makedirs(RESULTS, exist_ok=True)


def rf_scale(power, q_l=Q_L_1):
    """Field scale = sqrt(stored_energy / 1 J), stored_energy = 1e3·Q·P/(2π f_RF).

    Same expression as injector_sim.py: at P=0 (drift baseline) there is no cavity
    field, so the scale is 0.
    """
    if power <= 0:
        return 0.0
    return float(np.sqrt(1e3 * q_l * power / (2.0 * np.pi * F_RF)))


def rf_phase(phase, t_gap):
    """RF phase φ that puts the bunch-tail gap arrival at zero-crossing/crest.

    Mirrors injector_sim.py make_cavity() base term (phi_off omitted here — the
    waveform figure shows the zc/crest landing, not the GUI on-crest reference).
    """
    if phase == "crest":
        return -OMEGA * t_gap + np.pi
    return -OMEGA * t_gap + np.pi / 2.0


def load_cavity_axis():
    """On-axis Ez(z) of the raw 1-J Prebuncher-1 map, in LAB z."""
    s = io.Series(PREB1_FIELD, io.Access.read_only)
    E = s.iterations[0].meshes["E"]
    ez = E["z"].load_chunk()
    s.flush()
    ez = ez[0]                                   # (nr, nz), mode 0
    dz = E.grid_spacing[1]
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
    lam, _ = np.histogram(z, bins=edges, weights=w * Q_E)
    return float(lam.max() / dz * v_beam)


def analyse_case(path):
    """Per-snapshot metric arrays for one case directory (sorted by ⟨z⟩)."""
    ts = OpenPMDTimeSeries(os.path.join(path, "particles"))
    rec = dict(zmean=[], sigz=[], ke=[], dke=[], ipk=[], it=[])
    snaps = {}
    v_beam = None
    for i, it in enumerate(ts.iterations):
        z, ux, uy, uz, w = _retry_io(
            ts.get_particle,
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
    if not rec["zmean"]:                       # every dump empty (crashed/aborted run)
        return None
    order = np.argsort(rec["zmean"])
    for k in ("zmean", "sigz", "ke", "dke", "ipk"):
        rec[k] = np.asarray(rec[k])[order]
    rec["it"] = [rec["it"][i] for i in order]
    rec["v_beam"] = v_beam or 0.632 * c
    rec["snaps"] = snaps
    return rec


def case_label(name):
    """(power, phase) for a case dir name. The default chain run is 'main'."""
    if name == "main":
        return ("main", "main")
    if name == "P0_drift":
        return (0, "drift")
    m = re.match(r"P(\d+(?:\.\d+)?)_(zc|crest)$", name)   # accept fractional watts
    return (float(m.group(1)), m.group(2)) if m else None


def snapshot_picks(rec, base):
    """The three snapshots used by the phase-space and bunch-profile figures:
    injection / cavity exit / best bunching. Factored out so both figures sample
    the SAME iterations."""
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


def load_cavity_axis_file(field_path):
    """On-axis Ez(z) of a raw 1-J cavity map at ``field_path``, in LAB z."""
    s = io.Series(field_path, io.Access.read_only)
    E = s.iterations[0].meshes["E"]
    ez = E["z"].load_chunk()
    s.flush()
    ez = ez[0]
    dz = E.grid_spacing[1]
    z0 = E.grid_global_offset[1] if E.grid_global_offset else 0.0
    z_lab = z0 + np.arange(ez.shape[1]) * dz
    return z_lab, ez[0]


def cavity_phi(z_gap, v_at_gap, phi_off_deg, phase, rev_phase, t_offset=0.0):
    """Re-derive a cavity's drive phase φ exactly as injector_sim.make_cavity does
    (crest base for the GUI phi_off + the reversal phase), so the drawn waveform
    matches the run (single source of truth)."""
    t_gap = t_offset + (z_gap - Z_INJECT) / v_at_gap
    base = np.pi / 2.0 if phase == "zc" else np.pi
    return -OMEGA * t_gap + base + np.radians(phi_off_deg) + rev_phase, t_gap


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE A — injector_cavity.png : the RF DRIVE both prebunchers apply.
#
# Left  — on-axis Ez(z) of each scaled 1-J map placed at its lab gap (Preb 1 @ Z1,
#         Preb 2 @ Z2), so the two bunching lobes are visible along the line.
# Right — the temporal RF waveform E ∝ cos(ωt+φ) at each cavity's arrival, with the
#         bunch-centre marker, showing where the bunch lands on each cavity's phase.
# Preb-2's scale/phase are re-derived the SAME way the sim does (crest base + GUI
# phi_off + PREB2_REV_PHASE) so the figure can't drift from the run.
# ══════════════════════════════════════════════════════════════════════════════
def cavity_figure(name, rec, power, phase):
    from .injector_sim import (
        PREB1_KW, PREB1_Q, PREB1_PHI_OFF, PREB2_KW, PREB2_Q, PREB2_PHI_OFF,
        PREB2_REVERSED, PREB2_REV_PHASE, PREB1_FIELD, PREB2_FIELD)
    v_beam = rec["v_beam"]

    # Preb 1
    scale1 = rf_scale(PREB1_KW, PREB1_Q)
    phi1, t_gap1 = cavity_phi(Z_GAP_CENTER_1, v_beam, PREB1_PHI_OFF, phase, 0.0)
    z1_lab, ez1 = load_cavity_axis_file(PREB1_FIELD)
    # Preb 2 (reversed via PREB2_REV_PHASE; constant-v arrival uses v_beam)
    have2 = PREB2_KW > 0
    if have2:
        scale2 = rf_scale(PREB2_KW, PREB2_Q)
        rev = PREB2_REV_PHASE if PREB2_REVERSED else 0.0
        phi2, t_gap2 = cavity_phi(Z_GAP_CENTER_2, v_beam, PREB2_PHI_OFF, phase, rev)
        z2_lab, ez2 = load_cavity_axis_file(PREB2_FIELD)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.4), constrained_layout=True)

    # ── Left: spatial Ez(z) lobes of both cavities in the lab frame ─────────────
    a1.plot(z1_lab * 1e3, ez1 * scale1 / 1e6, color="C3",
            label=f"Preb 1 ({PREB1_KW:g} kW, $V_g$≈{scale1*V1J_KEV:.0f} kV)")
    if have2:
        a1.plot(z2_lab * 1e3, ez2 * scale2 / 1e6, color="C4",
                label=f"Preb 2 rev ({PREB2_KW:g} kW, $V_g$≈{scale2*V1J_KEV:.0f} kV)")
    a1.axhline(0, color="k", lw=0.6)
    a1.axvline(Z_GAP_CENTER_1 * 1e3, color="C3", ls=":", lw=0.8)
    if have2:
        a1.axvline(Z_GAP_CENTER_2 * 1e3, color="C4", ls=":", lw=0.8)
    a1.axvline(Z_INJECT * 1e3, color="C2", ls="--", label="injection")
    a1.set_xlabel("lab z  [mm]")
    a1.set_ylabel(r"on-axis $E_z \times$ scale  [MV/m]")
    a1.set_title(f"{name}: cavity fields (Preb 1 @ {Z_GAP_CENTER_1*1e3:.0f}, "
                 f"Preb 2 @ {Z_GAP_CENTER_2*1e3:.0f} mm)")
    a1.legend(fontsize=8)

    # ── Right: temporal RF waveform at each cavity's bunch arrival ──────────────
    T = 1.0 / F_RF
    tt = np.linspace(-1.0 * T, 1.0 * T, 800)
    a2.plot(tt * 1e9, np.cos(OMEGA * (t_gap1 + tt) + phi1), color="C3",
            label=r"Preb 1 $E\propto\cos$")
    if have2:
        a2.plot(tt * 1e9, np.cos(OMEGA * (t_gap2 + tt) + phi2), color="C4",
                label=r"Preb 2 $E\propto\cos$")
    a2.axhline(0, color="k", lw=0.6)
    a2.axvline(0.0, color="C1", lw=1.2, label="bunch centre @ gap")
    e1 = np.cos(OMEGA * t_gap1 + phi1)
    txt = f"Preb 1: E≈{e1:+.2f} of peak"
    if have2:
        e2 = np.cos(OMEGA * t_gap2 + phi2)
        txt += f"\nPreb 2: E≈{e2:+.2f} of peak"
    a2.annotate(txt, xy=(0.04, 0.06), xycoords="axes fraction", fontsize=8,
                bbox=dict(boxstyle="round", fc="white", alpha=0.8))
    a2.set_xlabel(r"$t - t_{gap}$  [ns]")
    a2.set_ylabel("normalised E field at gap arrival")
    a2.set_title(f"{name}: RF waveforms (φ1={phi1:.2f}"
                 + (f", φ2={phi2:.2f}" if have2 else "") + " rad)")
    a2.legend(fontsize=8, loc="upper right")

    path = f"{RESULTS}/injector_cavity.png"
    fig.savefig(path, dpi=140); plt.close(fig)
    print(f"wrote {path}")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE B — injector_bunch_profile.png : the REAL longitudinal bunch shape λ(z).
# ══════════════════════════════════════════════════════════════════════════════
def bunch_profile_figure(name, rec, base):
    snaps = rec["snaps"]
    picks, titles = snapshot_picks(rec, base)
    base_snaps = base["snaps"] if base is not None else None
    base_its = base["it"] if base is not None else None

    fig, axs = plt.subplots(1, 3, figsize=(13, 4.0), constrained_layout=True)
    for ax, it, ti in zip(axs, picks, titles):
        z, ke, w = snaps[it]
        zc = z - np.average(z, weights=w)
        sz = np.sqrt(np.average(zc**2, weights=w))
        span = max(4.0 * sz, 5e-4)
        edges = np.linspace(-span, span, 121)
        dzbin = edges[1] - edges[0]
        q_bin, _ = np.histogram(zc, bins=edges, weights=w * Q_E)
        lam = q_bin / dzbin * 1e9                   # C/m -> nC/m
        ctr = 0.5 * (edges[:-1] + edges[1:])
        ax.plot(ctr * 1e3, lam, color="C0", lw=1.4, label="cavity")
        ax.fill_between(ctr * 1e3, lam, color="C0", alpha=0.20)

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
    path = f"{RESULTS}/injector_bunch_profile.png"
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
    a1.axvline(Z_GAP_CENTER_1 * 1e3, color="C3", ls=":", label="Preb 1 gap")
    a1.axvline(Z_GAP_CENTER_2 * 1e3, color="C5", ls=":", label="Preb 2 gap")
    a1.set_xlabel("⟨z⟩  [mm]"); a1.set_ylabel("σ_z  [mm]")
    a1.set_title(f"{name}: bunch length vs. drift"); a1.legend(fontsize=8)

    a2b = a2.twinx()
    a2.plot(zmm, rec["ipk"], "o-", ms=3, color="C2")
    a2b.plot(zmm, rec["ke"], "s--", ms=3, color="C4")
    a2.axvline(Z_GAP_CENTER_1 * 1e3, color="C3", ls=":")
    a2.axvline(Z_GAP_CENTER_2 * 1e3, color="C5", ls=":")
    a2.set_xlabel("⟨z⟩  [mm]")
    a2.set_ylabel("peak current  [A]", color="C2")
    a2b.set_ylabel("mean KE  [keV]", color="C4")
    a2.set_title(f"{name}: peak current & mean energy")
    fig.savefig(f"{RESULTS}/injector_line.png", dpi=140); plt.close(fig)

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
    fig.savefig(f"{RESULTS}/injector_phasespace.png", dpi=140); plt.close(fig)

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
    """Plot injector case(s).

    With `cases=None`, plots `{DIAG_ROOT}/main` and every `{DIAG_ROOT}/P*` directory;
    pass a list of case dir names (e.g. ["main"]) to restrict.
    """
    if cases:
        dirs = [os.path.join(DIAG_ROOT, c) for c in cases]
    else:
        dirs = []
        if os.path.isdir(os.path.join(DIAG_ROOT, "main")):
            dirs.append(os.path.join(DIAG_ROOT, "main"))
        dirs += sorted(glob.glob(f"{DIAG_ROOT}/P*"))
    cases = [(d, case_label(os.path.basename(d))) for d in dirs if os.path.isdir(d)]
    for d, lab in cases:
        if lab is None:
            print(f"warning: skipping {os.path.basename(d)!r} — not a recognized case "
                  f"dir (expected 'main', P<power>_<zc|crest> or P0_drift)", flush=True)
    cases = [(d, lab) for d, lab in cases if lab]
    if not cases:
        print(f"No case directories found under {DIAG_ROOT}/ "
              f"(run injector_sim.py first).")
        return

    base = None
    for d, lab in cases:
        if lab[1] == "drift":
            print("analysing drift baseline …", flush=True)
            base = analyse_case(d)
            if base is None:
                print(f"  skipping {os.path.basename(d)}: no usable snapshots "
                      f"(empty/aborted run)", flush=True)

    powered = [(d, lab) for d, lab in cases if lab[1] != "drift"]
    if len(powered) > 1:
        print(f"note: {len(powered)} non-baseline cases present — the injector_*.png "
              f"per-case figures use fixed filenames and are OVERWRITTEN (last case "
              f"wins); see compare_power_phase.png for the cross-case scan.", flush=True)

    summary = []
    for d, (power, phase) in sorted(powered, key=lambda x: (str(x[1][1]), str(x[1][0]))):
        name = os.path.basename(d)
        print(f"analysing {name} …", flush=True)
        rec = analyse_case(d)
        if rec is None:
            print(f"  skipping {name}: no usable snapshots (empty/aborted run)",
                  flush=True)
            continue
        s = per_case_figure(name, rec, base)
        # For the default 'main' case the power/phase are read from the sim defaults.
        if phase == "main":
            from .injector_sim import PREB1_KW, PHASE
            cav_power, cav_phase = PREB1_KW, PHASE
        else:
            cav_power, cav_phase = power, phase
        cavity_figure(name, rec, cav_power, cav_phase)    # FIGURE A: the RF drive
        bunch_profile_figure(name, rec, base)             # FIGURE B: real λ(z) shape
        s.update(name=name, power=power, phase=phase, rec=rec)
        summary.append(s)

    # ── Headline comparison (only meaningful with several cases / a baseline) ──
    scan = [s for s in summary if s["phase"] in ("zc", "crest")]
    if (len(scan) > 1 or base is not None) and scan:
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
        if base is not None:
            a1.plot(base["zmean"] * 1e3, base["sigz"] * 1e3, "k--", lw=2,
                    label="drift (P=0)")
        for s in sorted([s for s in scan if s["phase"] == "zc"],
                        key=lambda s: s["power"]):
            r = s["rec"]
            a1.plot(r["zmean"] * 1e3, r["sigz"] * 1e3, "-", label=f"{s['power']:g} kW")
        a1.axvline(Z_GAP_CENTER * 1e3, color="C3", ls=":")
        a1.set_xlabel("⟨z⟩  [mm]"); a1.set_ylabel("σ_z  [mm]")
        a1.set_title("Bunch length: drift vs. zero-crossing cavity"); a1.legend(fontsize=8)

        for phase, col in (("zc", "C0"), ("crest", "C3")):
            pts = sorted([s for s in scan if s["phase"] == phase and "ratio" in s],
                         key=lambda s: s["power"])
            if pts:
                a2.plot([s["power"] for s in pts], [s["ratio"] for s in pts], "o-",
                        color=col, label="zero-crossing" if phase == "zc" else "on-crest")
        a2.axhline(1.0, color="k", lw=0.6)
        a2.set_xlabel("RF power P  [kW]")
        a2.set_ylabel("max bunching  σ_drift / σ_cavity")
        a2.set_title("Bunching vs. power"); a2.legend(fontsize=8)
        fig.savefig(f"{RESULTS}/compare_power_phase.png", dpi=140); plt.close(fig)
        print(f"wrote {RESULTS}/compare_power_phase.png")

    # ── Summary table ─────────────────────────────────────────────────────────
    print("\n" + "=" * 96)
    print(f"{'case':>10} {'P[kW]':>6} {'phase':>6} {'σz0[mm]':>8} {'σzmin[mm]':>10} "
          f"{'bunch×':>7} {'zbest[mm]':>9} {'Ipk[A]':>7} {'KEend[keV]':>11}")
    print("-" * 96)
    if base is not None:
        print(f"{'P0_drift':>10} {0:6d} {'drift':>6} {base['sigz'][0]*1e3:8.3f} "
              f"{base['sigz'].min()*1e3:10.3f} {'—':>7} {'—':>9} "
              f"{base['ipk'].max():7.2f} {base['ke'][-1]:11.1f}")
    for s in sorted(summary, key=lambda s: (str(s["phase"]), str(s["power"]))):
        ratio = f"{s['ratio']:7.2f}" if "ratio" in s else f"{'—':>7}"
        zbest = f"{s['zbest']*1e3:9.0f}" if "zbest" in s else f"{'—':>9}"
        pstr = f"{s['power']:>6g}" if not isinstance(s["power"], str) else f"{s['power']:>6}"
        print(f"{s['name']:>10} {pstr} {str(s['phase']):>6} {s['sz0']*1e3:8.3f} "
              f"{s['szmin']*1e3:10.3f} {ratio} {zbest} {s['ipk_max']:7.2f} "
              f"{s['ke_end']:11.1f}")
    print("=" * 96)


if __name__ == "__main__":
    main()
