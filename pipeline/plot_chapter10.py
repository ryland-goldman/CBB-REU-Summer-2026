"""
Chapter-10 (PARMELA) figure reproduction for the Cornell injector chain.

Reproduces the qualitative KE-vs-RF-phase views of Cheng-Yang Tan's dissertation
Chapter 10 (Tan Fig 10.2–10.5) from the EXISTING default pipeline diagnostics — no
pipeline re-run. Four 2-panel figures (KE-vs-phase scatter on top, peak-normalized
charge-vs-phase histogram on the bottom, sharing x) at four chain locations, plus a
comparison table against Tan's Table 10.2.

  * results/tan_fig10p2_at_gun_<case>.png        — Tan 10.2  (gun exit)
  * results/tan_fig10p3_before_preb2_<case>.png  — Tan 10.3  (injector, ~Z_PREB2)
  * results/tan_fig10p4_before_sec1_<case>.png   — Tan 10.4  (injector, ~Z_HANDOFF, pre-iris)
  * results/tan_fig10p5_after_sec1_<case>.png    — Tan 10.5  (linac_sec1 exit, captured core)
  * results/tan_comparison_<case>.md / .csv      — the comparison table

Run with (in the CBB env, from repo root):
    python pipeline/plot_chapter10.py

SCOPE (critical): the repo default operating point is NOT Tan condition (i). The default
is the LinacSim 8 kW / 10 kW two-cavity point (Preb1 −70°, Preb2 −45°), so these figures
are a QUALITATIVE comparison to Tan's shapes only — no claim of numeric agreement with
Table 10.2. Tan's published numbers appear as a labeled reference column with deltas. The
physics-correct Tan cond (i) is Deliverable B (an injector re-run) and is out of scope here.

PHYSICS / UNITS NOTES:
  - openPMD u-components (ux,uy,uz) are γβ already; γ = sqrt(1+ux²+uy²+uz²), KE = (γ−1)·MC2.
  - The RF-phase x-axis is the ARRIVAL TIME of each particle at the bunch's centroid plane,
    NOT a spatial z→φ map: t_i = −(z_i − z_ref)/v_z,i. A particle ahead (larger z) arrives
    EARLIER (t<0 ⇒ φ<0 = head). ONE v_z definition (`_vz`) feeds both φ and σ_t so a sign
    fix cannot desync them.
  - v_z guard: drop particles with v_z ≤ 1e-3·c (space-charge tails at the gun can have
    small/negative uz → |t|→∞); excluded from φ AND σ_t and counted in the caption.
  - σ_z@214 = 360·F_214·σ_t, σ_z@2856 = 360·F_2856·σ_t — both from one σ_t, so the repo
    @2856/@214 ratio is BY CONSTRUCTION F_2856/F_214 = 13.335 (an INTERNAL consistency
    check, NOT agreement with Tan's 13.52 which uses a different spatial-σ_z definition).
  - Histograms are PEAK-normalized (counts/counts.max()), NOT density=True: Tan's ~16 nC
    bunch vs the repo's ~0.8 nC means shapes compare, not absolute charge.
  - The linac_sec1 local-frame z-reset (linac_sec1_sim.py:178 `z = z - z.min() + Z_INJECT`)
    is offset-only and monotonic, so (z − z_ref) is frame-invariant — do NOT "fix" it.
"""

import os
import sys
import csv

# Put the repo root on sys.path so the plain `python pipeline/plot_chapter10.py` form
# works (not just `python -m pipeline.plot_chapter10`) — matches run_pipeline.py. Must
# precede the `pipeline`/`injector` imports below.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from openpmd_viewer import OpenPMDTimeSeries

from pipeline._runner import _raise_fd_limit
from injector.build_injector_field import (
    F_RF as F_214,                 # 214.1847857 MHz — import, don't re-derive
    Z_GAP_CENTER_2,                # 1.318 m — Prebuncher 2 gap center (mid-kick)
    MAP_HALF_Z,                    # 0.1524 m — half-length of the prebuncher map
    Z_HANDOFF,                     # 2.03 m — injector→linac handoff plane
)

# Tan's "before Preb2" plane is the cavity ENTRANCE (pre-kick), NOT the gap center: at the
# gap center the beam is mid-kick, inflating ⟨E⟩ (0.214 vs Tan's 0.139). The Preb2 map spans
# Z_GAP_CENTER_2 ± MAP_HALF_Z, so the entrance is one half-length upstream of the gap center.
Z_PREB2 = Z_GAP_CENTER_2 - MAP_HALF_Z   # ≈ 1.166 m (Preb2 entrance)

# F_2856 is the linac_sec1 S-band RF freq (= linac_sec1_sim.F_RF). Kept as a literal,
# NOT imported, because importing linac_sec1_sim loads pywarpx/PICMI which binds the
# global geometry at first .so load — unsafe at plot time (design memo FLAG D).
F_2856 = 2856.0e6                  # = linac_sec1_sim.F_RF (S-band)

MC2 = 0.51099895                   # electron rest energy [MeV]
C = 299792458.0
Q_E = 1.602176634e-19              # elementary charge [C]
RESULTS = "results"
SPECIES = "electrons"              # plural — the openPMD species key for all four reads
MIN_KE_MEV = 12.0                  # linac_sec1 captured-core / model-validity cut
MIN_LIVE = 50                      # < 50 live particles ⇒ skip dump (matches plot_chain.py)

# Tan Table 10.2 (condition i) published reference. Per location:
#   (Ebar_MeV, sigE_MeV, sigz_deg@214, sigz_deg@2856, cap_in_pct, cap_all_pct)
# Capture only at postSec1 (in-bucket 89.4 / all-buckets 96.8); upstream rows are 100%.
TAN_TABLE = {
    "at_gun":       dict(Ebar=0.150, sigE=0.000, sz214=31.7,  sz2856=428.8, cap_in=100.0, cap_all=100.0),
    "before_preb2": dict(Ebar=0.139, sigE=0.023, sz214=14.4,  sz2856=194.8, cap_in=100.0, cap_all=100.0),
    "before_sec1":  dict(Ebar=0.253, sigE=0.043, sz214=5.57,  sz2856=75.3,  cap_in=100.0, cap_all=100.0),
    "after_sec1":   dict(Ebar=27.2,  sigE=3.5,   sz214=None,  sz2856=11.1,  cap_in=89.4,  cap_all=96.8),
}


# ─────────────────────────────────────────────────────────────────────────────
# Snapshot reading + kinematics
# ─────────────────────────────────────────────────────────────────────────────
def _read_snapshot(series_path, iteration=None, target_z=None):
    """Return (z, ux, uy, uz, w, zbar) for one dump of an RZ stage.

    iteration='last' -> the last dump with >= MIN_LIVE live particles (walk reversed).
    target_z set     -> the dump whose charge-weighted <z> is nearest target_z.
    The ACTUAL zbar of the chosen dump is returned (for the caption).
    """
    ts = OpenPMDTimeSeries(series_path)
    its = list(ts.iterations)

    def _get(it):
        x, y, z, ux, uy, uz, w = ts.get_particle(
            ["x", "y", "z", "ux", "uy", "uz", "w"], species=SPECIES, iteration=it)
        return z, ux, uy, uz, w

    if target_z is not None:
        best = None
        for it in its:
            try:
                z, ux, uy, uz, w = _get(it)
            except Exception:
                continue
            if len(z) < MIN_LIVE:
                continue
            zbar = float(np.average(z, weights=w))
            d = abs(zbar - target_z)
            if best is None or d < best[0]:
                best = (d, it, z, ux, uy, uz, w, zbar)
        if best is None:
            raise RuntimeError(f"{series_path}: no dump with >= {MIN_LIVE} particles near z={target_z}")
        _, _, z, ux, uy, uz, w, zbar = best
        return z, ux, uy, uz, w, zbar

    # iteration == 'last': last dump with enough live particles
    for it in reversed(its):
        try:
            z, ux, uy, uz, w = _get(it)
        except Exception:
            continue
        if len(z) < MIN_LIVE:
            continue
        zbar = float(np.average(z, weights=w))
        return z, ux, uy, uz, w, zbar
    raise RuntimeError(f"{series_path}: no dump with >= {MIN_LIVE} particles")


def ke_mev(ux, uy, uz):
    gamma = np.sqrt(1.0 + ux**2 + uy**2 + uz**2)   # ux/uy/uz are γβ
    return (gamma - 1.0) * MC2


def _vz(ux, uy, uz):
    """Longitudinal velocity [m/s]. The ONE definition shared by φ and σ_t."""
    gamma = np.sqrt(1.0 + ux**2 + uy**2 + uz**2)
    return (uz / gamma) * C


def arrival_time(z, ux, uy, uz, w, z_ref=None):
    """t_i = -(z_i - z_ref)/v_z,i  (ahead -> earlier -> t<0). Returns (t, z_ref, good).

    GUARD: drop |v_z| below 1e-3·c (or non-positive) — small/negative uz in the gun's
    space-charge tail would give |t|->inf. Excluded particles carry t=nan and good=False;
    they are dropped from the phase axis AND σ_t and counted in the caption.

    z_ref defaults to the charge-weighted centroid. (z - z_ref) is frame-invariant under
    an offset-only monotonic z-reset, so the linac_sec1 local-frame reset is fine here.
    """
    if z_ref is None:
        z_ref = float(np.average(z, weights=w))
    vz = _vz(ux, uy, uz)
    good = vz > (1e-3 * C)
    t = np.full_like(z, np.nan, dtype=float)
    t[good] = -(z[good] - z_ref) / vz[good]
    return t, z_ref, good


def rf_phase_deg(t, f):
    """RF phase [deg] at frequency f. Head (t<0) -> phi<0."""
    return 360.0 * f * t


def sigma_t(t, w, good):
    """Charge-weighted RMS arrival-time spread [s] on the SAME (t, good) as the phase axis."""
    tt, ww = t[good], w[good]
    m = np.average(tt, weights=ww)
    return float(np.sqrt(np.average((tt - m) ** 2, weights=ww)))


def wstats(v, w):
    """Charge-weighted (mean, std) — the plot_chain.py idiom."""
    m = float(np.average(v, weights=w))
    return m, float(np.sqrt(np.average((v - m) ** 2, weights=w)))


# ─────────────────────────────────────────────────────────────────────────────
# Figure rendering
# ─────────────────────────────────────────────────────────────────────────────
def location_panel(phi, ke, w, title, caption, phi_label, out_path):
    """One 2-panel figure: KE-vs-phi scatter (top) + peak-normalized charge histogram (bottom)."""
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(7.0, 7.0), sharex=True, constrained_layout=True)

    ax_top.scatter(phi, ke, s=2, c="C3", alpha=0.35, edgecolors="none")
    ax_top.set_ylabel("Kinetic Energy (MeV)")
    ax_top.set_title(title, fontsize=10)
    ax_top.grid(alpha=0.25)

    counts, edges = np.histogram(phi, bins=120, weights=w)
    if counts.max() > 0:
        counts = counts / counts.max()
    centers = 0.5 * (edges[:-1] + edges[1:])
    ax_bot.step(centers, counts, where="mid", color="C0")
    ax_bot.fill_between(centers, counts, step="mid", color="C0", alpha=0.25)
    ax_bot.set_ylabel("normalized charge distribution")
    ax_bot.set_xlabel(f"RF phase φ (deg) {phi_label}")
    ax_bot.grid(alpha=0.25)

    fig.text(0.5, -0.02, caption, ha="center", va="top", fontsize=6.5, color="0.25", wrap=True)
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
def main(case="repo_default",
         gun_diag="gun/diags/particles",
         inj_diag="injector/diags/main/particles",
         sec1_diag="linac_sec1/diags/main/particles",
         sec1_summary="linac_sec1/diags/main/injection_summary.json"):
    _raise_fd_limit()
    os.makedirs(RESULTS, exist_ok=True)

    import json
    with open(sec1_summary) as fh:
        summ = json.load(fh)
    q_injected_C = summ["q_injected_C"]    # honest (pre-iris) denominator
    q_in_bore_C = summ["q_in_bore_C"]      # post-iris denominator
    q_in_domain_C = summ["q_in_domain_C"]

    ratio_expected = F_2856 / F_214        # 13.335

    # Operating-point label is case-aware: repo_default is a QUALITATIVE-only comparison
    # (NOT Tan cond i); cond_i IS the physics-correct Tan condition (i) and may be compared
    # numerically to Table 10.2. Add new cases here as they are run.
    if case == "cond_i":
        op_label = ("Tan condition (i): Preb1 50 kV @ phase-null (PREB1_PHI_OFF=-90, 5.83 kW), "
                    "Preb2 150 kV @ crest (PREB2_PHI_OFF=0, 36.61 kW) — achieved V_gap 50.0/150.0 kV; "
                    "numeric comparison to Tan Table 10.2 applies")
        title_tag = "Tan cond (i)"           # figure-title operating-point tag (case-aware)
    else:
        op_label = ("repo default operating point: 8 kW / 10 kW two-cavity "
                    "(Preb1 −70°, Preb2 −45°) — QUALITATIVE comparison only, NOT Tan cond (i)")
        title_tag = "repo default"
    sign_note = "sign: head = φ<0; φ=0 = charge-weighted centroid"
    iris_pct = 100.0 * q_in_bore_C / q_injected_C   # in-bore iris transmission for this case
    caveat = ("charge PEAK-normalized (Tan ~16 nC vs repo ~0.8 nC — shapes compare, not pC); "
              f"non-relativistic ES self-field (~γ² pessimism); real 9.547 mm iris scrape "
              f"(in-bore transmission {iris_pct:.0f}%)")
    # B/C bunching note: the repo default does NOT reproduce Tan's bunching (different
    # operating point), so a shape gap there is not a code error. cond_i IS Tan's point and
    # should bunch — the note flips to a numeric-comparison cue.
    if case == "cond_i":
        bunch_note = "NOTE: this is Tan cond (i) — expect Tan's bunching here; numeric comparison to Table 10.2 applies"
    else:
        bunch_note = ("NOTE: repo default does not bunch at this plane (operating-point divergence) "
                      "— a shape gap vs Tan is NOT a code error")

    table_rows = []   # (key, location, zbar, Ebar, sigE, sz214, sz2856, cap_in, cap_all)

    def _record(key, location, zbar, Ebar, sigE, sigt, cap_in, cap_all, sz214_blank=False):
        sz214 = None if sz214_blank else (360.0 * F_214 * sigt if sigt is not None else None)
        sz2856 = 360.0 * F_2856 * sigt if sigt is not None else None
        # internal consistency: both from one sigma_t
        if sigt is not None and not sz214_blank:
            assert abs(sz2856 / sz214 - ratio_expected) < 1e-6, "sigma_z internal ratio != 13.335"
        table_rows.append(dict(key=key, location=location, zbar=zbar, Ebar=Ebar, sigE=sigE,
                               sz214=sz214, sz2856=sz2856, cap_in=cap_in, cap_all=cap_all))
        return sz214, sz2856

    # ── A: at Gun (Tan 10.2) ──
    z, ux, uy, uz, w, zbar = _read_snapshot(gun_diag, iteration="last")
    ke = ke_mev(ux, uy, uz)
    t, z_ref, good = arrival_time(z, ux, uy, uz, w)
    phi = rf_phase_deg(t, F_214)
    Ebar, sigE = wstats(ke, w)
    sigt = sigma_t(t, w, good)
    ndrop = int((~good).sum())
    # SIGN ASSERT — gates all figures. Head (φ<0) must be at higher KE at the gun.
    cc = np.corrcoef(phi[good], ke[good])[0, 1]
    assert cc < 0, (f"gun corr(phi,KE)={cc:+.3f} >= 0 — sign is wrong; fix once in arrival_time, "
                    "not per-panel")
    print(f"[A gun]   zbar={zbar:.4f} m  <KE>={Ebar:.4f} MeV  sigE={sigE:.4f}  "
          f"corr(phi,KE)={cc:+.3f}  vz-dropped={ndrop}/{len(z)}")
    sz214_gun, _ = _record("at_gun", "A at Gun", zbar, Ebar, sigE, sigt, None, None)
    # Headline divergence: the repo gun bunch is ~270× shorter than Tan's 3.7 ns base, so
    # this is a narrow falling band, not Tan's full sine lobe (computed from sigt, not hard-coded).
    gun_len_note = (f"repo gun σ_t ≈ {sigt*1e12:.1f} ps (σ_z@214 = {sz214_gun:.2f}°) vs Tan's 3.7 ns base "
                    f"— bunch ~{3.7e-9/sigt:.0f}× shorter, so this is a narrow falling band, not Tan's full sine lobe")
    cap = (f"{op_label}\n{sign_note}; ref 214.18 MHz; <z>={zbar:.4f} m; "
           f"v_z-dropped={ndrop}/{len(z)}; {caveat}\n{gun_len_note}\n"
           f"Tan Fig 10.2 (expect sine KE(φ), high-KE on φ<0)")
    location_panel(phi[good], ke[good], w[good], f"Tan Fig 10.2 — at Gun ({title_tag})",
                   cap, "@214 MHz", f"{RESULTS}/tan_fig10p2_at_gun_{case}.png")

    # ── B: before Preb2 (Tan 10.3) ──
    z, ux, uy, uz, w, zbar = _read_snapshot(inj_diag, target_z=Z_PREB2)
    ke = ke_mev(ux, uy, uz)
    t, z_ref, good = arrival_time(z, ux, uy, uz, w)
    phi = rf_phase_deg(t, F_214)
    Ebar, sigE = wstats(ke, w)
    sigt = sigma_t(t, w, good)
    ndrop = int((~good).sum())
    print(f"[B preb2] zbar={zbar:.4f} m  <KE>={Ebar:.4f} MeV  sigE={sigE:.4f}  "
          f"vz-dropped={ndrop}/{len(z)}")
    _record("before_preb2", "B before Preb2", zbar, Ebar, sigE, sigt, None, None)
    cap = (f"{op_label}\n{sign_note}; ref 214.18 MHz; <z>={zbar:.4f} m "
           f"(target Preb2 ENTRANCE {Z_PREB2:.4f} m = gap {Z_GAP_CENTER_2}−half {MAP_HALF_Z}, pre-kick); "
           f"v_z-dropped={ndrop}/{len(z)}; {caveat}\n"
           f"{bunch_note}\n"
           f"Tan Fig 10.3 (expect monotonic S-curve; peakier histogram)")
    location_panel(phi[good], ke[good], w[good], f"Tan Fig 10.3 — before Preb2 ({title_tag})",
                   cap, "@214 MHz", f"{RESULTS}/tan_fig10p3_before_preb2_{case}.png")

    # ── C: before Sec1 (Tan 10.4) — PRE-iris full injector population ──
    z, ux, uy, uz, w, zbar = _read_snapshot(inj_diag, target_z=Z_HANDOFF)
    ke = ke_mev(ux, uy, uz)
    t, z_ref, good = arrival_time(z, ux, uy, uz, w)
    phi = rf_phase_deg(t, F_214)
    Ebar, sigE = wstats(ke, w)
    sigt = sigma_t(t, w, good)
    ndrop = int((~good).sum())
    print(f"[C sec1]  zbar={zbar:.4f} m  <KE>={Ebar:.4f} MeV  sigE={sigE:.4f}  "
          f"vz-dropped={ndrop}/{len(z)}")
    _record("before_sec1", "C before Sec1", zbar, Ebar, sigE, sigt, "n/a (pre-iris)", "n/a (pre-iris)")
    cap = (f"{op_label}\n{sign_note}; ref 214.18 MHz; <z>={zbar:.4f} m (target Z_HANDOFF={Z_HANDOFF}); "
           f"PRE-iris full population; v_z-dropped={ndrop}/{len(z)}; {caveat}\n"
           f"{bunch_note}\n"
           f"Tan Fig 10.4 (expect left-opening parabola/'C'; sharp spike near φ≈0 + tail)")
    location_panel(phi[good], ke[good], w[good], f"Tan Fig 10.4 — before Sec1 ({title_tag}, pre-iris)",
                   cap, "@214 MHz", f"{RESULTS}/tan_fig10p4_before_sec1_{case}.png")

    # ── D: after Sec1 (Tan 10.5) — captured core (KE >= 12 MeV), F_2856 ──
    z, ux, uy, uz, w, zbar = _read_snapshot(sec1_diag, iteration="last")
    ke = ke_mev(ux, uy, uz)
    core = ke >= MIN_KE_MEV
    zc, uxc, uyc, uzc, wc, kec = z[core], ux[core], uy[core], uz[core], w[core], ke[core]
    # phi=0 references the captured-core centroid
    t, z_ref, good = arrival_time(zc, uxc, uyc, uzc, wc)
    phi = rf_phase_deg(t, F_2856)
    Ebar, sigE = wstats(kec, wc)
    sigt = sigma_t(t, wc, good)
    ndrop = int((~good).sum())
    # capture bookkeeping
    q_exit = float(wc.sum()) * Q_E
    all_buckets_pct = 100.0 * q_exit / q_injected_C
    # in-bucket: within +-180 deg @2856 of the core centroid AND inside the KE window (already core).
    # The +-180 deg-of-centroid window is the INTENDED (non-wrapped) bucket definition: the
    # captured-core charge beyond +-180 deg is the 2856 MHz phase-wrap tail (faithful, NOT a
    # bug-drop). For cond_i the core spans ~2 RF periods at 2856 MHz, so σz@2856 (≈163 deg) is
    # the faithful consequence of that wrap, not an error.
    in_bucket = good & (np.abs(phi) <= 180.0)
    q_in_bucket = float(wc[in_bucket].sum()) * Q_E
    in_bucket_pct = 100.0 * q_in_bucket / q_injected_C
    in_bucket_over_bore_pct = 100.0 * q_in_bucket / q_in_bore_C
    print(f"[D sec1]  zbar={zbar:.4f} m  <KE_core>={Ebar:.4f} MeV  sigE={sigE:.4f}  "
          f"core={int(core.sum())}/{len(z)}  vz-dropped={ndrop}/{len(zc)}")
    print(f"          capture in-bucket={in_bucket_pct:.2f}%  all-buckets={all_buckets_pct:.2f}%  "
          f"in-bucket/in-bore={in_bucket_over_bore_pct:.2f}%")
    _record("after_sec1", "D after Sec1", zbar, Ebar, sigE, sigt, in_bucket_pct, all_buckets_pct,
            sz214_blank=True)
    cap = (f"{op_label}\n{sign_note} (captured-core centroid); ref 2856 MHz; <z>={zbar:.4f} m (local frame); "
           f"captured core KE>={MIN_KE_MEV} MeV ({int(core.sum())}/{len(z)} parts); v_z-dropped={ndrop}/{len(zc)}; "
           f"\nSec1 POWER_MW=11 → peak |Ez| 14.89 MV/m, on-crest-average ≈ 11.5 MV/m "
           f"(≈ Tan cond-i 11 MV/m section field); "
           f"\nrepo capture ≈ {in_bucket_pct:.1f}% in-bucket / ≈ {all_buckets_pct:.1f}% all-buckets (of q_injected_C); "
           f"iris transmission ≈ {100.0*q_in_bore_C/q_injected_C:.0f}% (q_in_bore/q_injected) reported separately "
           f"(far below Tan's 89.4/96.8 BY CONSTRUCTION: different denominator + γ² self-field + real iris); "
           f"from Sec1 exit dump (before the 12 MeV linac_rest handoff cut)\n"
           f"Tan Fig 10.5 (expect steeply falling KE(φ) near 25–30 MeV; narrow core spike)")
    location_panel(phi[good], kec[good], wc[good], f"Tan Fig 10.5 — after Sec1 ({title_tag}, captured core)",
                   cap, "@2856 MHz", f"{RESULTS}/tan_fig10p5_after_sec1_{case}.png")

    # ── comparison table (md + csv) ──
    _write_table(case, table_rows, q_injected_C, q_in_bore_C, q_in_domain_C, in_bucket_over_bore_pct)


def _fmt(v, nd=3):
    if v is None:
        return ""
    if isinstance(v, str):
        return v
    return f"{v:.{nd}f}"


def _write_table(case, rows, q_injected_C, q_in_bore_C, q_in_domain_C, in_bucket_over_bore_pct):
    """Write results/tan_comparison_<case>.{md,csv}: repo + Tan-published + delta rows."""
    cols = ["location", "zbar_m", "Ebar_MeV", "sigE_MeV", "sigz_deg@214",
            "sigz_deg@2856", "capture_in_bucket_pct", "capture_all_buckets_pct"]

    def repo_cells(r):
        return [r["location"], _fmt(r["zbar"], 4), _fmt(r["Ebar"]), _fmt(r["sigE"]),
                _fmt(r["sz214"], 2), _fmt(r["sz2856"], 2),
                r["cap_in"] if isinstance(r["cap_in"], str) else _fmt(r["cap_in"], 2),
                r["cap_all"] if isinstance(r["cap_all"], str) else _fmt(r["cap_all"], 2)]

    def tan_cells(key):
        t = TAN_TABLE[key]
        return ["Tan (published)", "", _fmt(t["Ebar"]), _fmt(t["sigE"]),
                _fmt(t["sz214"], 2), _fmt(t["sz2856"], 2),
                _fmt(t["cap_in"], 1), _fmt(t["cap_all"], 1)]

    def delta_cells(r, key):
        t = TAN_TABLE[key]
        d = []
        for repo_v, tan_v in [(r["Ebar"], t["Ebar"]), (r["sigE"], t["sigE"]),
                              (r["sz214"], t["sz214"]), (r["sz2856"], t["sz2856"])]:
            d.append(_fmt(repo_v - tan_v) if (repo_v is not None and tan_v is not None) else "")
        return ["Δ (repo−Tan)", "", d[0], d[1], d[2], d[3], "", ""]

    # CSV: flat — repo / Tan / delta interleaved per location, with a 'row_type' column.
    csv_path = f"{RESULTS}/tan_comparison_{case}.csv"
    with open(csv_path, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["row_type"] + cols)
        for r in rows:
            wr.writerow(["repo"] + repo_cells(r))
            wr.writerow(["tan_published"] + tan_cells(r["key"]))
            wr.writerow(["delta"] + delta_cells(r, r["key"]))
    print(f"wrote {csv_path}")

    # Markdown: a table block per location.
    md_path = f"{RESULTS}/tan_comparison_{case}.md"
    lines = []
    lines.append(f"# Chapter 10 comparison — {case}\n")
    if case == "cond_i":
        lines.append("**Operating point:** Tan condition (i) — Preb1 50 kV @ phase-null "
                     "(PREB1_PHI_OFF=−90, 5.83 kW), Preb2 150 kV @ crest (PREB2_PHI_OFF=0, 36.61 kW); "
                     "achieved V_gap 50.0 / 150.0 kV. **This IS Tan's operating point — the repo numbers "
                     "below may be compared NUMERICALLY to Table 10.2** (within the plan §6 tolerances; "
                     "σ_E and capture diverge in a documented direction).\n")
    else:
        lines.append("**Operating point:** repo default 8 kW / 10 kW two-cavity "
                     "(Preb1 −70°, Preb2 −45°). **QUALITATIVE comparison only — NOT Tan cond (i)** "
                     "(50 kV @ phase-null / 150 kV @ crest). Tan's published numbers are a reference column, "
                     "NOT a target.\n")
    lines.append("")
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    for r in rows:
        lines.append(f"### {r['location']}")
        lines.append(header)
        lines.append(sep)
        lines.append("| " + " | ".join(repo_cells(r)) + " |")
        lines.append("| " + " | ".join(tan_cells(r["key"])) + " |")
        lines.append("| " + " | ".join(delta_cells(r, r["key"])) + " |")
        lines.append("")
    lines.append("## Capture bookkeeping (denominators)\n")
    lines.append(f"- `q_injected_C` (honest, pre-iris) = {q_injected_C:.4e} C")
    lines.append(f"- `q_in_bore_C` (post-iris)         = {q_in_bore_C:.4e} C")
    lines.append(f"- `q_in_domain_C`                   = {q_in_domain_C:.4e} C")
    lines.append(f"- after-Sec1 in-bucket / q_in_bore  = {in_bucket_over_bore_pct:.2f}% "
                 "(separates iris loss from capture loss)\n")
    lines.append("## Footnotes\n")
    lines.append("- **σ_z internal ratio:** the repo @2856/@214 columns derive from ONE σ_t, so their "
                 "ratio is BY CONSTRUCTION F_2856/F_214 = 13.335 — an internal-consistency check only. Tan's "
                 "own ratio is 13.52 (a different spatial-σ_z definition with a per-location β); the 13.335 "
                 "check is NOT agreement with Tan.")
    # Capture/iris numbers computed from THIS case's data (not hard-coded — they differ per case).
    d_row = next((r for r in rows if r["key"] == "after_sec1"), None)
    cap_in_v = d_row["cap_in"] if d_row else None
    cap_all_v = d_row["cap_all"] if d_row else None
    iris_pct = 100.0 * q_in_bore_C / q_injected_C
    cap_in_s = f"{cap_in_v:.1f}%" if isinstance(cap_in_v, (int, float)) else str(cap_in_v)
    cap_all_s = f"{cap_all_v:.1f}%" if isinstance(cap_all_v, (int, float)) else str(cap_all_v)
    lines.append("- **Denominator mismatch:** Tan's 100% (upstream) and 89.4/96.8 (postSec1) reference the "
                 "gun-emitted bunch with no loss yet. The repo's `q_injected_C` already sits past the converging "
                 f"halo and includes the real 9.547 mm iris scrape. Repo capture = {cap_in_s} in-bucket / {cap_all_s} "
                 f"all-buckets (both vs `q_injected_C`); iris transmission = {iris_pct:.1f}% (`q_in_bore`/`q_injected`) "
                 "is a SEPARATE quantity (do NOT conflate it with capture). The repo numbers are far below Tan's "
                 "BY CONSTRUCTION (γ² self-field + real iris + operating point). Annotate, do NOT 'fix'.")
    lines.append("- **Location C capture = n/a (pre-iris):** C is the pre-scrape injector population; its charge "
                 "differs from D's iris survivors, so the C→D ratio is not a clean repo capture fraction.")
    lines.append("- **after-Sec1 σ_z@214 blank:** Tan leaves it blank; the captured-core σ_z is reported @2856 only.")
    with open(md_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"wrote {md_path}")


if __name__ == "__main__":
    main()
