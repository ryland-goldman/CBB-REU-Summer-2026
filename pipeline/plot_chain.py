"""
Cross-stage beam-evolution figures for the whole Cornell Linac chain.

In-process (no pywarpx): reads each stage's existing openPMD particle series, builds
ONE per-dump moment table per stage (placed at each dump's own lab ⟨z⟩), and renders
four figures — all views of that single table:

  * results/chain_evolution.png      — 3×2 panels vs lab ⟨z⟩ across cathode→gun→
                                       injector→linac: ⟨KE⟩ (log-y), ε_n,x, σ_x/σ_r,
                                       σ_z (log-y), within-stage charge fraction, I_peak.
  * results/emittance_budget.png     — ε_n,x entry vs exit per stage (waterfall).
  * results/transmission_waterfall.png — the end-to-end charge chain with the TWO
                                       distinct loss stages (enters bore vs captured).
  * results/chain_scorecard.png      — per-stage entry/exit table (also printed to stdout).

Output lands in the repo-root `results/` (git-ignored by the existing `results/` line;
PNGs are committed with `git add -f results/*.png`).

Run with:
    conda run -n CBB python -c "import pipeline; pipeline.plot_chain()"
or it is called automatically at the end of pipeline/run_pipeline.py:main().

PHYSICS / UNITS NOTES (reviewer-flagged):
  - Transverse ε_n,x = sqrt(⟨x²⟩⟨ux²⟩ − ⟨x·ux⟩²)·1e6 [mm·mrad]. openPMD u = γβ already
    carries γ — NO extra γ multiply. The ×1e6 is the m·rad → mm·mrad scaling.
  - Longitudinal ε_n,z = sqrt(⟨z²⟩⟨uz²⟩ − ⟨z·uz⟩²) is mm·(dimensionless), NOT mm·mrad —
    so it gets a DIFFERENT scaling (×1e3 for the z[m]→mm only, uz is dimensionless γβ_z).
    Labelled as the z–(γβ_z) longitudinal emittance in mm.
  - The cathode is 2D x–z (no y/uy); its ε_n,x is the slab x-emittance. The cathode→gun
    seam on the ε_n panel is a 2D→RZ DEFINITIONAL discontinuity, NOT physical growth —
    annotated as such (downstream RZ ε_n,x is the projected emittance from reconstructed
    Cartesian x,y).
  - Capture is reported vs the TRUE injected charge (linac injection_summary.json q_inj),
    never the post-scrape first-dump charge. Absolute charge resets at each handoff
    (gun/injector downsample-reweight), so panel (5) shows WITHIN-stage transmission and the
    waterfall stitches the true end-to-end chain via the recorded denominators.
  - Lab-frame ES self-field overestimates transverse space charge by ~γ² (≈1.66× at β≈0.7);
    the injector/linac σ_r and capture numbers are conservative (pessimistic) — noted on the
    panels.
"""

import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from openpmd_viewer import OpenPMDTimeSeries

from pipeline.beam_metrics import rms_emit   # shared normalized-rms-emittance helper

c = 299792458.0
MC2_KEV = 0.51099895e3           # electron rest energy [keV]
Q_E = 1.602176634e-19           # elementary charge [C]
RESULTS = "results"             # repo-root results/ (git-ignored; git add -f)

# Each stage's openPMD particle series, in chain order. z0 is an escape-hatch lab-z
# shift (0 everywhere: the openPMD z is already lab-frame and the segments abut
# naturally — cathode ~0, gun ~0, injector 0.06–2.1, linac entrance ~2.03→3.5 m). The
# cathode is 2D x–z (no y/uy); the rest are RZ.
STAGES = [
    {"name": "cathode",  "path": "cathode/diags/particles",         "z0": 0.0, "geom": "2d", "color": "C0"},
    {"name": "gun",      "path": "gun/diags/particles",             "z0": 0.0, "geom": "rz", "color": "C1"},
    {"name": "injector", "path": "injector/diags/main/particles",   "z0": 0.0, "geom": "rz", "color": "C2"},
    {"name": "linac",    "path": "linac_sec1/diags/main/particles", "z0": 0.0, "geom": "rz", "color": "C3"},
]
LINAC_INJ_SUMMARY = "linac_sec1/diags/main/injection_summary.json"
Z_HANDOFF = 2.03                # [m] linac entrance (Z_acc_1); the injector→linac handoff plane


def _exit_row(name, rows):
    """The row representing a stage's EXIT.

    For the injector, the domain extends past the 2.03 m handoff to ZMAX=2.10 m and the run
    stops while the bunch is partially draining through the absorbing exit, so the
    largest-⟨z⟩ dump (rows[-1], ~2.076 m) is depleted and NOT what the linac ingests. The
    linac reader selects the dump nearest the 2.03 m handoff, so the cross-stage figures
    must use the SAME plane — otherwise the 'injector exit' charge disagrees with the linac
    input. Other stages end at their physical exit, so rows[-1] is correct there.
    """
    if name == "injector":
        return min(rows, key=lambda r: abs(r["z_mean"] - Z_HANDOFF))
    return rows[-1]


def _peak_current(z, w, v_beam, nbins=400):
    """Peak current = max line-charge density × beam velocity [A]."""
    if len(z) < 2 or z.max() <= z.min():
        return 0.0
    edges = np.linspace(z.min(), z.max(), nbins + 1)
    dz = edges[1] - edges[0]
    lam, _ = np.histogram(z, bins=edges, weights=w * Q_E)
    return float(lam.max() / dz * v_beam)


def build_moment_table(stage):
    """Per-dump beam-moment rows for one stage, sorted by ⟨z⟩.

    Returns a list of dicts (one per well-populated dump). Returns [] gracefully when the
    series is missing/empty (so plot_chain works when only some stages have run — matches
    the "run one stage off existing upstream output" workflow). The cathode is 2D: request
    only x/ux; the rest are RZ (x,y reconstructed from r → projected ε_n,x).
    """
    path = stage["path"]
    if not os.path.isdir(path):
        return []
    try:
        ts = OpenPMDTimeSeries(path)
    except Exception:
        return []
    rows = []
    is_rz = stage["geom"] == "rz"
    for it in ts.iterations:
        try:
            if is_rz:
                x, y, z, ux, uy, uz, w = ts.get_particle(
                    ["x", "y", "z", "ux", "uy", "uz", "w"], species="electrons", iteration=it)
            else:
                x, z, ux, uz, w = ts.get_particle(
                    ["x", "z", "ux", "uz", "w"], species="electrons", iteration=it)
                y = uy = None
        except Exception:
            continue
        if len(z) < 50:                       # skip near-empty (boundary) dumps
            continue
        gb = np.sqrt(1.0 + ux**2 + uy**2 + uz**2) if is_rz else np.sqrt(1.0 + ux**2 + uz**2)
        ke = (gb - 1.0) * MC2_KEV
        zm = float(np.average(z, weights=w))
        v_beam = float(np.average(uz / gb, weights=w) * c)
        sig_x = float(np.sqrt(np.average((x - np.average(x, weights=w)) ** 2, weights=w)))
        sig_z = float(np.sqrt(np.average((z - zm) ** 2, weights=w)))
        rows.append(dict(
            z_mean=zm + stage["z0"],
            ke_mean=float(np.average(ke, weights=w)),
            ke_std=float(np.sqrt(np.average((ke - np.average(ke, weights=w)) ** 2, weights=w))),
            emit_nx=rms_emit(x, ux, w) * 1e6,             # mm·mrad (transverse)
            emit_nz=rms_emit(z, uz, w) * 1e3,             # mm·(γβ_z) — NOT mm·mrad
            sig_x=sig_x, sig_z=sig_z,
            q=float(w.sum()) * Q_E,                       # absolute charge in THIS series [C]
            i_peak=_peak_current(z, w, v_beam),
        ))
    rows.sort(key=lambda r: r["z_mean"])
    return rows


def _arr(rows, key):
    return np.array([r[key] for r in rows])


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — chain_evolution.png : 3×2 panels of the whole chain vs lab ⟨z⟩.
# ══════════════════════════════════════════════════════════════════════════════
def render_chain_evolution(tables):
    fig, axs = plt.subplots(3, 2, figsize=(13, 11), constrained_layout=True)
    (a_ke, a_ex), (a_sx, a_sz), (a_q, a_ip) = axs

    # cathode→gun seam (first gun ⟨z⟩) for the 2D→RZ ε_n annotation.
    seam_z = tables["gun"][0]["z_mean"] * 1e3 if tables.get("gun") else None

    for st in STAGES:
        rows = tables.get(st["name"]) or []
        if not rows:
            continue
        z = _arr(rows, "z_mean") * 1e3
        col, nm = st["color"], st["name"]
        ke = _arr(rows, "ke_mean"); dke = _arr(rows, "ke_std")
        a_ke.plot(z, ke, "-", color=col, label=nm)
        a_ke.fill_between(z, np.maximum(ke - dke, 1e-3), ke + dke, color=col, alpha=0.18)
        a_ex.plot(z, _arr(rows, "emit_nx"), "-", color=col, label=nm)
        a_sx.plot(z, _arr(rows, "sig_x") * 1e3, "-", color=col, label=nm)
        a_sz.plot(z, np.maximum(_arr(rows, "sig_z") * 1e3, 1e-3), "-", color=col, label=nm)
        # within-stage charge fraction + I_peak: EXCLUDE the cathode. The cathode is an
        # emitter (q grows over time → a within-stage "transmission" rises >1, misleading)
        # and its 2D-slab pre-renorm I_peak is ~10 kA (non-physical), which on a linear axis
        # flatlines every real stage. Both panels start at the gun.
        if nm != "cathode":
            q = _arr(rows, "q")
            a_q.plot(z, q / q[0] if q[0] > 0 else q, "-", color=col, label=nm)
            a_ip.plot(z, _arr(rows, "i_peak"), "-", color=col, label=nm)

    a_ke.set_yscale("log"); a_ke.set_ylabel("⟨KE⟩  [keV]  (±σ band)")
    a_ke.set_title("Mean kinetic energy"); a_ke.legend(fontsize=8)
    a_ex.set_ylabel(r"$\varepsilon_{n,x}$  [mm·mrad]")
    a_ex.set_title("Transverse normalized emittance")
    if seam_z is not None:
        a_ex.axvline(seam_z, color="0.5", ls=":", lw=1)
        a_ex.annotate("cathode→gun: 2D→RZ\ndefinitional step (not physical)",
                      xy=(seam_z, a_ex.get_ylim()[1]), xytext=(0.30, 0.92),
                      textcoords="axes fraction", fontsize=7, color="0.3",
                      arrowprops=dict(arrowstyle="->", color="0.5", lw=0.8))
    a_sx.set_ylabel("σ_x / σ_r  [mm]"); a_sx.set_title("Transverse size")
    a_sx.annotate("σ_r conservative:\nES omits 1/γ² pinch (~γ²≈1.7×)",
                  xy=(0.50, 0.92), xycoords="axes fraction", fontsize=7, color="0.3")
    a_sz.set_yscale("log"); a_sz.set_ylabel("σ_z  [mm]"); a_sz.set_title("Bunch length")
    a_q.set_ylabel("q(z) / q(stage entry)")
    a_q.set_title("Within-stage charge fraction (gun→linac; cathode emitter excluded; q resets each handoff)")
    a_ip.set_ylabel("I_peak  [A]")
    a_ip.set_title("Peak current (gun→linac; cathode 2D-slab I_peak non-physical, excluded)")
    for ax in (a_ke, a_ex, a_sx, a_sz, a_q, a_ip):
        ax.set_xlabel("lab ⟨z⟩  [mm]"); ax.grid(alpha=0.25)

    fig.suptitle("Cornell Linac chain: beam evolution  (cathode → gun → injector → linac)",
                 fontsize=13)
    os.makedirs(RESULTS, exist_ok=True)
    path = f"{RESULTS}/chain_evolution.png"
    fig.savefig(path, dpi=140); plt.close(fig)
    print(f"wrote {path}")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — emittance_budget.png : ε_n,x entry vs exit per stage (waterfall).
# ══════════════════════════════════════════════════════════════════════════════
def render_emittance_budget(tables):
    names, e_in, e_out = [], [], []
    for st in STAGES:
        rows = tables.get(st["name"]) or []
        if not rows:
            continue
        names.append(st["name"]); e_in.append(rows[0]["emit_nx"])
        e_out.append(_exit_row(st["name"], rows)["emit_nx"])
    if not names:
        return
    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    x = np.arange(len(names)); wbar = 0.38
    ax.bar(x - wbar / 2, e_in, wbar, label="entry", color="C0")
    ax.bar(x + wbar / 2, e_out, wbar, label="exit", color="C3")
    for i, (ei, eo) in enumerate(zip(e_in, e_out)):
        ax.annotate(f"{ei:.2f}", (i - wbar / 2, ei), ha="center", va="bottom", fontsize=8)
        ax.annotate(f"{eo:.2f}", (i + wbar / 2, eo), ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(names)
    ax.set_ylabel(r"$\varepsilon_{n,x}$  [mm·mrad]")
    ax.set_title("Transverse emittance budget: entry vs exit per stage")
    ax.legend()
    ax.annotate("Footnotes: (1) the cathode→gun jump is a 2D→RZ DEFINITIONAL change (slab "
                "x-emittance → projected RZ), not physical emittance growth. (2) the injector "
                "ε_n growth is space-charge + solenoid-aberration dominated over the 2 m "
                "low-energy drift; the γ²≈1.7× ES transverse-SC overestimate makes it an UPPER "
                "bound (real growth is somewhat less — opposite sense to the capture lower bound).",
                xy=(0.0, -0.15), xycoords="axes fraction", fontsize=7, color="0.3")
    os.makedirs(RESULTS, exist_ok=True)
    path = f"{RESULTS}/emittance_budget.png"
    fig.savefig(path, dpi=140); plt.close(fig)
    print(f"wrote {path}")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — transmission_waterfall.png : the end-to-end charge chain, two loss stages.
# ══════════════════════════════════════════════════════════════════════════════
def render_transmission_waterfall(tables, linac_inj):
    """Charge milestones along the chain. Absolute charge resets at each downsample/
    handoff, so we show the within-stage-transmission-scaled end-to-end picture and,
    crucially, separate the injector's TWO loss stages: 'enters bore' (the 9.547 mm iris
    scrape) and 'captured' (the linac RF-bucket loss) — the separation that motivates the
    solenoids. Capture denominator = TRUE injected charge (linac injection_summary.json).

    The waterfall starts at GUN EXIT, where the beam carries its physical charge (the gun
    renormalizes to the 1 nC bunch). The CATHODE dump's weight sum (~82 nC of raw
    macroparticle weight, pre-renorm) is NOT a physical charge and is excluded — plotting it
    would dwarf the physical ≤1 nC bars on one axis; it's noted in the caption instead.

    The 'injector exit' bar uses the dump at the 2.03 m HANDOFF plane (via _exit_row), NOT
    the largest-⟨z⟩ dump (which is partially drained through the 2.10 m absorbing exit and
    under-counts what the linac actually ingests)."""
    bars, vals = [], []
    gun = tables.get("gun") or []
    inj = tables.get("injector") or []
    if gun:
        bars.append("gun exit\n(renorm ~1 nC)"); vals.append(_exit_row("gun", gun)["q"] * 1e9)
    if inj:
        bars.append("injector exit\n(@2.03m handoff)"); vals.append(_exit_row("injector", inj)["q"] * 1e9)
    # The two distinct downstream losses, anchored on the linac's recorded true-injected
    # breakdown (q_injected_C at the handoff, q_in_bore_C through the 9.547 mm iris) and the
    # captured charge from the last linac dump.
    if linac_inj:
        bars.append("enters bore\n(9.547mm iris)"); vals.append(linac_inj["q_in_bore_C"] * 1e9)
        lin = tables.get("linac") or []
        if lin:
            bars.append("captured\n(~26 MeV)"); vals.append(lin[-1]["q"] * 1e9)
    if not bars:
        return
    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    cols = ["C1", "C2", "C4", "C3"][:len(bars)]   # gun / injector / bore / captured
    ax.bar(range(len(bars)), vals, color=cols)
    for i, v in enumerate(vals):
        ax.annotate(f"{v:.4f} nC", (i, v), ha="center", va="bottom", fontsize=8)
    # Annotate the two loss steps that matter most (bore scrape, capture) vs true injected.
    if linac_inj:
        qinj = linac_inj["q_injected_C"] * 1e9
        ax.axhline(qinj, color="0.6", ls="--", lw=0.8)
        ax.annotate(f"true injected at handoff = {qinj:.3f} nC (capture denominator)",
                    xy=(0.02, qinj), fontsize=7, color="0.3", va="bottom")
    ax.set_xticks(range(len(bars))); ax.set_xticklabels(bars, fontsize=8)
    ax.set_ylabel("charge  [nC]")
    ax.set_title("End-to-end charge / transmission waterfall\n"
                 "(from gun exit; bore-scrape and capture are SEPARATE losses)")
    ax.annotate("Starts at gun exit (physical ~1 nC renorm); cathode dump weight (~82 nC, "
                "pre-renorm, not physical) excluded. 'injector exit' is the 2.03 m handoff dump. "
                "Capture vs TRUE injected; γ²≈1.7× ES self-field overestimate ⇒ a conservative LOWER bound.",
                xy=(0.0, -0.18), xycoords="axes fraction", fontsize=7, color="0.3")
    os.makedirs(RESULTS, exist_ok=True)
    path = f"{RESULTS}/transmission_waterfall.png"
    fig.savefig(path, dpi=140); plt.close(fig)
    print(f"wrote {path}")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 / #10 — chain_scorecard.png + stdout : per-stage entry/exit table.
# ══════════════════════════════════════════════════════════════════════════════
def render_scorecard(tables, linac_inj):
    cols = ["stage", "⟨KE⟩in", "⟨KE⟩out", "σ_KEout", "ε_nx,in", "ε_nx,out", "ε_nz,out[mm]",
            "σ_x,out[mm]", "σ_z,out[mm]", "q_out[nC]"]
    table_rows = []
    for st in STAGES:
        rows = tables.get(st["name"]) or []
        if not rows:
            continue
        a, b = rows[0], _exit_row(st["name"], rows)   # injector exit = 2.03 m handoff dump, not drained tail
        table_rows.append([
            st["name"], f"{a['ke_mean']:.1f}", f"{b['ke_mean']:.1f}", f"{b['ke_std']:.2f}",
            f"{a['emit_nx']:.2f}", f"{b['emit_nx']:.2f}", f"{b['emit_nz']:.2f}",
            f"{b['sig_x']*1e3:.2f}", f"{b['sig_z']*1e3:.2f}", f"{b['q']*1e9:.4f}",
        ])
    # Capture line vs true injected (the legible end-to-end number).
    cap_note = ""
    if linac_inj and (tables.get("linac")):
        qinj = linac_inj["q_injected_C"]; qcap = tables["linac"][-1]["q"]
        cap_note = (f"linac capture = {qcap/qinj*100:.2f}% of true injected "
                    f"({qcap*1e9:.4f}/{qinj*1e9:.4f} nC); "
                    f"iris transmission = {linac_inj['q_in_bore_C']/qinj*100:.1f}%. "
                    f"γ²≈1.7× → capture is a conservative LOWER bound. σ_KE charge-conditional.")
    # Two reader notes (physics-flagged) so adjacent-dump and emittance effects aren't misread:
    note_handoff = ("the injector-exit row is the dump at the 2.03 m handoff plane (same dump the "
                    "linac reader ingests), not the drained tail at the 2.10 m absorbing exit; any "
                    "small ⟨KE⟩ difference vs the linac-entry row is dump spacing, not a discontinuity.")
    note_emit = ("injector ε_n,x growth is space-charge + solenoid-aberration dominated over the "
                 "2 m low-energy drift; the γ²≈1.7× ES transverse-SC overestimate makes this an "
                 "UPPER bound on emittance growth (opposite direction to the capture lower bound).")

    # ── stdout ──
    print("\n" + "=" * 100)
    print("CHAIN SCORECARD  (KE in keV; entry=first dump, exit=last dump)")
    print("-" * 100)
    print("  ".join(f"{c:>11}" for c in cols))
    for r in table_rows:
        print("  ".join(f"{v:>11}" for v in r))
    if cap_note:
        print("-" * 100); print("  " + cap_note)
    print("  " + note_handoff)
    print("  " + note_emit)
    print("=" * 100)

    # ── figure (table render) ──
    fig, ax = plt.subplots(figsize=(13, 1.6 + 0.5 * (len(table_rows) + 1)),
                           constrained_layout=True)
    ax.axis("off")
    tbl = ax.table(cellText=table_rows, colLabels=cols, loc="center", cellLoc="center")
    tbl.auto_set_font_size(False); tbl.set_fontsize(9); tbl.scale(1, 1.5)
    ttl = "Cornell Linac chain scorecard"
    if cap_note:
        ttl += "\n" + cap_note
    ax.set_title(ttl, fontsize=10)
    # The two reader notes as a footnote below the table.
    ax.annotate(note_handoff + "\n" + note_emit, xy=(0.5, 0.0), xycoords="axes fraction",
                ha="center", va="top", fontsize=7, color="0.3")
    os.makedirs(RESULTS, exist_ok=True)
    path = f"{RESULTS}/chain_scorecard.png"
    fig.savefig(path, dpi=140); plt.close(fig)
    print(f"wrote {path}")


def main():
    """Build the moment table once per stage, then render all cross-stage figures."""
    tables = {st["name"]: build_moment_table(st) for st in STAGES}
    present = [n for n, r in tables.items() if r]
    if not present:
        print("plot_chain: no stage particle series found — run the pipeline first.")
        return
    print(f"plot_chain: building cross-stage figures from stages {present}")
    linac_inj = None
    if os.path.isfile(LINAC_INJ_SUMMARY):
        with open(LINAC_INJ_SUMMARY) as fh:
            linac_inj = json.load(fh)
    render_chain_evolution(tables)
    render_emittance_budget(tables)
    render_transmission_waterfall(tables, linac_inj)
    render_scorecard(tables, linac_inj)


if __name__ == "__main__":
    main()
