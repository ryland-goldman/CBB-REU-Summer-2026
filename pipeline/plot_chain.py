"""
Cross-stage beam-evolution figures for the whole Cornell Linac chain.

In-process (no pywarpx): reads each stage's existing openPMD particle series, builds
ONE per-dump moment table per stage, and renders four figures into the repo-root
results/ (chain_evolution, emittance_budget, transmission_waterfall, chain_scorecard).

See pipeline/README.md for the figures, units conventions, and z0/frame handling.
"""

import os
import json
import textwrap
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from openpmd_viewer import OpenPMDTimeSeries

from pipeline.beam_metrics import rms_emit   # shared normalized-rms-emittance helper
from pipeline.constants import C_LIGHT as c, MC2_EV

MC2_KEV = MC2_EV / 1e3           # electron rest energy [keV]
Q_E = 1.602176634e-19           # elementary charge [C]
RESULTS = "results"             # repo-root results/ (git-ignored; git add -f)

# z0 = lab-z shift per dump. cathode/gun/injector are already lab-frame (z0=0); linac
# stages are sim-LOCAL frame and get z0 filled at runtime (see _apply_linac*_z0).
STAGES = [
    {"name": "cathode",    "path": "cathode/diags/particles",         "z0": 0.0, "geom": "2d", "color": "C0"},
    {"name": "gun",        "path": "gun/diags/particles",             "z0": 0.0, "geom": "rz", "color": "C1"},
    {"name": "injector",   "path": "injector/diags/main/particles",   "z0": 0.0, "geom": "rz", "color": "C2"},
    {"name": "linac",      "path": "linac_sec1/diags/main/particles", "z0": 0.0, "geom": "rz", "color": "C3"},
    # linac_rest (sections 2–8, Impact-T): local-frame z0 filled at runtime. Color C5 NOT
    # C4 — C4 is the iris bar in the transmission waterfall.
    {"name": "linac_rest", "path": "linac_rest/diags/main/particles", "z0": 0.0, "geom": "rz", "color": "C5"},
]
LINAC_INJ_SUMMARY = "linac_sec1/diags/main/injection_summary.json"
LINAC_REST_INJ_SUMMARY = "linac_rest/diags/main/injection_summary.json"
Z_HANDOFF = 2.03                # [m] linac entrance (Z_acc_1); the injector→linac handoff plane


def _apply_linac_z0(linac_inj):
    """Set the linac stage's lab-z offset (local frame → lab).

    The linac sim resets the imported beam to a linac-LOCAL frame, so its openPMD z is NOT
    lab-frame. Offset = z_handoff_m − z_inject_mean_m (both recorded). Falls back to
    Z_HANDOFF (older runs lacking the fields).
    """
    linac = next(st for st in STAGES if st["name"] == "linac")
    if linac_inj and "z_handoff_m" in linac_inj and "z_inject_mean_m" in linac_inj:
        linac["z0"] = linac_inj["z_handoff_m"] - linac_inj["z_inject_mean_m"]
    else:
        linac["z0"] = Z_HANDOFF


def _apply_linac_rest_z0(rest_inj, tables=None):
    """Set the linac_rest (sections 2–8, Impact-T) stage's lab-z offset (local → lab).

    Offset = z_inject_lab_m − z_inject_local_m (recorded; z_inject_local_m defaults to 0
    since the sim zeroes z at injection). Fallback (no summary): linac_sec1's exit lab-z
    from the tables (its z0 + EXIT dump local ⟨z⟩), then the linac stage's z0.
    """
    rest = next(st for st in STAGES if st["name"] == "linac_rest")
    if rest_inj and "z_inject_lab_m" in rest_inj:
        rest["z0"] = rest_inj["z_inject_lab_m"] - rest_inj.get("z_inject_local_m", 0.0)
        return
    linac_rows = (tables or {}).get("linac") or []
    if linac_rows:
        rest["z0"] = linac_rows[-1]["z_mean"]    # z_mean already carries the linac z0
    else:
        linac = next(st for st in STAGES if st["name"] == "linac")
        rest["z0"] = linac["z0"]


def _exit_row(name, rows):
    """The row representing a stage's EXIT.

    Injector exit = dump nearest the 2.03 m handoff, NOT rows[-1]: the domain extends to
    ZMAX=2.10 m and the run stops mid-drain, so the largest-⟨z⟩ dump is depleted and the
    linac reader (which selects the 2.03 m dump) ingests a different one. Other stages end
    at their physical exit, so rows[-1] is correct.
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
        # σ_z + I_peak exclude linac_rest (only 2 Impact-T dumps → meaningless 2-point line).
        if nm != "linac_rest":
            a_sz.plot(z, np.maximum(_arr(rows, "sig_z") * 1e3, 1e-3), "-", color=col, label=nm)
        # charge-fraction + I_peak exclude the cathode (emitter q grows >1; 2D-slab I_peak non-physical).
        if nm != "cathode":
            q = _arr(rows, "q")
            a_q.plot(z, q / q[0] if q[0] > 0 else q, "-", color=col, label=nm)
            if nm != "linac_rest":
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
    a_sx.set_ylabel("σ_x  [mm]"); a_sx.set_title("Transverse size (per-plane RMS)")
    a_sx.annotate("gun+injector: EMS (γ² pinch in);\nlinac_sec1 ES but SC small @25 MeV",
                  xy=(0.50, 0.92), xycoords="axes fraction", fontsize=7, color="0.3")
    a_sz.set_yscale("log"); a_sz.set_ylabel("σ_z  [mm]")
    a_sz.set_title("Bunch length (linac_rest excluded: only 2 Impact-T dumps)")
    a_q.set_ylabel("q(z) / q(stage entry)")
    a_q.set_title("Within-stage charge fraction (gun→linac; cathode emitter excluded; q resets each handoff)")
    a_ip.set_ylabel("I_peak  [A]")
    a_ip.set_title("Peak current (gun→linac; cathode 2D-slab I_peak non-physical and "
                   "2-dump linac_rest excluded)")
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
    # Wrap to the figure width — a single long line makes constrained_layout collapse the
    # axes ("axes sizes collapsed to zero") and the text run off-canvas.
    footnote = (
        "Footnotes: (1) the cathode→gun jump is a 2D→RZ DEFINITIONAL change: the slab "
        "x-emittance (uniform x∈[−R,R], ⟨x²⟩=R²/3) becomes the gun's projected RZ "
        "emittance after the r-importance resample builds a uniform DISC (⟨x²⟩=R²/4), so "
        "ε_n,x drops ×√(3/4)≈0.87 (~2.3→~2.0 mm·mrad) — a geometry correction (the disc is "
        "more physical), NOT physical growth. (2) the injector ε_n growth is space-charge + "
        "solenoid-aberration dominated over the 2 m low-energy drift, now with the relativistic "
        "EMS self-field (the 1/γ² magnetic pinch is included, so it is no longer a γ²-inflated "
        "upper bound). (3) the injector-exit bar is the "
        "UN-collimated 2.03 m handoff beam (no iris mask); the iris-survivor beam "
        "linac_sec1 actually receives is ~13% lower (≈326 vs ≈375 mm·mrad).")
    ax.annotate(textwrap.fill(footnote, width=150),
                xy=(0.0, -0.12), xycoords="axes fraction", va="top",
                fontsize=7, color="0.3")
    os.makedirs(RESULTS, exist_ok=True)
    path = f"{RESULTS}/emittance_budget.png"
    fig.savefig(path, dpi=140); plt.close(fig)
    print(f"wrote {path}")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — transmission_waterfall.png : the end-to-end charge chain, two loss stages.
# ══════════════════════════════════════════════════════════════════════════════
def render_transmission_waterfall(tables, linac_inj):
    """Charge milestones along the chain, separating the two loss stages (iris scrape vs
    captured). Capture denominator = TRUE injected charge. Starts at gun exit (physical
    charge); cathode pre-renorm weight is excluded. Injector-exit bar uses the recorded
    handoff q_injected_C (fallback: the 2.03 m _exit_row, NOT the drained largest-⟨z⟩ dump)."""
    bars, vals = [], []
    gun = tables.get("gun") or []
    inj = tables.get("injector") or []
    if gun:
        bars.append("gun exit\n(renorm ~1 nC)"); vals.append(_exit_row("gun", gun)["q"] * 1e9)
    if inj:
        # Prefer the recorded handoff q_injected_C (one source w/ the iris/injected bars;
        # _exit_row's nearest-2.03m selector differs from the reader's population-gated pick).
        bars.append("injector exit\n(@2.03m handoff)")
        if linac_inj and "q_injected_C" in linac_inj:
            vals.append(linac_inj["q_injected_C"] * 1e9)
        else:
            vals.append(_exit_row("injector", inj)["q"] * 1e9)
    # Iris bar uses q_in_domain_C (9.547 mm multi-plane iris scrape), NOT q_in_bore_C (9.55 mm RF bore).
    if linac_inj:
        bars.append("passes iris\n(9.547mm)"); vals.append(linac_inj["q_in_domain_C"] * 1e9)
        lin = tables.get("linac") or []
        if lin:
            bars.append("captured\n(~26 MeV)"); vals.append(lin[-1]["q"] * 1e9)
    if not bars:
        return
    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    cols = ["C1", "C2", "C4", "C3"][:len(bars)]   # gun / injector / iris / captured
    ax.bar(range(len(bars)), vals, color=cols)
    for i, v in enumerate(vals):
        ax.annotate(f"{v:.4f} nC", (i, v), ha="center", va="bottom", fontsize=8)
    # Annotate the two loss steps that matter most (iris scrape, capture) vs true injected.
    if linac_inj:
        qinj = linac_inj["q_injected_C"] * 1e9
        ax.axhline(qinj, color="0.6", ls="--", lw=0.8)
        ax.annotate(f"true injected at handoff = {qinj:.3f} nC (capture denominator)",
                    xy=(len(bars) - 0.55, qinj), fontsize=7, color="0.3",
                    va="bottom", ha="right")
    ax.set_xticks(range(len(bars))); ax.set_xticklabels(bars, fontsize=8)
    ax.set_ylabel("charge  [nC]")
    ax.set_title("End-to-end charge / transmission waterfall\n"
                 "(from gun exit; bore-scrape and capture are SEPARATE losses)")
    # Wrapped for constrained_layout — see the emittance-budget footnote note.
    footnote = (
        "Starts at gun exit (physical ~1 nC renorm); cathode dump weight (~82 nC, "
        "pre-renorm, not physical) excluded. 'injector exit' is the recorded 2.03 m "
        "handoff charge (q_injected_C; dump fallback). "
        "Capture vs TRUE injected (injector now relativistic EMS; linac_sec1 SC small at 25 MeV).")
    ax.annotate(textwrap.fill(footnote, width=165),
                xy=(0.0, -0.13), xycoords="axes fraction", va="top",
                fontsize=7, color="0.3")
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
                    f"iris transmission = {linac_inj['q_in_domain_C']/qinj*100:.1f}% "
                    f"(multi-plane 9.547 mm scrape). "
                    f"capture vs TRUE injected (injector now EMS; linac_sec1 SC small @25 MeV). "
                    f"σ_KE charge-conditional.")
    # Two reader notes (physics-flagged) so adjacent-dump and emittance effects aren't misread:
    note_handoff = ("the injector-exit row is the dump at the 2.03 m handoff plane, not the drained "
                    "tail at the 2.10 m absorbing exit. NOTE: linac_sec1's selector applies an "
                    "n≥0.8·nmax population gate, so under heavy radial loss it can fall back to an "
                    "earlier off-plane dump rather than this 2.03 m one (the deferred re-tune); any "
                    "small ⟨KE⟩ difference vs the linac-entry row is dump spacing, not a discontinuity.")
    note_emit = ("injector ε_n,x growth is space-charge + solenoid-aberration dominated over the "
                 "2 m low-energy drift, now with the relativistic EMS self-field (the 1/γ² magnetic "
                 "pinch is included, so it is no longer a γ²-inflated upper bound).")

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
    # Raise the fd limit: openpmd-viewer leaks an fd per get_particle and the per-stage
    # loops exhaust macOS's 256-fd default (standalone entry point — no in-process raise).
    from pipeline._runner import _raise_fd_limit
    _raise_fd_limit()
    linac_inj = None
    if os.path.isfile(LINAC_INJ_SUMMARY):
        with open(LINAC_INJ_SUMMARY) as fh:
            linac_inj = json.load(fh)
    rest_inj = None
    if os.path.isfile(LINAC_REST_INJ_SUMMARY):
        with open(LINAC_REST_INJ_SUMMARY) as fh:
            rest_inj = json.load(fh)
    _apply_linac_z0(linac_inj)
    # Build all but linac_rest first: its lab z0 (fallback) derives from the linac table, so
    # resolve z0 then build its rows ONCE (so they carry the offset).
    tables = {st["name"]: build_moment_table(st)
              for st in STAGES if st["name"] != "linac_rest"}
    _apply_linac_rest_z0(rest_inj, tables)
    rest_stage = next(st for st in STAGES if st["name"] == "linac_rest")
    tables["linac_rest"] = build_moment_table(rest_stage)
    present = [n for n, r in tables.items() if r]
    if not present:
        print("plot_chain: no stage particle series found — run the pipeline first.")
        return
    print(f"plot_chain: building cross-stage figures from stages {present}")
    render_chain_evolution(tables)
    render_emittance_budget(tables)
    render_transmission_waterfall(tables, linac_inj)
    render_scorecard(tables, linac_inj)


if __name__ == "__main__":
    main()
