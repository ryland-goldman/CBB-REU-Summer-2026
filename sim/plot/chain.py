"""Cross-stage beam-evolution figures for the whole chain.

Reads each stage's openPMD particle series + injection_summary.json sidecars, builds one
per-dump moment table per stage, and renders four figures into logs/plots/chain/
(chain_evolution, emittance_budget, transmission_waterfall, chain_scorecard).
"""

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

import json
import textwrap

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from openpmd_viewer import OpenPMDTimeSeries

from sim.helpers.tools import C_LIGHT as c, MC2_KEV, E_CHARGE as Q_E, prepare_env
from sim.helpers.metrics import rms_emit

RESULTS = "logs/plots/chain"
Z_HANDOFF = 2.03                 # [m] injector->linac handoff plane
_SECTION_PITCH = 3.2             # nominal SLAC section+drift pitch for z0 fallbacks

# z0 = lab-z shift per dump. cathode/gun/injector are lab-frame (z0=0); the linac sections
# (sim-local frame) get z0 filled at runtime from their injection_summary.json.
STAGES = [
    {"name": "cathode",  "path": "logs/diags/cathode/particles",            "z0": 0.0, "geom": "2d", "color": "C0"},
    {"name": "gun",      "path": "logs/diags/gun/particles",                "z0": 0.0, "geom": "rz", "color": "C1"},
    {"name": "injector", "path": "logs/diags/injector/main/particles",      "z0": 0.0, "geom": "rz", "color": "C2"},
    {"name": "linac1",   "path": "logs/diags/linac1-3/sec1/main/particles", "z0": 0.0, "geom": "rz", "color": "C3"},
    {"name": "linac2",   "path": "logs/diags/linac1-3/sec2/main/particles", "z0": 0.0, "geom": "rz", "color": "C6"},
    {"name": "linac3",   "path": "logs/diags/linac1-3/sec3/main/particles", "z0": 0.0, "geom": "rz", "color": "C7"},
    {"name": "linac4-8", "path": "logs/diags/linac4-8/main/particles",      "z0": 0.0, "geom": "rz", "color": "C5"},
]
LINAC_SECTION_SUMMARIES = {
    "linac1": "logs/diags/linac1-3/sec1/main/injection_summary.json",
    "linac2": "logs/diags/linac1-3/sec2/main/injection_summary.json",
    "linac3": "logs/diags/linac1-3/sec3/main/injection_summary.json",
}
LINAC_INJ_SUMMARY = LINAC_SECTION_SUMMARIES["linac1"]            # the iris/capture summary (sec1)
LINAC_REST_INJ_SUMMARY = "logs/diags/linac4-8/main/injection_summary.json"


def _section_z0(summary):
    """local->lab z offset for a linac section from its summary: (z_inject_lab_m | z_handoff_m)
    - z_inject_mean_m. None if unavailable."""
    if not summary:
        return None
    z_lab = summary.get("z_inject_lab_m", summary.get("z_handoff_m"))
    if z_lab is None or "z_inject_mean_m" not in summary:
        return None
    return float(z_lab) - float(summary["z_inject_mean_m"])


def _apply_linac_z0(summaries):
    """Set the lab-z offset for each WarpX linac section (sec1->Z_HANDOFF fallback, then +pitch)."""
    fallback = Z_HANDOFF
    for name in ("linac1", "linac2", "linac3"):
        st = next(s for s in STAGES if s["name"] == name)
        z0 = _section_z0(summaries.get(name))
        st["z0"] = z0 if z0 is not None else fallback
        fallback = st["z0"] + _SECTION_PITCH


def _apply_linac_rest_z0(rest_inj, tables=None):
    """Set the linac4-8 (Impact-T) stage lab-z offset (z_inject_lab_m - z_inject_local_m)."""
    rest = next(st for st in STAGES if st["name"] == "linac4-8")
    if rest_inj and "z_inject_lab_m" in rest_inj:
        rest["z0"] = rest_inj["z_inject_lab_m"] - rest_inj.get("z_inject_local_m", 0.0)
        return
    sec3_rows = (tables or {}).get("linac3") or []
    if sec3_rows:
        rest["z0"] = sec3_rows[-1]["z_mean"]
    else:
        rest["z0"] = next(st for st in STAGES if st["name"] == "linac3")["z0"]


def _exit_row(name, rows):
    """The row representing a stage's EXIT. Injector exit = dump nearest the 2.03 m handoff
    (the run drains past it to ZMAX); other stages end at their physical exit (rows[-1])."""
    if name == "injector":
        return min(rows, key=lambda r: abs(r["z_mean"] - Z_HANDOFF))
    return rows[-1]


def _peak_current(z, w, v_beam, nbins=400):
    """Peak current = max line-charge density * beam velocity [A]."""
    if len(z) < 2 or z.max() <= z.min():
        return 0.0
    edges = np.linspace(z.min(), z.max(), nbins + 1)
    dz = edges[1] - edges[0]
    lam, _ = np.histogram(z, bins=edges, weights=w * Q_E)
    return float(lam.max() / dz * v_beam)


def build_moment_table(stage):
    """Per-dump beam-moment rows for one stage, sorted by <z>. [] if the series is missing."""
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
        if len(z) < 50:
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
            emit_nx=rms_emit(x, ux, w) * 1e6,             # mm.mrad (transverse)
            emit_nz=rms_emit(z, uz, w) * 1e3,             # mm.(gamma*beta_z)
            sig_x=sig_x, sig_z=sig_z,
            q=float(w.sum()) * Q_E,
            i_peak=_peak_current(z, w, v_beam),
        ))
    rows.sort(key=lambda r: r["z_mean"])
    return rows


def _arr(rows, key):
    return np.array([r[key] for r in rows])


# ── FIGURE 1: chain_evolution.png (3x2 panels vs lab <z>) ────────────────────────
def render_chain_evolution(tables):
    fig, axs = plt.subplots(3, 2, figsize=(13, 11), constrained_layout=True)
    (a_ke, a_ex), (a_sx, a_sz), (a_q, a_ip) = axs
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
        if nm != "linac4-8":           # sig_z + i_peak excluded (only 2 Impact-T dumps)
            a_sz.plot(z, np.maximum(_arr(rows, "sig_z") * 1e3, 1e-3), "-", color=col, label=nm)
        if nm != "cathode":            # cathode emitter q grows >1; 2D-slab I_peak non-physical
            q = _arr(rows, "q")
            a_q.plot(z, q / q[0] if q[0] > 0 else q, "-", color=col, label=nm)
            if nm != "linac4-8":
                a_ip.plot(z, _arr(rows, "i_peak"), "-", color=col, label=nm)

    a_ke.set_yscale("log"); a_ke.set_ylabel("<KE>  [keV]  (+-sigma band)")
    a_ke.set_title("Mean kinetic energy"); a_ke.legend(fontsize=8)
    a_ex.set_ylabel(r"$\varepsilon_{n,x}$  [mm.mrad]")
    a_ex.set_title("Transverse normalized emittance")
    if seam_z is not None:
        a_ex.axvline(seam_z, color="0.5", ls=":", lw=1)
        a_ex.annotate("cathode->gun: 2D->RZ\ndefinitional step (not physical)",
                      xy=(seam_z, a_ex.get_ylim()[1]), xytext=(0.30, 0.92),
                      textcoords="axes fraction", fontsize=7, color="0.3",
                      arrowprops=dict(arrowstyle="->", color="0.5", lw=0.8))
    a_sx.set_ylabel("sigma_x  [mm]"); a_sx.set_title("Transverse size (per-plane RMS)")
    a_sz.set_yscale("log"); a_sz.set_ylabel("sigma_z  [mm]")
    a_sz.set_title("Bunch length (linac4-8 excluded: only 2 Impact-T dumps)")
    a_q.set_ylabel("q(z) / q(stage entry)")
    a_q.set_title("Within-stage charge fraction (gun->linac; q resets each handoff)")
    a_ip.set_ylabel("I_peak  [A]")
    a_ip.set_title("Peak current (gun->linac; cathode + 2-dump linac4-8 excluded)")
    for ax in (a_ke, a_ex, a_sx, a_sz, a_q, a_ip):
        ax.set_xlabel("lab <z>  [mm]"); ax.grid(alpha=0.25)

    fig.suptitle("Cornell Linac chain: beam evolution  (cathode -> gun -> injector -> "
                 "linac1/2/3 -> linac4-8)", fontsize=13)
    os.makedirs(RESULTS, exist_ok=True)
    path = f"{RESULTS}/chain_evolution.png"
    fig.savefig(path, dpi=140); plt.close(fig)
    print(f"wrote {path}")


# ── FIGURE 2: emittance_budget.png (entry vs exit per stage) ─────────────────────
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
    ax.set_ylabel(r"$\varepsilon_{n,x}$  [mm.mrad]")
    ax.set_title("Transverse emittance budget: entry vs exit per stage")
    ax.legend()
    footnote = (
        "(1) cathode->gun is a 2D->RZ DEFINITIONAL change (slab x-emittance -> projected RZ "
        "emittance after the r-importance resample; ~x0.87), not physical growth. (2) injector "
        "growth is space-charge + solenoid-aberration over the low-energy drift (relativistic EMS "
        "self-field). (3) the injector-exit bar is the UN-collimated 2.03 m handoff (no iris mask); "
        "the iris-survivor beam linac1 receives is lower.")
    ax.annotate(textwrap.fill(footnote, width=150),
                xy=(0.0, -0.12), xycoords="axes fraction", va="top", fontsize=7, color="0.3")
    os.makedirs(RESULTS, exist_ok=True)
    path = f"{RESULTS}/emittance_budget.png"
    fig.savefig(path, dpi=140); plt.close(fig)
    print(f"wrote {path}")


# ── FIGURE 3: transmission_waterfall.png (the charge chain, two loss stages) ─────
def render_transmission_waterfall(tables, linac_inj):
    bars, vals = [], []
    gun = tables.get("gun") or []
    inj = tables.get("injector") or []
    if gun:
        bars.append("gun exit\n(renorm ~1 nC)"); vals.append(_exit_row("gun", gun)["q"] * 1e9)
    if inj:
        bars.append("injector exit\n(@2.03m handoff)")
        if linac_inj and "q_injected_C" in linac_inj:
            vals.append(linac_inj["q_injected_C"] * 1e9)
        else:
            vals.append(_exit_row("injector", inj)["q"] * 1e9)
    if linac_inj and "q_in_domain_C" in linac_inj:
        bars.append("passes iris\n(9.547mm)"); vals.append(linac_inj["q_in_domain_C"] * 1e9)
        lin = tables.get("linac1") or []
        if lin:
            bars.append(f"captured\n(~{lin[-1]['ke_mean']/1e3:.0f} MeV)"); vals.append(lin[-1]["q"] * 1e9)
    if not bars:
        return
    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    cols = ["C1", "C2", "C4", "C3"][:len(bars)]
    ax.bar(range(len(bars)), vals, color=cols)
    for i, v in enumerate(vals):
        ax.annotate(f"{v:.4f} nC", (i, v), ha="center", va="bottom", fontsize=8)
    if linac_inj and "q_injected_C" in linac_inj:
        qinj = linac_inj["q_injected_C"] * 1e9
        ax.axhline(qinj, color="0.6", ls="--", lw=0.8)
        ax.annotate(f"true injected at handoff = {qinj:.3f} nC (capture denominator)",
                    xy=(len(bars) - 0.55, qinj), fontsize=7, color="0.3", va="bottom", ha="right")
    ax.set_xticks(range(len(bars))); ax.set_xticklabels(bars, fontsize=8)
    ax.set_ylabel("charge  [nC]")
    ax.set_title("End-to-end charge / transmission waterfall\n"
                 "(from gun exit; bore-scrape and capture are SEPARATE losses)")
    footnote = (
        "Starts at gun exit (physical ~1 nC renorm); cathode dump weight (pre-renorm) excluded. "
        "'injector exit' is the recorded 2.03 m handoff charge. Capture vs TRUE injected.")
    ax.annotate(textwrap.fill(footnote, width=165),
                xy=(0.0, -0.13), xycoords="axes fraction", va="top", fontsize=7, color="0.3")
    os.makedirs(RESULTS, exist_ok=True)
    path = f"{RESULTS}/transmission_waterfall.png"
    fig.savefig(path, dpi=140); plt.close(fig)
    print(f"wrote {path}")


# ── FIGURE 4: chain_scorecard.png + stdout (per-stage entry/exit table) ──────────
def render_scorecard(tables, linac_inj):
    cols = ["stage", "<KE>in", "<KE>out", "sig_KEout", "enx,in", "enx,out", "enz,out[mm]",
            "sx,out[mm]", "sz,out[mm]", "q_out[nC]"]
    table_rows = []
    for st in STAGES:
        rows = tables.get(st["name"]) or []
        if not rows:
            continue
        a, b = rows[0], _exit_row(st["name"], rows)
        table_rows.append([
            st["name"], f"{a['ke_mean']:.1f}", f"{b['ke_mean']:.1f}", f"{b['ke_std']:.2f}",
            f"{a['emit_nx']:.2f}", f"{b['emit_nx']:.2f}", f"{b['emit_nz']:.2f}",
            f"{b['sig_x']*1e3:.2f}", f"{b['sig_z']*1e3:.2f}", f"{b['q']*1e9:.4f}",
        ])
    cap_note = ""
    if linac_inj and tables.get("linac1") and "q_injected_C" in linac_inj:
        qinj = linac_inj["q_injected_C"]; qcap = tables["linac1"][-1]["q"]
        cap_note = (f"linac capture = {qcap/qinj*100:.2f}% of true injected "
                    f"({qcap*1e9:.4f}/{qinj*1e9:.4f} nC); "
                    f"iris transmission = {linac_inj['q_in_domain_C']/qinj*100:.1f}% "
                    f"(multi-plane 9.547 mm scrape).")

    print("\n" + "=" * 100)
    print("CHAIN SCORECARD  (KE in keV; entry=first dump, exit=last dump)")
    print("-" * 100)
    print("  ".join(f"{c:>11}" for c in cols))
    for r in table_rows:
        print("  ".join(f"{v:>11}" for v in r))
    if cap_note:
        print("-" * 100); print("  " + cap_note)
    print("=" * 100)

    fig, ax = plt.subplots(figsize=(13, 1.6 + 0.5 * (len(table_rows) + 1)), constrained_layout=True)
    ax.axis("off")
    tbl = ax.table(cellText=table_rows, colLabels=cols, loc="center", cellLoc="center")
    tbl.auto_set_font_size(False); tbl.set_fontsize(9); tbl.scale(1, 1.5)
    ttl = "Cornell Linac chain scorecard"
    if cap_note:
        ttl += "\n" + cap_note
    ax.set_title(ttl, fontsize=10)
    os.makedirs(RESULTS, exist_ok=True)
    path = f"{RESULTS}/chain_scorecard.png"
    fig.savefig(path, dpi=140); plt.close(fig)
    print(f"wrote {path}")


def main():
    prepare_env()                       # repo-root chdir + fd-limit raise (openpmd-viewer leaks fds)
    summaries = {}
    for name, path in LINAC_SECTION_SUMMARIES.items():
        if os.path.isfile(path):
            with open(path) as fh:
                summaries[name] = json.load(fh)
    linac_inj = summaries.get("linac1")
    rest_inj = None
    if os.path.isfile(LINAC_REST_INJ_SUMMARY):
        with open(LINAC_REST_INJ_SUMMARY) as fh:
            rest_inj = json.load(fh)
    _apply_linac_z0(summaries)
    tables = {st["name"]: build_moment_table(st)
              for st in STAGES if st["name"] != "linac4-8"}
    _apply_linac_rest_z0(rest_inj, tables)
    tables["linac4-8"] = build_moment_table(next(st for st in STAGES if st["name"] == "linac4-8"))
    present = [n for n, r in tables.items() if r]
    if not present:
        print("chain: no stage particle series found -- run the pipeline first.")
        return
    print(f"chain: building cross-stage figures from stages {present}")
    render_chain_evolution(tables)
    render_emittance_budget(tables)
    render_transmission_waterfall(tables, linac_inj)
    render_scorecard(tables, linac_inj)


if __name__ == "__main__":
    main()
