"""
Figures for the SLAC Linac Section 1 stage (linac_sec1_sim.py).

Reads the traveling-wave field maps (linac_sec1/linac_sec1_field/linac_rf{1,2}.h5)
and the openPMD beam diagnostics of the run (linac_sec1/diags/main). Writes five
figures to linac_sec1/results/:

  1. linac_field.png         — on-axis traveling-wave |Ez|(z) envelope (× scale) and a
                               fixed-t snapshot of Ez(z,t): the accelerating structure.
  2. energy_gain.png         — ⟨KE⟩ and max KE vs ⟨z⟩ (~220 keV → ~26 MeV captured mean,
                               ~32 MeV max at the default 11 MW) with β → 1.
  3. long_phase_space.png    — (z − ⟨z⟩) vs KE at injection / mid / exit: RF capture.
  4. beam_envelope.png       — σ_x and surviving charge vs ⟨z⟩: focusing + adiabatic damping.
  5. exit_spectrum_capture.png — exit energy spectrum and the captured-charge fraction.

Run with:
    conda run -n CBB python -c "import linac_sec1; linac_sec1.plot()"
"""

import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import openpmd_api as io
from openpmd_viewer import OpenPMDTimeSeries

from .build_linac_sec1_field import Z_STRUCT, RMAX, BORE_R   # geometry, kept in sync
from . import DEFAULT_OUTDIR                                  # default diags dir for run()

MC2 = 0.51099895                 # electron rest energy [MeV]
Q_E = 1.602176634e-19
RF_NORM_MW = 0.001
POWER_MW = 11.0                  # config(POWER_MW=...) updates this too (mirrors the sim default)
L_STRUCT = 3.016                 # structure length [m]

RF1 = "linac_sec1/linac_sec1_field/linac_rf1.h5"
RF2 = "linac_sec1/linac_sec1_field/linac_rf2.h5"
OUTDIR = None                    # config(OUTDIR=...) sets this; None → DEFAULT_OUTDIR
RESULTS = "linac_sec1/results"


def gamma_of(ux, uy, uz):
    return np.sqrt(1.0 + ux**2 + uy**2 + uz**2)


def wstat(a, w):
    """Weighted mean and standard deviation."""
    m = np.average(a, weights=w)
    return m, np.sqrt(np.average((a - m) ** 2, weights=w))


def on_axis_ez(path):
    """Return (z [m], Ez on axis [V/m]) of an RF quadrature map."""
    s = io.Series(path, io.Access.read_only)
    E = s.iterations[0].meshes["E"]
    ez = E["z"].load_chunk()
    s.flush()
    ez = ez[0][0]                                # mode 0, r = 0 row -> (nz,)
    dz, off = E.grid_spacing[1], E.grid_global_offset[1]
    z = off + np.arange(ez.size) * dz
    del s
    return z, ez


def beam_track(diag):
    """Per-snapshot beam metrics for the run directory; None if unreadable/empty."""
    pdir = os.path.join(diag, "particles")
    if not os.path.isdir(pdir):
        return None
    ts = OpenPMDTimeSeries(pdir)
    its = list(ts.iterations)
    if not its:
        return None
    rec = dict(z=[], ke=[], kemax=[], beta=[], sigx=[], q=[])
    snaps = {}                                    # cache raw (z, ke, w) for the figures
    q_entered = None                              # charge in the FIRST dump (already post-scrape)
    for it in its:
        x, y, z, ux, uy, uz, w = ts.get_particle(
            ["x", "y", "z", "ux", "uy", "uz", "w"], species="electrons", iteration=it)
        if q_entered is None:
            q_entered = w.sum()
        if len(z) < 5:
            continue
        g = gamma_of(ux, uy, uz)
        ke = (g - 1.0) * MC2
        zm, _ = wstat(z, w)
        km, _ = wstat(ke, w)
        rec["z"].append(zm); rec["ke"].append(km); rec["kemax"].append(ke.max())
        rec["beta"].append(np.average(uz / g, weights=w))
        rec["sigx"].append(np.sqrt(np.average(x**2, weights=w)))
        rec["q"].append(w.sum())
        snaps[it] = (z, ke, w)
    for k in rec:
        rec[k] = np.asarray(rec[k])
    if not rec["z"].size:                         # every snapshot near-empty -> no usable beam
        return None
    if not q_entered:                             # degenerate injection (zero baseline charge)
        return None
    # The TRUE injected charge is recorded by the sim (injection_summary.json), because WarpX
    # drops r>RMAX particles before the first dump — so q_entered already hides the injection
    # loss. Report capture against the injected charge; fall back to q_entered if the sidecar
    # is missing (e.g. an old run), in which case the injection loss is simply not shown.
    summ_path = os.path.join(diag, "injection_summary.json")
    inj = None
    if os.path.isfile(summ_path):
        with open(summ_path) as fh:
            inj = json.load(fh)
    q_inj = (inj["q_injected_C"] / Q_E) if inj else q_entered   # in macroparticle-weight units
    rec["q_entered"] = q_entered                  # entered the domain (first dump)
    rec["q0"] = q_inj                             # capture denominator = true injected charge
    rec["inj"] = inj
    rec["snaps"] = snaps
    return rec


def main():
    os.makedirs(RESULTS, exist_ok=True)
    scale = float(np.sqrt(POWER_MW / RF_NORM_MW))

    # ══ Fig 1: applied traveling-wave field ════════════════════════════════════
    if not (os.path.exists(RF1) and os.path.exists(RF2)):
        print(f"no RF maps in {os.path.dirname(RF1)}; run build first. Skipping field figure.")
    else:
        z, ez1 = on_axis_ez(RF1)
        _, ez2 = on_axis_ez(RF2)
        amp = np.sqrt(ez1**2 + ez2**2)                         # traveling-wave amplitude
        env = amp * scale
        snap = (ez1 * np.cos(0.0) - ez2 * np.sin(0.0)) * scale  # Ez(z, t) at one instant
        vgain = np.trapezoid(amp, z) * scale
        fig, (a1, a2) = plt.subplots(2, 1, figsize=(9.2, 6.2), constrained_layout=True, sharex=True)
        a1.plot(z, env / 1e6, color="C3")
        a1.fill_between(z, env / 1e6, alpha=0.12, color="C3")
        a1.set_ylabel(r"$|E_z|$ amplitude  [MV/m]")
        a1.set_title(f"SLAC Section 1 traveling wave at P = {POWER_MW:g} MW "
                     f"(peak {env.max()/1e6:.1f} MV/m, ∫|Ez|dz = {vgain/1e6:.1f} MV)")
        a2.plot(z, snap / 1e6, color="C0", lw=0.7)
        a2.axhline(0, color="k", lw=0.5)
        a2.set_xlabel("z  [m]"); a2.set_ylabel(r"$E_z(z, t_0)$  [MV/m]")
        a2.set_title("On-axis field snapshot (2π/3 traveling-wave structure)")
        a2.set_xlim(Z_STRUCT, Z_STRUCT + 0.4)                   # zoom to show the cell structure
        fig.savefig(f"{RESULTS}/linac_field.png", dpi=140)
        print(f"wrote {RESULTS}/linac_field.png")

    main_diag = OUTDIR or DEFAULT_OUTDIR           # honour config(OUTDIR=...) overrides
    rec = beam_track(main_diag)
    if rec is None:
        print(f"no beam diagnostics in {main_diag}; run the sim first. Skipping beam figures.")
        return

    # ══ Fig 2: energy gain — KE with both the Lorentz factor γ and β ═══════════
    zmm = rec["z"] * 1e3
    gamma = 1.0 + rec["ke"] / MC2                       # γ = 1 + KE/mc²
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    fig.subplots_adjust(left=0.08, right=0.79, bottom=0.13, top=0.91)
    h_struct = ax.axvspan(Z_STRUCT * 1e3, (Z_STRUCT + L_STRUCT) * 1e3, color="0.92",
                          zorder=0, label="structure")
    hmean, = ax.plot(zmm, rec["ke"], "o-", color="C2", ms=3, label="mean KE")
    hmax, = ax.plot(zmm, rec["kemax"], "^--", color="C1", ms=3, label="max KE")
    ax.set_xlabel("mean beam position  ⟨z⟩  [mm]")
    ax.set_ylabel("kinetic energy  [MeV]")
    ax.set_title("Beam energy gain through SLAC Section 1")
    # γ (energy, → ~70) on the inner right axis; β (velocity, → 1) on an offset right axis.
    axg = ax.twinx()
    hg, = axg.plot(zmm, gamma, "-.", color="C5", lw=1.6, label=r"$\gamma$ (Lorentz factor)")
    axg.set_ylabel(r"Lorentz factor  $\gamma$", color="C5")
    axg.tick_params(axis="y", labelcolor="C5"); axg.set_ylim(0, gamma.max() * 1.08)
    axb = ax.twinx()
    axb.spines["right"].set_position(("axes", 1.14))
    hb, = axb.plot(zmm, rec["beta"], ":", color="C4", lw=1.8, label=r"$\beta = v/c$")
    axb.set_ylabel(r"$\beta = v/c$", color="C4")
    axb.tick_params(axis="y", labelcolor="C4"); axb.set_ylim(0.5, 1.02)
    ax.legend(handles=[hmean, hmax, hg, hb, h_struct], loc="center right", fontsize=8)
    fig.savefig(f"{RESULTS}/energy_gain.png", dpi=140)
    print(f"wrote {RESULTS}/energy_gain.png")

    # ══ Fig 3: longitudinal phase space at injection / mid / exit ═══════════════
    snaps = rec["snaps"]
    its = list(snaps)
    # Pick the mid panel by beam position (⟨z⟩ nearest the capture region, ≈Z_STRUCT
    # + 0.2 m), not the middle iteration index — the latter lands ~1.65 m downstream,
    # well past where the RF bucket forms.
    zmeans = {it: np.average(snaps[it][0], weights=snaps[it][2]) for it in its}
    mid = min(its, key=lambda it: abs(zmeans[it] - (Z_STRUCT + 0.2)))
    picks = [its[0], mid, its[-1]]
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True, squeeze=False)
    for ax, it in zip(axs[0], picks):
        z, ke, w = snaps[it]                          # cached from beam_track
        zm = np.average(z, weights=w)
        ax.scatter((z - zm) * 1e3, ke, s=2, alpha=0.25, color="C0")
        ax.set_xlabel(r"$z - \langle z\rangle$  [mm]"); ax.set_ylabel("KE  [MeV]")
        ax.set_title(f"⟨z⟩ = {zm*1e3:.0f} mm   (N = {len(z)})")
    fig.suptitle("Longitudinal phase space: capture into the RF bucket", fontsize=12)
    fig.savefig(f"{RESULTS}/long_phase_space.png", dpi=140)
    print(f"wrote {RESULTS}/long_phase_space.png")

    # ══ Fig 4: transverse envelope + survival ═══════════════════════════════════
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.2, 6.4), constrained_layout=True, sharex=True)
    a1.plot(rec["z"] * 1e3, rec["sigx"] * 1e3, "o-", color="C0", ms=3)
    a1.axhline(BORE_R * 1e3, color="k", ls=":", lw=1, label="structure bore")
    a1.axhline(RMAX * 1e3, color="0.5", ls=":", lw=1, label="domain wall")
    a1.set_ylabel(r"RMS size  $\sigma_x$  [mm]")
    a1.set_title("Transverse envelope and beam survival")
    a1.legend(loc="upper right", fontsize=8)
    # Normalised to the TRUE injected charge (q0). The first tracked dump already sits below 1
    # because WarpX scrapes the r>RMAX particles at injection; prepend the injection point
    # (q/q0 = 1 at ⟨z⟩_inject) so that step-0 radial loss is visible rather than hidden.
    # z_inject_mean_m is the FULL injected beam's ⟨z⟩, so the marker leads the first dump by the
    # scraped (large-r) population's z-offset (a few mm — negligible on the ~3.5 m axis).
    qfrac = rec["q"] / rec["q0"]
    zmm_q = rec["z"] * 1e3
    if rec.get("inj"):
        z_inj_mm = rec["inj"]["z_inject_mean_m"] * 1e3
        zmm_q = np.concatenate([[z_inj_mm], zmm_q])
        qfrac = np.concatenate([[1.0], qfrac])
        a2.annotate("injection scraping\n(r > domain wall)",
                    xy=(rec["z"][0] * 1e3, rec["q"][0] / rec["q0"]),
                    xytext=(0.30, 0.55), textcoords="axes fraction", fontsize=8,
                    arrowprops=dict(arrowstyle="->", color="C3", lw=1), color="C3")
    a2.plot(zmm_q, qfrac, "o-", color="C0", ms=3)
    a2.set_xlabel("mean beam position  ⟨z⟩  [mm]")
    a2.set_ylabel("surviving charge  q / q$_{inj}$")
    a2.set_ylim(-0.03, 1.05)
    fig.savefig(f"{RESULTS}/beam_envelope.png", dpi=140)
    print(f"wrote {RESULTS}/beam_envelope.png")

    # ══ Fig 5: exit energy spectrum + capture fraction ══════════════════════════
    z, ke, w = snaps[its[-1]]                         # cached exit snapshot
    km, sk = wstat(ke, w)
    cap = rec["q"][-1] / rec["q0"]                    # captured / TRUE injected
    q_cap_pC = rec["q"][-1] * Q_E * 1e12
    q_inj_pC = rec["q0"] * Q_E * 1e12
    fig, ax = plt.subplots(figsize=(7.8, 4.8), constrained_layout=True)
    cnt, edges, _ = ax.hist(ke, bins=60, weights=w * Q_E * 1e12, color="C3", alpha=0.85)
    ax.axvline(km, color="k", ls="--", label=f"⟨KE⟩ = {km:.1f} ± {sk:.1f} MeV")
    ax.set_xlabel("KE  [MeV]"); ax.set_ylabel("charge per bin  [pC]")
    ax.set_title(f"Exit energy spectrum — captured {q_cap_pC:.1f} pC "
                 f"= {cap*100:.1f}% of {q_inj_pC:.1f} pC injected")
    ax.legend(loc="upper left")
    # Make the injection loss explicit: how much charge ever entered the domain, and the
    # capture fraction relative to that in-domain charge (so both denominators are visible).
    # Use the sidecar's exact step-0 in-domain charge (q_in_domain_C) rather than the first-dump
    # charge: the two are equal only if nothing scrapes between step 0 and the first dump, so the
    # sidecar value is the correct baseline and keeps "scraped at injection" from absorbing any
    # early-transit loss.
    if rec.get("inj"):
        q_dom = rec["inj"]["q_in_domain_C"] / Q_E       # weight units; exact step-0 baseline
        q_dom_pC = q_dom * Q_E * 1e12
        ax.text(0.985, 0.97,
                f"{q_dom_pC:.0f} pC entered the {RMAX*1e3:.0f} mm domain\n"
                f"({rec['q'][-1]/q_dom*100:.0f}% of those captured;\n"
                f"{(1-q_dom/rec['q0'])*100:.0f}% scraped at injection)",
                transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color="0.3")
    # Inset: zoom into the low-energy tail (phase-slipped / off-crest captured particles),
    # which the dominant captured-energy peak otherwise hides. Same bins; y-axis fit to the tail.
    centers = 0.5 * (edges[:-1] + edges[1:])
    cut = 0.85 * km
    tail = cnt[centers < cut]
    if tail.size and tail.max() > 0:
        axin = ax.inset_axes([0.30, 0.36, 0.50, 0.56])
        axin.hist(ke, bins=edges, weights=w * Q_E * 1e12, color="C3", alpha=0.85)
        axin.set_xlim(max(0.0, ke.min() - 1.0), cut)
        axin.set_ylim(0, tail.max() * 1.35)
        axin.set_title("low-energy tail (zoom)", fontsize=8)
        axin.set_xlabel("KE  [MeV]", fontsize=7); axin.set_ylabel("pC/bin", fontsize=7)
        axin.tick_params(labelsize=7)
    fig.savefig(f"{RESULTS}/exit_spectrum_capture.png", dpi=140)
    print(f"wrote {RESULTS}/exit_spectrum_capture.png")


if __name__ == "__main__":
    main()
