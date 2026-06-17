"""
Figures for the SLAC Linac Section 1 stage: reads the RF field maps and the run's
openPMD diagnostics, writes five PNGs to linac_sec1/results/.

See linac_sec1/README.md for physics, parameters, outputs, and gotchas.
"""

import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import openpmd_api as io
from openpmd_viewer import OpenPMDTimeSeries

from .build_linac_sec1_field import Z_STRUCT, RMAX, BORE_R
from . import DEFAULT_OUTDIR

MC2 = 0.51099895                 # electron rest energy [MeV]
Q_E = 1.602176634e-19
RF_NORM_MW = 0.001
POWER_MW = 11.0                  # config()-overridable; mirrors the sim default
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
    ez = ez[0][0]                                # mode 0, r = 0 row -> (nz,)  [thetaMode axis order]
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
    snaps = {}
    q_entered = None                              # charge in the FIRST dump (already post-scrape)
    for it in its:
        x, y, z, ux, uy, uz, w = ts.get_particle(   # species plural "electrons"
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
        rec["sigx"].append(wstat(x, w)[1])        # centered weighted RMS
        rec["q"].append(w.sum())
        snaps[it] = (z, ke, w)
    for k in rec:
        rec[k] = np.asarray(rec[k])
    if not rec["z"].size:
        return None
    if not q_entered:
        return None
    # Capture denominator = TRUE injected charge from injection_summary.json; the first dump
    # is already post-collimation so q_entered hides the injection loss. Fall back to
    # q_entered if the sidecar is absent (old run).
    summ_path = os.path.join(diag, "injection_summary.json")
    inj = None
    if os.path.isfile(summ_path):
        with open(summ_path) as fh:
            inj = json.load(fh)
    q_inj = (inj["q_injected_C"] / Q_E) if inj else q_entered   # macroparticle-weight units
    rec["q_entered"] = q_entered
    rec["q0"] = q_inj
    rec["inj"] = inj
    rec["snaps"] = snaps
    return rec


def main():
    os.makedirs(RESULTS, exist_ok=True)
    scale = float(np.sqrt(POWER_MW / RF_NORM_MW))

    # Fig 1: applied traveling-wave field
    if not (os.path.exists(RF1) and os.path.exists(RF2)):
        print(f"no RF maps in {os.path.dirname(RF1)}; run build first. Skipping field figure.")
    else:
        z, ez1 = on_axis_ez(RF1)
        _, ez2 = on_axis_ez(RF2)
        amp = np.sqrt(ez1**2 + ez2**2)                         # traveling-wave amplitude
        env = amp * scale
        snap = (ez1 * np.cos(0.0) - ez2 * np.sin(0.0)) * scale  # Ez(z, t0)
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
        a2.set_xlim(Z_STRUCT, Z_STRUCT + 0.4)                   # zoom to the cell structure
        fig.savefig(f"{RESULTS}/linac_field.png", dpi=140)
        print(f"wrote {RESULTS}/linac_field.png")

    main_diag = OUTDIR or DEFAULT_OUTDIR
    rec = beam_track(main_diag)
    if rec is None:
        print(f"no beam diagnostics in {main_diag}; run the sim first. Skipping beam figures.")
        return

    # Fig 2: energy gain — KE with the Lorentz factor γ and β
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

    # Fig 3: longitudinal phase space at injection / mid / exit
    snaps = rec["snaps"]
    its = list(snaps)
    # Mid panel = ⟨z⟩ nearest the capture region (≈Z_STRUCT+0.2 m), NOT the middle
    # iteration index (which lands well past where the RF bucket forms).
    zmeans = {it: np.average(snaps[it][0], weights=snaps[it][2]) for it in its}
    mid = min(its, key=lambda it: abs(zmeans[it] - (Z_STRUCT + 0.2)))
    if mid == its[0] and len(its) > 1:            # keep mid ≠ injection
        mid = its[1]
    picks = [its[0], mid, its[-1]]
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True, squeeze=False)
    for ax, it in zip(axs[0], picks):
        z, ke, w = snaps[it]
        zm = np.average(z, weights=w)
        ax.scatter((z - zm) * 1e3, ke, s=2, alpha=0.25, color="C0")
        ax.set_xlabel(r"$z - \langle z\rangle$  [mm]"); ax.set_ylabel("KE  [MeV]")
        ax.set_title(f"⟨z⟩ = {zm*1e3:.0f} mm   (N = {len(z)})")
    fig.suptitle("Longitudinal phase space: capture into the RF bucket", fontsize=12)
    fig.savefig(f"{RESULTS}/long_phase_space.png", dpi=140)
    print(f"wrote {RESULTS}/long_phase_space.png")

    # Fig 4: transverse envelope + survival
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.2, 6.4), constrained_layout=True, sharex=True)
    a1.plot(rec["z"] * 1e3, rec["sigx"] * 1e3, "o-", color="C0", ms=3)
    a1.axhline(BORE_R * 1e3, color="k", ls=":", lw=1, label="structure bore")
    a1.axhline(RMAX * 1e3, color="0.5", ls=":", lw=1, label="domain wall")
    a1.set_ylabel(r"RMS size  $\sigma_x$  [mm]")
    a1.set_title("Transverse envelope and beam survival")
    a1.legend(loc="upper right", fontsize=8)
    # Normalised to the TRUE injected charge (q0); prepend (q/q0=1 at ⟨z⟩_inject) so the
    # step-0 injection scraping is visible rather than hidden in the first post-scrape dump.
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

    # Fig 5: exit energy spectrum + capture fraction
    z, ke, w = snaps[its[-1]]
    km, sk = wstat(ke, w)
    cap = rec["q"][-1] / rec["q0"]                    # captured / TRUE injected
    q_cap_pC = rec["q"][-1] * Q_E * 1e12
    q_inj_pC = rec["q0"] * Q_E * 1e12
    fig, ax = plt.subplots(figsize=(7.8, 4.8), constrained_layout=True)
    # √N bin count (clamped): the captured core is too few macroparticles for a fixed 60 bins
    # (would render as a spiky comb).
    nbins = int(np.clip(round(np.sqrt(ke.size)), 12, 60))
    cnt, edges, _ = ax.hist(ke, bins=nbins, weights=w * Q_E * 1e12, color="C3", alpha=0.85)
    ax.axvline(km, color="k", ls="--", label=f"⟨KE⟩ = {km:.1f} ± {sk:.1f} MeV")
    ax.set_xlabel("KE  [MeV]"); ax.set_ylabel("charge per bin  [pC]")
    ax.set_title(f"Exit energy spectrum — captured {q_cap_pC:.1f} pC "
                 f"= {cap*100:.1f}% of {q_inj_pC:.1f} pC injected")
    ax.legend(loc="upper left")
    # Use the sidecar's exact step-0 in-domain charge (q_in_domain_C), NOT the first-dump
    # charge: the two differ if anything scrapes between step 0 and the first dump, so the
    # sidecar keeps "scraped at injection" from absorbing early-transit loss.
    if rec.get("inj"):
        q_dom = rec["inj"]["q_in_domain_C"] / Q_E       # weight units; exact step-0 baseline
        q_dom_pC = q_dom * Q_E * 1e12
        ax.text(0.985, 0.97,
                f"{q_dom_pC:.0f} pC entered the {RMAX*1e3:.0f} mm domain\n"
                f"({rec['q'][-1]/q_dom*100:.0f}% of those captured;\n"
                f"{(1-q_dom/rec['q0'])*100:.0f}% scraped at injection)",
                transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color="0.3")
    # Inset: zoom into the low-energy (phase-slipped) tail the captured-energy peak hides.
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
