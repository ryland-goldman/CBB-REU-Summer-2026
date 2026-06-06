"""
Figures for the SLAC Linac Section 1 stage (linac_sec1_sim.py).

Reads the traveling-wave field maps (linac_sec1/linac_sec1_field/linac_rf{1,2}.h5)
and the openPMD beam diagnostics of the headline run (linac_sec1/diags/main) plus,
when present, the focus-off comparison (linac_sec1/diags/focusoff) and the RF-phase
scan cases (linac_sec1/diags/scan_phi*). Writes six figures to linac_sec1/results/:

  1. linac_field.png         — on-axis traveling-wave |Ez|(z) envelope (× scale) and a
                               fixed-t snapshot of Ez(z,t): the accelerating structure.
  2. energy_gain.png         — ⟨KE⟩ and max KE vs ⟨z⟩ (148 keV → ~37 MeV) with β → 1.
  3. long_phase_space.png    — (z − ⟨z⟩) vs KE at injection / mid / exit: RF capture.
  4. beam_envelope.png       — σ_r and surviving charge vs ⟨z⟩, focusing ON vs OFF,
                               with the structure bore: why the solenoid is needed.
  5. exit_spectrum_capture.png — exit energy spectrum and the captured-charge fraction.
  6. phase_acceptance.png    — energy gain and capture fraction vs injection RF phase.

Run with:
    conda run -n CBB python -c "import linac_sec1; linac_sec1.plot()"
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

from .build_linac_sec1_field import Z_STRUCT, RMAX            # geometry, kept in sync

MC2 = 0.51099895                 # electron rest energy [MeV]
c = 299792458.0
Q_E = 1.602176634e-19
F_RF = 2856.0e6                  # must match linac_sec1_sim.py
RF_NORM_MW = 0.001
POWER_MW = 15.0                  # config(POWER_MW=...) updates this too (mirrors the sim)
BORE_R = 0.00955                 # structure bore radius [m] (SLAC map r-extent)
L_STRUCT = 3.016                 # structure length [m]

RF1 = "linac_sec1/linac_sec1_field/linac_rf1.h5"
RF2 = "linac_sec1/linac_sec1_field/linac_rf2.h5"
DIAG_ROOT = "linac_sec1/diags"
MAIN = f"{DIAG_ROOT}/main"
FOCUSOFF = f"{DIAG_ROOT}/focusoff"
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
    """Per-snapshot beam metrics for one case directory; None if unreadable/empty."""
    pdir = os.path.join(diag, "particles")
    if not os.path.isdir(pdir):
        return None
    ts = OpenPMDTimeSeries(pdir)
    its = list(ts.iterations)
    if not its:
        return None
    rec = dict(z=[], ke=[], kemax=[], ske=[], beta=[], sigr=[], sigz=[], q=[], n=[])
    q0 = None
    for it in its:
        x, y, z, ux, uy, uz, w = ts.get_particle(
            ["x", "y", "z", "ux", "uy", "uz", "w"], species="electrons", iteration=it)
        if q0 is None:
            q0 = w.sum()
        if len(z) < 5:
            continue
        g = gamma_of(ux, uy, uz)
        ke = (g - 1.0) * MC2
        zm, sz = wstat(z, w)
        km, sk = wstat(ke, w)
        rec["z"].append(zm); rec["ke"].append(km); rec["kemax"].append(ke.max())
        rec["ske"].append(sk); rec["beta"].append(np.average(uz / g, weights=w))
        rec["sigr"].append(np.sqrt(np.average(x**2, weights=w)))
        rec["sigz"].append(sz); rec["q"].append(w.sum()); rec["n"].append(len(z))
    for k in rec:
        rec[k] = np.asarray(rec[k])
    rec["q0"] = q0
    rec["ts"] = ts
    return rec


def main():
    os.makedirs(RESULTS, exist_ok=True)
    scale = float(np.sqrt(POWER_MW / RF_NORM_MW))

    # ══ Fig 1: applied traveling-wave field ════════════════════════════════════
    z, ez1 = on_axis_ez(RF1)
    _, ez2 = on_axis_ez(RF2)
    env = np.sqrt(ez1**2 + ez2**2) * scale                 # traveling-wave amplitude
    snap = (ez1 * np.cos(0.0) - ez2 * np.sin(0.0)) * scale  # Ez(z, t) at one instant
    vgain = np.trapezoid(np.sqrt(ez1**2 + ez2**2), z) * scale
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

    main_rec = beam_track(MAIN)
    if main_rec is None:
        print(f"no beam diagnostics in {MAIN}; run the sim first. Skipping beam figures.")
        return

    # ══ Fig 2: energy gain + β ══════════════════════════════════════════════════
    zmm = main_rec["z"] * 1e3
    fig, ax = plt.subplots(figsize=(7.8, 4.8), constrained_layout=True)
    ax.plot(zmm, main_rec["ke"], "o-", color="C2", ms=3, label="mean KE")
    ax.plot(zmm, main_rec["kemax"], "^--", color="C1", ms=3, label="max KE")
    ax.axvspan(Z_STRUCT * 1e3, (Z_STRUCT + L_STRUCT) * 1e3, color="0.9", zorder=0,
               label="structure")
    ax.set_xlabel("mean beam position  ⟨z⟩  [mm]")
    ax.set_ylabel("kinetic energy  [MeV]")
    ax.set_title("Beam energy gain through SLAC Section 1")
    axb = ax.twinx()
    axb.plot(zmm, main_rec["beta"], ":", color="C4", label=r"$\beta$")
    axb.set_ylabel(r"$\beta = v/c$", color="C4"); axb.tick_params(axis="y", labelcolor="C4")
    axb.set_ylim(0.5, 1.02)
    ax.legend(loc="lower right")
    fig.savefig(f"{RESULTS}/energy_gain.png", dpi=140)
    print(f"wrote {RESULTS}/energy_gain.png")

    # ══ Fig 3: longitudinal phase space at injection / mid / exit ═══════════════
    ts = main_rec["ts"]
    its = list(ts.iterations)
    picks = [its[0], its[len(its) // 2], its[-1]]
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True, squeeze=False)
    for ax, it in zip(axs[0], picks):
        z, ux, uy, uz, w = ts.get_particle(
            ["z", "ux", "uy", "uz", "w"], species="electrons", iteration=it)
        ke = (gamma_of(ux, uy, uz) - 1.0) * MC2
        zm = np.average(z, weights=w)
        ax.scatter((z - zm) * 1e3, ke, s=2, alpha=0.25, color="C0")
        ax.set_xlabel(r"$z - \langle z\rangle$  [mm]"); ax.set_ylabel("KE  [MeV]")
        ax.set_title(f"⟨z⟩ = {zm*1e3:.0f} mm   (N = {len(z)})")
    fig.suptitle("Longitudinal phase space: capture into the RF bucket", fontsize=12)
    fig.savefig(f"{RESULTS}/long_phase_space.png", dpi=140)
    print(f"wrote {RESULTS}/long_phase_space.png")

    # ══ Fig 4: transverse envelope + survival, focus ON vs OFF ══════════════════
    off_rec = beam_track(FOCUSOFF)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.2, 6.4), constrained_layout=True, sharex=True)
    a1.plot(main_rec["z"] * 1e3, main_rec["sigr"] * 1e3, "o-", color="C0", ms=3,
            label="focus ON")
    a2.plot(main_rec["z"] * 1e3, main_rec["q"] / main_rec["q0"], "o-", color="C0", ms=3,
            label="focus ON")
    if off_rec is not None:
        a1.plot(off_rec["z"] * 1e3, off_rec["sigr"] * 1e3, "s--", color="C3", ms=3,
                label="focus OFF (I=0)")
        a2.plot(off_rec["z"] * 1e3, off_rec["q"] / off_rec["q0"], "s--", color="C3", ms=3,
                label="focus OFF (I=0)")
    a1.axhline(BORE_R * 1e3, color="k", ls=":", lw=1, label="structure bore")
    a1.axhline(RMAX * 1e3, color="0.5", ls=":", lw=1, label="domain wall")
    a1.set_ylabel(r"RMS size  $\sigma_x$  [mm]")
    a1.set_title("Transverse envelope and beam survival (solenoid focusing on/off)")
    a1.legend(loc="upper right", fontsize=8)
    a2.set_xlabel("mean beam position  ⟨z⟩  [mm]")
    a2.set_ylabel("surviving charge  q/q₀"); a2.set_ylim(0, 1.05)
    a2.legend(loc="upper right", fontsize=8)
    fig.savefig(f"{RESULTS}/beam_envelope.png", dpi=140)
    print(f"wrote {RESULTS}/beam_envelope.png")

    # ══ Fig 5: exit energy spectrum + capture fraction ══════════════════════════
    it_exit = its[-1]
    z, ux, uy, uz, w = ts.get_particle(
        ["z", "ux", "uy", "uz", "w"], species="electrons", iteration=it_exit)
    ke = (gamma_of(ux, uy, uz) - 1.0) * MC2
    km, sk = wstat(ke, w)
    cap = main_rec["q"][-1] / main_rec["q0"]
    fig, ax = plt.subplots(figsize=(7.6, 4.6), constrained_layout=True)
    ax.hist(ke, bins=60, weights=w * Q_E * 1e12, color="C3", alpha=0.85)
    ax.axvline(km, color="k", ls="--", label=f"⟨KE⟩ = {km:.1f} ± {sk:.1f} MeV")
    ax.set_xlabel("KE  [MeV]"); ax.set_ylabel("charge per bin  [pC]")
    ax.set_title(f"Exit energy spectrum — captured fraction {cap*100:.0f}% "
                 f"({main_rec['q'][-1]*Q_E*1e12:.1f} pC of {main_rec['q0']*Q_E*1e12:.1f} pC)")
    ax.legend()
    fig.savefig(f"{RESULTS}/exit_spectrum_capture.png", dpi=140)
    print(f"wrote {RESULTS}/exit_spectrum_capture.png")

    # ══ Fig 6: RF-phase acceptance scan ═════════════════════════════════════════
    scan_dirs = sorted(glob.glob(f"{DIAG_ROOT}/scan_phi*"))
    pts = []
    for d in scan_dirs:
        m = re.search(r"scan_phi(-?\d+)", os.path.basename(d))
        rec = beam_track(d) if m else None
        if rec is None or rec["ke"].size == 0:
            continue
        pts.append((float(m.group(1)), rec["ke"][-1], rec["q"][-1] / rec["q0"]))
    if len(pts) >= 3:
        pts.sort()
        phi = np.array([p[0] for p in pts])
        kef = np.array([p[1] for p in pts])
        capf = np.array([p[2] for p in pts])
        fig, ax = plt.subplots(figsize=(7.8, 4.8), constrained_layout=True)
        l1, = ax.plot(phi, kef, "o-", color="C2", ms=4, label="final ⟨KE⟩")
        ax.set_xlabel("injection RF phase offset  [deg]")
        ax.set_ylabel("final mean KE  [MeV]", color="C2")
        ax.tick_params(axis="y", labelcolor="C2")
        ax.set_title("RF-phase acceptance: energy gain and capture vs injection phase")
        axc = ax.twinx()
        l2, = axc.plot(phi, capf * 100, "s--", color="C0", ms=4, label="capture fraction")
        axc.set_ylabel("captured charge  [%]", color="C0")
        axc.tick_params(axis="y", labelcolor="C0")
        ax.legend(handles=[l1, l2], loc="best")
        fig.savefig(f"{RESULTS}/phase_acceptance.png", dpi=140)
        print(f"wrote {RESULTS}/phase_acceptance.png")
    else:
        print(f"phase_acceptance.png skipped (need ≥3 scan_phi* cases, found {len(pts)})")


if __name__ == "__main__":
    main()
