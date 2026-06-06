"""
Figures for the SLAC Linac Section 1 stage (linac_sec1_sim.py).

Reads the traveling-wave field maps (linac_sec1/linac_sec1_field/linac_rf{1,2}.h5)
and the openPMD beam diagnostics of the run (linac_sec1/diags/main). Writes five
figures to linac_sec1/results/:

  1. linac_field.png         — on-axis traveling-wave |Ez|(z) envelope (× scale) and a
                               fixed-t snapshot of Ez(z,t): the accelerating structure.
  2. energy_gain.png         — ⟨KE⟩ and max KE vs ⟨z⟩ (148 keV → ~37 MeV) with β → 1.
  3. long_phase_space.png    — (z − ⟨z⟩) vs KE at injection / mid / exit: RF capture.
  4. beam_envelope.png       — σ_r and surviving charge vs ⟨z⟩: focusing + adiabatic damping.
  5. exit_spectrum_capture.png — exit energy spectrum and the captured-charge fraction.

Run with:
    conda run -n CBB python -c "import linac_sec1; linac_sec1.plot()"
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import openpmd_api as io
from openpmd_viewer import OpenPMDTimeSeries

from .build_linac_sec1_field import Z_STRUCT, RMAX            # geometry, kept in sync

MC2 = 0.51099895                 # electron rest energy [MeV]
Q_E = 1.602176634e-19
RF_NORM_MW = 0.001
POWER_MW = 15.0                  # config(POWER_MW=...) updates this too (mirrors the sim)
BORE_R = 0.00955                 # structure bore radius [m] (SLAC map r-extent)
L_STRUCT = 3.016                 # structure length [m]

RF1 = "linac_sec1/linac_sec1_field/linac_rf1.h5"
RF2 = "linac_sec1/linac_sec1_field/linac_rf2.h5"
MAIN = "linac_sec1/diags/main"
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

    rec = beam_track(MAIN)
    if rec is None:
        print(f"no beam diagnostics in {MAIN}; run the sim first. Skipping beam figures.")
        return

    # ══ Fig 2: energy gain + β ══════════════════════════════════════════════════
    zmm = rec["z"] * 1e3
    fig, ax = plt.subplots(figsize=(7.8, 4.8), constrained_layout=True)
    ax.plot(zmm, rec["ke"], "o-", color="C2", ms=3, label="mean KE")
    ax.plot(zmm, rec["kemax"], "^--", color="C1", ms=3, label="max KE")
    ax.axvspan(Z_STRUCT * 1e3, (Z_STRUCT + L_STRUCT) * 1e3, color="0.9", zorder=0,
               label="structure")
    ax.set_xlabel("mean beam position  ⟨z⟩  [mm]")
    ax.set_ylabel("kinetic energy  [MeV]")
    ax.set_title("Beam energy gain through SLAC Section 1")
    axb = ax.twinx()
    axb.plot(zmm, rec["beta"], ":", color="C4", label=r"$\beta$")
    axb.set_ylabel(r"$\beta = v/c$", color="C4"); axb.tick_params(axis="y", labelcolor="C4")
    axb.set_ylim(0.5, 1.02)
    ax.legend(loc="lower right")
    fig.savefig(f"{RESULTS}/energy_gain.png", dpi=140)
    print(f"wrote {RESULTS}/energy_gain.png")

    # ══ Fig 3: longitudinal phase space at injection / mid / exit ═══════════════
    ts = rec["ts"]
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

    # ══ Fig 4: transverse envelope + survival ═══════════════════════════════════
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.2, 6.4), constrained_layout=True, sharex=True)
    a1.plot(rec["z"] * 1e3, rec["sigr"] * 1e3, "o-", color="C0", ms=3)
    a1.axhline(BORE_R * 1e3, color="k", ls=":", lw=1, label="structure bore")
    a1.axhline(RMAX * 1e3, color="0.5", ls=":", lw=1, label="domain wall")
    a1.set_ylabel(r"RMS size  $\sigma_x$  [mm]")
    a1.set_title("Transverse envelope and beam survival (solenoid focused, on crest)")
    a1.legend(loc="upper right", fontsize=8)
    a2.plot(rec["z"] * 1e3, rec["q"] / rec["q0"], "o-", color="C0", ms=3)
    a2.set_xlabel("mean beam position  ⟨z⟩  [mm]")
    a2.set_ylabel("surviving charge  q/q₀"); a2.set_ylim(0, 1.05)
    fig.savefig(f"{RESULTS}/beam_envelope.png", dpi=140)
    print(f"wrote {RESULTS}/beam_envelope.png")

    # ══ Fig 5: exit energy spectrum + capture fraction ══════════════════════════
    it_exit = its[-1]
    z, ux, uy, uz, w = ts.get_particle(
        ["z", "ux", "uy", "uz", "w"], species="electrons", iteration=it_exit)
    ke = (gamma_of(ux, uy, uz) - 1.0) * MC2
    km, sk = wstat(ke, w)
    cap = rec["q"][-1] / rec["q0"]
    fig, ax = plt.subplots(figsize=(7.6, 4.6), constrained_layout=True)
    ax.hist(ke, bins=60, weights=w * Q_E * 1e12, color="C3", alpha=0.85)
    ax.axvline(km, color="k", ls="--", label=f"⟨KE⟩ = {km:.1f} ± {sk:.1f} MeV")
    ax.set_xlabel("KE  [MeV]"); ax.set_ylabel("charge per bin  [pC]")
    ax.set_title(f"Exit energy spectrum — captured fraction {cap*100:.0f}% "
                 f"({rec['q'][-1]*Q_E*1e12:.1f} pC of {rec['q0']*Q_E*1e12:.1f} pC)")
    ax.legend()
    fig.savefig(f"{RESULTS}/exit_spectrum_capture.png", dpi=140)
    print(f"wrote {RESULTS}/exit_spectrum_capture.png")


if __name__ == "__main__":
    main()
