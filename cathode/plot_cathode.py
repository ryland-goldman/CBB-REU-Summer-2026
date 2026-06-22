"""
Figures for the finite-cathode space-charge-limited (Child–Langmuir) diode (cathode_diode.py)
over cathode/diags/. Writes PNGs to cathode/results/.

Two layers: generic phase-space / field figures via lume-warpx's helpers and the shared
`pipeline.plot_extras` beam figures, plus the stage-specific rich figures that validate the
emission physics — the on-axis potential and field against the Child–Langmuir law, the
transmitted current saturating at J_CL despite 2× over-injection, and the source's intrinsic
thermal transverse phase space. See cathode/README.md for the physics each figure shows.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from openpmd_viewer import OpenPMDTimeSeries

from pipeline.constants import C_LIGHT, MC2_EV
from pipeline.emission import child_langmuir_current_density, thermal_velocity_sigma
from pipeline import plot_extras as px

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG = os.path.join(HERE, "cathode.yaml")
RESULTS = "cathode/results"
FIELDS = "cathode/diags/fields"
PARTICLES = "cathode/diags/particles"

MC2_KEV = MC2_EV / 1e3


def _last_populated(diag, species="electrons"):
    ts = OpenPMDTimeSeries(diag)
    for it in reversed([int(i) for i in ts.iterations]):
        try:
            x, = ts.get_particle(["x"], species=species, iteration=it)
        except Exception:
            continue
        if len(x):
            return it
    return int(ts.iterations[-1])


def _save(fig, name):
    fig.savefig(f"{RESULTS}/{name}.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {RESULTS}/{name}.png")


def child_langmuir_figure(V_anode, gap_d):
    """On-axis potential and field vs the planar Child–Langmuir and vacuum references."""
    ts = OpenPMDTimeSeries(FIELDS)
    it = ts.iterations[-1]                               # steady-state snapshot
    phi, meta = ts.get_field("phi", iteration=it)
    ez, _ = ts.get_field("E", "z", iteration=it)
    ix0 = np.argmin(np.abs(meta.x))                      # column nearest the axis
    z = meta.z
    phi_axis, ez_axis = phi[:, ix0], ez[:, ix0]

    s = z / gap_d
    phi_cl, ez_cl = V_anode * s ** (4 / 3), -(4 / 3) * (V_anode / gap_d) * s ** (1 / 3)
    phi_vac, ez_vac = V_anode * s, -(V_anode / gap_d) * np.ones_like(s)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    a1.plot(z * 1e3, phi_axis, "o-", color="C0", ms=3, label="WarpX (on axis)")
    a1.plot(z * 1e3, phi_cl, "k--", label=r"Child–Langmuir  $V(z/d)^{4/3}$")
    a1.plot(z * 1e3, phi_vac, ":", color="gray", label="vacuum (no space charge)")
    a1.set_xlabel("z  [mm]"); a1.set_ylabel("φ  [V]")
    a1.set_title("On-axis potential"); a1.legend()

    a2.plot(z * 1e3, ez_axis / 1e3, "o-", color="C3", ms=3, label="WarpX (on axis)")
    a2.plot(z * 1e3, ez_cl / 1e3, "k--", label=r"Child–Langmuir  $-\frac{4V}{3d}(z/d)^{1/3}$")
    a2.plot(z * 1e3, ez_vac / 1e3, ":", color="gray", label="vacuum")
    a2.set_xlabel("z  [mm]"); a2.set_ylabel(r"$E_z$  [kV/m]")
    a2.set_title("On-axis longitudinal field"); a2.legend()
    fig.suptitle("Space-charge depression of the field follows the Child–Langmuir law", fontsize=12)
    _save(fig, "child_langmuir")


def current_saturation_figure(V_anode, gap_d, R_cathode, J_CL, over_inject):
    """Transmitted current density at the anode saturating at J_CL despite over-injection."""
    ts = OpenPMDTimeSeries(FIELDS)
    dx = None
    times, j_trans = [], []
    for i, it in enumerate(ts.iterations):
        jz, meta = ts.get_field("j", "z", iteration=it)
        if dx is None:
            dx = meta.x[1] - meta.x[0]
        line_current = np.abs(jz[-2, :].sum() * dx)      # ∫|jz| dx just inside the anode [A/m depth]
        j_trans.append(line_current / (2.0 * R_cathode)) # referenced to the cathode width
        times.append(ts.t[i] * 1e9)                      # ns

    fig, ax = plt.subplots(figsize=(7, 4.8), constrained_layout=True)
    ax.plot(times, j_trans, "o-", color="C2", label="WarpX transmitted current")
    ax.axhline(J_CL, color="k", ls="--", label=r"Child–Langmuir limit $J_{CL}$")
    ax.axhline(over_inject * J_CL, color="r", ls=":",
               label=f"injected ({over_inject:.0f}× $J_{{CL}}$)")
    ax.set_xlabel("time  [ns]"); ax.set_ylabel(r"current density at anode  $|J_z|$  [A/m²]")
    ax.set_title("Space charge limits the transmitted current toward $J_{CL}$")
    ax.set_xlim(0, 0.15); ax.set_ylim(0, over_inject * J_CL * 1.1); ax.legend()
    _save(fig, "current_saturation")


def emission_phase_space_figure(T_cathode):
    """Intrinsic thermal transverse phase space and emittance of the source beam."""
    ts = OpenPMDTimeSeries(PARTICLES)
    x, ux, w = ts.get_particle(["x", "ux", "w"], iteration=ts.iterations[-1])
    xbar, uxbar = np.average(x, weights=w), np.average(ux, weights=w)
    x2 = np.average((x - xbar) ** 2, weights=w)
    ux2 = np.average((ux - uxbar) ** 2, weights=w)
    xux = np.average((x - xbar) * (ux - uxbar), weights=w)
    emit_mm_mrad = np.sqrt(max(x2 * ux2 - xux ** 2, 0.0)) * 1e6   # m·(γβ) → mm·mrad

    px_kev = ux * MC2_KEV                                 # p_x·c = γβ_x·m_ec² [keV]
    u_thermal = thermal_velocity_sigma(T_cathode) / C_LIGHT      # dimensionless γβ rms = √(kT/mₑc²)
    p_thermal = u_thermal * MC2_KEV                              # √(kT·mₑc²) [keV/c]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    hb = a1.hexbin(x * 1e3, px_kev, gridsize=70, cmap="viridis", mincnt=1)
    fig.colorbar(hb, ax=a1, label="macroparticles / bin")
    a1.set_xlabel("x  [mm]"); a1.set_ylabel(r"$p_x$  [keV/$c$]")
    a1.set_title("Transverse phase space at the cathode")
    a1.text(0.03, 0.97,
            rf"$\varepsilon_{{n,x}} = {emit_mm_mrad:.3f}$ mm·mrad"
            "\n" rf"$\sqrt{{\langle x^2\rangle}} = {np.sqrt(x2)*1e3:.2f}$ mm"
            "\n" rf"gun disc receives $\times\sqrt{{3/4}}\approx{emit_mm_mrad*np.sqrt(0.75):.2f}$ mm·mrad",
            transform=a1.transAxes, va="top", fontsize=9,
            bbox=dict(boxstyle="round", fc="white", alpha=0.85))

    a2.hist(px_kev, bins=120, color="C0", alpha=0.8, density=True)
    for sign, lbl in ((+1, r"$\pm\sqrt{kT\,m_ec^2}$"), (-1, None)):
        a2.axvline(sign * p_thermal, color="k", ls="--", lw=1.2, label=lbl)
    a2.set_xlabel(r"$p_x$  [keV/$c$]"); a2.set_ylabel("probability density")
    a2.set_title(f"Thermal transverse momentum spread ({T_cathode:.0f} K)")
    a2.legend(loc="upper right", fontsize=9)
    fig.suptitle("Intrinsic thermal emittance of the cathode — the source quality the gun inherits",
                 fontsize=12)
    _save(fig, "emission_phase_space")


def main():
    from warpx import WarpX
    os.makedirs(RESULTS, exist_ok=True)

    w = WarpX(input_file=CONFIG, path="cathode")
    w.load_output()                                      # cathode/diags/{fields, particles}
    it = _last_populated(PARTICLES)
    pg = w._particle_group(iteration=it)

    # Generic phase-space / field figures (lume-warpx helpers + shared plot_extras).
    for name, fig in [
        ("phase_space_z_KE", w.plot2D("z", "kinetic_energy", iteration=it)),
        ("transverse_x_px",  w.plot2D("x", "px", iteration=it)),
        ("potential_xz",     w.plot_fields("phi", "x", "z")),
        ("charge_density_xz", w.plot_fields("rho", "x", "z")),
        ("centroid_vs_t",    w.plot1D("t", "mean_z")),
        ("charge_vs_t",      w.plot1D("t", "charge")),
        ("energy_spectrum",  px.energy_spectrum(pg)),
        ("current_profile",  px.current_profile(pg)),
    ]:
        _save(fig, name)

    # Stage-specific rich figures (raw openPMD; emission-physics validation).
    V_anode = w.get("grid/warpx_potential_hi_z")
    gap_d = w.get("grid/upper_bound")[1]
    R_cathode = w.get("species")[0]["upper_bound"][0]
    p = w.get("params")
    J_CL = float(child_langmuir_current_density(V_anode, gap_d))

    child_langmuir_figure(V_anode, gap_d)
    current_saturation_figure(V_anode, gap_d, R_cathode, J_CL, p["over_inject"])
    emission_phase_space_figure(p["T_cathode"])


if __name__ == "__main__":
    main()
