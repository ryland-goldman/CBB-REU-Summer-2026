"""
Figures for the cathode stage (sim/cathode.py) from logs/diags/cathode/. Writes PNGs to
logs/plots/cathode/. main() runs ONLY plotting (the sim must have been run first). See
docs/cathode.md for what each figure shows.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from openpmd_viewer import OpenPMDTimeSeries

from sim.helpers.tools import C_LIGHT, MC2_KEV, child_langmuir_current_density, thermal_velocity_sigma
from sim.helpers.loadparticles import anode_beam_mask
from sim.plot import common

CONFIG = "config/cathode.yaml"
RESULTS = "logs/plots/cathode"
DIAG_DIR = "logs/diags/cathode"
FIELDS = os.path.join(DIAG_DIR, "fields")
PARTICLES = os.path.join(DIAG_DIR, "particles")


def _load_series(w):
    """Populate w._outputs with the cathode's two openPMD series (fields, particles).

    sim/cathode.py overrides the diagnostics' write_dir to logs/diags/cathode, so the two
    named series sit directly under DIAG_DIR rather than under <path>/diags — load them by
    hand the way WarpX.load_output() would, so plot2D/plot1D/plot_fields/_particle_group work.
    """
    w._outputs = {}
    for name, sub in (("fields", FIELDS), ("particles", PARTICLES)):
        if os.path.isdir(sub):
            try:
                w._outputs[name] = OpenPMDTimeSeries(sub)
            except Exception:
                continue
    if not w._outputs:
        raise FileNotFoundError(f"No readable openPMD diagnostics found under {DIAG_DIR}")
    w._diag_dir = DIAG_DIR
    w._output = next(iter(w._outputs.values()))
    return w._outputs


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


def _pulse_voltage(t, p):
    """Grid bias V(t) = V_OFF + V_PULSE·tent(t) evaluated at times t [s] (mirrors the parser tent)."""
    t_rise = p["V_PULSE"] / p["V_SLOPE"]
    tent = np.clip(np.minimum(t / t_rise, (2.0 * t_rise - t) / t_rise), 0.0, None)
    return p["V_OFF"] + p["V_PULSE"] * np.clip(tent, None, 1.0)


def _peak_field_iteration(ts, p):
    """The field dump nearest the pulse crest (t = t_rise) — the SCL snapshot for the CL comparison."""
    t_rise = p["V_PULSE"] / p["V_SLOPE"]
    return ts.iterations[int(np.argmin(np.abs(np.asarray(ts.t) - t_rise)))]


def _crest_iteration(ts, summary, fallback):
    """The particle dump nearest the pulse crest — the gap-full beam. WarpX force-writes a drained
    final-step dump, so iterations[-1] is the wrong template; select by crest_time_s if available."""
    if summary and summary.get("crest_time_s") is not None:
        return ts.iterations[int(np.argmin(np.abs(np.asarray(ts.t) - summary["crest_time_s"])))]
    return fallback


def child_langmuir_figure(V_anode, gap_d, it):
    """On-axis potential and field vs the planar Child–Langmuir and vacuum references at the crest."""
    ts = OpenPMDTimeSeries(FIELDS)
    phi, meta = ts.get_field("phi", iteration=it)
    ez, _ = ts.get_field("E", "z", iteration=it)
    ix0 = np.argmin(np.abs(meta.x))
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


def grid_pulse_figure(gap_d, R_cathode, p, summary):
    """Grid V(t) vs the transmitted current; valid since transit time ≪ pulse duration."""
    ts = OpenPMDTimeSeries(FIELDS)
    times = np.asarray(ts.t)
    j_trans = []
    for it in ts.iterations:
        jz, meta = ts.get_field("j", "z", iteration=it)
        nz = jz.shape[0]
        xm = np.abs(meta.x) <= R_cathode
        j_trans.append(np.abs(jz[nz // 2, xm]).mean())
    j_trans = np.asarray(j_trans)

    v_t = _pulse_voltage(times, p)
    j_cl = np.array([float(child_langmuir_current_density(max(v, 0.0), gap_d)) for v in v_t])
    over = p["over_inject"] * float(child_langmuir_current_density(p["V_OFF"] + p["V_PULSE"], gap_d))

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 6.5), sharex=True, constrained_layout=True)
    a1.plot(times * 1e9, v_t, color="C0")
    a1.axhline(0, color="gray", lw=0.8, ls=":")
    a1.fill_between(times * 1e9, 0, v_t, where=v_t > 0, color="C0", alpha=0.15)
    a1.set_ylabel("grid bias  V(t)  [V]")
    a1.set_title(f"Pulsed grid: peak {p['V_OFF']+p['V_PULSE']:.0f} V, "
                 f"slope {p['V_SLOPE']/1e9:.0f} V/ns, FWHM {p['PULSE_WIDTH']*1e9:.1f} ns")

    a2.plot(times * 1e9, j_trans, "o-", color="C2", ms=3, label="WarpX transmitted (mid-gap)")
    a2.plot(times * 1e9, j_cl, "k--", label=r"$J_{CL}(V(t))$ (quasi-static)")
    a2.axhline(over, color="r", ls=":", label=f"injected ({p['over_inject']:.0f}× peak $J_{{CL}}$)")
    a2.set_xlabel("time  [ns]"); a2.set_ylabel(r"current density  $|J_z|$  [A/m²]")
    a2.set_ylim(0, over * 1.1); a2.legend(loc="upper right", fontsize=9)
    if summary:
        a2.text(0.03, 0.95,
                rf"measured $Q = {summary['q_emit_C']*1e9:.3f}$ nC"
                "\n" rf"(pre-grid {summary['q_pre_grid_C']*1e9:.3f} nC, "
                rf"grid {summary['grid_trans']*100:.0f}%)",
                transform=a2.transAxes, va="top", fontsize=9,
                bbox=dict(boxstyle="round", fc="white", alpha=0.85))
    fig.suptitle("Transmitted current tracks the pulsed Child–Langmuir limit", fontsize=12)
    _save(fig, "grid_pulse")


def anode_spectrum_figure(gap_d, anode_frac, it):
    """Unlike energy_spectrum (whole-gap snapshot), z-cut to the anode slab — the delivered flux."""
    ts = OpenPMDTimeSeries(PARTICLES)
    z, ux, uy, uz, w = ts.get_particle(["z", "ux", "uy", "uz", "w"], iteration=it)
    m = anode_beam_mask(z, uz, gap_d, anode_frac)
    ke = common.ke_kev_from_u(ux[m], uy[m], uz[m]) * 1e3     # keV → eV
    wq = w[m]
    frac = wq.sum() / w.sum()
    mean = np.average(ke, weights=wq)
    std = np.sqrt(np.average((ke - mean) ** 2, weights=wq))

    # Normalize to a fraction: 2D-slab weights are per-unit-out-of-plane-length, not physical
    # charge, so an absolute pC axis would be meaningless — only the spectral shape is.
    fig, ax = plt.subplots(figsize=(7, 4.8), constrained_layout=True)
    ax.hist(ke, bins=60, weights=wq / wq.sum(), color="C2", alpha=0.85)
    ax.axvline(mean, color="k", ls="--", lw=1, label=f"<KE> = {mean:.1f} eV\nσ = {std:.1f} eV")
    ax.set_xlabel("kinetic energy  [eV]"); ax.set_ylabel("fraction of delivered charge / bin")
    ax.set_title(f"Delivered beam at the anode (top {anode_frac*100:.0f}% of gap, forward-moving)\n"
                 f"the flux that seeds the gun — {frac*100:.1f}% of the gap charge")
    ax.legend(loc="upper right", fontsize=9)
    _save(fig, "anode_spectrum")


def emission_phase_space_figure(T_cathode, it):
    """Intrinsic thermal transverse phase space and emittance of the source beam (crest dump)."""
    ts = OpenPMDTimeSeries(PARTICLES)
    x, ux, w = ts.get_particle(["x", "ux", "w"], iteration=it)
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

    w = WarpX(input_file=CONFIG)
    _load_series(w)

    summary = None
    summary_path = os.path.join(DIAG_DIR, "injection_summary.json")
    if os.path.isfile(summary_path):
        with open(summary_path) as f:
            summary = json.load(f)

    # Snapshot figures use the CREST dump (gap full of transiting electrons), not iterations[-1]
    # (the drained final-step dump WarpX force-writes after the grid pulses off).
    ts = OpenPMDTimeSeries(PARTICLES)
    it = _crest_iteration(ts, summary, _last_populated(PARTICLES))
    pg = w._particle_group(iteration=it)

    for name, fig in [
        ("phase_space_z_KE", w.plot2D("z", "kinetic_energy", iteration=it)),
        ("potential_xz",     w.plot_fields("phi", "x", "z")),
        ("charge_density_xz", w.plot_fields("rho", "x", "z")),
        ("energy_spectrum",  common.energy_spectrum(pg, e_unit="eV")),
    ]:
        _save(fig, name)

    # Rising-edge dumps only up to the crest: post-crest dumps are drained and would skew the
    # near-cathode screens with a cold layer.
    crest_t = summary["crest_time_s"] if summary and summary.get("crest_time_s") else ts.t[-1]
    rising = [i for i, t in zip(ts.iterations, ts.t) if t <= crest_t * 1.001]
    pool = common.pool_trajectories(ts, rising, with_y=False)          # 2D slab: no y
    z_m, ke, emit, sigma, _q_pc = common.evolution_screens(pool)   # 2D slab: charge non-physical
    _save(common.evolution_vs_z(z_m, ke, emit, sigma,
                                title="Cathode beam evolution across the gap"),
          "evolution_vs_z")

    # The grid bias is a V(t) parser string, so the SCL reference voltage is the pulse PEAK from
    # params, not w.get().
    gap_d = w.get("grid/upper_bound")[1]
    R_cathode = w.get("species")[0]["upper_bound"][0]
    p = w.get("params")
    v_peak = p["V_OFF"] + p["V_PULSE"]

    fts = OpenPMDTimeSeries(FIELDS)
    child_langmuir_figure(v_peak, gap_d, _peak_field_iteration(fts, p))
    grid_pulse_figure(gap_d, R_cathode, p, summary)
    anode_spectrum_figure(gap_d, p["ANODE_FRAC"], it)
    emission_phase_space_figure(p["T_cathode"], it)


if __name__ == "__main__":
    main()
