"""
Figures for the finite-cathode space-charge-limited diode (cathode_diode.py):
reads cathode/diags/ openPMD output and writes six PNGs to cathode/results/.

See cathode/README.md for the physics, parameters, and a description of each figure.
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from openpmd_viewer import OpenPMDTimeSeries
from scipy.constants import e as q_e, epsilon_0, m_e
from scipy.constants import k as k_B, c

# Display/overlay constants for theory curves and titles; kept in sync with
# cathode_diode.py by hand (not read by the sim).
V_anode   = 60.0
gap_d     = 200.0e-6
R_cathode = 8.0e-3
over_inject = 2.0
T_cathode = 1425.0

RESULTS = "cathode/results"


def main():
    os.makedirs(RESULTS, exist_ok=True)

    # Computed here (not at module top) so a config() override applied via setattr
    # after import is reflected, matching cathode_diode.main().
    J_CL = (4.0 / 9.0) * epsilon_0 * np.sqrt(2.0 * q_e / m_e) * V_anode**1.5 / gap_d**2

    ts  = OpenPMDTimeSeries("cathode/diags/fields")
    it  = ts.iterations[-1]            # final (steady-state) snapshot

    # openPMD field arrays are (nz, nx); meta.z, meta.x give the axes.
    phi, meta = ts.get_field("phi", iteration=it)
    ex,  _    = ts.get_field("E", "x", iteration=it)
    ez,  _    = ts.get_field("E", "z", iteration=it)
    rho, _    = ts.get_field("rho", iteration=it)
    z, x = meta.z, meta.x
    ix0 = np.argmin(np.abs(x))          # column nearest the axis x = 0
    Emag = np.sqrt(ex**2 + ez**2)

    extent = [x.min() * 1e3, x.max() * 1e3, z.min() * 1e3, z.max() * 1e3]  # mm

    # Figure 1 — 2D maps
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True)

    # |ρ| (electrons make ρ < 0) on a √ scale to expose the bulk gradient.
    absrho = np.abs(rho) * 1e6
    rho_norm = PowerNorm(gamma=0.5, vmin=0.0, vmax=absrho.max())

    panels = [
        (absrho, "|Charge density|  |ρ|  [µC/m³]  (√)", "viridis", rho_norm),
        (phi,       "Potential  φ  [V]",          "plasma",  None),
        (Emag * 1e-3, "Field magnitude  |E|  [kV/m]", "inferno", None),
    ]
    for ax, (data, title, cmap, norm) in zip(axs, panels):
        im = ax.imshow(data, extent=extent, origin="lower", aspect="auto",
                       cmap=cmap, norm=norm)
        fig.colorbar(im, ax=ax, shrink=0.9)
        # mark the emitting cathode patch (z = 0, |x| < R)
        ax.plot([-R_cathode * 1e3, R_cathode * 1e3], [0, 0], "w-", lw=3,
                solid_capstyle="butt")
        ax.set_title(title)
        ax.set_xlabel("x  [mm]")
    axs[0].set_ylabel("z  [mm]   (cathode → anode)")
    fig.suptitle("Finite thermionic cathode in WarpX — emission from |x| < "
                 f"{R_cathode*1e3:.0f} mm (white bar); note field transition at edges",
                 fontsize=12)
    fig.savefig(f"{RESULTS}/cathode_2d.png", dpi=140)
    print(f"wrote {RESULTS}/cathode_2d.png")

    # Figure 2 — on-axis profiles vs. Child–Langmuir theory
    phi_axis = phi[:, ix0]
    ez_axis  = ez[:, ix0]

    # Child–Langmuir (1D, planar) reference for the same V and gap:
    zt = z
    phi_cl = V_anode * (zt / gap_d) ** (4.0 / 3.0)
    ez_cl  = -(4.0 / 3.0) * (V_anode / gap_d) * (zt / gap_d) ** (1.0 / 3.0)
    # vacuum (no space charge) reference: linear potential, uniform field
    phi_vac = V_anode * (zt / gap_d)
    ez_vac  = -(V_anode / gap_d) * np.ones_like(zt)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)

    a1.plot(zt * 1e3, phi_axis, "o-", color="C0", ms=3, label="WarpX (on axis)")
    a1.plot(zt * 1e3, phi_cl, "k--", label=r"Child–Langmuir  $V(z/d)^{4/3}$")
    a1.plot(zt * 1e3, phi_vac, ":", color="gray", label="vacuum (no space charge)")
    a1.set_xlabel("z  [mm]"); a1.set_ylabel("φ  [V]")
    a1.set_title("On-axis potential"); a1.legend()

    a2.plot(zt * 1e3, ez_axis * 1e-3, "o-", color="C3", ms=3, label="WarpX (on axis)")
    a2.plot(zt * 1e3, ez_cl * 1e-3, "k--",
            label=r"Child–Langmuir  $-\frac{4V}{3d}(z/d)^{1/3}$")
    a2.plot(zt * 1e3, ez_vac * 1e-3, ":", color="gray", label="vacuum")
    a2.set_xlabel("z  [mm]"); a2.set_ylabel("$E_z$  [kV/m]")
    a2.set_title("On-axis longitudinal field"); a2.legend()

    fig.suptitle("Space-charge depression of the field at the cathode follows "
                 "the Child–Langmuir law", fontsize=12)
    fig.savefig(f"{RESULTS}/child_langmuir.png", dpi=140)
    print(f"wrote {RESULTS}/child_langmuir.png")

    # Figure 3 — transmitted current saturates at J_CL despite 2× over-injection
    # Integrate ∫jz dx across the full domain (robust to transverse spreading) and
    # reference to the cathode width 2R. jz < 0 (electrons +z carry -q); use |jz|.
    iz_anode = -2                       # row just inside the anode
    dx       = x[1] - x[0]

    times, J_trans, rho_zt = [], [], []
    for i, itr in enumerate(ts.iterations):
        jz, _ = ts.get_field("j", "z", iteration=itr)
        I_line = np.abs(jz[iz_anode, :].sum() * dx)     # total current [A/m depth]
        J_trans.append(I_line / (2.0 * R_cathode))      # referenced to cathode width
        rho_it, _ = ts.get_field("rho", iteration=itr)
        rho_zt.append(rho_it[:, ix0])                   # on-axis charge density column
        times.append(ts.t[i])
    times = np.array(times) * 1e9       # ns
    J_trans = np.array(J_trans)
    rho_zt = np.array(rho_zt).T         # shape (nz, n_times): rho(z, t) on axis

    fig, ax = plt.subplots(figsize=(7, 4.8), constrained_layout=True)
    ax.plot(times, J_trans, "o-", color="C2", label="WarpX transmitted current")
    ax.axhline(J_CL, color="k", ls="--", label=r"Child–Langmuir limit $J_{CL}$")
    ax.axhline(over_inject * J_CL, color="r", ls=":",
               label=f"injected current ({over_inject:.0f}× $J_{{CL}}$)")
    ax.set_xlabel("time  [ns]")
    ax.set_ylabel(r"current density at anode  $|J_z|$  [A/m²]")
    ax.set_title("Emission self-limits to the Child–Langmuir value")
    ax.set_xlim(0, 0.15)
    # headroom for the 2×J_CL injected-current reference line (the over-injection story)
    ax.set_ylim(0, over_inject * J_CL * 1.1)
    ax.legend()
    fig.savefig(f"{RESULTS}/current_saturation.png", dpi=140)
    print(f"wrote {RESULTS}/current_saturation.png")

    # Figure 4 — on-axis charge density ρ(z, t): build-up of the space-charge cloud
    fig, ax = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
    abs_rho_zt = np.abs(rho_zt) * 1e6                   # |ρ| in µC/m³
    # pcolormesh on the true (non-uniform) time coords; imshow would distort the axis.
    im = ax.pcolormesh(
        times, z * 1e3, abs_rho_zt,
        shading="nearest", cmap="viridis",
        norm=PowerNorm(gamma=0.5, vmin=0.0, vmax=abs_rho_zt.max()),
    )
    fig.colorbar(im, ax=ax, label="|charge density|  |ρ|  [µC/m³]  (√)")
    ax.set_xlim(0, 0.15)            # densely-sampled transient; cloud is steady after
    ax.set_xlabel("time  [ns]")
    ax.set_ylabel("z  [mm]   (cathode → anode)")
    ax.set_title("On-axis charge density vs. time — space-charge cloud fills the gap")
    fig.savefig(f"{RESULTS}/rho_z_time.png", dpi=140)
    print(f"wrote {RESULTS}/rho_z_time.png")

    # Figure 5 — equipotentials + E-field streamlines: the 2D cathode-EDGE transition
    # streamplot needs 1-D coord vectors matching the array axes: fields are (nz, nx),
    # so call streamplot(x_mm, z_mm, Ex, Ez) with len(x_mm)=nx, len(z_mm)=nz.
    x_mm = x * 1e3                       # nx, transverse [mm]
    z_mm = z * 1e3                       # nz, longitudinal [mm]
    Ex_k, Ez_k = ex * 1e-3, ez * 1e-3   # kV/m; sign/orientation untouched
    speed = np.sqrt(Ex_k**2 + Ez_k**2)   # |E| [kV/m], colours the streamlines
    phi_levels = np.linspace(phi.min(), V_anode, 12)   # ~12 equipotential lines

    fig, (p1, p2) = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)

    def draw_field(ax, lw_speed=True):
        """φ contours + E streamlines on the (x, z) gap; shared by both panels."""
        cs = ax.contour(x_mm, z_mm, phi, levels=phi_levels,
                        colors="0.25", linewidths=0.8)
        ax.clabel(cs, inline=True, fontsize=6, fmt="%.0f V")
        lw = (1.4 * speed / speed.max() + 0.3) if lw_speed else 1.0
        strm = ax.streamplot(x_mm, z_mm, Ex_k, Ez_k, color=speed, cmap="inferno",
                             density=1.4, linewidth=lw, arrowsize=0.8)
        # emitting patch (z = 0, |x| < R) as a white bar; cathode edges as dotted lines
        ax.plot([-R_cathode * 1e3, R_cathode * 1e3], [0, 0], "w-", lw=4,
                solid_capstyle="butt", zorder=5)
        for xe in (-R_cathode * 1e3, R_cathode * 1e3):
            ax.axvline(xe, color="cyan", ls=":", lw=1.2, zorder=4)
        ax.set_xlabel("x  [mm]")
        return strm

    # Panel 1 — full gap
    strm = draw_field(p1)
    p1.set_ylabel("z  [mm]   (cathode → anode)")
    p1.set_xlim(x_mm.min(), x_mm.max())
    p1.set_ylim(0, z_mm.max())
    p1.set_title("Full gap: equipotentials + E streamlines")

    # Panel 2 — zoom on the +x edge
    Redge = R_cathode * 1e3
    draw_field(p2)
    p2.set_xlim(max(0.0, Redge - 5.0), Redge + 1.0)
    p2.set_ylim(0, z_mm.max())
    p2.set_xlabel("x  [mm]")
    p2.set_title(f"Zoom on +x cathode edge (x = +{Redge:.0f} mm)")

    fig.colorbar(strm.lines, ax=[p1, p2], shrink=0.85, label="|E|  [kV/m]")
    fig.suptitle(f"2D field at the finite cathode — equipotentials crowd and field "
                 f"lines splay at the edges x = ±{Redge:.0f} mm (the field transition at the "
                 "emission edge)", fontsize=12)
    fig.savefig(f"{RESULTS}/field_lines.png", dpi=140)
    print(f"wrote {RESULTS}/field_lines.png")
    plt.close(fig)

    # Figure 6 — intrinsic thermal transverse phase space + emittance of the source
    # openPMD/WarpX store ux = γβ_x, so εn,x = sqrt(⟨x²⟩⟨ux²⟩ − ⟨x·ux⟩²) directly.
    tsp = OpenPMDTimeSeries("cathode/diags/particles")
    itp = tsp.iterations[-1]                 # last particle snapshot
    xp, uxp, uzp, wp = tsp.get_particle(["x", "ux", "uz", "w"], iteration=itp)

    # Weighted central moments (correct for any weighting).
    xbar  = np.average(xp,  weights=wp)
    uxbar = np.average(uxp, weights=wp)
    x2  = np.average((xp - xbar) ** 2, weights=wp)
    ux2 = np.average((uxp - uxbar) ** 2, weights=wp)
    xux = np.average((xp - xbar) * (uxp - uxbar), weights=wp)
    emit_n = np.sqrt(max(x2 * ux2 - xux ** 2, 0.0))      # [m·rad]; ux is γβ_x
    emit_n_mm_mrad = emit_n * 1e6                         # m·rad → mm·mrad

    # p_x = γβ_x · m_ec² in keV/c (p_x·c = γβ_x·m_ec² is an energy).
    MC2_keV = m_e * c ** 2 / q_e / 1e3                   # electron rest energy [keV]
    px = uxp * MC2_keV                                   # transverse momentum [keV/c]
    u_th  = np.sqrt(k_B * T_cathode / (m_e * c ** 2))    # rms of γβ_x (dimensionless)
    p_th  = u_th * MC2_keV                               # rms thermal momentum [keV/c]

    fig, (b1, b2) = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)

    # Panel 1 — transverse phase space x [mm] vs p_x [keV/c], density via hexbin.
    hb = b1.hexbin(xp * 1e3, px, gridsize=70, cmap="viridis", mincnt=1)
    fig.colorbar(hb, ax=b1, label="macroparticles / bin")
    b1.set_xlabel("x  [mm]")
    b1.set_ylabel(r"$p_x$  [keV/$c$]")
    b1.set_title("Transverse phase space at the cathode")
    # 2D-slab εn,x; the gun's RZ remap (uniform disc) receives ×√(3/4) of it.
    b1.text(0.03, 0.97,
            rf"$\varepsilon_{{n,x}} = {emit_n_mm_mrad:.3f}$ mm·mrad"
            "\n" rf"$\sqrt{{\langle x^2\rangle}} = {np.sqrt(x2)*1e3:.2f}$ mm"
            "\n" rf"gun receives disc value $\times\sqrt{{3/4}}"
            rf" \approx {emit_n_mm_mrad*np.sqrt(0.75):.2f}$ mm·mrad"
            "\n(full-gap snapshot; drift correlation negligible)",
            transform=b1.transAxes, va="top", ha="left", fontsize=9,
            bbox=dict(boxstyle="round", fc="white", alpha=0.85))

    # Panel 2 — histogram of p_x with the expected ±√(kT·m_ec²) scale overlaid.
    b2.hist(px, bins=120, color="C0", alpha=0.8, density=True)
    for s, lbl in ((+1, r"$\pm\sqrt{kT\,m_ec^2}$"), (-1, None)):
        b2.axvline(s * p_th, color="k", ls="--", lw=1.2, label=lbl)
    b2.set_xlabel(r"$p_x$  [keV/$c$]")
    b2.set_ylabel("probability density")
    b2.set_title(f"Thermal transverse momentum spread ({T_cathode:.0f} K)")
    b2.legend(loc="upper right", fontsize=9)
    b2.text(0.03, 0.97,
            rf"rms $p_x = {np.sqrt(ux2)*MC2_keV:.3f}$ keV/$c$"
            "\n" rf"$\sqrt{{kT\,m_ec^2}} = {p_th:.3f}$ keV/$c$",
            transform=b2.transAxes, va="top", ha="left", fontsize=9,
            bbox=dict(boxstyle="round", fc="white", alpha=0.85))

    fig.suptitle("Intrinsic thermal emittance of the cathode — the source beam "
                 "quality inherited by the downstream gun", fontsize=12)
    fig.savefig(f"{RESULTS}/emission_phase_space.png", dpi=140)
    print(f"wrote {RESULTS}/emission_phase_space.png")
    plt.close(fig)

    # Quantitative sanity check
    J_final = J_trans[-1]
    print(f"\nJ_CL (theory)          = {J_CL:8.1f} A/m²")
    print(f"injected               = {over_inject*J_CL:8.1f} A/m²")
    print(f"transmitted (steady)   = {J_final:8.1f} A/m²  "
          f"({100*J_final/J_CL:.0f}% of J_CL)")
    print(f"thermal emittance εn,x = {emit_n_mm_mrad:8.4f} mm·mrad  "
          f"(rms u_x = {np.sqrt(ux2)*1e3:.3f}e-3 vs √(kT/mc²) = {u_th*1e3:.3f}e-3)")


if __name__ == "__main__":
    main()
