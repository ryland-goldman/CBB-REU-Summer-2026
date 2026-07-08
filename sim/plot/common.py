"""Reusable custom beam figures shared by the cathode/gun/injector plot modules.

Each helper takes a `pmd_beamphysics.ParticleGroup` (one diagnostic dump) and returns a
matplotlib Figure -- the physics overlays lume-warpx's generic plot2D/plot1D/plot_fields
do not cover (energy spectrum, longitudinal current, transverse spot, slice chirp).
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from sim.helpers.tools import C_LIGHT, E_CHARGE, MC2_KEV
from sim.helpers.metrics import screen_profile


def gamma_from_u(ux, uy, uz):
    """Lorentz gamma from openPMD normalized momenta u = gamma*beta."""
    return np.sqrt(1.0 + ux ** 2 + uy ** 2 + uz ** 2)


def ke_kev_from_u(ux, uy, uz):
    """Kinetic energy [keV] from u = gamma*beta."""
    return (gamma_from_u(ux, uy, uz) - 1.0) * MC2_KEV


def _wstd(x, w):
    mean = np.average(x, weights=w)
    return float(np.sqrt(np.average((x - mean) ** 2, weights=w)))


_E_UNIT_PER_EV = {"eV": 1.0, "keV": 1e3, "MeV": 1e6}   # divide eV by this to get e_unit


def energy_spectrum(pg, n_bins=80, use_ke=True, e_unit="keV"):
    """Charge-weighted histogram of (kinetic) energy -- the beam's energy spread.

    `e_unit` ("eV"/"keV"/"MeV") scales the energy axis.
    """
    key = "kinetic_energy" if use_ke else "energy"
    div = _E_UNIT_PER_EV[e_unit]
    e = pg[key] / div
    w = pg["weight"] * 1e12                             # C -> pC
    mean, std = np.average(e, weights=pg["weight"]), _wstd(e, pg["weight"])

    fig, ax = plt.subplots()
    ax.hist(e, bins=n_bins, weights=w, color="C0", alpha=0.85)
    ax.axvline(mean, color="k", ls="--", lw=1,
               label=f"<E> = {mean:.3g} {e_unit}\nsigma_E = {std:.3g} {e_unit}  ({100*std/mean:.2f}%)")
    ax.set_xlabel(("kinetic energy" if use_ke else "energy") + f" [{e_unit}]")
    ax.set_ylabel("charge / bin [pC]")
    ax.set_title("Energy spectrum")
    ax.legend(loc="upper left", fontsize=9)
    return fig


def current_profile(pg, n_bins=80):
    """Longitudinal beam current I(z) ~ <v_z>*dQ/dz -- the bunching profile.

    Keyed on z (a WarpX diagnostic is an instantaneous snapshot, sigma_t ~ 0): each particle
    contributes weight*v_z/dz to the local current (per-particle v_z).
    """
    z0 = np.average(pg["z"], weights=pg["weight"])
    z = (pg["z"] - z0) * 1e3                                           # m -> mm, centred
    qv = pg["weight"] * pg["beta_z"] * C_LIGHT                         # C*(m/s) per particle
    counts, edges = np.histogram(z, bins=n_bins, weights=qv)
    dz = (edges[1] - edges[0]) * 1e-3                                  # mm -> m
    centres = 0.5 * (edges[:-1] + edges[1:])
    current = counts / dz                                              # I = sum(w*v_z)/dz [A]

    fig, ax = plt.subplots()
    ax.fill_between(centres, current, step="mid", color="C3", alpha=0.85)
    ax.set_xlabel("z - <z> [mm]")
    ax.set_ylabel("current [A]")
    ax.set_title(f"Longitudinal current profile (peak {current.max():.1f} A)")
    return fig


def beam_spot(pg, bins=160):
    """Transverse x-y spot (RZ stages reconstruct y, so this checks beam roundness)."""
    fig = pg.plot("x", "y", bins=bins, return_figure=True)
    fig.axes[0].set_aspect("equal")
    return fig


def energy_chirp(pg, n_slice=60):
    """Slice <KE> vs z -- the longitudinal energy chirp that drives velocity bunching."""
    fig = pg.slice_plot("mean_kinetic_energy", slice_key="z", n_slice=n_slice,
                        return_figure=True)
    fig.axes[0].set_title("Slice mean energy vs z (longitudinal chirp)")
    return fig


def pool_trajectories(ts, iters, species="electrons", with_y=True):
    """Concatenate (id, z, x, y, ux, uy, uz, ke, w) over every dump for virtual screens.

    `with_y=False` for a 2D (x, z) slab that has no y coordinate (y is then filled with zeros).
    """
    keys = ["id", "z", "x"] + (["y"] if with_y else []) + ["ux", "uy", "uz", "w"]
    cols = {k: [] for k in ("id", "z", "x", "y", "ux", "uy", "uz", "ke", "w")}
    for it in iters:
        got = dict(zip(keys, ts.get_particle(keys, species=species, iteration=it)))
        z = got["z"]
        if not len(z):
            continue
        for k in ("id", "z", "x", "ux", "uy", "uz", "w"):
            cols[k].append(got[k])
        cols["y"].append(got["y"] if with_y else np.zeros_like(z))
        cols["ke"].append(ke_kev_from_u(got["ux"], got["uy"], got["uz"]))
    return {k: (np.concatenate(v) if v else np.array([])) for k, v in cols.items()}


def evolution_screens(pool, n_screen=80):
    """(z [m], <KE> [keV], eps_n,x [mm.mrad], sigma_x [mm], charge [pC]) on fixed-z screens.

    Charge is only physical for the RZ stages -- the 2D cathode slab weight is per-unit-length,
    not Coulombs, so skip its charge panel there.
    """
    screens, prof = screen_profile(
        pool["id"], pool["z"], pool["w"],
        {"x": pool["x"], "ux": pool["ux"], "ke": pool["ke"]},
        emit_pairs=[("x", "ux")], n_screen=n_screen)
    return (screens, prof["mean"]["ke"],
            prof["emit"][("x", "ux")] * 1e6,            # m*(gamma*beta) -> mm.mrad
            prof["rms"]["x"] * 1e3,                      # m -> mm
            prof["charge"] * E_CHARGE * 1e12)            # real-particle weight -> pC


def evolution_vs_z(z_m, ke, emit, sigma, charge_pc=None, ke_unit="keV", title="", notes=None):
    """Beam evolution along z: mean KE, normalized emittance eps_n,x, RMS spot size, and -- when
    `charge_pc` is given -- the charge [pC] crossing each screen (a 4th panel showing loss vs z).

    Arrays are aligned on z_m [m]. `notes` is an optional {panel: str} (panel in
    {"ke","emit","sigma","charge"}) for per-panel caveats.
    """
    notes = notes or {}
    z_mm = np.asarray(z_m) * 1e3
    panels = [
        (ke, "C2", f"mean KE  [{ke_unit}]", "ke"),
        (emit, "C3", r"$\varepsilon_{n,x}$  [mm$\cdot$mrad]", "emit"),
        (sigma, "C0", r"RMS size  $\sigma_x$  [mm]", "sigma"),
    ]
    if charge_pc is not None:
        panels.append((charge_pc, "C4", "charge  [pC]", "charge"))
    n = len(panels)
    fig, axs = plt.subplots(n, 1, figsize=(7.4, 2.8 * n), constrained_layout=True, sharex=True)
    axs = np.atleast_1d(axs)
    for ax, (y, color, ylab, key) in zip(axs, panels):
        y = np.asarray(y, float)
        ok = np.isfinite(y)
        ax.plot(z_mm[ok], y[ok], "o-", color=color, ms=3)
        ax.set_ylabel(ylab)
        ax.grid(alpha=0.3)
        # Charge loss can span decades (aperture scrape): log-y when it does, else linear.
        if key == "charge":
            pos = y[ok][y[ok] > 0]
            if pos.size and pos.max() / pos.min() > 20:
                ax.set_yscale("log")
        if key in notes:
            ax.set_title(notes[key], fontsize=8)
    axs[-1].set_xlabel("beam position  z  [mm]")
    if title:
        fig.suptitle(title, fontsize=12)
    return fig


def transverse_rpr(x, y, ux, uy, w=None, title="Transverse phase space  (r, p_r)",
                   p_unit="keV"):
    """Transverse r-p_r phase space: r = hypot(x, y) [mm], p_r = (x.ux + y.uy)/r . MC2 [p_unit/c].

    `p_unit` ("keV"/"MeV") scales the momentum axis.
    """
    x, y, ux, uy = (np.asarray(a, float) for a in (x, y, ux, uy))
    r = np.hypot(x, y)
    safe = r > 0
    pr_u = np.zeros_like(r)
    pr_u[safe] = (x[safe] * ux[safe] + y[safe] * uy[safe]) / r[safe]
    r_mm = r * 1e3
    pr_p = pr_u * MC2_KEV * 1e3 / _E_UNIT_PER_EV[p_unit]   # u*MC2[eV/c] -> p_unit/c

    fig, ax = plt.subplots(figsize=(7.0, 4.8), constrained_layout=True)
    hb = ax.hexbin(r_mm, pr_p, gridsize=70, cmap="viridis", mincnt=1)
    fig.colorbar(hb, ax=ax, label="macroparticles / bin")
    ax.set_xlabel("r  [mm]"); ax.set_ylabel(rf"$p_r$  [{p_unit}/$c$]")
    ax.set_title(title)
    return fig
