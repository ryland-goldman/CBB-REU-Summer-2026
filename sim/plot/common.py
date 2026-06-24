"""Reusable custom beam figures shared by the cathode/gun/injector plot modules.

Each helper takes a `pmd_beamphysics.ParticleGroup` (one diagnostic dump) and returns a
matplotlib Figure -- the physics overlays lume-warpx's generic plot2D/plot1D/plot_fields
do not cover (energy spectrum, longitudinal current, transverse spot, slice chirp).
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from sim.helpers.tools import C_LIGHT, MC2_KEV


def gamma_from_u(ux, uy, uz):
    """Lorentz gamma from openPMD normalized momenta u = gamma*beta."""
    return np.sqrt(1.0 + ux ** 2 + uy ** 2 + uz ** 2)


def ke_kev_from_u(ux, uy, uz):
    """Kinetic energy [keV] from u = gamma*beta."""
    return (gamma_from_u(ux, uy, uz) - 1.0) * MC2_KEV


def _wstd(x, w):
    mean = np.average(x, weights=w)
    return float(np.sqrt(np.average((x - mean) ** 2, weights=w)))


def energy_spectrum(pg, n_bins=80, use_ke=True):
    """Charge-weighted histogram of (kinetic) energy -- the beam's energy spread."""
    key = "kinetic_energy" if use_ke else "energy"
    e = pg[key] / 1e3                                   # eV -> keV
    w = pg["weight"] * 1e12                             # C -> pC
    mean, std = np.average(e, weights=pg["weight"]), _wstd(e, pg["weight"])

    fig, ax = plt.subplots()
    ax.hist(e, bins=n_bins, weights=w, color="C0", alpha=0.85)
    ax.axvline(mean, color="k", ls="--", lw=1,
               label=f"<E> = {mean:.1f} keV\nsigma_E = {std*1e3:.0f} eV  ({100*std/mean:.2f}%)")
    ax.set_xlabel(("kinetic energy" if use_ke else "energy") + " [keV]")
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
