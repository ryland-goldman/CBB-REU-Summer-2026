"""Thermionic-emission physics shared by the cathode and gun stages.

Child-Langmuir space-charge-limited current and the cathode thermal momentum
spread, single-sourced so the sim and its plot overlays cannot drift. See
cathode/README.md for the emission model and its provenance.
"""

import numpy as np

from pipeline.constants import C_LIGHT, E_CHARGE, EPS0, K_B_EV, M_E, MC2_EV


def child_langmuir_current_density(voltage, gap):
    """Space-charge-limited J [A/m^2] across a planar gap; 0 below grid cutoff (V<=0)."""
    v = np.maximum(np.asarray(voltage, dtype=float), 0.0)
    return (4.0 / 9.0) * EPS0 * np.sqrt(2.0 * E_CHARGE / M_E) * v**1.5 / gap**2


def thermal_velocity_sigma(t_k):
    """RMS thermal velocity per Cartesian component [m/s] for a Maxwellian cathode.

    Thermal energies are non-relativistic so u~v: sigma = sqrt(kT/m_e)
    = sqrt(kT[eV]/mc2[eV]) * c. WarpX consumes this directly as `rms_velocity`;
    the plotter consumes the dimensionless sqrt(kT[eV]/mc2[eV]) = sigma / c.
    """
    kt_ev = K_B_EV * t_k
    return np.sqrt(kt_ev / MC2_EV) * C_LIGHT
