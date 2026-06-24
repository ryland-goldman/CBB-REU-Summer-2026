"""Shared physics constants, emission physics, RF drive strings, and a little runtime
plumbing (environment setup, file-path module loading).

Everything here is stage-agnostic so the stage drivers read as physics. Physical constants
come from scipy so no stage carries a divergent literal (u = gamma*beta; momentum in eV/c;
rest energy in eV).
"""

import os
import resource
import sys

import numpy as np
import scipy.constants as _sc

# ── Physical constants (SI + eV conventions) ─────────────────────────────────────
C_LIGHT = _sc.c                            # m/s
E_CHARGE = _sc.e                           # C (elementary charge, positive)
M_E = _sc.m_e                              # kg (electron mass)
EPS0 = _sc.epsilon_0                       # F/m
K_B = _sc.k                                # J/K
K_B_EV = _sc.k / _sc.e                     # eV/K
MC2_EV = _sc.m_e * _sc.c ** 2 / _sc.e      # electron rest energy [eV] ~= 510998.95
MC2_KEV = MC2_EV / 1e3

# Repo root: this file is <root>/sim/helpers/tools.py.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# ── Emission physics (cathode + gun) ─────────────────────────────────────────────
def child_langmuir_current_density(voltage, gap):
    """Space-charge-limited current density J [A/m^2] across a planar gap; 0 for V<=0."""
    v = np.maximum(np.asarray(voltage, dtype=float), 0.0)
    return (4.0 / 9.0) * EPS0 * np.sqrt(2.0 * E_CHARGE / M_E) * v ** 1.5 / gap ** 2


def thermal_velocity_sigma(t_k):
    """RMS thermal velocity per Cartesian component [m/s] for a Maxwellian cathode.

    Thermal energies are non-relativistic (u ~ v): sigma = sqrt(kT/m_e)
    = sqrt(kT[eV]/mc2[eV]) * c. WarpX consumes this as `rms_velocity`.
    """
    return np.sqrt(K_B_EV * t_k / MC2_EV) * C_LIGHT


# ── RF cavity drive ──────────────────────────────────────────────────────────────
def rf_time_functions(scale, omega, phi, amp_prec=10, phase_prec=10):
    """(E, B) `warpx_*_time_function` strings for a standing-wave TM cavity drive:
    E ~ scale*cos(omega*t + phi), B ~ scale*sin(omega*t + phi). omega keeps .10e --
    its truncation accumulates over the ~ns transit.
    """
    e = f"{scale:.{amp_prec}e}*cos({omega:.10e}*t + ({phi:.{phase_prec}e}))"
    b = f"{scale:.{amp_prec}e}*sin({omega:.10e}*t + ({phi:.{phase_prec}e}))"
    return e, b


# ── Runtime plumbing ─────────────────────────────────────────────────────────────
def prepare_env():
    """Set the env WarpX/HDF5 latch at import time, raise the fd limit, and run from the
    repo root (every stage uses repo-relative paths). Call before importing warpx.
    """
    # OMP_NUM_THREADS must be set before pywarpx loads OpenMP; OMP_THREADS overrides.
    os.environ["OMP_NUM_THREADS"] = os.environ.get(
        "OMP_THREADS", os.environ.get("OMP_NUM_THREADS", "1"))
    os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    try:                                   # openpmd-viewer leaks one fd per get_particle()
        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        if soft < 16384:
            resource.setrlimit(resource.RLIMIT_NOFILE, (min(hard, 16384), hard))
    except (ValueError, OSError):
        pass
    if os.getcwd() != REPO_ROOT:
        os.chdir(REPO_ROOT)
    if REPO_ROOT not in sys.path:
        sys.path.insert(0, REPO_ROOT)
