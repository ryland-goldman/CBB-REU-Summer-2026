"""Shared beam-moment helpers (pure numpy; no pywarpx).

Home for the small, stage-independent beam statistics that several plotters compute
identically — so the math lives in one place and the stage plotters can't drift apart.
(Deliberately NOT in `pipeline/_runner.py`, which is the pywarpx/subprocess/tqdm
machinery — wrong home for pure-numpy moment math.)
"""

import numpy as np


def rms_emit(q, uq, w):
    """Normalized rms emittance sqrt(⟨q²⟩⟨uq²⟩ − ⟨q·uq⟩²) for one phase plane.

    ``q`` is a position [m], ``uq`` the conjugate normalized momentum γβ (openPMD ``u``,
    which already includes γ — do NOT multiply by γ again). ``w`` are the macroparticle
    weights. Returns the RAW emittance in m·(γβ); the caller applies the plane-specific
    unit scaling: ×1e6 for the transverse plane (→ mm·mrad), or ×1e3 for the longitudinal
    z–(γβ_z) plane (→ mm·dimensionless, NOT mm·mrad — the longitudinal plane is not an
    angle, so it gets no mrad factor).

    All moments are charge-weighted by ``w``. (Note: this differs from the older inline
    gun emittance, which used unweighted np.mean — see gun/plot_gun.py; the weighted form
    is correct for the downsampled+reweighted snapshots.)
    """
    w = np.asarray(w, dtype=float)
    sw = w.sum()
    if sw <= 0 or len(q) < 2:
        return 0.0
    qm = np.average(q, weights=w)
    um = np.average(uq, weights=w)
    q2 = np.average((q - qm) ** 2, weights=w)
    u2 = np.average((uq - um) ** 2, weights=w)
    qu = np.average((q - qm) * (uq - um), weights=w)
    return float(np.sqrt(max(q2 * u2 - qu * qu, 0.0)))
