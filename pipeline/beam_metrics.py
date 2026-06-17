"""Shared beam-moment helpers (pure numpy; no pywarpx).

See pipeline/README.md for the emittance unit conventions and which plotters use these.
"""

import numpy as np


def rms_emit(q, uq, w):
    """Charge-weighted normalized rms emittance sqrt(⟨q²⟩⟨uq²⟩ − ⟨q·uq⟩²) for one phase plane.

    ``uq`` is the openPMD ``u`` = γβ (already includes γ — do NOT multiply by γ again).
    Returns the RAW emittance in m·(γβ); the caller applies the unit scaling.
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
