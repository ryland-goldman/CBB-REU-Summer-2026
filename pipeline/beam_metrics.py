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


def screen_profile(idP, zP, wP, quantities, emit_pairs=(), nscreen=80,
                   min_cross=20, zlo=None, zhi=None):
    """Charge-weighted beam moments on fixed-z virtual SCREENS (station diagnostic).

    Pools every ``(id, dump)`` row and treats each macroparticle's id-trajectory as a
    path that crosses every z-plane once (valid where motion is forward / monotonic in
    z — accelerating WarpX stages). For each screen each crossing particle's listed
    quantities are linearly interpolated along its own trajectory to the plane, then
    charge-weighted moments are accumulated per screen. This is a true LOCAL-in-z phase
    space: every screen sees each particle exactly once — no z-binning and no
    pooling-stream assumption (a z-histogram of the pooled quasi-DC stream groups
    unrelated particles and makes a difference-of-large-moments quantity like εn,x
    jitter bin to bin). See gun/README.md → ``beam_envelope.png``.

    Parameters
    ----------
    idP, zP, wP : (N,) arrays pooled over all dumps (particle id, z [m], weight).
    quantities  : dict name -> (N,) array; per-screen mean/rms/max are returned for each.
    emit_pairs  : iterable of (qname, uname) (both must be in ``quantities``); returns
                  the raw rms emittance sqrt(⟨q²⟩⟨u²⟩−⟨qu⟩²) per screen for each pair.

    Returns ``(screens, out)``: ``screens`` are the plane positions [m]; ``out`` has
    ``count`` and dicts ``mean``/``rms``/``max`` (keyed by quantity name) and ``emit``
    (keyed by the ``(qname, uname)`` tuple). All moments are RAW (caller scales units)
    and NaN on screens crossed by < ``min_cross`` particles.
    """
    idP = np.asarray(idP); zP = np.asarray(zP, float); wP = np.asarray(wP, float)
    names = list(quantities)
    Q = {n: np.asarray(quantities[n], float) for n in names}
    nan = np.full(nscreen, np.nan)
    empty = dict(count=np.zeros(nscreen),
                 mean={n: nan.copy() for n in names},
                 rms={n: nan.copy() for n in names},
                 max={n: nan.copy() for n in names},
                 emit={p: nan.copy() for p in emit_pairs})
    if zP.size == 0:
        return np.linspace(0.0, 1.0, nscreen), empty
    lo = zP.min() if zlo is None else zlo
    hi = zP.max() if zhi is None else zhi
    screens = np.linspace(lo, hi, nscreen)

    W = np.zeros(nscreen); Cnt = np.zeros(nscreen)
    S1 = {n: np.zeros(nscreen) for n in names}      # Σ w·v
    S2 = {n: np.zeros(nscreen) for n in names}      # Σ w·v²
    Mx = {n: np.full(nscreen, -np.inf) for n in names}
    SP = {p: np.zeros(nscreen) for p in emit_pairs}  # Σ w·q·u (cross term)

    order = np.argsort(idP, kind="stable")
    ids = idP[order]
    zs = zP[order]; ws = wP[order]
    Qs = {n: Q[n][order] for n in names}
    bnd = np.flatnonzero(np.r_[True, ids[1:] != ids[:-1], True])  # id-group boundaries
    for g in range(len(bnd) - 1):
        a, b = bnd[g], bnd[g + 1]
        if b - a < 2:                               # need ≥2 dumps to interpolate
            continue
        o = a + np.argsort(zs[a:b])                 # np.interp needs increasing z
        zi = zs[o]
        m = (screens >= zi[0]) & (screens <= zi[-1])
        if not m.any():
            continue
        si = np.flatnonzero(m); sc = screens[si]
        wi = ws[a]                                  # weight is constant per id
        vals = {n: np.interp(sc, zi, Qs[n][o]) for n in names}
        W[si] += wi; Cnt[si] += 1.0
        for n in names:
            v = vals[n]
            S1[n][si] += wi * v; S2[n][si] += wi * v * v
            np.maximum.at(Mx[n], si, v)
        for (qn, un) in emit_pairs:
            SP[(qn, un)][si] += wi * vals[qn] * vals[un]

    ok = Cnt >= min_cross
    Wm = np.where(ok, W, np.nan)
    out = dict(count=Cnt, mean={}, rms={}, max={}, emit={})
    for n in names:
        mu = S1[n] / Wm
        out["mean"][n] = mu
        out["rms"][n] = np.sqrt(np.clip(S2[n] / Wm - mu * mu, 0.0, None))
        out["max"][n] = np.where(ok, Mx[n], np.nan)
    for (qn, un) in emit_pairs:
        mq = S1[qn] / Wm; mu = S1[un] / Wm
        vq = S2[qn] / Wm - mq * mq
        vu = S2[un] / Wm - mu * mu
        cqu = SP[(qn, un)] / Wm - mq * mu
        out["emit"][(qn, un)] = np.sqrt(np.clip(vq * vu - cqu * cqu, 0.0, None))
    return screens, out
