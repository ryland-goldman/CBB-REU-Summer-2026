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


def _group_bounds(sorted_ids):
    """Start/stop indices of each run of equal ids in an id-sorted array."""
    edges = np.flatnonzero(np.r_[True, sorted_ids[1:] != sorted_ids[:-1], True])
    return zip(edges[:-1], edges[1:])


def screen_profile(ids, z, weight, quantities, emit_pairs=(),
                   n_screen=80, min_cross=20, z_range=None):
    """Charge-weighted beam moments on fixed-z virtual screens (a station diagnostic).

    Pool every (id, dump) row from a volumetric run and treat each macroparticle's
    id-trajectory as a path that crosses every z-plane once (valid for forward,
    monotonic-in-z motion — the accelerating WarpX stages). For each screen, every
    crossing particle's listed quantities are linearly interpolated along its own
    trajectory to that plane, then charge-weighted moments accumulate. The result is a
    true local-in-z phase space — each screen sees a particle exactly once, with no
    z-binning and no quasi-DC pooling assumption (which would make a difference-of-large-
    moments quantity like εn,x jitter plane to plane). See gun/README.md → beam_envelope.

    Parameters
    ----------
    ids, z, weight : (N,) arrays pooled over all dumps (particle id, z [m], weight).
    quantities     : {name: (N,) array}; per-screen mean/rms/max returned for each.
    emit_pairs     : iterable of (q_name, u_name); per-screen raw rms emittance for each.

    Weight is assumed constant per id (set once at injection). Interpolation is linear in
    z, so an oscillation under-resolved by the dump cadence biases ``emit`` low. A particle
    contributes only to screens inside its own z-span, so downstream of a loss (e.g. an
    iris) the screens describe the surviving subset.

    Returns ``(screens, out)``: ``screens`` are plane positions [m]; ``out`` has ``count``
    and dicts ``mean``/``rms``/``max`` (by quantity) and ``emit`` (by pair). Moments are
    RAW (caller scales units) and NaN on screens with fewer than ``min_cross`` crossings.
    """
    ids = np.asarray(ids)
    z = np.asarray(z, float)
    weight = np.asarray(weight, float)
    names = list(quantities)
    quantities = {n: np.asarray(quantities[n], float) for n in names}

    nan = np.full(n_screen, np.nan)
    if z.size == 0:
        empty = dict(count=np.zeros(n_screen),
                     mean={n: nan.copy() for n in names},
                     rms={n: nan.copy() for n in names},
                     max={n: nan.copy() for n in names},
                     emit={p: nan.copy() for p in emit_pairs})
        return np.linspace(0.0, 1.0, n_screen), empty

    lo, hi = z_range if z_range else (z.min(), z.max())
    screens = np.linspace(lo, hi, n_screen)

    # Charge-weighted accumulators over all id-trajectories that cross each screen.
    count = np.zeros(n_screen)
    sum_w = np.zeros(n_screen)
    sum_wv = {n: np.zeros(n_screen) for n in names}          # Σ w·v
    sum_wvv = {n: np.zeros(n_screen) for n in names}         # Σ w·v²
    peak = {n: np.full(n_screen, -np.inf) for n in names}
    sum_wqu = {p: np.zeros(n_screen) for p in emit_pairs}    # Σ w·q·u (emittance cross term)

    order = np.argsort(ids, kind="stable")
    ids, z, weight = ids[order], z[order], weight[order]
    quantities = {n: q[order] for n, q in quantities.items()}

    for a, b in _group_bounds(ids):
        if b - a < 2:                                        # need ≥2 dumps to interpolate
            continue
        traj = a + np.argsort(z[a:b])                        # np.interp needs increasing z
        z_traj = z[traj]
        on = (screens >= z_traj[0]) & (screens <= z_traj[-1])
        if not on.any():
            continue
        w_id = weight[a]                                     # constant per id
        vals = {n: np.interp(screens[on], z_traj, quantities[n][traj]) for n in names}
        count[on] += 1.0
        sum_w[on] += w_id
        for n in names:
            sum_wv[n][on] += w_id * vals[n]
            sum_wvv[n][on] += w_id * vals[n] ** 2
            np.maximum.at(peak[n], np.flatnonzero(on), vals[n])
        for q_name, u_name in emit_pairs:
            sum_wqu[(q_name, u_name)][on] += w_id * vals[q_name] * vals[u_name]

    ok = count >= min_cross
    norm = np.where(ok, sum_w, np.nan)                       # NaN-out sparse screens
    out = dict(count=count, mean={}, rms={}, max={}, emit={})
    for n in names:
        mean = sum_wv[n] / norm
        out["mean"][n] = mean
        out["rms"][n] = np.sqrt(np.clip(sum_wvv[n] / norm - mean ** 2, 0.0, None))
        out["max"][n] = np.where(ok, peak[n], np.nan)
    for q_name, u_name in emit_pairs:
        mean_q = sum_wv[q_name] / norm
        mean_u = sum_wv[u_name] / norm
        var_q = sum_wvv[q_name] / norm - mean_q ** 2
        var_u = sum_wvv[u_name] / norm - mean_u ** 2
        cov = sum_wqu[(q_name, u_name)] / norm - mean_q * mean_u
        out["emit"][(q_name, u_name)] = np.sqrt(np.clip(var_q * var_u - cov ** 2, 0.0, None))
    return screens, out
