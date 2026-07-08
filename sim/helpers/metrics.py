"""Beam-moment helpers (pure numpy). Emittance conventions: uq is the openPMD
u = gamma*beta (already includes gamma); moments are returned RAW and the caller
applies unit scaling.
"""

import numpy as np

# Exit beam-quality keys (NaN-filled on the no-survivor case so the summary schema is stable).
BEAM_QUALITY_KEYS = ("eps_n_x_m", "eps_n_y_m", "sigma_E_mev", "sigma_E_rel",
                     "sigma_x_mm", "sigma_y_mm")


def beam_quality(pg):
    """sigma_E_rel = std(KE)/mean(KE); std(KE)==std(total energy) since rest mass is constant.
    All-NaN when `pg` is None or empty so callers never KeyError.
    """
    if pg is None or getattr(pg, "n_particle", 0) == 0:
        return {k: float("nan") for k in BEAM_QUALITY_KEYS}
    return dict(
        eps_n_x_m=float(pg.norm_emit_x),
        eps_n_y_m=float(pg.norm_emit_y),
        sigma_E_mev=float(pg["sigma_energy"]) / 1e6,
        sigma_E_rel=float(pg["sigma_energy"]) / float(pg["mean_kinetic_energy"]),
        sigma_x_mm=float(pg["sigma_x"]) * 1e3,
        sigma_y_mm=float(pg["sigma_y"]) * 1e3,
    )


def _group_bounds(sorted_ids):
    edges = np.flatnonzero(np.r_[True, sorted_ids[1:] != sorted_ids[:-1], True])
    return zip(edges[:-1], edges[1:])


def screen_profile(ids, z, weight, quantities, emit_pairs=(),
                   n_screen=80, min_cross=20, z_range=None):
    """Charge-weighted beam moments on fixed-z virtual screens, interpolated per-id along z
    (assumes forward, monotonic-in-z motion); NaN on screens with < min_cross crossings.
    """
    ids = np.asarray(ids)
    z = np.asarray(z, float)
    weight = np.asarray(weight, float)
    names = list(quantities)
    quantities = {n: np.asarray(quantities[n], float) for n in names}

    nan = np.full(n_screen, np.nan)
    if z.size == 0:
        empty = dict(count=np.zeros(n_screen),
                     charge=nan.copy(),
                     mean={n: nan.copy() for n in names},
                     rms={n: nan.copy() for n in names},
                     max={n: nan.copy() for n in names},
                     emit={p: nan.copy() for p in emit_pairs})
        return np.linspace(0.0, 1.0, n_screen), empty

    lo, hi = z_range if z_range else (z.min(), z.max())
    screens = np.linspace(lo, hi, n_screen)

    count = np.zeros(n_screen)
    sum_w = np.zeros(n_screen)
    sum_wv = {n: np.zeros(n_screen) for n in names}
    sum_wvv = {n: np.zeros(n_screen) for n in names}
    peak = {n: np.full(n_screen, -np.inf) for n in names}
    sum_wqu = {p: np.zeros(n_screen) for p in emit_pairs}

    order = np.argsort(ids, kind="stable")
    ids, z, weight = ids[order], z[order], weight[order]
    quantities = {n: q[order] for n, q in quantities.items()}

    for a, b in _group_bounds(ids):
        if b - a < 2:                                        # need >=2 dumps to interpolate
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
    out = dict(count=count, charge=np.where(ok, sum_w, np.nan),
               mean={}, rms={}, max={}, emit={})
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
