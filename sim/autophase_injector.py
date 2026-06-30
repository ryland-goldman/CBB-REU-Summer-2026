"""
Auto-phase the injector prebuncher cavities and rewrite config/injector.yaml.

The analog of sim/autophase.py (WarpX linac) / sim/autophase_impact.py (Impact-T) for the injector
RF -- but the prebunchers are VELOCITY BUNCHERS, not accelerators, so the objective is the opposite:
not the crest (max bunch energy) but the phase that MINIMISES the bunch length at the linac handoff
plane. The driver lands each cavity centroid on the RF zero-crossing analytically (`PHASE: zc`, net
mean kick 0); the autophased knobs are the per-cavity offsets from that base (PREB1_PHI_OFF /
PREB2_PHI_OFF), which set where on the bunching slope the cavity sits and thus where the velocity
focus lands. Re-run it whenever the gun beam shifts and the hardcoded offsets go stale.

Method: import the gun exit beam exactly as sim/injector.py does (same downsample, same
tail->Z_INJECT shift, same centroid kinematics), reproduce the driver's per-cavity phase
construction byte-for-byte by calling sim/injector.cavity_drive (scale = sqrt(1e3 Q_L P /
2 pi f), the centroid-arrival t_gap, the +pi reversed-Preb-2 phase, and Preb-2's two-segment arrival
through the analytic post-Preb-1 speed), and RK4-integrate the WHOLE bunch longitudinally through the
two on-axis prebuncher gaps to the INJ_Z_HANDOFF plane. The on-axis field is the same 1-J map the
driver loads (Bphi vanishes on axis, so only Ez acts longitudinally). The objective is the
weighted-RMS arrival-time spread sigma_t = sigma_z / v_beam at the handoff snapshot; the optimum
offset minimises it. Two cavities are optimised by coordinate descent (each a coarse -> fine ->
parabolic min, two passes); Preb-2's arrival depends on Preb-1's offset, so the phases are rebuilt at
every evaluation. This is a 1D longitudinal model (no transverse / space-charge back-reaction) -- the
focus position the buncher exists to set is a longitudinal quantity, so it captures the tuning knob.

  python sim/autophase_injector.py            # phase both prebunchers, rewrite the YAML
  python sim/autophase_injector.py 1          # only Prebuncher 1's offset
  python sim/autophase_injector.py --dry-run  # scan + report, write nothing

No WarpX is launched; only the existing logs/diags/gun dump + the prebuncher GDF map are read.
"""

import os
import re
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np
import yaml

from sim.helpers.tools import C_LIGHT as c, E_CHARGE as q_e, M_E as m_e, MC2_KEV, prepare_env
from sim.helpers.buildfields import (
    Z_GAP_CENTER_1, Z_GAP_CENTER_2, V1J_KEV, INJ_Z_HANDOFF as Z_HANDOFF, _load_prebuncher_map)
from sim import injector as drv          # driver loaders + cavity_drive (warpx-free at import)

CONFIG = "config/injector.yaml"
# Which cavity each CLI section selects, and its offset knob in the YAML.
PHASE_KEY = {1: "PREB1_PHI_OFF", 2: "PREB2_PHI_OFF"}

SUBSAMPLE = 1536            # macroparticles integrated per phase (bunch-length objective converges)
COARSE_HALF_DEG = 60.0     # +-window about the current offset for the wide scan
COARSE_STEP_DEG = 10.0
FINE_HALF_DEG = 10.0       # +-window about the coarse min
FINE_STEP_DEG = 1.0
DESCENT_PASSES = 2         # coordinate-descent sweeps over (off1, off2)


def _load_bunch(p):
    """Gun exit beam as sim/injector.py injects it, subsampled for the scan. Returns (z [m],
    uz = longitudinal gamma*beta per particle, uperp2 = ux^2+uy^2 (frozen transverse gamma*beta^2),
    w, z_centroid [m], v_beam [m/s], ke_mean [keV]). uz and uperp are carried separately so the Ez
    kick acts only on uz while gamma keeps the transverse momentum -- matching the uz/gamma v_beam in
    helpers.loadparticles.beam_kinematics that converts sigma_z to sigma_t."""
    bunch, v_beam, ke_mean, z_centroid = drv.load_gun_bunch(
        p["MAX_PART"], p["RNG_SEED"], p["Z_INJECT"])
    z, w = np.asarray(bunch["z"], float), np.asarray(bunch["w"], float)
    uz = np.asarray(bunch["uz"], float)                                   # longitudinal gamma*beta
    uperp2 = np.asarray(bunch["ux"], float) ** 2 + np.asarray(bunch["uy"], float) ** 2
    if z.size > SUBSAMPLE:                                   # weighted resample -> equal weight
        rng = np.random.default_rng(p["RNG_SEED"])
        sel = rng.choice(z.size, SUBSAMPLE, replace=False, p=w / w.sum())
        z, uz, uperp2, w = z[sel], uz[sel], uperp2[sel], np.ones(SUBSAMPLE)
    return z, uz, uperp2, w, z_centroid, v_beam, ke_mean


def _cavity_phases(p, v_beam, ke_mean, z_centroid, off1, off2, omega):
    """Reproduce sim/injector.py's two cavity_drive calls for given offsets, returning
    (scale1, phi1, scale2, phi2). Preb-2's arrival time runs through the analytic post-Preb-1 speed
    (its t_gap, hence phi2, depends on off1 via Preb-1's mean kick), exactly as the driver bakes it."""
    PHASE, F_RF = p["PHASE"], p["F_RF"]
    base = np.pi / 2.0 if PHASE == "zc" else np.pi
    _, _, scale1, phi1, t_gap1 = drv.cavity_drive(
        p["PREB1_KW"], p["Q_L_1"], F_RF, Z_GAP_CENTER_1, v_beam, off1, PHASE, omega,
        z_ref=z_centroid)
    kick1 = -np.cos(base + np.radians(off1)) * scale1 * V1J_KEV
    ke_after1 = max(ke_mean + (kick1 if p["PREB1_KW"] > 0 else 0.0), 1.0)
    v_after_preb1 = c * np.sqrt(1.0 - 1.0 / (1.0 + ke_after1 / MC2_KEV) ** 2)
    rev_phase = p["PREB2_REV_PHASE"] if p["PREB2_REVERSED"] else 0.0
    _, _, scale2, phi2, _ = drv.cavity_drive(
        p["PREB2_KW"], p["Q_L_2"], F_RF, Z_GAP_CENTER_2, v_after_preb1, off2, PHASE, omega,
        t_offset=t_gap1, z_ref=Z_GAP_CENTER_1, rev_phase=rev_phase)
    return scale1, phi1, scale2, phi2


def _sigma_t_ps(z0, uz0, uperp2, w, scale1, phi1, scale2, phi2, ez_axis, zmap1, zmap2, omega,
                preb1_on, preb2_on, v_beam):
    """Weighted-RMS arrival-time spread [ps] at the Z_HANDOFF plane for one (phi1, phi2) pair.

    RK4-integrates the WHOLE bunch (dz/dt = c*uz/gamma, duz/dt = -(e/m_e c) Ez) through both on-axis
    prebuncher gaps until the centroid reaches Z_HANDOFF, then returns sigma_z/v_beam at that
    snapshot. The longitudinal Ez kicks only uz; gamma = sqrt(1 + uz^2 + uperp2) carries the frozen
    transverse momentum, so v_z = c*uz/gamma matches the uz/gamma convention of beam_kinematics.
    Ez = scale1 Ez_map(z-gap1) cos(wt+phi1) + scale2 Ez_map(z-gap2) cos(wt+phi2); the field is exactly
    zero outside each map (np.interp left/right=0) so the bunch coasts between cavities and after the
    second gap. Bphi vanishes on axis, so Ez is the only longitudinal field."""
    k = q_e / (m_e * c)
    z, uz = z0.copy(), uz0.copy()
    dt = (2.0 * np.pi / omega) / 100.0                        # resolve the RF cycle
    n_max = int(1.6 * (Z_HANDOFF - z.mean()) / v_beam / dt) + 1

    def ez(zz, t):
        e = np.zeros_like(zz)
        if preb1_on:
            e += scale1 * np.interp(zz, zmap1, ez_axis, left=0.0, right=0.0) * np.cos(omega * t + phi1)
        if preb2_on:
            e += scale2 * np.interp(zz, zmap2, ez_axis, left=0.0, right=0.0) * np.cos(omega * t + phi2)
        return e

    def dstate(zz, uuz, t):
        g = np.sqrt(1.0 + uuz * uuz + uperp2)                 # gamma keeps the transverse momentum
        return c * uuz / g, -k * ez(zz, t)                    # Ez kicks only the longitudinal uz

    t = 0.0
    for _ in range(n_max):
        dz1, du1 = dstate(z, uz, t)
        dz2, du2 = dstate(z + 0.5 * dt * dz1, uz + 0.5 * dt * du1, t + 0.5 * dt)
        dz3, du3 = dstate(z + 0.5 * dt * dz2, uz + 0.5 * dt * du2, t + 0.5 * dt)
        dz4, du4 = dstate(z + dt * dz3, uz + dt * du3, t + dt)
        z = z + (dt / 6.0) * (dz1 + 2 * dz2 + 2 * dz3 + dz4)
        uz = uz + (dt / 6.0) * (du1 + 2 * du2 + 2 * du3 + du4)
        t += dt
        if np.average(z, weights=w) >= Z_HANDOFF:
            break

    zbar = np.average(z, weights=w)
    sigma_z = float(np.sqrt(np.average((z - zbar) ** 2, weights=w)))
    return sigma_z / v_beam * 1e12


def _min_offset(objective, x0):
    """Offset [deg] minimising `objective` near x0: a coarse scan, a fine scan about the min, then a
    parabolic refine of the bottom 3 (the min mirror of sim/autophase.find_crest's argmax)."""
    coarse = x0 + np.arange(-COARSE_HALF_DEG, COARSE_HALF_DEG + 1e-9, COARSE_STEP_DEG)
    c0 = coarse[int(np.argmin([objective(x) for x in coarse]))]
    fine = c0 + np.arange(-FINE_HALF_DEG, FINE_HALF_DEG + 1e-9, FINE_STEP_DEG)
    g = np.array([objective(x) for x in fine])
    i = int(np.argmin(g))
    if 0 < i < len(fine) - 1:                                 # parabolic vertex of the bottom 3
        y0, y1, y2 = g[i - 1], g[i], g[i + 1]
        denom = y0 - 2 * y1 + y2
        d = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
        return float(fine[i] + d * (fine[1] - fine[0]))
    return float(fine[i])


def set_yaml_param(path, key, value_str):
    """Replace the value of `key` in a YAML file, preserving its inline comment (shared with the other
    autophase scripts' writeback). Returns the old value string."""
    with open(path) as fh:
        txt = fh.read()
    pat = re.compile(rf"^(\s*{re.escape(key)}:\s*)(\S+)(.*)$", re.M)
    m = pat.search(txt)
    if not m:
        raise KeyError(f"{key} not found in {path}")
    new = pat.sub(lambda mm: mm.group(1) + value_str + mm.group(3), txt, count=1)
    with open(path, "w") as fh:
        fh.write(new)
    return m.group(2)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    dry = "--dry-run" in sys.argv or "-n" in sys.argv
    cavities = sorted({int(a) for a in args}) if args else [1, 2]
    if any(N not in (1, 2) for N in cavities):
        sys.exit("usage: python sim/autophase_injector.py [1] [2] [--dry-run]")

    prepare_env()                                             # chdir out_root (else repo root); no WarpX import
    with open(CONFIG) as fh:
        p = yaml.safe_load(fh)["params"]
    omega = 2.0 * np.pi * p["F_RF"]
    preb1_on = p["PREB1_KW"] > 0
    preb2_on = p["PREB2_KW"] > 0

    # On-axis Ez (r=0 row of the 1-J map) placed at each gap, exactly as the driver loads it.
    r, z_native, _Er, Ez, _Bphi = _load_prebuncher_map()
    ez_axis = Ez[0]
    zmap1 = Z_GAP_CENTER_1 + z_native
    zmap2 = Z_GAP_CENTER_2 + z_native
    print(f"Prebuncher on-axis map: {ez_axis.size} samples, gaps at z = {Z_GAP_CENTER_1*1e3:.0f}, "
          f"{Z_GAP_CENTER_2*1e3:.0f} mm; handoff plane z = {Z_HANDOFF*1e3:.0f} mm\n")

    try:
        z0, uz0, uperp2, w, z_centroid, v_beam, ke_mean = _load_bunch(p)
    except Exception as e:
        sys.exit(f"cannot read the gun exit beam ({drv.GUN_DIAG}) -- run sim/gun.py first: {e}")

    sigma_t0 = float(np.sqrt(np.average((z0 - np.average(z0, weights=w)) ** 2, weights=w)))
    print(f"Gun beam: {z0.size} parts (subsampled), ⟨KE⟩ {ke_mean:.1f} keV, β {v_beam/c:.4f}, "
          f"injected σ_z {sigma_t0*1e3:.2f} mm\n")

    def sigma_t(off1, off2):
        s1, ph1, s2, ph2 = _cavity_phases(p, v_beam, ke_mean, z_centroid, off1, off2, omega)
        return _sigma_t_ps(z0, uz0, uperp2, w, s1, ph1, s2, ph2, ez_axis, zmap1, zmap2, omega,
                           preb1_on, preb2_on, v_beam)

    off1 = float(p["PREB1_PHI_OFF"])
    off2 = float(p["PREB2_PHI_OFF"])
    st_in = sigma_t(off1, off2)
    print(f"Current offsets: PREB1_PHI_OFF={off1:.4f}°, PREB2_PHI_OFF={off2:.4f}°  → "
          f"σ_t at handoff {st_in:.3f} ps\n")

    # Coordinate descent: Preb-2's arrival depends on off1, so re-optimise each axis with the other
    # pinned to its current best, sweeping twice for the (mild) cross-coupling to settle.
    do1, do2 = (1 in cavities) and preb1_on, (2 in cavities) and preb2_on
    for _ in range(DESCENT_PASSES):
        if do1:
            off1 = _min_offset(lambda x: sigma_t(x, off2), off1)
        if do2:
            off2 = _min_offset(lambda x: sigma_t(off1, x), off2)
        if not (do1 and do2):                                # single free axis -> one pass suffices
            break

    st_out = sigma_t(off1, off2)
    print(f"── Prebunchers ─────────────────────────────────────────────")
    results = []
    for N, (on, val, key) in {1: (preb1_on, off1, "PREB1_PHI_OFF"),
                              2: (preb2_on, off2, "PREB2_PHI_OFF")}.items():
        if N not in cavities:
            continue
        if not on:
            print(f"  Prebuncher {N}: off (PREB{N}_KW=0) — offset not phased")
            continue
        old = float(p[key])
        print(f"  {key} = {val:.4f}°  (was {old:.4f}°, Δ {val - old:+.4f}°)")
        results.append((key, val))
    print(f"  σ_t at handoff {st_in:.3f} → {st_out:.3f} ps "
          f"({100*(st_in-st_out)/st_in:+.1f}% bunching)\n", flush=True)

    if dry:
        print("--dry-run: YAML unchanged.")
        return
    for key, val in results:
        old = set_yaml_param(CONFIG, key, f"{val:.4f}")
        print(f"Wrote {CONFIG}: {key} {old} → {val:.4f}")
    if results:
        print("\nSetpoints updated. Re-run sim/injector.py to propagate the new offsets.")


if __name__ == "__main__":
    main()
