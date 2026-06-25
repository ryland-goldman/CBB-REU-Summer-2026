"""
Auto-phase the linac 1-4 RF cavities and rewrite the section YAMLs.

`sim/main.py` runs this **before each WarpX linac stage** (one section per call), so the chain
re-derives the crest each run; it also works standalone. It re-derives the frozen RF crest phase
(CREST_PHASE_DEG for sections 2/3/4, PHASE_DEG for section 1) the linac driver consumes, then
writes it back into config/linacN.yaml (comments preserved). Run it whenever an upstream change
shifts the beam and the hardcoded setpoints go stale.

Method: read the section's upstream exit beam exactly as sim/linac1-4.py does (same kinematics,
same iris scrape for section 1), reproduce the driver's arrival-referenced phase convention
phi = -omega*(Z_STRUCT - z_center)/v_beam + base_deg, and RK4-integrate the WHOLE bunch
longitudinally through the real on-axis SLAC quadrature field over a phase scan. The crest is the
base phase maximising the bunch-averaged exit energy (parabolic-refined). Integrating the whole
bunch — not a centroid proxy — is essential: the captured core spans ~140 deg of RF, so its
phase-averaged crest sits ~70 deg from the single-particle crest (validated against a WarpX phase
scan). This is a 1D longitudinal model (no transverse / space-charge back-reaction): exact for the
relativistic sections 2/3/4; for the 150 keV capture section 1 it is the max-energy phase, which a
deliberately off-crest bunching setpoint would differ from.

  python sim/autophase.py            # phase sections 1 2 3 4, rewrite the YAMLs
  python sim/autophase.py 2 3 4      # only the relativistic sections
  python sim/autophase.py --dry-run  # scan + report, write nothing

No WarpX is launched; only the existing logs/diags/ dumps + the GDF field maps are read.
"""

import os
import re
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import importlib.util

import numpy as np
import yaml

from sim.helpers.tools import (
    C_LIGHT as c, E_CHARGE as q_e, M_E as m_e, MC2_KEV, prepare_env)
from sim.helpers.buildfields import onaxis_quadrature_ez, Z_STRUCT, V1KW_KEV

CONFIG = {N: f"config/linac{N}.yaml" for N in (1, 2, 3, 4)}
# The crest the driver applies is base_deg; for section 1 that knob is PHASE_DEG, for 2/3/4 it is
# CREST_PHASE_DEG (with PHASE_DEG kept as the detune-from-crest = 0).
PHASE_KEY = {1: "PHASE_DEG", 2: "CREST_PHASE_DEG", 3: "CREST_PHASE_DEG", 4: "CREST_PHASE_DEG"}


def _num(v):
    """Coerce a YAML scalar to int/float (YAML 1.1 leaves unsigned-exponent forms like 2856.0e6
    as strings); leave non-numeric strings untouched."""
    if isinstance(v, str):
        try:
            return int(v)
        except ValueError:
            try:
                return float(v)
            except ValueError:
                return v
    return v


def _load_driver():
    """Import sim/linac1-4.py (hyphenated, not a normal module name) for its beam loaders so the
    centroid kinematics are byte-identical to what the driver injects."""
    path = os.path.join(os.path.dirname(__file__), "linac1-4.py")
    spec = importlib.util.spec_from_file_location("linac13_driver", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SUBSAMPLE = 2048            # macroparticles integrated per phase (bunch-averaged crest converges)


def _load_bunch(drv, N, p):
    """Section N's injected beam (matching sim/linac1-4.py: section 1 reads the iris-scraped
    injector handoff; 2/3/4 read the previous section's captured-core exit), subsampled for the
    scan. Returns (z [m], u = |gamma*beta| per particle, w, z_center [m], v_beam [m/s],
    ke_mean [keV], scale). `scale` is the field-map amplitude the driver would apply."""
    resample_n = p.get("RESAMPLE_N", 0)                     # match the driver's injected beam size
    if N == 1:
        bunch, v_beam, ke_mean, _ = drv.load_injector_bunch(
            p["MAX_PART"], p["RNG_SEED"], p["Z_INJECT"], p["Z_HANDOFF"], p["COLLIM_Z"],
            resample_n=resample_n)
        scale = float(np.sqrt(p["POWER_MW"] / p["RF_NORM_MW"]))
    else:
        bunch, v_beam, ke_mean, _ = drv.load_warpx_exit_bunch(
            drv.PREV_PARTICLES[N], drv.PREV_LABEL[N], p["MAX_PART"], p["RNG_SEED"], p["Z_INJECT"],
            resample_n=resample_n)
        scale = float(p["FIELD_SCALE"])
    z, w = np.asarray(bunch["z"], float), np.asarray(bunch["w"], float)
    u = np.sqrt(bunch["ux"] ** 2 + bunch["uy"] ** 2 + bunch["uz"] ** 2)   # |gamma*beta|
    if z.size > SUBSAMPLE:                                   # weighted resample → equal weight
        rng = np.random.default_rng(p["RNG_SEED"])
        sel = rng.choice(z.size, SUBSAMPLE, replace=False, p=w / w.sum())
        z, u, w = z[sel], u[sel], np.ones(SUBSAMPLE)
    z_center = float(np.average(z, weights=w))
    return z, u, w, z_center, v_beam, ke_mean, scale


PROBE_LEN = 1.0            # integrate this far into the structure to fix the crest (TW gradient is
                          # phase-uniform; a WarpX scan confirmed the crest at a 0.6 m probe)


def _mean_exit_ke(base_deg, z0, u0, w, z_center, v_beam, scale, omega, zmap, ez1, ez2, z_probe):
    """Bunch-averaged KE [keV] at the `z_probe` plane for a single base phase.

    Integrates the WHOLE bunch (not a centroid proxy — the captured core spans ~140 deg of RF, so
    its phase-averaged crest is far from the single-particle crest). Reproduces the driver phase
    reference phi = -omega*t_in + base (t_in from the centroid, as the driver sets it) and the
    90 deg quadrature sum Ez = scale*[Ez1 cos(wt+phi) + Ez2 cos(wt+phi+pi/2)]; RK4-integrates
    dz/dt = c*u/gamma, du/dt = -(e/m_e c) Ez until the centroid passes z_probe. Field is exactly
    zero outside the map (np.interp left/right=0), so exited particles coast.
    """
    phi = -omega * (Z_STRUCT - z_center) / v_beam + np.deg2rad(base_deg)
    k = q_e / (m_e * c)
    z, u = z0.copy(), u0.copy()
    dt = (2.0 * np.pi / omega) / 100.0                        # resolve the RF cycle
    n_max = int(1.5 * (z_probe - z_center) / v_beam / dt) + 1

    def dstate(zz, uu, t):
        e1 = np.interp(zz, zmap, ez1, left=0.0, right=0.0)
        e2 = np.interp(zz, zmap, ez2, left=0.0, right=0.0)
        ez = scale * (e1 * np.cos(omega * t + phi) + e2 * np.cos(omega * t + phi + np.pi / 2.0))
        g = np.sqrt(1.0 + uu * uu)
        return c * uu / g, -k * ez

    t = 0.0
    for _ in range(n_max):
        dz1, du1 = dstate(z, u, t)
        dz2, du2 = dstate(z + 0.5 * dt * dz1, u + 0.5 * dt * du1, t + 0.5 * dt)
        dz3, du3 = dstate(z + 0.5 * dt * dz2, u + 0.5 * dt * du2, t + 0.5 * dt)
        dz4, du4 = dstate(z + dt * dz3, u + dt * du3, t + dt)
        z = z + (dt / 6.0) * (dz1 + 2 * dz2 + 2 * dz3 + dz4)
        u = u + (dt / 6.0) * (du1 + 2 * du2 + 2 * du3 + du4)
        t += dt
        if z.mean() >= z_probe:
            break

    ke = (np.sqrt(1.0 + u * u) - 1.0) * MC2_KEV
    return float(np.average(ke, weights=w))


def find_crest(z0, u0, w, z_center, v_beam, ke_mean, scale, omega, zmap, ez1, ez2):
    """Crest base phase [deg, in [0,360)] and its bunch-averaged gain [keV] at the probe plane: a
    2 deg coarse scan (light subsample), a 0.1 deg fine scan about the peak, then a parabolic
    refine. The probe plane caps the integration short of the full structure for speed."""
    z_probe = min(zmap[-1], Z_STRUCT + PROBE_LEN)
    coarse_sub = slice(0, min(512, z0.size))                  # lighter bunch for the wide scan

    def ke_at(phases, zz, uu, ww):
        return np.array([_mean_exit_ke(b, zz, uu, ww, z_center, v_beam, scale, omega,
                                       zmap, ez1, ez2, z_probe) for b in np.atleast_1d(phases)])

    coarse = np.arange(0.0, 360.0, 2.0)
    c0 = coarse[int(np.argmax(ke_at(coarse, z0[coarse_sub], u0[coarse_sub], w[coarse_sub])))]
    fine = c0 + np.arange(-3.0, 3.0 + 1e-9, 0.1)
    g = ke_at(fine, z0, u0, w)
    i = int(np.argmax(g))
    if 0 < i < len(fine) - 1:                                 # parabolic vertex of the top 3
        y0, y1, y2 = g[i - 1], g[i], g[i + 1]
        denom = y0 - 2 * y1 + y2
        d = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
        crest = fine[i] + d * (fine[1] - fine[0])
    else:
        crest = fine[i]
    # Full-structure ΔE at the crest (one integration to the map exit) for an accurate report.
    full = _mean_exit_ke(crest, z0, u0, w, z_center, v_beam, scale, omega, zmap, ez1, ez2,
                         zmap[-1] + 0.05)
    return float(crest % 360.0), float(full - ke_mean)


def set_yaml_param(path, key, value_str):
    """Replace the value of `key` in a YAML file, preserving its inline comment. Returns the old
    value string."""
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
    sections = [int(a) for a in args] if args else [1, 2, 3, 4]
    if any(N not in (1, 2, 3, 4) for N in sections):
        sys.exit("usage: python sim/autophase.py [1] [2] [3] [4] [--dry-run]")

    prepare_env()                                             # chdir repo root; no WarpX import
    drv = _load_driver()
    zmap_local, ez1, ez2 = onaxis_quadrature_ez()
    zmap = Z_STRUCT + zmap_local                              # lab z of the on-axis field samples

    print(f"On-axis SLAC map: {zmap.size} samples, entrance z={Z_STRUCT*1e3:.0f} mm, "
          f"exit z={zmap[-1]*1e3:.0f} mm, 1-kW int|Ez|dz={V1KW_KEV:.1f} keV\n")

    results = []
    for N in sorted(set(sections)):
        with open(CONFIG[N]) as fh:
            p = yaml.safe_load(fh)["params"]
        p = {k: _num(v) for k, v in p.items()}              # YAML 1.1 leaves 2856.0e6 a string
        omega = 2.0 * np.pi * p["F_RF"]
        print(f"── Section {N} ─────────────────────────────────────────────")
        try:
            z0, u0, w, z_center, v_beam, ke_mean, scale = _load_bunch(drv, N, p)
        except Exception as e:
            print(f"  SKIP: cannot read section {N}'s upstream beam — {e}\n", flush=True)
            continue

        crest, gain = find_crest(z0, u0, w, z_center, v_beam, ke_mean, scale, omega, zmap, ez1, ez2)
        key = PHASE_KEY[N]
        old = float(p.get(key, 0.0))
        print(f"  scale={scale:.4g}, ⟨KE⟩_in={ke_mean/1e3:.3f} MeV, β_in={v_beam/c:.4f}")
        print(f"  crest {key} = {crest:.4f}°  (was {old:.4f}°, Δ {crest - old:+.4f}°); "
              f"on-crest gain {gain/1e3:.3f} MeV", flush=True)
        if N != 1 and p.get("DE_TARGET_MEV"):                # info only — field scale is a separate knob
            print(f"  note: hits {gain/1e3:.2f} MeV vs DE_TARGET {p['DE_TARGET_MEV']} MeV — "
                  f"re-derive FIELD_SCALE separately if ΔE is off.")
        results.append((N, key, crest))
        print()

    if dry:
        print("--dry-run: YAMLs unchanged.")
        return
    for N, key, crest in results:
        old = set_yaml_param(CONFIG[N], key, f"{crest:.4f}")
        print(f"Wrote {CONFIG[N]}: {key} {old} → {crest:.4f}")
    if results:
        print("\nSetpoints updated. Re-run the affected sections to propagate the new crest.")


if __name__ == "__main__":
    main()
