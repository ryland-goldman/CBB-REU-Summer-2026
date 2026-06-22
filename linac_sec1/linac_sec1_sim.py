"""
SLAC Linac Section 1 in WarpX (RZ): capture the injector beam in a 3 m 2π/3
traveling-wave accelerating structure with self-consistent space charge. Reads the
injector handoff beam, applies the iris scrape, tracks through the two quadrature RF
maps, and writes openPMD diagnostics.

Drives lume-warpx from linac_sec1/linac_sec1.yaml (which holds every constant); this module
reads those back, imports the captured beam via WarpX(initial_particles=...), and overrides only
runtime-computed values (the two RF time functions, step count, dt, diagnostic period). See
linac_sec1/README.md for physics, parameters, and gotchas.
"""

import os
import json
import shutil
import numpy as np
from openpmd_viewer import OpenPMDTimeSeries

from pipeline.collimator import pipe_violator_ids, survivor_mask
from pipeline.constants import C_LIGHT as c, E_CHARGE as q_e, MC2_EV
from .build_linac_sec1_field import Z_STRUCT, RMAX, BORE_R, V1KW_KEV
from . import DEFAULT_OUTDIR

MC2_KEV = MC2_EV / 1e3                  # electron rest energy [keV] ≈ 511
HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG = os.path.join(HERE, "linac_sec1.yaml")
INJECTOR_DIAG = "injector/diags/main/particles"


def load_injector_bunch(max_part, rng_seed, z_inject, z_handoff, collim_z):
    """Import the injector beam at the z ≈ z_handoff plane and shift it to entry.

    Selects the dump whose bunch ⟨z⟩ is nearest z_handoff, applies the multi-plane iris scrape,
    and returns (dict [γβ momenta], v_beam, mean KE [keV], inj summary).
    """
    ts = OpenPMDTimeSeries(INJECTOR_DIAG)
    if len(ts.iterations) == 0:
        raise RuntimeError(f"{INJECTOR_DIAG} has no iterations — did the injector stage run?")

    # Find the well-populated dump nearest the handoff plane (the n ≥ 0.8·nmax gate avoids a
    # depleted late dump that happens to sit near 2.03 m).
    min_count = max(50, max_part // 50)
    recs = []
    for it in ts.iterations:
        z, w = ts.get_particle(["z", "w"], species="electrons", iteration=it)
        if len(z) < min_count:
            continue
        recs.append((it, len(z), float(np.average(z, weights=w))))
    if not recs:
        raise RuntimeError(f"{INJECTOR_DIAG}: no snapshot with ≥{min_count} macroparticles")
    nmax = max(n for _, n, _ in recs)
    cands = [(it, zm) for it, n, zm in recs if n >= 0.8 * nmax]
    it_handoff, zm_handoff = min(cands, key=lambda t: abs(t[1] - z_handoff))
    if abs(zm_handoff - z_handoff) > 0.02:
        print(f"  WARNING: nearest injector dump to the {z_handoff*1e3:.0f} mm handoff is at "
              f"⟨z⟩={zm_handoff*1e3:.0f} mm — the handoff diagnostic may be too coarse.", flush=True)

    idh, x, y, z, ux, uy, uz, w = ts.get_particle(
        ["id", "x", "y", "z", "ux", "uy", "uz", "w"], species="electrons", iteration=it_handoff)
    q_inj = float(w.sum()) * q_e                        # TRUE injected charge (all r), capture denominator

    # Iris/pipe collimation: multi-plane id scrape, NOT a single 2.03 m cut. See pipeline.collimator.
    scan_iters = [it for it, _n, zm in recs if (collim_z - 0.05) <= zm <= (z_handoff + 0.03)]
    violators = pipe_violator_ids(ts, scan_iters, RMAX, collim_z)
    keep = survivor_mask(idh, violators)
    idh, x, y, z, ux, uy, uz, w = (a[keep] for a in (idh, x, y, z, ux, uy, uz, w))
    r = np.hypot(x, y)
    q_dom = float(w.sum()) * q_e                         # in-iris survivors
    q_bore = float(w[r <= BORE_R].sum()) * q_e           # of those, within the RF bore

    if z.size > max_part:
        rng = np.random.default_rng(rng_seed)
        sel = rng.choice(z.size, max_part, replace=False)
        scale_w = z.size / max_part
        x, y, z, ux, uy, uz, w = (a[sel] for a in (x, y, z, ux, uy, uz, w))
        w = w * scale_w
    z = z - z.min() + z_inject                          # bunch tail (smallest z) → z_inject

    gb = np.sqrt(1.0 + ux**2 + uy**2 + uz**2)           # γ (ux/uy/uz are γβ)
    v_beam = float(np.average(uz / gb, weights=w) * c)
    ke_mean = float(np.average(gb - 1.0, weights=w) * MC2_KEV)
    sz = float(np.sqrt(np.average((z - np.average(z, weights=w)) ** 2, weights=w)))
    rmax = float(np.hypot(x, y).max())
    inj = dict(it_handoff=int(it_handoff), z_handoff_m=float(zm_handoff),
               n_injected=int(z.size), q_injected_C=q_inj, q_in_domain_C=q_dom,
               q_in_bore_C=q_bore, z_inject_mean_m=float(np.average(z, weights=w)),
               rmax_m=rmax, sigma_z_m=sz, ke_mean_keV=ke_mean)
    print(f"Injected {z.size} macroparticles (iris survivors) from injector handoff "
          f"(iter {it_handoff}, ⟨z⟩={zm_handoff*1e3:.0f} mm); ⟨KE⟩ {ke_mean:.1f} keV, "
          f"σ_z {sz*1e3:.2f} mm, r_max {rmax*1e3:.2f} mm, true-injected q {q_inj*1e9:.4f} nC", flush=True)
    print(f"  multi-plane iris scrape at the {RMAX*1e3:.3f} mm aperture (z≥{collim_z*1e3:.0f} mm, "
          f"{len(scan_iters)} planes): {q_dom/q_inj*100:.1f}% survives the pipe, "
          f"{q_bore/q_inj*100:.1f}% within the {BORE_R*1e3:.2f} mm RF bore.", flush=True)
    return dict(x=x, y=y, z=z, ux=ux, uy=uy, uz=uz, w=w), v_beam, ke_mean, inj


def main():
    from warpx import WarpX
    from pmd_beamphysics import ParticleGroup

    w = WarpX(input_file=CONFIG, path="linac_sec1")
    NR, NZ = w.get("grid/number_of_cells")
    _, ZMAX = w.get("grid/upper_bound")
    outdir = w.get("diagnostics/0/write_dir") or DEFAULT_OUTDIR
    p = w.get("params")
    omega = 2.0 * np.pi * p["F_RF"]

    # The last applied field must load_E or picmi forces the global E_ext style to "none" (RF dark).
    # Both RF maps load_E here, so this is always satisfied — kept so a future reorder fails loudly.
    assert w.get("fields")[-1].get("load_E"), \
        "linac_sec1.yaml: last applied field must have load_E:true"

    if os.path.isdir(outdir):                           # fresh diags (WarpX appends per dump)
        shutil.rmtree(outdir)

    bunch, v_beam, ke_mean, inj = load_injector_bunch(
        p["MAX_PART"], p["RNG_SEED"], p["Z_INJECT"], p["Z_HANDOFF"], p["COLLIM_Z"])
    z_center = float(np.average(bunch["z"], weights=bunch["w"]))
    # Persist the true injected charge so the plotter reports capture against it. WarpX makes
    # outdir at first diag; make it now for this sidecar.
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "injection_summary.json"), "w") as fh:
        json.dump(inj, fh, indent=2)

    pg = ParticleGroup(data=dict(
        x=bunch["x"], y=bunch["y"], z=bunch["z"],
        px=bunch["ux"] * MC2_EV, py=bunch["uy"] * MC2_EV, pz=bunch["uz"] * MC2_EV,
        t=np.zeros(bunch["x"].size), weight=bunch["w"] * q_e,
        status=np.ones(bunch["x"].size, dtype=np.int64), species="electron"))
    w.initial_particles = pg

    # ── RF amplitude + phase (field2 is the 90° quadrature half) ──────────────
    scale = float(np.sqrt(p["POWER_MW"] / p["RF_NORM_MW"]))
    t_in = (Z_STRUCT - z_center) / v_beam
    phi = -omega * t_in + np.deg2rad(p["PHASE_DEG"])
    phi2 = phi + np.pi / 2.0
    e1 = f"{scale:.8e}*cos({omega:.10e}*t + ({phi:.8e}))"
    b1 = f"{scale:.8e}*sin({omega:.10e}*t + ({phi:.8e}))"
    e2 = f"{scale:.8e}*cos({omega:.10e}*t + ({phi2:.8e}))"
    b2 = f"{scale:.8e}*sin({omega:.10e}*t + ({phi2:.8e}))"
    print(f"Case: P={p['POWER_MW']:g} MW (scale={scale:.1f}), phase_off={p['PHASE_DEG']:g}°, "
          f"f_RF={p['F_RF']/1e6:.0f} MHz → {outdir}/", flush=True)

    # ── Time step / duration: segmented transit to a stop plane short of ZMAX so the bunch
    # finishes in the field-free exit drift, NOT at the absorbing wall (empty domain aborts MLMG).
    dt = p["CFL"] * (ZMAX / NZ) / v_beam
    beta_in = v_beam / c
    gain_keV = scale * V1KW_KEV                         # ≈ on-crest gain [keV] (1-kW V1kW from the maps)
    beta_hi = float(np.sqrt(1.0 - 1.0 / (1.0 + (ke_mean + gain_keV) / MC2_KEV) ** 2))
    z_end = ZMAX - 0.20                                 # stop in the drift, clear of the wall
    L_cap = 0.40                                        # capture length over which β: in→hi
    beta_cap = 0.5 * (beta_in + beta_hi)
    transit = ((Z_STRUCT - z_center) / v_beam
               + min(L_cap, max(0.0, z_end - Z_STRUCT)) / (beta_cap * c)
               + max(0.0, z_end - Z_STRUCT - L_cap) / (beta_hi * c))
    n_steps = p["MAX_STEPS"] if p.get("MAX_STEPS") else int(p["TRANSIT_MARGIN"] * transit / dt)
    period = max(1, n_steps // p["N_DIAGS"])
    print(f"dt = {dt:.3e} s, max_steps = {n_steps}, β_in={beta_in:.3f}, "
          f"est. on-crest gain {gain_keV/1e3:.1f} MeV", flush=True)

    w.update({
        "simulation/max_steps": n_steps,
        "simulation/time_step_size": dt,
        "diagnostics/0/period": period,
        "diagnostics/0/write_dir": outdir,
        "fields/0/warpx_E_time_function": e1, "fields/0/warpx_B_time_function": b1,
        "fields/1/warpx_E_time_function": e2, "fields/1/warpx_B_time_function": b2,
    })

    print(f"\nRunning {n_steps} steps (diag every {period}) → {outdir}/")
    w.run(progress="linac_sec1")
    print("\nDone.")


if __name__ == "__main__":
    main()
