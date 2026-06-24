"""
SLAC / Cornell Linac sections 1-3 in WarpX (RZ), merged into ONE parametrized driver.

  Section 1 (capture): import the injector handoff beam at the z ≈ Z_HANDOFF plane, apply the
    multi-plane 9.547 mm iris scrape, and capture it in the 3 m 2π/3 traveling-wave SLAC
    structure with self-consistent space charge. RF amplitude = sqrt(POWER_MW/RF_NORM_MW),
    phase referenced to bunch arrival (PHASE_DEG offset, default on-crest for the slipping beam).
  Sections 2, 3 (accelerate): import the previous section's captured-core exit beam and accelerate
    the relativistic core through the SAME reused SLAC quadrature maps, scaled by a FROZEN
    FIELD_SCALE and phased to a FROZEN CREST_PHASE_DEG (the old runtime crest-finding +
    ΔE-target field-scale loop is dropped — the setpoints were derived once and hardcoded in
    the section yaml).

Run as:  python sim/linac1-3.py <N>   with N in {1, 2, 3}.

Drives lume-warpx from config/linacN.yaml (which holds every constant); this module reads those
back, imports the upstream beam via WarpX(initial_particles=...), and overrides only the
runtime-computed values (the two quadrature RF time functions, step count, dt, diagnostic
period). The two SLAC quadrature maps are shared across all three sections. See docs/linac1-3.md
for physics, the captured-core cut, the frozen setpoints, lab-z chaining, and gotchas.

main() runs ONLY the simulation; sim/plot/linac1-3.py produces the figures.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import json
import shutil

import numpy as np

from sim.helpers.tools import (
    C_LIGHT as c, E_CHARGE as q_e, MC2_KEV, prepare_env, rf_time_functions)
from sim.helpers.loadparticles import (
    open_particle_series, make_particle_group, downsample, beam_kinematics,
    load_warpx_exit_bunch, upstream_exit_lab_z, pipe_violator_ids, survivor_mask)
from sim.helpers.buildfields import build_linac_slac_fields, Z_STRUCT, RMAX, BORE_R, V1KW_KEV

# Section 1 reads the injector handoff; sections 2/3 read the previous section's exit. All
# paths are repo-root-relative (prepare_env() chdir's to the repo root).
INJECTOR_DIAG = "logs/diags/injector/main/particles"
PREV_PARTICLES = {2: "logs/diags/linac1-3/sec1/main/particles",
                  3: "logs/diags/linac1-3/sec2/main/particles"}
PREV_SUMMARY = {2: "logs/diags/linac1-3/sec1/main/injection_summary.json",
                3: "logs/diags/linac1-3/sec2/main/injection_summary.json"}
PREV_LABEL = {2: "linac1-3/sec1", 3: "linac1-3/sec2"}


def load_injector_bunch(max_part, rng_seed, z_inject, z_handoff, collim_z):
    """Import the injector beam at the z ≈ z_handoff plane and shift it to entry (section 1).

    Selects the dump whose bunch ⟨z⟩ is nearest z_handoff, applies the multi-plane iris scrape,
    and returns (dict [γβ momenta], v_beam, mean KE [keV], inj summary).
    """
    ts = open_particle_series(INJECTOR_DIAG, "injector")

    # Find the well-populated dump nearest the handoff plane (the n ≥ 0.8·nmax gate avoids a
    # depleted late dump that happens to sit near z_handoff).
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
    q_inj = float(w.sum()) * q_e                          # TRUE injected charge (all r), capture denominator

    # Iris/pipe collimation: multi-plane id scrape, NOT a single z_handoff cut. See loadparticles.
    scan_iters = [it for it, _n, zm in recs if (collim_z - 0.05) <= zm <= (z_handoff + 0.03)]
    violators = pipe_violator_ids(ts, scan_iters, RMAX, collim_z)
    keep = survivor_mask(idh, violators)
    idh, x, y, z, ux, uy, uz, w = (a[keep] for a in (idh, x, y, z, ux, uy, uz, w))
    r = np.hypot(x, y)
    q_dom = float(w.sum()) * q_e                          # in-iris survivors
    q_bore = float(w[r <= BORE_R].sum()) * q_e            # of those, within the RF bore

    (x, y, z, ux, uy, uz), w = downsample(
        (x, y, z, ux, uy, uz), w, max_part, np.random.default_rng(rng_seed))
    z = z - z.min() + z_inject                            # bunch tail (smallest z) → z_inject

    v_beam, ke_mean = beam_kinematics(ux, uy, uz, w)
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
    if len(sys.argv) < 2 or sys.argv[1] not in ("1", "2", "3"):
        sys.exit("usage: python sim/linac1-3.py <N>   with N in {1, 2, 3}")
    N = int(sys.argv[1])

    prepare_env()
    build_linac_slac_fields()                            # idempotent; shared RF maps for all 3 sections
    from warpx import WarpX

    config = f"config/linac{N}.yaml"
    w = WarpX(input_file=config, path=f"logs/diags/linac1-3/sec{N}")
    NR, NZ = w.get("grid/number_of_cells")
    _, ZMAX = w.get("grid/upper_bound")
    outdir = w.get("diagnostics/0/write_dir")
    p = w.get("params")
    omega = 2.0 * np.pi * p["F_RF"]

    # The last applied field must load_E or picmi forces the global E_ext style to "none" (RF dark).
    # Both RF maps load_E here, so this is always satisfied — kept so a future reorder fails loudly.
    assert w.get("fields")[-1].get("load_E"), \
        f"config/linac{N}.yaml: last applied field must have load_E:true"

    if os.path.isdir(outdir):                            # fresh diags (WarpX appends per dump)
        shutil.rmtree(outdir)

    # ── Input beam + frozen RF setpoints (branch on section) ──────────────────────────────────
    if N == 1:
        bunch, v_beam, ke_mean, inj = load_injector_bunch(
            p["MAX_PART"], p["RNG_SEED"], p["Z_INJECT"], p["Z_HANDOFF"], p["COLLIM_Z"])
        summary = inj                                   # records z_handoff_m for downstream lab-z chaining
        # Capture beam (slipping ~150 keV): scale from input power, phase referenced to arrival.
        scale = float(np.sqrt(p["POWER_MW"] / p["RF_NORM_MW"]))
        base_deg = p["PHASE_DEG"]
        z_span = 0.0                                    # section 1 stops the centroid (old sec1 behaviour)
    else:
        bunch, v_beam, ke_mean, info = load_warpx_exit_bunch(
            PREV_PARTICLES[N], PREV_LABEL[N], p["MAX_PART"], p["RNG_SEED"], p["Z_INJECT"])
        # Lab-z of this section's injection: chain the upstream local→lab offset + its exit ⟨z⟩.
        z_inject_lab = upstream_exit_lab_z(PREV_SUMMARY[N], info["exit_zmean_local_m"])
        # FROZEN setpoints (no runtime crest-finding / ΔE-target loop): read from the yaml.
        scale = float(p["FIELD_SCALE"])
        base_deg = p["CREST_PHASE_DEG"] + p["PHASE_DEG"]   # PHASE_DEG = detune from the frozen crest (0)
        # Stop the bunch HEAD (not just the centroid) clear of the absorbing wall: the captured core
        # is ~1 RF wavelength long, so subtract its z-extent or the run drains the domain (MLMG abort).
        z_span = float(bunch["z"].max() - bunch["z"].min())
        summary = dict(n_injected=info["n_injected"], q_injected_C=info["q_injected_C"],
                       z_inject_lab_m=float(z_inject_lab), z_inject_mean_m=info["z_inject_mean_m"],
                       ke_in_mev=float(ke_mean / 1e3), de_target_mev=float(p["DE_TARGET_MEV"]),
                       power_equiv_mw=float(p["RF_NORM_MW"] * scale ** 2),
                       sigma_z_m=info["sigma_z_m"], rmax_m=info["rmax_m"])

    z_center = float(np.average(bunch["z"], weights=bunch["w"]))

    # ── RF amplitude + phase (uniform across sections; field2 is the 90° quadrature half) ──────
    # phi is referenced to the bunch arrival time at the structure entrance; base_deg is either the
    # capture phase offset (sec 1) or the frozen relativistic crest + detune (sec 2/3).
    t_in = (Z_STRUCT - z_center) / v_beam
    phi = -omega * t_in + np.deg2rad(base_deg)
    phi2 = phi + np.pi / 2.0
    e1, b1 = rf_time_functions(scale, omega, phi, amp_prec=8, phase_prec=8)
    e2, b2 = rf_time_functions(scale, omega, phi2, amp_prec=8, phase_prec=8)
    gain_keV = scale * V1KW_KEV                          # ≈ on-crest mean gain [keV] (1-kW V1kW from the maps)
    print(f"Section {N}: scale={scale:.1f}, base_phase={base_deg:g}°, f_RF={p['F_RF']/1e6:.0f} MHz, "
          f"⟨KE⟩_in={ke_mean/1e3:.2f} MeV → {outdir}/", flush=True)

    # ── Time step / duration: segmented transit to a stop plane short of ZMAX so the bunch
    # finishes in the field-free exit drift, NOT at the absorbing wall (empty domain aborts MLMG).
    dt = p["CFL"] * (ZMAX / NZ) / v_beam
    beta_in = v_beam / c
    beta_hi = float(np.sqrt(1.0 - 1.0 / (1.0 + (ke_mean + gain_keV) / MC2_KEV) ** 2))
    z_end = ZMAX - 0.20 - z_span                         # centroid stop ⇒ head ends short of the wall
    L_cap = 0.40                                         # length over which β: in→hi
    beta_cap = 0.5 * (beta_in + beta_hi)
    transit = ((Z_STRUCT - z_center) / v_beam
               + min(L_cap, max(0.0, z_end - Z_STRUCT)) / (beta_cap * c)
               + max(0.0, z_end - Z_STRUCT - L_cap) / (beta_hi * c))
    n_steps = p["MAX_STEPS"] if p.get("MAX_STEPS") else int(p["TRANSIT_MARGIN"] * transit / dt)
    period = max(1, n_steps // p["N_DIAGS"])
    print(f"dt = {dt:.3e} s, max_steps = {n_steps}, β_in={beta_in:.3f}, "
          f"est. on-crest gain {gain_keV/1e3:.1f} MeV → ⟨KE⟩_out≈{(ke_mean+gain_keV)/1e3:.1f} MeV",
          flush=True)

    # Persist the injection bookkeeping (charge in + local→lab z chain) for the chain plotter.
    # WarpX makes outdir at first diag; make it now for this sidecar.
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "injection_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    pg = make_particle_group(bunch["x"], bunch["y"], bunch["z"],
                             bunch["ux"], bunch["uy"], bunch["uz"], bunch["w"])
    w.initial_particles = pg

    w.update({
        "simulation/max_steps": n_steps,
        "simulation/time_step_size": dt,
        "diagnostics/0/period": period,
        "diagnostics/0/write_dir": outdir,
        "fields/0/warpx_E_time_function": e1, "fields/0/warpx_B_time_function": b1,
        "fields/1/warpx_E_time_function": e2, "fields/1/warpx_B_time_function": b2,
    })

    print(f"\nRunning {n_steps} steps (diag every {period}) → {outdir}/")
    w.run(progress=f"linac{N}")
    print("\nDone.")


if __name__ == "__main__":
    main()
