"""CESR injector in WarpX (RZ): the full injector subsection in one self-consistent
space-charge run -- two 214 MHz prebuncher cavities (Preb 2 reversed) and six solenoid
lenses (Lens 0A / Sol 0 / Lens 0E carry current at the default tune) -- reading the gun
exit beam and handing a focused, velocity-bunched beam to linac_sec1 at z ~= 2.03 m.

Drives lume-warpx from config/injector.yaml (which holds every constant); this module reads
those back, imports the gun handoff via WarpX(initial_particles=...), and overrides only
runtime-computed values (the per-field RF/solenoid time functions, step count, dt, diag period).
See docs/injector.md for physics, parameters, field maps, and gotchas.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# Must precede `import numpy`: OpenMP latches OMP_NUM_THREADS at load, so prepare_env()'s later
# set is ignored and the tiny grid oversubscribes.
os.environ.setdefault("OMP_NUM_THREADS", os.environ.get("OMP_THREADS", "1"))

import json
import shutil

import numpy as np

from sim.helpers.tools import C_LIGHT as c, E_CHARGE as q_e, MC2_KEV, prepare_env, rf_time_functions
from sim.helpers.loadparticles import (
    make_particle_group, downsample, beam_kinematics, open_particle_series,
    pipe_violator_ids, survivor_mask)
from sim.helpers.buildfields import (
    build_injector_fields, Z_GAP_CENTER_1, Z_GAP_CENTER_2, V1J_KEV, SOL_FILES,
    INJ_Z_HANDOFF as Z_HANDOFF, RMAX as IRIS_R)   # IRIS_R: the one iris radius (= the linac scrape RMAX)

CONFIG = "config/injector.yaml"
OUTDIR = "logs/diags/injector/main"
# Prefer the gun's reconstructed time-release exit beam when present, else the legacy snapshot.
GUN_DIAG = ("logs/diags/gun/handoff" if os.path.isdir("logs/diags/gun/handoff")
            else "logs/diags/gun/particles")


def load_gun_bunch(max_part, rng_seed, z_inject):
    """Import the gun's last snapshot (already RZ) and shift it to the entrance.

    Returns (dict [gamma*beta momenta], v_beam, mean KE [keV], z_centroid). The cavities are
    phased to put the CENTROID (not the tail) at the zero-crossing. See docs -> RF drive.
    """
    ts = open_particle_series(GUN_DIAG, "gun")
    it = ts.iterations[-1]
    x, y, z, ux, uy, uz, w = ts.get_particle(
        ["x", "y", "z", "ux", "uy", "uz", "w"], species="electrons", iteration=it)
    (x, y, z, ux, uy, uz), w = downsample(
        (x, y, z, ux, uy, uz), w, max_part, np.random.default_rng(rng_seed))
    z = z - z.min() + z_inject                         # bunch tail (smallest z) -> z_inject

    v_beam, ke_mean = beam_kinematics(ux, uy, uz, w)
    z_centroid = float(np.average(z, weights=w))
    print(f"Imported {z.size} macroparticles from gun (iter {it}); "
          f"z {z.min()*1e3:.1f}-{z.max()*1e3:.1f} mm, <z> {z_centroid*1e3:.1f} mm, "
          f"<KE> {ke_mean:.1f} keV, v_beam {v_beam:.3e} m/s, q {w.sum()*q_e*1e9:.3f} nC", flush=True)
    return dict(x=x, y=y, z=z, ux=ux, uy=uy, uz=uz, w=w), v_beam, ke_mean, z_centroid


def cavity_drive(power, q_l, f_rf, z_gap, v_at_gap, phi_off_deg, phase, omega,
                 t_offset=0.0, rev_phase=0.0, z_ref=0.0):
    """Time-function strings (warpx_E/B_time_function) for one prebuncher cavity.

    Drives the 1-J map as a standing-wave TM mode: E ~ scale*cos(wt+phi), B ~ scale*sin(wt+phi).
    The zc base lands the bunch centroid on the RF zero-crossing (net mean kick 0). Keep .10e
    precision -- w*t truncation accumulates over the ~5 ns transit. See docs -> RF drive.
    Returns (e_time, b_time, scale, phi, t_gap).
    """
    scale = float(np.sqrt(1e3 * q_l * power / (2.0 * np.pi * f_rf)))
    t_gap = t_offset + (z_gap - z_ref) / v_at_gap
    base = np.pi / 2.0 if phase == "zc" else np.pi
    phi = -omega * t_gap + base + np.radians(phi_off_deg) + rev_phase
    e_time, b_time = rf_time_functions(scale, omega, phi)   # amp/phase precision .10e
    return e_time, b_time, scale, phi, t_gap


def _report_collimated_handoff(outdir, collim_r, collim_z):
    """Report and return the multi-plane-collimated handoff charge at the ~Z_HANDOFF plane
    (sanity log; the physical cut is the linac reader's at injection). Returns a dict of
    metrics (or {} if no populated snapshot near the plane). See docs -> The 9.547 mm collimator.
    """
    from openpmd_viewer import OpenPMDTimeSeries
    ts = OpenPMDTimeSeries(os.path.join(outdir, "particles"))
    recs = []
    for it in ts.iterations:
        z, w = ts.get_particle(["z", "w"], species="electrons", iteration=it)
        if len(z) < 50:
            continue
        recs.append((it, float(np.average(z, weights=w))))
    if not recs:
        print("  collimated handoff: no populated snapshot near the plane", flush=True)
        return {}
    it_h, zm_h = min(recs, key=lambda t: abs(t[1] - Z_HANDOFF))
    idh, x, y, z, w = ts.get_particle(["id", "x", "y", "z", "w"],
                                      species="electrons", iteration=it_h)
    q_dom = float(w.sum()) * q_e
    scan_iters = [it for it, zm in recs if (collim_z - 0.05) <= zm <= (Z_HANDOFF + 0.03)]
    violators = pipe_violator_ids(ts, scan_iters, collim_r, collim_z)
    q_coll = float(w[survivor_mask(idh, violators)].sum()) * q_e
    print(f"  COLLIMATED handoff (<z>={zm_h*1e3:.1f} mm, iris {collim_r*1e3:.3f} mm, "
          f"multi-plane {len(scan_iters)} planes): {q_coll*1e9:.3f} nC survives the pipe / "
          f"{q_dom*1e9:.3f} nC in-domain = {100*q_coll/q_dom:.0f}% through the aperture", flush=True)
    return dict(it_handoff=int(it_h), z_handoff_mean_m=zm_h, collim_r_m=collim_r, collim_z_m=collim_z,
                n_scan_planes=len(scan_iters), q_in_domain_C=q_dom, q_collimated_C=q_coll,
                transmission=(q_coll / q_dom if q_dom else 0.0))


def main():
    prepare_env()
    from warpx import WarpX

    build_injector_fields()                            # idempotent; (re)build the openPMD maps

    w = WarpX(input_file=CONFIG, path="logs/diags/injector")
    NR, NZ = w.get("grid/number_of_cells")
    RMAX, ZMAX = w.get("grid/upper_bound")
    outdir = w.get("diagnostics/0/write_dir") or OUTDIR
    p = w.get("params")
    F_RF = p["F_RF"]
    omega = 2.0 * np.pi * F_RF
    PHASE = p["PHASE"]
    base = np.pi / 2.0 if PHASE == "zc" else np.pi

    # The last applied field MUST load_E or picmi forces the global E_ext style to "none" and the
    # RF cavities go dark -- solenoids (load_E:false) must precede them in the YAML fields list.
    assert w.get("fields")[-1].get("load_E"), \
        "config/injector.yaml: last applied field must have load_E:true (RF cavity), solenoids first"

    if os.path.isdir(outdir):                           # fresh diags (WarpX appends per dump)
        shutil.rmtree(outdir)

    bunch, v_beam, ke_mean, z_centroid = load_gun_bunch(p["MAX_PART"], p["RNG_SEED"], p["Z_INJECT"])
    pg = make_particle_group(bunch["x"], bunch["y"], bunch["z"],
                             bunch["ux"], bunch["uy"], bunch["uz"], bunch["w"])
    w.initial_particles = pg                           # imported beam for FromInitialParticles

    for nm, cur in (("0A", p["I_LENS0A"]), ("0B", p["I_LENS0B"]), ("0C", p["I_LENS0C"]),
                    ("0D", p["I_LENS0D"]), ("Sol0", p["I_SOL0"]), ("0E", p["I_LENS0E"])):
        if cur != 0.0:
            print(f"Solenoid {nm}: I={cur:g} A", flush=True)

    # `V_gap~` (= scale*V1J_KEV) is a transit-time-free upper bound for diagnostics/step sizing only,
    # not a physics input -- WarpX integrates the real time-varying field independently.
    # -- Prebuncher 1 (forward map): centroid arrival uses v_beam over z_centroid->gap --
    e1, b1, scale1, phi1, t_gap1 = cavity_drive(
        p["PREB1_KW"], p["Q_L_1"], F_RF, Z_GAP_CENTER_1, v_beam, p["PREB1_PHI_OFF"],
        PHASE, omega, z_ref=z_centroid)
    print(f"Preb 1: P={p['PREB1_KW']:g} kW, scale={scale1:.3f}, V_gap~={scale1*V1J_KEV:.1f} kV, "
          f"phi={phi1:.3f} rad, t_gap={t_gap1*1e9:.3f} ns", flush=True)

    # -- Prebuncher 2 (reversed install): two-segment arrival accounts for Preb-1's kick --
    # Baked BEFORE WarpX integrates Preb 1, so estimate post-Preb-1 speed analytically (mean kick).
    kick1 = -np.cos(base + np.radians(p["PREB1_PHI_OFF"])) * scale1 * V1J_KEV
    ke_after1 = max(ke_mean + (kick1 if p["PREB1_KW"] > 0 else 0.0), 1.0)
    v_after_preb1 = c * np.sqrt(1.0 - 1.0 / (1.0 + ke_after1 / MC2_KEV) ** 2)
    rev_phase = p["PREB2_REV_PHASE"] if p["PREB2_REVERSED"] else 0.0
    e2, b2, scale2, phi2, t_gap2 = cavity_drive(
        p["PREB2_KW"], p["Q_L_2"], F_RF, Z_GAP_CENTER_2, v_after_preb1, p["PREB2_PHI_OFF"],
        PHASE, omega, t_offset=t_gap1, z_ref=Z_GAP_CENTER_1, rev_phase=rev_phase)
    if p["PREB2_KW"] > 0:
        print(f"Preb 2 (reversed): P={p['PREB2_KW']:g} kW, scale={scale2:.3f}, "
              f"V_gap~={scale2*V1J_KEV:.1f} kV, phi={phi2:.3f} rad, t_gap={t_gap2*1e9:.3f} ns "
              f"(v_after_preb1={v_after_preb1:.3e} m/s from +{kick1:.1f} keV Preb-1 kick)", flush=True)

    # -- Time step / duration: 3-leg transit estimate with the real per-leg speed --
    dt = p["CFL"] * (ZMAX / NZ) / v_beam
    kick_frac1 = -np.cos(base + np.radians(p["PREB1_PHI_OFF"]))
    ke_after1 = max(ke_mean + (kick_frac1 * scale1 * V1J_KEV if p["PREB1_KW"] > 0 else 0.0), 1.0)
    v_after1 = c * np.sqrt(1.0 - 1.0 / (1.0 + ke_after1 / MC2_KEV) ** 2)
    if p["PREB2_KW"] > 0:
        kick_frac2 = -np.cos(base + np.radians(p["PREB2_PHI_OFF"]) + rev_phase)
        ke_after2 = max(ke_after1 + kick_frac2 * scale2 * V1J_KEV, 1.0)
        v_after2 = c * np.sqrt(1.0 - 1.0 / (1.0 + ke_after2 / MC2_KEV) ** 2)
    else:
        ke_after2, v_after2 = ke_after1, v_after1
    transit = ((Z_GAP_CENTER_1 - z_centroid) / v_beam
               + (Z_GAP_CENTER_2 - Z_GAP_CENTER_1) / v_after1
               + (ZMAX - Z_GAP_CENTER_2) / v_after2)
    n_steps = p["MAX_STEPS"] if p.get("MAX_STEPS") else int(p["TRANSIT_MARGIN"] * transit / dt)
    # Size `period` so dump spacing near the handoff is <= HANDOFF_DZ (post-Preb-2 speed v_after2).
    period_handoff = max(1, int(p["HANDOFF_DZ"] / (v_after2 * dt)))
    period = min(max(1, n_steps // p["N_DIAGS"]), period_handoff)
    print(f"  <KE> after Preb-1 ~= {ke_after1:.1f} keV, after Preb-2 ~= {ke_after2:.1f} keV", flush=True)
    print(f"dt = {dt:.3e} s, max_steps = {n_steps}, diag period {period} steps "
          f"(~{period*v_after2*dt*1e3:.1f} mm near handoff)", flush=True)

    w.update({
        "simulation/max_steps": n_steps,
        "simulation/time_step_size": dt,
        "diagnostics/0/period": period,
        "diagnostics/0/write_dir": outdir,
        "fields/0/warpx_B_time_function": f"{p['I_LENS0A']:.8e}",
        "fields/1/warpx_B_time_function": f"{p['I_LENS0B']:.8e}",
        "fields/2/warpx_B_time_function": f"{p['I_LENS0C']:.8e}",
        "fields/3/warpx_B_time_function": f"{p['I_LENS0D']:.8e}",
        "fields/4/warpx_B_time_function": f"{p['I_SOL0']:.8e}",
        "fields/5/warpx_B_time_function": f"{p['I_LENS0E']:.8e}",
        "fields/6/warpx_E_time_function": e1 if p["PREB1_KW"] > 0 else "0.0",
        "fields/6/warpx_B_time_function": b1 if p["PREB1_KW"] > 0 else "0.0",
        "fields/7/warpx_E_time_function": e2 if p["PREB2_KW"] > 0 else "0.0",
        "fields/7/warpx_B_time_function": b2 if p["PREB2_KW"] > 0 else "0.0",
    })

    print(f"\nRunning {n_steps} steps (diag every {period}) -> {outdir}/")
    w.run(progress="injector")
    print("\nDone.")

    summary = dict(v_beam_mps=v_beam, ke_mean_keV=ke_mean, z_centroid_m=z_centroid,
                   n_steps=int(n_steps), dt_s=float(dt), period=int(period))
    if p["COLLIMATE"]:
        try:
            summary.update(_report_collimated_handoff(outdir, IRIS_R, p["COLLIM_Z"]))
        except Exception as e:
            print(f"  (collimated-handoff report failed: {e}; writing base summary)", flush=True)
    with open(os.path.join(outdir, "injection_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    main()
