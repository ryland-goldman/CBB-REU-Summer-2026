"""
CESR injector in WarpX (RZ): the full LinacSim injector subsection in one
self-consistent space-charge run — two 214 MHz prebuncher cavities (Preb 2 reversed)
and three solenoid lenses (Lens 0A / Sol 0 / Lens 0E) — reading the gun exit beam and
handing a focused, velocity-bunched beam to linac_sec1 at z ≈ 2.03 m.

Drives lume-warpx from injector/injector.yaml (which holds every constant); this module
reads those back, imports the gun handoff via WarpX(initial_particles=...), and overrides only
runtime-computed values (the per-field RF/solenoid time functions, step count, dt, diag period).
See injector/README.md for physics, parameters, field maps, and gotchas.
"""

import os
import shutil
import time

import numpy as np
import openpmd_api as io
from openpmd_viewer import OpenPMDTimeSeries

from pipeline.constants import C_LIGHT as c, M_E as m_e, E_CHARGE as q_e, MC2_EV
from .build_injector_field import (
    Z_GAP_CENTER_1, Z_GAP_CENTER_2, V1J_KEV, F_RF, Q_L_1, Q_L_2, SOL_FILES, Z_HANDOFF)
from . import DEFAULT_OUTDIR

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG = os.path.join(HERE, "injector.yaml")
# Prefer the gun's reconstructed time-release exit beam when present, else the legacy snapshot.
GUN_DIAG = ("gun/diags/handoff" if os.path.isdir("gun/diags/handoff")
            else "gun/diags/particles")


def _retry_io(fn, *args, tries=6, base=0.25, **kwargs):
    """Call an openPMD read, retrying a transient HDF5 "Inaccessible" open error (backstop only;
    the production failure is fd exhaustion, fixed by raising RLIMIT_NOFILE)."""
    for i in range(tries):
        try:
            return fn(*args, **kwargs)
        except io.Error:
            if i == tries - 1:
                raise
            time.sleep(base * 2 ** i)


def load_gun_bunch(max_part, rng_seed, z_inject):
    """Import the gun's last snapshot (already RZ) and shift it to the entrance.

    Returns (dict [γβ momenta], v_beam, mean KE [keV], z_centroid). The cavities are phased to
    put the CENTROID (not the tail) at the zero-crossing. See README -> RF drive.
    """
    ts = OpenPMDTimeSeries(GUN_DIAG)
    if len(ts.iterations) == 0:
        raise RuntimeError(f"{GUN_DIAG} has no iterations — did the gun stage run?")
    it = ts.iterations[-1]
    x, y, z, ux, uy, uz, w = _retry_io(
        ts.get_particle, ["x", "y", "z", "ux", "uy", "uz", "w"], species="electrons", iteration=it)
    if z.size > max_part:
        rng = np.random.default_rng(rng_seed)
        sel = rng.choice(z.size, max_part, replace=False)
        scale_w = z.size / max_part
        x, y, z, ux, uy, uz, w = (a[sel] for a in (x, y, z, ux, uy, uz, w))
        w = w * scale_w
    z = z - z.min() + z_inject                         # bunch tail (smallest z) → z_inject

    gb = np.sqrt(1.0 + ux**2 + uy**2 + uz**2)          # γ (ux/uy/uz are γβ)
    v_beam = float(np.average(uz / gb, weights=w) * c)
    ke_mean = float(np.average(gb - 1.0, weights=w) * m_e * c**2 / q_e / 1e3)
    z_centroid = float(np.average(z, weights=w))
    print(f"Imported {z.size} macroparticles from gun (iter {it}); "
          f"z {z.min()*1e3:.1f}–{z.max()*1e3:.1f} mm, ⟨z⟩ {z_centroid*1e3:.1f} mm, "
          f"⟨KE⟩ {ke_mean:.1f} keV, v_beam {v_beam:.3e} m/s, q {w.sum()*q_e*1e9:.3f} nC", flush=True)
    return dict(x=x, y=y, z=z, ux=ux, uy=uy, uz=uz, w=w), v_beam, ke_mean, z_centroid


def cavity_drive(power, q_l, z_gap, v_at_gap, phi_off_deg, phase, omega,
                 t_offset=0.0, rev_phase=0.0, z_ref=0.0):
    """Time-function strings (warpx_E/B_time_function) for one prebuncher cavity.

    Drives the 1-J map as a standing-wave TM mode: E ∝ scale·cos(ωt+φ), B ∝ scale·sin(ωt+φ).
    The zc base lands the bunch centroid on the RF zero-crossing (net mean kick 0). Keep .10e
    precision — ω·t truncation accumulates over the ~5 ns transit. See README -> RF drive.
    Returns (e_time, b_time, scale, phi, t_gap).
    """
    scale = float(np.sqrt(1e3 * q_l * power / (2.0 * np.pi * F_RF)))
    t_gap = t_offset + (z_gap - z_ref) / v_at_gap
    base = np.pi / 2.0 if phase == "zc" else np.pi
    phi = -omega * t_gap + base + np.radians(phi_off_deg) + rev_phase
    e_time = f"{scale:.10e}*cos({omega:.10e}*t + ({phi:.10e}))"
    b_time = f"{scale:.10e}*sin({omega:.10e}*t + ({phi:.10e}))"
    return e_time, b_time, scale, phi, t_gap


def _report_collimated_handoff(outdir, collim_r, collim_z):
    """Report the multi-plane-collimated handoff charge at the ~Z_HANDOFF plane (sanity log only;
    the physical cut is the linac reader's at injection). See README -> The 9.547 mm collimator."""
    try:
        from pipeline.collimator import pipe_violator_ids, survivor_mask
        ts = OpenPMDTimeSeries(os.path.join(outdir, "particles"))
        recs = []
        for it in ts.iterations:
            z, w = _retry_io(ts.get_particle, ["z", "w"], species="electrons", iteration=it)
            if len(z) < 50:
                continue
            recs.append((it, float(np.average(z, weights=w))))
        if not recs:
            print("  collimated handoff: no populated snapshot near the plane", flush=True)
            return
        it_h, zm_h = min(recs, key=lambda t: abs(t[1] - Z_HANDOFF))
        idh, x, y, z, w = _retry_io(ts.get_particle, ["id", "x", "y", "z", "w"],
                                    species="electrons", iteration=it_h)
        q_dom = float(w.sum()) * q_e
        scan_iters = [it for it, zm in recs if (collim_z - 0.05) <= zm <= (Z_HANDOFF + 0.03)]
        violators = pipe_violator_ids(ts, scan_iters, collim_r, collim_z)
        q_coll = float(w[survivor_mask(idh, violators)].sum()) * q_e
        print(f"  COLLIMATED handoff (⟨z⟩={zm_h*1e3:.1f} mm, iris {collim_r*1e3:.3f} mm, "
              f"multi-plane {len(scan_iters)} planes): {q_coll*1e9:.3f} nC survives the pipe / "
              f"{q_dom*1e9:.3f} nC in-domain = {100*q_coll/q_dom:.0f}% through the aperture", flush=True)
    except Exception as e:
        print(f"  collimated-handoff report unavailable: {e}", flush=True)


def main():
    from warpx import WarpX
    from pmd_beamphysics import ParticleGroup

    w = WarpX(input_file=CONFIG, path="injector")
    NR, NZ = w.get("grid/number_of_cells")
    RMAX, ZMAX = w.get("grid/upper_bound")
    outdir = w.get("diagnostics/0/write_dir") or DEFAULT_OUTDIR
    p = w.get("params")
    omega = 2.0 * np.pi * F_RF                          # in main() so a config(F_RF=...) is honoured
    PHASE = p["PHASE"]
    base = np.pi / 2.0 if PHASE == "zc" else np.pi
    MC2_KEV = MC2_EV / 1e3

    # The last applied field MUST load_E or picmi forces the global E_ext style to "none" and the
    # RF cavities go dark — solenoids (load_E:false) must precede them in the YAML fields list.
    assert w.get("fields")[-1].get("load_E"), \
        "injector.yaml: last applied field must have load_E:true (RF cavity), solenoids first"

    if os.path.isdir(outdir):                           # fresh diags (WarpX appends per dump)
        shutil.rmtree(outdir)

    bunch, v_beam, ke_mean, z_centroid = load_gun_bunch(p["MAX_PART"], p["RNG_SEED"], p["Z_INJECT"])
    pg = ParticleGroup(data=dict(
        x=bunch["x"], y=bunch["y"], z=bunch["z"],
        px=bunch["ux"] * MC2_EV, py=bunch["uy"] * MC2_EV, pz=bunch["uz"] * MC2_EV,
        t=np.zeros(bunch["x"].size), weight=bunch["w"] * q_e,
        status=np.ones(bunch["x"].size, dtype=np.int64), species="electron"))
    w.initial_particles = pg                           # imported beam for FromInitialParticles

    for nm, cur in (("0A", p["I_LENS0A"]), ("0B", p["I_LENS0B"]), ("0C", p["I_LENS0C"]),
                    ("0D", p["I_LENS0D"]), ("Sol0", p["I_SOL0"]), ("0E", p["I_LENS0E"])):
        if cur != 0.0:
            print(f"Solenoid {nm}: I={cur:g} A", flush=True)

    # ── Prebuncher 1 (forward map): centroid arrival uses v_beam over z_centroid→gap ──
    e1, b1, scale1, phi1, t_gap1 = cavity_drive(
        p["PREB1_KW"], Q_L_1, Z_GAP_CENTER_1, v_beam, p["PREB1_PHI_OFF"], PHASE, omega, z_ref=z_centroid)
    print(f"Preb 1: P={p['PREB1_KW']:g} kW, scale={scale1:.3f}, V_gap≈{scale1*V1J_KEV:.1f} kV, "
          f"φ={phi1:.3f} rad, t_gap={t_gap1*1e9:.3f} ns", flush=True)

    # ── Prebuncher 2 (reversed install): two-segment arrival accounts for Preb-1's kick ──
    # Baked BEFORE WarpX integrates Preb 1, so estimate post-Preb-1 speed analytically (mean kick).
    kick1 = -np.cos(base + np.radians(p["PREB1_PHI_OFF"])) * scale1 * V1J_KEV
    ke_after1 = max(ke_mean + (kick1 if p["PREB1_KW"] > 0 else 0.0), 1.0)
    v_after_preb1 = c * np.sqrt(1.0 - 1.0 / (1.0 + ke_after1 / MC2_KEV) ** 2)
    rev_phase = p["PREB2_REV_PHASE"] if p["PREB2_REVERSED"] else 0.0
    e2, b2, scale2, phi2, t_gap2 = cavity_drive(
        p["PREB2_KW"], Q_L_2, Z_GAP_CENTER_2, v_after_preb1, p["PREB2_PHI_OFF"], PHASE, omega,
        t_offset=t_gap1, z_ref=Z_GAP_CENTER_1, rev_phase=rev_phase)
    if p["PREB2_KW"] > 0:
        print(f"Preb 2 (reversed): P={p['PREB2_KW']:g} kW, scale={scale2:.3f}, "
              f"V_gap≈{scale2*V1J_KEV:.1f} kV, φ={phi2:.3f} rad, t_gap={t_gap2*1e9:.3f} ns "
              f"(v_after_preb1={v_after_preb1:.3e} m/s from +{kick1:.1f} keV Preb-1 kick)", flush=True)

    # ── Time step / duration: 3-leg transit estimate with the real per-leg speed ──
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
    # Size `period` so dump spacing near the handoff is ≤ HANDOFF_DZ (post-Preb-2 speed v_after2).
    period_handoff = max(1, int(p["HANDOFF_DZ"] / (v_after2 * dt)))
    period = min(max(1, n_steps // p["N_DIAGS"]), period_handoff)
    print(f"  ⟨KE⟩ after Preb-1 ≈ {ke_after1:.1f} keV, after Preb-2 ≈ {ke_after2:.1f} keV", flush=True)
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

    if p["COLLIMATE"]:
        _report_collimated_handoff(outdir, p["COLLIM_R"], p["COLLIM_Z"])


if __name__ == "__main__":
    main()
