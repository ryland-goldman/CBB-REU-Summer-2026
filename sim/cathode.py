"""
Finite thermionic cathode emitting through a pulsed grid — WarpX 2D.

Drives lume-warpx from config/cathode.yaml (which holds every constant); this module reads
those back and overrides the runtime-computed values (flux, thermal velocity, dt, the diagnostic
periods, and the V(t) grid-pulse boundary). The grid bias is pulsed (LinacSim CESR operating
point): the anode rides V(t) = V_OFF + V_PULSE·tent(t), and over-injecting at 2× the peak-voltage
J_CL lets WarpX's self-consistent fields self-limit the transmitted current to J_CL(V(t)) through
the pulse — so the EMITTED CHARGE is measured (∫J_z over the disc × grid transmission), not imposed.
See docs/cathode.md for physics, parameters, and gotchas.

main() runs the simulation and writes logs/diags/cathode/injection_summary.json (the measured
bunch charge the gun renormalizes to); sim/plot/cathode.py produces the figures.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# OpenMP latches OMP_NUM_THREADS when its runtime loads (at the numpy/h5py import below);
# prepare_env()'s later set is ignored, so a standalone run would oversubscribe this tiny grid
# (slower). Pin it here, first. OMP_THREADS overrides; main.py sets it in the child env.
os.environ.setdefault("OMP_NUM_THREADS", os.environ.get("OMP_THREADS", "1"))

import glob
import json
import shutil

import h5py
import numpy as np

from sim.helpers.tools import (
    E_CHARGE as q_e,
    M_E as m_e,
    child_langmuir_current_density,
    thermal_velocity_sigma,
    prepare_env,
)

CONFIG = "config/cathode.yaml"
DIAG_DIR = "logs/diags/cathode"


def _pulse_string(v_off, v_pulse, v_slope):
    """WarpX parser expression for V(t): a rounded triangle rising/falling at v_slope, peak swing
    v_pulse, peaking at t = t_rise. tent(t) ∈ [0,1] is 0 outside [0, 2·t_rise]."""
    t_rise = v_pulse / v_slope
    return (f"{v_off!r} + {v_pulse!r}*max(0.0, "
            f"min(t/{t_rise!r}, ({2.0 * t_rise!r}-t)/{t_rise!r}))"), t_rise


def _measure_emitted_charge(emit_r, grid_trans):
    """Integrate the transmitted current density over the pulse to get the physical emitted charge.

    Q = π·emit_r²·∫⟨J_z⟩_midgap(t) dt · grid_trans. J_z [A/m²] is a real local current density even
    in 2D (the planar diode is locally 1D), so this is the physical disc charge — the naive Σ(weight)
    is not (2D weights are per-unit-out-of-plane-length). Mid-gap row, not the anode, avoids the
    ~14% collection-edge inflation seen at the absorbing boundary."""
    files = sorted(glob.glob(os.path.join(DIAG_DIR, "fields", "openpmd_*.h5")),
                   key=lambda f: int(f.split("_")[-1].split(".")[0]))
    times, jbar = [], []
    for fn in files:
        with h5py.File(fn, "r") as f:
            it = f["data"][list(f["data"].keys())[0]]
            jg = it["fields"]["j"]
            jz = np.array(jg["z"]) * jg["z"].attrs["unitSI"]      # A/m², (nz, nx)
            gs, off = jg.attrs["gridSpacing"], jg.attrs["gridGlobalOffset"]
            nx = jz.shape[1]
            x = off[1] + (np.arange(nx) + 0.5) * gs[1]
            xm = np.abs(x) <= emit_r
            times.append(float(it.attrs["time"]))
            jbar.append(float(np.abs(jz[jz.shape[0] // 2, xm]).mean()))
    times, jbar = np.asarray(times), np.asarray(jbar)
    q_area = float(np.trapezoid(jbar, times))                    # C/m² transmitted over the pulse
    area = np.pi * emit_r**2
    q_pre_grid = q_area * area
    return dict(q_emit_C=q_pre_grid * grid_trans, q_pre_grid_C=q_pre_grid,
                peak_current_A=float(jbar.max() * area), grid_trans=grid_trans,
                emit_radius_m=emit_r, pulse_charge_density_C_per_m2=q_area)


def main():
    prepare_env()
    from warpx import WarpX

    w = WarpX(input_file=CONFIG, path=DIAG_DIR)
    gap_d = w.get("grid/upper_bound")[1]
    nz = w.get("grid/number_of_cells")[1]
    p = w.get("params")

    v_peak = p["V_OFF"] + p["V_PULSE"]                  # peak grid bias [V] (the SCL operating point)
    pulse_str, t_rise = _pulse_string(p["V_OFF"], p["V_PULSE"], p["V_SLOPE"])

    # Diode sized at the PEAK voltage: 2×J_CL(v_peak) constant flux stays over-injected (SCL) at every
    # instant of the pulse since J_CL(V(t)) ≤ J_CL(v_peak), so the transmitted current self-limits.
    J_CL = float(child_langmuir_current_density(v_peak, gap_d))
    flux = p["over_inject"] * J_CL / q_e
    v_th = thermal_velocity_sigma(p["T_cathode"])
    v_final = np.sqrt(2.0 * q_e * v_peak / m_e)         # cold final velocity through the peak bias
    dt = p["CFL"] * (gap_d / nz) / v_final

    # Run spans the full pulse base (2·t_rise) plus drift for the last-emitted slug to clear the gap.
    drift = 4.0 * gap_d / v_final
    max_steps = int((2.0 * t_rise + drift) / dt) + 1
    peak_step = int(t_rise / dt)                        # tent peaks at t_rise → representative template

    print(f"Diode : pulsed grid, peak V = {v_peak:.0f} V (off {p['V_OFF']:.0f} V), gap d = {gap_d*1e3:.1f} mm")
    print(f"Pulse : V_PULSE {p['V_PULSE']:.0f} V, slope {p['V_SLOPE']/1e9:.0f} V/ns, "
          f"t_rise {t_rise*1e9:.2f} ns, base {2*t_rise*1e9:.1f} ns, grid trans {p['GRID_TRANS']*100:.0f}%")
    print(f"Child-Langmuir J_CL(peak) = {J_CL:.1f} A/m^2  "
          f"(injecting {p['over_inject']:.0f}x = {p['over_inject']*J_CL:.1f} A/m^2)")
    print(f"v_th = {v_th:.2e} m/s, v_final = {v_final:.2e} m/s")

    if w.get("species/0/warpx_do_not_deposit"):
        print("WARNING: cathode warpx_do_not_deposit=true — beam self-field is OFF, so the "
              "space-charge-limited (Child-Langmuir) mechanism is disabled and the validation "
              "figures are NOT valid. Forces-off diagnostic only.", flush=True)

    # Fresh diags: the h5 backend appends one file per dump, so stale files would corrupt plots.
    if os.path.isdir(DIAG_DIR):
        shutil.rmtree(DIAG_DIR)

    # Field period: uniform sampling resolving the pulse (~200 dumps) for the ∫J_z charge integral.
    # Particle period: a short series whose LAST dump lands on peak_step, so the gun's iterations[-1]
    # grabs the peak-emission template (the gap full of transiting electrons), not the drained tail.
    field_period = str(p["DIAG_PERIOD"]) if p["DIAG_PERIOD"] else str(max(1, max_steps // 200))
    part_stride = max(1, peak_step // 4)
    w.update({
        "grid/warpx_potential_hi_z": pulse_str,
        "simulation/max_steps": max_steps,
        "simulation/time_step_size": dt,
        "species/0/flux": flux,
        "species/0/rms_velocity": [v_th, v_th, v_th],
        "diagnostics/0/period": field_period,
        "diagnostics/1/period": f"0:{peak_step}:{part_stride}, {peak_step}:{peak_step+1}:1",
    })

    print(f"\nRunning {max_steps} steps  dt = {dt:.3e} s  "
          f"(pulse base ≈ {int(2*t_rise/dt)} steps, peak template at step {peak_step})")
    w.run(progress="cathode")

    info = _measure_emitted_charge(p["R_CATHODE"], p["GRID_TRANS"])
    # The crest time tags which particle dump is the representative beam: WarpX force-writes a
    # diagnostic at the final step (gap drained, grid off), so iterations[-1] is NOT the beam.
    # gap_d_m + anode_frac define the anode handoff slab (the delivered flux the gun seeds from).
    info["crest_time_s"] = t_rise
    info["gap_d_m"] = gap_d
    info["anode_frac"] = p["ANODE_FRAC"]
    with open(os.path.join(DIAG_DIR, "injection_summary.json"), "w") as f:
        json.dump(info, f, indent=2)
    print(f"\nEmitted charge (measured) = {info['q_emit_C']*1e9:.3f} nC "
          f"(pre-grid {info['q_pre_grid_C']*1e9:.3f} nC, peak I = {info['peak_current_A']:.3f} A)")
    print(f"Done. openPMD output → {DIAG_DIR}/{{fields,particles}}/ + injection_summary.json")


if __name__ == "__main__":
    main()
