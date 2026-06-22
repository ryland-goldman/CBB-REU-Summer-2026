"""
Finite thermionic cathode at the space-charge (Child–Langmuir) limit — WarpX 2D.

Drives lume-warpx from cathode/cathode.yaml (which holds every constant); this module
reads those back and overrides only the runtime-computed values (flux, thermal velocity,
dt, diagnostic periods). Over-injects at 2×J_CL and lets WarpX's self-consistent fields
build a virtual cathode that self-limits the transmitted current to J_CL. See
cathode/README.md for physics, parameters, and gotchas.

Run from the repo root: python -c "import cathode; cathode.run()".
"""

import os
import shutil

import numpy as np

from pipeline.constants import E_CHARGE as q_e, M_E as m_e
from pipeline.emission import child_langmuir_current_density, thermal_velocity_sigma

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG = os.path.join(HERE, "cathode.yaml")
DIAG_DIR = "cathode/diags"


def main():
    from warpx import WarpX

    w = WarpX(input_file=CONFIG, path="cathode")
    V_anode = w.get("grid/warpx_potential_hi_z")
    gap_d = w.get("grid/upper_bound")[1]
    nz = w.get("grid/number_of_cells")[1]
    max_steps = w.get("simulation/max_steps")
    p = w.get("params")

    J_CL = float(child_langmuir_current_density(V_anode, gap_d))
    flux = p["over_inject"] * J_CL / q_e
    v_th = thermal_velocity_sigma(p["T_cathode"])
    v_final = np.sqrt(2.0 * q_e * V_anode / m_e)        # cold final velocity through full bias
    dt = p["CFL"] * (gap_d / nz) / v_final

    print(f"Diode : V = {V_anode:.0f} V, gap d = {gap_d*1e3:.1f} mm")
    print(f"Child–Langmuir J_CL = {J_CL:.1f} A/m^2  "
          f"(injecting {p['over_inject']:.0f}× = {p['over_inject']*J_CL:.1f} A/m^2)")
    print(f"v_th = {v_th:.2e} m/s, v_final = {v_final:.2e} m/s")

    if w.get("species/0/warpx_do_not_deposit"):
        print("WARNING: cathode warpx_do_not_deposit=true — beam self-field is OFF, so the "
              "space-charge-limited (Child–Langmuir) mechanism is disabled and the validation "
              "figures are NOT valid. Forces-off diagnostic only.", flush=True)

    # Fresh diags: the h5 backend appends one file per dump, so stale files would corrupt plots.
    if os.path.isdir(DIAG_DIR):
        shutil.rmtree(DIAG_DIR)

    # Union slice: dense through the gap-fill transient (every 5 steps to 470), sparse after.
    # max(max_steps, 471) guards the second slice from inverting if max_steps ≤ 470.
    field_period = (str(p["DIAG_PERIOD"]) if p["DIAG_PERIOD"]
                    else f"0:470:5, 470:{max(max_steps, 471)}:80")
    w.update({
        "simulation/time_step_size": dt,
        "species/0/flux": flux,
        "species/0/rms_velocity": [v_th, v_th, v_th],
        "diagnostics/0/period": field_period,
        "diagnostics/1/period": p["DIAG_PERIOD"] or 200,
    })

    print(f"\nRunning {max_steps} steps  dt = {dt:.3e} s  "
          f"(gap-fill ≈ {int(3*gap_d/v_final/dt)} steps)")
    w.run(progress="cathode")
    print(f"\nDone. openPMD output → {DIAG_DIR}/{{fields,particles}}/")


if __name__ == "__main__":
    main()
