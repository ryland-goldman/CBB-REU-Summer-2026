"""
Figures for the WarpX RZ CESR injector (injector_sim.py), generated entirely with
lume-warpx's plotting helpers (WarpX.plot2D / plot1D) over injector/diags/main/. Writes
PNGs to injector/results/. (No field diagnostic is dumped, so no plot_fields here.)

See injector/README.md for the physics each figure shows (velocity bunching to the
σ_z waist at the 2.03 m handoff).
"""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from openpmd_viewer import OpenPMDTimeSeries

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG = os.path.join(HERE, "injector.yaml")
RESULTS = "injector/results"
DIAG = "injector/diags/main/particles"


def _last_populated(diag, species="electrons"):
    ts = OpenPMDTimeSeries(diag)
    for it in reversed([int(i) for i in ts.iterations]):
        try:
            x, = ts.get_particle(["x"], species=species, iteration=it)
        except Exception:
            continue
        if len(x):
            return it
    return int(ts.iterations[-1])


def main():
    from warpx import WarpX
    os.makedirs(RESULTS, exist_ok=True)

    w = WarpX(input_file=CONFIG, path="injector")
    w.load_output(diag_dir=DIAG)
    it = _last_populated(DIAG)

    figs = [
        ("phase_space_z_KE", w.plot2D("z", "kinetic_energy", iteration=it)),   # energy-flat ~150 keV
        ("transverse_x_px",  w.plot2D("x", "px", iteration=it)),
        ("centroid_vs_t",    w.plot1D("t", "mean_z")),
        ("bunch_length_vs_t", w.plot1D("t", "sigma_z")),                       # velocity bunching → waist
        ("emittance_vs_t",   w.plot1D("t", "norm_emit_x")),
    ]
    for name, fig in figs:
        fig.savefig(f"{RESULTS}/injector_{name}.png", dpi=140, bbox_inches="tight")
        print(f"wrote {RESULTS}/injector_{name}.png")
    plt.close("all")


if __name__ == "__main__":
    main()
