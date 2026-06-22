"""
Figures for the WarpX RZ SLAC Linac Section 1 (linac_sec1_sim.py), generated entirely
with lume-warpx's plotting helpers (WarpX.plot2D / plot1D) over linac_sec1/diags/main/.
Writes PNGs to linac_sec1/results/. (No field diagnostic is dumped, so no plot_fields.)

See linac_sec1/README.md for the physics each figure shows (capture + acceleration to
~21 MeV).
"""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from openpmd_viewer import OpenPMDTimeSeries

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG = os.path.join(HERE, "linac_sec1.yaml")
RESULTS = "linac_sec1/results"
DIAG = "linac_sec1/diags/main/particles"


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

    w = WarpX(input_file=CONFIG, path="linac_sec1")
    w.load_output(diag_dir=DIAG)
    it = _last_populated(DIAG)

    figs = [
        ("phase_space_z_KE", w.plot2D("z", "kinetic_energy", iteration=it)),   # capture to ~21 MeV
        ("transverse_x_px",  w.plot2D("x", "px", iteration=it)),
        ("centroid_vs_t",    w.plot1D("t", "mean_z")),
        ("emittance_vs_t",   w.plot1D("t", "norm_emit_x")),
    ]
    for name, fig in figs:
        fig.savefig(f"{RESULTS}/linac_sec1_{name}.png", dpi=140, bbox_inches="tight")
        print(f"wrote {RESULTS}/linac_sec1_{name}.png")
    plt.close("all")


if __name__ == "__main__":
    main()
