"""
Figures for the WarpX RZ CESR-gun simulation, generated entirely with lume-warpx's
plotting helpers (WarpX.plot2D / plot_fields / plot1D) over gun/diags/. Writes PNGs to
gun/results/.

See gun/README.md for the physics each figure shows.
"""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from openpmd_viewer import OpenPMDTimeSeries

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG = os.path.join(HERE, "gun.yaml")
RESULTS = "gun/results"


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

    w = WarpX(input_file=CONFIG, path="gun")
    it = _last_populated("gun/diags/particles")

    # Load the particles and fields series explicitly (gun/diags also holds a `handoff` series).
    w.load_output(diag_dir="gun/diags/particles")
    figs = [
        ("phase_space_z_KE", w.plot2D("z", "kinetic_energy", iteration=it)),   # gain to ~150 keV
        ("transverse_x_px",  w.plot2D("x", "px", iteration=it)),
        ("centroid_vs_t",    w.plot1D("t", "mean_z")),
        ("emittance_vs_t",   w.plot1D("t", "norm_emit_x")),
    ]
    w.load_output(diag_dir="gun/diags/fields")
    figs += [
        ("Ez_rz",          w.plot_fields("E", "z", "r")),                       # |E| in the gun gap
        ("self_charge_rz", w.plot_fields("rho", "z", "r")),                     # beam self charge
    ]
    for name, fig in figs:
        fig.savefig(f"{RESULTS}/{name}.png", dpi=140, bbox_inches="tight")
        print(f"wrote {RESULTS}/{name}.png")
    plt.close("all")


if __name__ == "__main__":
    main()
