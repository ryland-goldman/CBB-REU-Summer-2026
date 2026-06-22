"""
Figures for the finite-cathode space-charge-limited diode (cathode_diode.py),
generated entirely with lume-warpx's plotting helpers (WarpX.plot2D / plot_fields /
plot1D) over cathode/diags/. Writes PNGs to cathode/results/.

See cathode/README.md for the physics each figure shows.
"""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from openpmd_viewer import OpenPMDTimeSeries

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG = os.path.join(HERE, "cathode.yaml")
RESULTS = "cathode/results"


def _last_populated(diag, species="electrons"):
    """Last diagnostic iteration that actually has particles (the exit beam can be sparse)."""
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

    w = WarpX(input_file=CONFIG, path="cathode")
    w.load_output()                                    # cathode/diags/{fields, particles}
    it = _last_populated("cathode/diags/particles")

    figs = [
        ("phase_space_z_KE", w.plot2D("z", "kinetic_energy", iteration=it)),
        ("transverse_x_px",  w.plot2D("x", "px", iteration=it)),
        ("potential_xz",     w.plot_fields("phi", "x", "z")),
        ("charge_density_xz", w.plot_fields("rho", "x", "z")),
        ("centroid_vs_t",    w.plot1D("t", "mean_z")),
        ("charge_vs_t",      w.plot1D("t", "charge")),
    ]
    for name, fig in figs:
        fig.savefig(f"{RESULTS}/{name}.png", dpi=140, bbox_inches="tight")
        print(f"wrote {RESULTS}/{name}.png")
    plt.close("all")


if __name__ == "__main__":
    main()
