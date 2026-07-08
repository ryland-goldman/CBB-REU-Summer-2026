"""
Figures for WarpX linac sections 1-4 (sim/linac1-4.py) from logs/diags/linac1-4/secN/main/.
Writes PNGs to logs/plots/linac1-4/. Run as: python sim/plot/linac1-4.py <N>  (N in 1-4).
See docs/linac1-4.md.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from openpmd_viewer import OpenPMDTimeSeries

from sim.plot import common as px

RESULTS = "logs/plots/linac1-4"


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
    if len(sys.argv) < 2 or sys.argv[1] not in ("1", "2", "3", "4"):
        sys.exit("usage: python sim/plot/linac1-4.py <N>   with N in {1, 2, 3, 4}")
    N = int(sys.argv[1])

    from warpx import WarpX
    os.makedirs(RESULTS, exist_ok=True)

    config = f"config/linac{N}.yaml"
    diag_dir = f"logs/diags/linac1-4/sec{N}/main"
    particles = os.path.join(diag_dir, "particles")

    # w has no path= set, so WarpX.load_output can't find the series; populate w._outputs by hand instead.
    ts = OpenPMDTimeSeries(particles)
    w = WarpX(input_file=config)
    w._outputs = {"particles": ts}
    w._diag_dir = diag_dir
    w._output = ts
    it = _last_populated(particles)
    pg = w._particle_group(iteration=it)

    x, y, ux, uy, wgt = ts.get_particle(["x", "y", "ux", "uy", "w"],
                                        species="electrons", iteration=it)
    z_m, ke, emit, sigma, q_pc = px.evolution_screens(px.pool_trajectories(ts, ts.iterations))

    figs = [
        ("phase_space_z_KE", w.plot2D("z", "kinetic_energy", iteration=it)),
        ("transverse_r_pr",  px.transverse_rpr(x, y, ux, uy, wgt, p_unit="MeV",
                             title=f"Section {N} exit transverse phase space  (r, p_r)")),
        ("energy_spectrum",  px.energy_spectrum(pg, e_unit="MeV")),
        ("evolution_vs_z",   px.evolution_vs_z(z_m, ke / 1e3, emit, sigma, charge_pc=q_pc,
                             ke_unit="MeV",
                             title=f"Section {N} beam evolution  (fixed-z virtual screens)")),
    ]
    for name, fig in figs:
        out = f"{RESULTS}/sec{N}_{name}.png"
        fig.savefig(out, dpi=140, bbox_inches="tight")
        print(f"wrote {out}")
    plt.close("all")


if __name__ == "__main__":
    main()
