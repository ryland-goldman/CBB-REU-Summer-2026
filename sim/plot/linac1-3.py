"""
Figures for the WarpX RZ Cornell Linac sections 1-3 (sim/linac1-3.py), generated entirely with
lume-warpx's plotting helpers (WarpX.plot2D / plot1D) over logs/diags/linac1-3/secN/main/.
Writes PNGs to logs/plots/linac1-3/. (No field diagnostic is dumped, so no plot_fields.)

Run as:  python sim/plot/linac1-3.py <N>   with N in {1, 2, 3}.

The four figures per section show the capture / acceleration physics:
  phase_space_z_KE — the captured/accelerated slice (sec 1 ~21 MeV, sec 2 ~46 MeV, sec 3 ~78 MeV),
  transverse_x_px  — the exit transverse phase space within the bore,
  centroid_vs_t    — the bunch crossing the 3 m structure + drift,
  emittance_vs_t   — transverse emittance over the run.

main() runs ONLY plotting (the section sim must have been run first).
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from openpmd_viewer import OpenPMDTimeSeries

RESULTS = "logs/plots/linac1-3"


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
    if len(sys.argv) < 2 or sys.argv[1] not in ("1", "2", "3"):
        sys.exit("usage: python sim/plot/linac1-3.py <N>   with N in {1, 2, 3}")
    N = int(sys.argv[1])

    from warpx import WarpX
    os.makedirs(RESULTS, exist_ok=True)

    config = f"config/linac{N}.yaml"
    diag_dir = f"logs/diags/linac1-3/sec{N}/main"
    particles = os.path.join(diag_dir, "particles")

    # The sim overrides write_dir to logs/diags/linac1-3/secN/main, so the particles series sits
    # directly there (not under <path>/diags). Populate w._outputs by hand the way the cathode
    # plotter does, so plot2D/plot1D work without WarpX.load_output (which needs a path= on w).
    w = WarpX(input_file=config)
    w._outputs = {"particles": OpenPMDTimeSeries(particles)}
    w._diag_dir = diag_dir
    w._output = w._outputs["particles"]
    it = _last_populated(particles)

    figs = [
        ("phase_space_z_KE", w.plot2D("z", "kinetic_energy", iteration=it)),
        ("transverse_x_px",  w.plot2D("x", "px", iteration=it)),
        ("centroid_vs_t",    w.plot1D("t", "mean_z")),
        ("emittance_vs_t",   w.plot1D("t", "norm_emit_x")),
    ]
    for name, fig in figs:
        out = f"{RESULTS}/sec{N}_{name}.png"
        fig.savefig(out, dpi=140, bbox_inches="tight")
        print(f"wrote {out}")
    plt.close("all")


if __name__ == "__main__":
    main()
