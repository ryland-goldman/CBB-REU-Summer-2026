"""Cornell Linac — Section 1 stage (WarpX RZ SLAC 3 m traveling-wave structure).

First of the downstream accelerating sections (later sections → `linac_sec2`, …).
Reads the prebuncher's bunched beam at its focus, focuses it with a solenoid, and
accelerates it through the SLAC Section-1 traveling wave to ~37 MeV.

Usage:
    import linac_sec1
    linac_sec1.config(POWER_MW=15, PHASE_DEG=0, I_SOL=400)   # optional overrides
    linac_sec1.run()                                          # build field + sim + plots
    linac_sec1.run(plots=False)                               # build + sim only
    linac_sec1.plot()                                         # plots only

Parameter names match the module-level constants in `linac_sec1/build_linac_sec1_field.py`
and `linac_sec1/linac_sec1_sim.py`.
"""

import os

from pipeline._runner import Stage

# Default case directory for a bare run() (the headline operating point). The scan
# driver in run_pipeline.py sets OUTDIR explicitly per case (scan_phi*, focusoff).
DEFAULT_OUTDIR = "linac_sec1/diags/main"

_stage = Stage(
    name="linac_sec1",
    build_module="linac_sec1.build_linac_sec1_field",
    sim_module="linac_sec1.linac_sec1_sim",
    plot_module="linac_sec1.plot_linac_sec1",
)
config = _stage.config
run = _stage.run
plot = _stage.plot


def resolve_outdir():
    """Return the diags dir the next run() will write to (OUTDIR override or default).

    Used by `pipeline/run_pipeline.py` so the final-beam summary reads the same
    directory the sim wrote, without importing the pywarpx-laden sim module.
    """
    return _stage._params.get("OUTDIR") or DEFAULT_OUTDIR


def _case_metrics(diag):
    """(capture fraction, mean exit KE [MeV]) at a case's last snapshot, or (0, 0)."""
    import numpy as np
    from openpmd_viewer import OpenPMDTimeSeries
    pdir = os.path.join(diag, "particles")
    if not os.path.isdir(pdir):
        return 0.0, 0.0
    ts = OpenPMDTimeSeries(pdir)
    its = list(ts.iterations)
    if not its:
        return 0.0, 0.0
    _, _, _, _, _, _, w0 = ts.get_particle(
        ["x", "y", "z", "ux", "uy", "uz", "w"], species="electrons", iteration=its[0])
    _, _, _, ux, uy, uz, w = ts.get_particle(
        ["x", "y", "z", "ux", "uy", "uz", "w"], species="electrons", iteration=its[-1])
    if len(w) < 5:
        return 0.0, 0.0
    ke = float(np.average((np.sqrt(1 + ux**2 + uy**2 + uz**2) - 1) * 0.51099895, weights=w))
    return float(w.sum() / w0.sum()), ke


def demo(phases=None, scan_nz=1024, full_nz=1664, i_sol=1000.0):
    """Run the full Section-1 demonstration and write all six figures.

    1. RF-phase acceptance scan (focus on, reduced resolution) → pick the crest;
    2. headline run at the crest (full resolution, focus on) → ``diags/main``;
    3. focus-off comparison at the same phase (I=0)          → ``diags/focusoff``;
    4. aggregate plots.

    Scan cases land in ``diags/scan_phi<deg>``. POWER_MW (and any other config set
    before calling) is preserved. Tune the grid via ``phases`` / resolution via
    ``scan_nz``/``full_nz``.
    """
    if phases is None:
        phases = list(range(-150, 181, 30))
    print(f"\n[linac_sec1.demo] RF-phase acceptance scan over {phases} deg "
          f"(I_sol={i_sol:g} A)")
    scan = []
    for ph in phases:
        config(PHASE_DEG=ph, I_SOL=i_sol, NZ=scan_nz, N_DIAGS=10,
               OUTDIR=f"linac_sec1/diags/scan_phi{ph:d}")
        run(plots=False)
        cap, ke = _case_metrics(f"linac_sec1/diags/scan_phi{ph:d}")
        scan.append((ph, cap, ke))
        print(f"  scan phase {ph:+4d}deg: capture {cap*100:4.0f}%  mean KE {ke:5.1f} MeV")
    best = max(scan, key=lambda t: t[2])[0]          # crest = max mean energy
    print(f"[linac_sec1.demo] crest phase ≈ {best:+d} deg; headline + focus-off runs")

    config(PHASE_DEG=best, I_SOL=i_sol, NZ=full_nz, N_DIAGS=60,
           OUTDIR="linac_sec1/diags/main")
    run(plots=False)
    config(I_SOL=0.0, OUTDIR="linac_sec1/diags/focusoff")
    run(plots=False)

    # restore the headline operating point so resolve_outdir()/plot() point at main
    config(PHASE_DEG=best, I_SOL=i_sol, NZ=full_nz, OUTDIR="linac_sec1/diags/main")
    plot()
