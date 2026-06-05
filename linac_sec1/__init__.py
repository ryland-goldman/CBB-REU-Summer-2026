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
