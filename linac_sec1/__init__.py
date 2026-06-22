"""Cornell Linac Section 1 stage facade: run/plot over the WarpX RZ SLAC 3 m
traveling-wave structure (driven via lume-warpx). Constants live in
linac_sec1/linac_sec1.yaml (edit it to retune); config() is not the knob API here.

See linac_sec1/README.md for physics, parameters, and gotchas.
"""

from pipeline._runner import Stage

# Diags dir for run() (the single operating point).
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
