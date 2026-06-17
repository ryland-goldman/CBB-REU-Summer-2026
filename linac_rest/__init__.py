"""Cornell Linac sections 2–8 stage facade (Impact-T, in-process).

Exposes config/run/plot over an ImpactStage. config() overrides match the
module-level constants in build_linac_rest_lattice.py and linac_rest_sim.py.

See linac_rest/README.md for physics, parameters, and gotchas.
"""

from pipeline._impact_runner import ImpactStage

DEFAULT_OUTDIR = "linac_rest/diags/main"

_stage = ImpactStage(
    name="linac_rest",
    build_module="linac_rest.build_linac_rest_lattice",
    sim_module="linac_rest.linac_rest_sim",
    plot_module="linac_rest.plot_linac_rest",
)
config = _stage.config
run = _stage.run
plot = _stage.plot


def resolve_outdir():
    """Return the diags dir the next run() will write to (OUTDIR override or default).

    Lets run_pipeline.py find the sim's output dir without importing the
    lume-impact-laden sim module.
    """
    return _stage._params.get("OUTDIR") or DEFAULT_OUTDIR
