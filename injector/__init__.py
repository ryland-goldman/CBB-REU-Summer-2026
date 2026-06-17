"""Cornell Linac injector stage facade: config()/run()/plot().

See injector/README.md for physics, parameters, and gotchas.
"""

from pipeline._runner import Stage

# Single source of truth for the default diags dir; imported by injector_sim.py so the
# parent can resolve the path without importing the pywarpx-laden sim module.
DEFAULT_OUTDIR = "injector/diags/main"

_stage = Stage(
    name="injector",
    build_module="injector.build_injector_field",
    sim_module="injector.injector_sim",
    plot_module="injector.plot_injector",
)
config = _stage.config
run = _stage.run
plot = _stage.plot


def resolve_outdir():
    """Return the diags dir the next run() will write to (OUTDIR override or default)."""
    return _stage._params.get("OUTDIR") or DEFAULT_OUTDIR
