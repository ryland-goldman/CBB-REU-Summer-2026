"""Cornell Linac injector stage — WarpX RZ two-prebuncher + solenoid drift.

Collapses the full LinacSim injector subsection (Lens 0A → Prebuncher 1 →
Prebuncher 2 → Sol 0 / Lens 0E) into one self-consistent RZ space-charge run,
handing a focused, velocity-bunched beam to `linac_sec1` at z ≈ 2.03 m.

Usage:
    import injector
    injector.config(PREB1_KW=8, PREB2_KW=10)                 # optional overrides
    injector.run()                                            # build field + sim + plots
    injector.run(plots=False)                                 # build + sim only
    injector.plot()                                           # plots only

    # Optional power/phase scan (writes one OUTDIR per case, default is diags/main):
    injector.config(PREB1_KW=300, OUTDIR="injector/diags/P300_zc")

Parameter names match the module-level constants in `injector/build_injector_field.py`
and `injector/injector_sim.py`.
"""

from pipeline._runner import Stage

# Diags dir for the default chain run (the single operating point). A scan overrides
# OUTDIR per case via config(). Single source of truth — imported by injector_sim.py
# so the parent can resolve the path without importing the pywarpx-laden sim module.
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
    """Return the diags dir the next run() will write to (OUTDIR override or default).

    OUTDIR explicitly set via config() wins; otherwise the default `injector/diags/main`.
    Used by `pipeline/run_pipeline.py` so the final-beam summary reads the same
    directory the sim wrote, without importing the pywarpx-laden sim module.
    """
    return _stage._params.get("OUTDIR") or DEFAULT_OUTDIR
