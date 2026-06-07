"""Cornell Linac — Section 1 stage (WarpX RZ SLAC 3 m traveling-wave structure).

First of the downstream accelerating sections (later sections → `linac_sec2`, …).
Reads the injector's focused, velocity-bunched beam at the z ≈ 2.03 m handoff plane
(already collimated to the 9.547 mm iris) and accelerates it through the SLAC
Section-1 traveling wave. Transverse focusing lives upstream in the injector now —
the linac carries no solenoid.

Usage:
    import linac_sec1
    linac_sec1.config(POWER_MW=11, PHASE_DEG=0)             # optional overrides
    linac_sec1.run()                                        # build field + sim + plots
    linac_sec1.run(plots=False)                             # build + sim only
    linac_sec1.plot()                                       # plots only

The default operating point (top of `linac_sec1_sim.py`) is the original LinacSim
section-1 setting (`PHASE_DEG=0`, `POWER_MW=11`). Capture is reported against the true
injected charge through the 9.547 mm iris (a conservative lower bound — the lab-frame ES
self-field overestimates transverse space charge by ~γ²).

Parameter names match the module-level constants in `linac_sec1/build_linac_sec1_field.py`
and `linac_sec1/linac_sec1_sim.py`.
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
