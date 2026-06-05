"""Cornell Linac prebuncher stage — WarpX RZ standing-wave TM cavity.

Usage:
    import prebuncher
    prebuncher.config(POWER_W=800, PHASE="zc",
                      OUTDIR="prebuncher/diags/P800_zc")    # optional overrides
    prebuncher.run()                                         # build field + sim + plots
    prebuncher.run(plots=False)                              # build + sim only
    prebuncher.plot()                                        # plots only

Parameter names match the module-level constants in `prebuncher/build_prebuncher_field.py` and
`prebuncher/prebuncher_sim.py`.
"""

from pipeline._runner import Stage

# Defaults mirrored by prebuncher_sim.py — single source of truth so the parent
# can resolve OUTDIR without importing the sim module (which would pull in pywarpx
# and defeat the per-stage subprocess isolation).
DEFAULT_POWER_W = 800.0
DEFAULT_PHASE = "zc"


def _derive_outdir(power_w, phase):
    # `:g` keeps whole watts compact (800.0 -> "P800_zc", matching the shipped
    # config) while preserving fractional operating points (160.5 -> "P160.5_zc")
    # so a power scan over non-integer watts can't collide distinct cases into one
    # OUTDIR. Integer-valued int() would truncate 160.5 and 160.9 to the same dir.
    return f"prebuncher/diags/P{power_w:g}_{phase}"


_stage = Stage(
    name="prebuncher",
    build_module="prebuncher.build_prebuncher_field",
    sim_module="prebuncher.prebuncher_sim",
    plot_module="prebuncher.plot_prebuncher",
)
config = _stage.config
run = _stage.run
plot = _stage.plot


def resolve_outdir():
    """Return the OUTDIR the next run() will write to, given the current config().

    OUTDIR explicitly set via config() wins; otherwise derive from POWER_W/PHASE
    (falling back to the module defaults). Used by `pipeline/run_pipeline.py` so
    the final-beam summary reads the same directory the sim wrote.
    """
    p = _stage._params
    if p.get("OUTDIR"):
        return p["OUTDIR"]
    return _derive_outdir(p.get("POWER_W", DEFAULT_POWER_W),
                          p.get("PHASE", DEFAULT_PHASE))
