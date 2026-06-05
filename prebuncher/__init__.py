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

_stage = Stage(
    name="prebuncher",
    build_module="prebuncher.build_prebuncher_field",
    sim_module="prebuncher.prebuncher_sim",
    plot_module="prebuncher.plot_prebuncher",
)
config = _stage.config
run = _stage.run
plot = _stage.plot
