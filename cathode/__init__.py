"""Cornell Linac cathode stage — WarpX 2D Child–Langmuir diode.

Usage:
    import cathode
    cathode.config(V_anode=50.0, gap_d=100e-6)   # optional overrides
    cathode.run()                                  # sim + plots
    cathode.run(plots=False)                       # sim only
    cathode.plot()                                 # plots only

Parameter names match the module-level constants at the top of `cathode/cathode_diode.py`.
"""

from pipeline._runner import Stage

_stage = Stage(
    name="cathode",
    sim_module="cathode.cathode_diode",
    plot_module="cathode.plot_cathode",
)
config = _stage.config
run = _stage.run
plot = _stage.plot
