"""Cornell Linac cathode stage facade — WarpX 2D Child–Langmuir diode.

Exposes config()/run()/plot() over the Stage runner.
See cathode/README.md for physics, parameters, and gotchas.
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
