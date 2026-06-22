"""Cornell Linac gun stage facade — WarpX RZ CESR gun (driven via lume-warpx).

Exposes run()/plot() over the build+sim+plot scripts. Constants live in gun/gun.yaml
(edit it to retune); config() is not the knob API for this WarpX stage.
See gun/README.md for physics, parameters, and gotchas.
"""

from pipeline._runner import Stage

_stage = Stage(
    name="gun",
    build_module="gun.build_gun_field",
    sim_module="gun.gun_sim",
    plot_module="gun.plot_gun",
)
config = _stage.config
run = _stage.run
plot = _stage.plot
