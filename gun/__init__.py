"""Cornell Linac gun stage facade — WarpX RZ CESR gun.

Exposes config()/run()/plot() over the build+sim+plot scripts. config() keys
match the module-level constants in gun/build_gun_field.py and gun/gun_sim.py.
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
