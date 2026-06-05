"""Cornell Linac gun stage — WarpX RZ CESR gun.

Usage:
    import gun
    gun.config(GUN_VOLTAGE=150e3, BUNCH_CHARGE=0.1e-9)   # optional overrides
    gun.run()                                             # build field + sim + plots
    gun.run(plots=False)                                  # build + sim only
    gun.plot()                                            # plots only

Parameter names match the module-level constants in `gun/build_gun_field.py` and `gun/gun_sim.py`.
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
