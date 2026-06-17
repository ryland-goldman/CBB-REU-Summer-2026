"""Injector->linac iris/pipe collimation as a multi-plane particle-id scrape.

Post-hoc collimation on the openPMD dumps (this pywarpx RZ build cannot scrape in-run:
the position SoA accessor raises "Component x does not exist"). A particle outside the
9.547 mm pipe at ANY plane from z=1.922 m on is scraped, by id-tracking across dumps —
NOT a single 2.03 m cut, since the beam converges across the 1.922->2.03 m tail.
Shared by ``linac_sec1.load_injector_bunch`` and ``injector._report_collimated_handoff``.

See injector/README.md -> "The 9.547 mm collimator" for physics and measured numbers.
"""

import numpy as np


def pipe_violator_ids(ts, scan_iterations, collim_r, z_iris, species="electrons"):
    """Union of ids scraped by the 9.547 mm pipe over ``scan_iterations``.

    A particle is a violator if its own z >= ``z_iris`` and r = hypot(x, y) > ``collim_r``
    in ANY scanned dump. No upper z bound needed (the domain absorbs anything past ZMAX).
    """
    violators = set()
    for it in scan_iterations:
        idv, xv, yv, zv = ts.get_particle(
            ["id", "x", "y", "z"], species=species, iteration=it)
        r = np.hypot(xv, yv)
        bad = (zv >= z_iris) & (r > collim_r)
        if bad.any():
            violators.update(idv[bad].tolist())
    return violators


def survivor_mask(ids, violator_ids):
    """Boolean mask over ``ids`` (True = survives): id not in ``violator_ids``."""
    ids = np.asarray(ids)
    if not violator_ids:
        return np.ones(ids.shape, dtype=bool)
    return ~np.isin(ids, np.fromiter(violator_ids, dtype=ids.dtype))
