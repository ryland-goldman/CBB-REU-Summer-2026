"""Injector→linac iris/pipe collimation as a multi-plane particle-id scrape.

The LinacSim injector→linac aperture is a 9.547 mm `scatteriris` at z = 1.922 m
followed by a 9.547 mm beam pipe to 2.1 m (gpt_master.in). A particle is scraped if
its radius exceeds the iris anywhere from 1.922 m onward. This pywarpx RZ build cannot
scrape particles in-run (the position SoA accessor raises "Component x does not exist"),
so collimation is applied **post-hoc** on the openPMD dumps.

A SINGLE radial cut at the 2.03 m handoff plane — what the code used to do — is WRONG.
The Sol 0 / Lens 0E matching telescope focuses the beam **hard** across the 1.922→2.03 m
tail (measured on the faithful 6/40/10 A tune: in-iris 38 % @1.92 m → 93 % @2.03 m;
σ_r 12.4 → 4.9 mm). So a particle OUTSIDE the iris at 1.922 m — which the real machine's
aperture scrapes — can converge back INSIDE it by 2.03 m and be wrongly counted as
transmitted. The physical aperture is at the iris ENTRANCE (and continuously down the
pipe), not at the focus.

We therefore emulate the continuous pipe by tracking particle ``id`` across every dump
that samples the pipe region: a particle whose own z ≥ z_iris with r > collim_r in ANY
dump hit the wall and is removed. This is exact in the dense-dump limit and reduces to
the entrance-plane cut when the envelope is monotone. WarpX always writes the ``id``
record in its openPMD particle output, so the tracking needs no extra ``data_list`` entry.

Used by ``linac_sec1.load_injector_bunch`` (the physical cut: only survivors are injected)
and ``injector._report_collimated_handoff`` (the sanity-log transmission number), so the
two stages cannot drift on how the iris is applied.
"""

import numpy as np


def pipe_violator_ids(ts, scan_iterations, collim_r, z_iris, species="electrons"):
    """Set of particle ids scraped by the 9.547 mm pipe over ``scan_iterations``.

    For each dump in ``scan_iterations`` (the caller selects the dumps whose ⟨z⟩ lies in
    the pipe region, so we only read the relevant snapshots), a particle is a violator if
    its OWN z ≥ ``z_iris`` and its radius r = √(x²+y²) > ``collim_r``. Returns the union of
    violator ids across all scanned planes — a particle outside the aperture at ANY plane
    in the pipe hit the wall. (No upper z bound is needed: the domain absorbs anything past
    ZMAX, and the pipe runs from ``z_iris`` to the absorbing wall.)
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
    """Boolean mask over ``ids`` (a handoff-dump id array) — True = passes the iris.

    A particle survives iff its id is not in ``violator_ids`` (the pipe scrape set from
    ``pipe_violator_ids``). Empty violator set → everything survives.
    """
    ids = np.asarray(ids)
    if not violator_ids:
        return np.ones(ids.shape, dtype=bool)
    return ~np.isin(ids, np.fromiter(violator_ids, dtype=ids.dtype))
