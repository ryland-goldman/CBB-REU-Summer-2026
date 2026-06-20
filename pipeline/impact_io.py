"""Impact-T <-> WarpX-openPMD particle adapters for the `linac_rest/` stage.

`read_warpx_dump` reads linac_sec1's exit dump into a ParticleGroup for Impact-T;
`write_openpmd_particles` writes Impact-T output back in WarpX's exact openPMD
layout so cross-stage readers ingest it unchanged.
See pipeline/README.md and .claude/CLAUDE.md for the linac_rest handoff contract.

Hand-rolled writer (not ParticleGroup.write): the latter emits openPMD 2.0 with a
STRING openPMDextension that openpmd-viewer rejects; we replicate WarpX's openpmd-api
byte-layout (openPMD 1.1.0, integer ED-PIC ext). Output species is "electrons"
(PLURAL) though ParticleGroup.species is "electron" (SINGULAR) — readers key plural.
"""

import numpy as np
import openpmd_api as io

# Physical constants — single-sourced from pipeline.constants (scipy) so no stage
# carries a divergent literal. C/M_E/Q_E kept as local aliases for the existing call sites.
from pipeline.constants import C_LIGHT as C, E_CHARGE as Q_E, M_E, MC2_EV


def _is_electron_species(pg):
    # ParticleGroup uses SINGULAR 'electron'; tolerate either spelling.
    return str(getattr(pg, "species", "electron")).lower().startswith("electron")


def write_openpmd_particles(pg, out_dir, iteration=0, time=0.0):
    """Write a `ParticleGroup` to `out_dir` as a WarpX-style openPMD particle dump.

    Emits ``<out_dir>/openpmd_%06T.h5`` with species "electrons"; records
    position [m], momentum [gamma*beta*m_e*c, kg*m/s], weighting [count], charge,
    mass, id. `pg.weight` is CHARGE [C] (converted to count); `pg.px/py/pz` are
    [eV/c]. Returns the written file path.
    See pipeline/README.md for the cross-stage contract.
    """
    import os
    os.makedirs(out_dir, exist_ok=True)
    n = pg.n_particle
    if n == 0:
        raise ValueError("write_openpmd_particles: ParticleGroup is empty")

    # px[eV/c]/mass[eV] = gamma*beta, then *m_e*c -> kg*m/s so viewer's
    # momentum/(mass*c) recovers gamma*beta (the WarpX u convention).
    gbx = np.asarray(pg.px, dtype=np.float64) / MC2_EV
    gby = np.asarray(pg.py, dtype=np.float64) / MC2_EV
    gbz = np.asarray(pg.pz, dtype=np.float64) / MC2_EV
    mom_scale = M_E * C
    px = gbx * mom_scale
    py = gby * mom_scale
    pz = gbz * mom_scale

    x = np.asarray(pg.x, dtype=np.float64)
    y = np.asarray(pg.y, dtype=np.float64)
    z = np.asarray(pg.z, dtype=np.float64)

    # weighting = macroparticle COUNT; PG.weight is per-particle CHARGE [C], so /q_e.
    w = np.asarray(pg.weight, dtype=np.float64) / Q_E

    series = io.Series(os.path.join(out_dir, "openpmd_%06T.h5"),
                       io.Access.create)
    # openpmd-viewer accepts 1.x base + INTEGER ED-PIC ext; rejects 2.0/STRING-ext.
    series.set_openPMD("1.1.0")
    series.set_openPMD_extension(1)            # ED-PIC (integer, NOT a string)
    series.set_software("linac_rest.impact_io")
    series.set_particles_path("particles")

    it = series.iterations[int(iteration)]
    it.set_time(float(time)).set_dt(1.0).set_time_unit_SI(1.0)

    sp = it.particles["electrons"]             # PLURAL — the cross-stage contract

    dset_f = io.Dataset(np.dtype("float64"), [n])
    dset_i = io.Dataset(np.dtype("int64"), [n])

    def _tag(record, weighting_power, macro_weighted=0):
        # ED-PIC per-record attrs openpmd-viewer REQUIRES (raises "Error:
        # macroWeighted" if absent). macroWeighted=1 only on weighting;
        # weightingPower=1 on the extensive records (momentum/charge/mass/weighting).
        record.set_attribute("macroWeighted", np.int32(macro_weighted))
        record.set_attribute("weightingPower", float(weighting_power))

    # position [m] (unitSI = 1)
    pos = sp["position"]
    pos.set_unit_dimension({io.Unit_Dimension.L: 1})
    _tag(pos, weighting_power=0)
    for comp, arr in (("x", x), ("y", y), ("z", z)):
        pos[comp].reset_dataset(dset_f)
        pos[comp].store_chunk(np.ascontiguousarray(arr))
        pos[comp].unit_SI = 1.0

    # positionOffset = 0 (WarpX writes it; some readers add position+offset).
    off = sp["positionOffset"]
    off.set_unit_dimension({io.Unit_Dimension.L: 1})
    _tag(off, weighting_power=0)
    zeros = np.zeros(n, dtype=np.float64)
    for comp in ("x", "y", "z"):
        off[comp].reset_dataset(dset_f)
        off[comp].store_chunk(np.ascontiguousarray(zeros))
        off[comp].unit_SI = 1.0

    # momentum = γβ·m_e·c  [kg·m/s] (unitDimension M·L·T⁻¹). viewer: u = mom/(mass·c).
    mom = sp["momentum"]
    mom.set_unit_dimension({io.Unit_Dimension.M: 1,
                            io.Unit_Dimension.L: 1,
                            io.Unit_Dimension.T: -1})
    _tag(mom, weighting_power=1)
    for comp, arr in (("x", px), ("y", py), ("z", pz)):
        mom[comp].reset_dataset(dset_f)
        mom[comp].store_chunk(np.ascontiguousarray(arr))
        mom[comp].unit_SI = 1.0

    # weighting = macroparticle count (dimensionless Scalar).
    wt = sp["weighting"][io.Mesh_Record_Component.SCALAR]
    sp["weighting"].set_unit_dimension({})
    _tag(sp["weighting"], weighting_power=1, macro_weighted=1)
    wt.reset_dataset(dset_f)
    wt.store_chunk(np.ascontiguousarray(w))
    wt.unit_SI = 1.0

    # charge [C] and mass [kg] per macroparticle (WarpX writes both as Scalars).
    ch = sp["charge"][io.Mesh_Record_Component.SCALAR]
    sp["charge"].set_unit_dimension({io.Unit_Dimension.T: 1, io.Unit_Dimension.I: 1})
    _tag(sp["charge"], weighting_power=1)
    ch.reset_dataset(dset_f)
    ch.store_chunk(np.ascontiguousarray(np.full(n, -Q_E, dtype=np.float64)))
    ch.unit_SI = 1.0

    ms = sp["mass"][io.Mesh_Record_Component.SCALAR]
    sp["mass"].set_unit_dimension({io.Unit_Dimension.M: 1})
    _tag(sp["mass"], weighting_power=1)
    ms.reset_dataset(dset_f)
    ms.store_chunk(np.ascontiguousarray(np.full(n, M_E, dtype=np.float64)))
    ms.unit_SI = 1.0

    # id: carry the ParticleGroup id if present, else 1..n (WarpX writes a Scalar id).
    if "id" in pg:
        ids = np.asarray(pg["id"], dtype=np.int64)
    else:
        ids = np.arange(1, n + 1, dtype=np.int64)
    sp["id"].set_unit_dimension({})
    _tag(sp["id"], weighting_power=0)
    idc = sp["id"][io.Mesh_Record_Component.SCALAR]
    idc.reset_dataset(dset_i)
    idc.store_chunk(np.ascontiguousarray(ids))
    idc.unit_SI = 1.0

    series.flush()
    del series                                 # close (openpmd-api flushes on destruct)
    return os.path.join(out_dir, "openpmd_%06T.h5".replace("%06T", f"{int(iteration):06d}"))


def read_warpx_dump(particles_dir, iteration=None, species="electrons"):
    """Read a WarpX-style openPMD particle dump into a `ParticleGroup`.

    Handoff-IN reader for `linac_rest`: returns the chosen iteration (default last,
    linac_sec1's exit dump) as a ParticleGroup ready for `I.initial_particles`. The
    caller does any `drift_to_t()` / z-zeroing (Impact-T wants z==0 at injection).
    See pipeline/README.md.
    """
    from pmd_beamphysics import ParticleGroup
    from openpmd_viewer import OpenPMDTimeSeries

    ts = OpenPMDTimeSeries(particles_dir)
    if len(ts.iterations) == 0:
        raise RuntimeError(f"{particles_dir} has no iterations")
    it = ts.iterations[-1] if iteration is None else iteration

    x, y, z, ux, uy, uz, w = ts.get_particle(
        ["x", "y", "z", "ux", "uy", "uz", "w"], species=species, iteration=it)

    # ux/uy/uz are gamma*beta; PG wants px/py/pz [eV/c] = gamma*beta*MC2_EV.
    # weight is per-macro CHARGE [C] = w*q_e.
    data = dict(
        x=x, y=y, z=z,
        px=ux * MC2_EV, py=uy * MC2_EV, pz=uz * MC2_EV,
        t=np.zeros_like(x),
        weight=w * Q_E,
        status=np.ones_like(x, dtype=int),
        species="electron",
    )
    return ParticleGroup(data=data)
