"""Impact-T ↔ WarpX-openPMD adapters for the `linac_rest/` stage.

The rest-of-linac stage is an **Impact-T** run driven through lume-impact, not a
pywarpx run, so its beam crosses two format boundaries that the WarpX stages don't:

  * **Handoff IN**  — read `linac_sec1`'s captured exit beam (a WarpX-written
    openPMD series) into a `pmd_beamphysics.ParticleGroup` that Impact-T ingests
    via `I.initial_particles`.  (`read_warpx_dump` here.)
  * **Handoff OUT** — write Impact-T's output `ParticleGroup` back into the *exact*
    WarpX openPMD layout the cross-stage tools expect, so `plot_chain` /
    `_beam_summary` / `build_moment_table` read `linac_rest/diags/main/particles`
    with the same `get_particle([...], species="electrons")` call they use on every
    WarpX stage.  (`write_openpmd_particles` here.)

Why a hand-rolled writer and not `ParticleGroup.write()`:
    `ParticleGroup.write()` emits openPMD **2.0** with a STRING `openPMDextension`
    attribute. `openpmd-viewer` (what `plot_chain`/`_beam_summary` use) only accepts
    the ED-PIC extension as an **integer** (1) and the 1.x base standard, so it
    rejects a `ParticleGroup.write()` file outright. We replicate WarpX's own
    openpmd-api output byte-layout instead (verified against
    injector/linac_sec1 dumps): openPMD "1.1.0", openPMDextension=1 (int), records
    position/momentum/weighting (+ charge/mass/id), groupBased encoding.

THE CONTRACT (matched to WarpX exactly — every value verified against a WarpX dump):
  * species name is literally **"electrons"** (PLURAL) in the openPMD output. NOTE
    the asymmetry: `ParticleGroup.species` is "electron" (SINGULAR); the WarpX
    stages and every cross-stage reader key on "electrons". The writer translates.
  * record `position`  : x,y,z components in METERS (unitSI=1).
  * record `momentum`  : x,y,z = γβ·m_e·c  [kg·m/s] (unitSI=1). openpmd-viewer
    divides momentum by (mass-record × c) to recover ``ux=uy=uz=γβ`` (dimensionless),
    which is the WarpX `u` convention every downstream metric assumes — do NOT write
    bare γβ here or every ⟨KE⟩/emittance downstream is wrong by a factor m_e·c.
  * record `weighting` : macroparticle COUNT (physical electrons per macroparticle),
    NOT charge. `get_particle(["w"])` returns this; `q = w·q_e`. `ParticleGroup.weight`
    is a CHARGE [C], so the writer converts w = weight / e.
"""

import numpy as np
import openpmd_api as io

# Physical constants (SI), matched to the values the WarpX stages / plot_chain use.
M_E = 9.1093837015e-31          # electron rest mass [kg]
C = 299792458.0                 # speed of light [m/s]
Q_E = 1.602176634e-19           # elementary charge [C]
MC2_EV = 510998.95069           # electron rest energy [eV] (ParticleGroup mass unit)


def _is_electron_species(pg):
    """ParticleGroup uses the SINGULAR 'electron'; tolerate either spelling."""
    return str(getattr(pg, "species", "electron")).lower().startswith("electron")


def write_openpmd_particles(pg, out_dir, iteration=0, time=0.0):
    """Write a `ParticleGroup` to `out_dir` as a WarpX-style openPMD particle dump.

    Emits ``<out_dir>/openpmd_%06T.h5`` (groupBased, same filename pattern WarpX
    uses) with species ``"electrons"`` and records position [m], momentum
    [γβ·m_e·c, kg·m/s], weighting [count], plus charge/mass/id — readable by
    BOTH ``OpenPMDTimeSeries(out_dir).get_particle([...], species="electrons")``
    and the pipeline's `_beam_summary`/`build_moment_table`.

    `pg`        : `pmd_beamphysics.ParticleGroup` (Impact-T `final_particles`, or
                  any along-z slice). Its `.weight` is CHARGE [C]; converted to a
                  count here. `.px/.py/.pz` are [eV/c]; converted to γβ·m_e·c.
    `iteration` : openPMD iteration index (use a per-slice index when dumping
                  several z-slices so they form a series).
    `time`      : iteration time [s] (Impact-T output z is local-frame; the
                  cross-stage z0 shift is applied in plot_chain, not here).

    Returns the written file path.
    """
    import os
    os.makedirs(out_dir, exist_ok=True)
    n = pg.n_particle
    if n == 0:
        raise ValueError("write_openpmd_particles: ParticleGroup is empty")

    # γβ per component: px[eV/c] / mass[eV] = γβ_x. Then ×m_e·c → kg·m/s so the
    # viewer's momentum/(mass·c) recovers γβ. (Equivalent: momentum = px[eV/c]
    # converted to SI; we go via γβ to keep the mass/charge constants in one place.)
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

    # weighting = macroparticle COUNT (physical electrons per macro). PG.weight is
    # a per-particle CHARGE [C]; q/e = count. (Sum of w · q_e == total charge.)
    w = np.asarray(pg.weight, dtype=np.float64) / Q_E

    series = io.Series(os.path.join(out_dir, "openpmd_%06T.h5"),
                       io.Access.create)
    # Match WarpX's standard/extension exactly: openpmd-viewer accepts the 1.x base
    # standard + INTEGER ED-PIC extension; a ParticleGroup.write() 2.0/STRING-ext
    # file is rejected. openpmd-api sets `openPMD`/`openPMDextension` from these.
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
        """Set the ED-PIC per-record attributes openpmd-viewer REQUIRES.

        The viewer's particle reader reads `macroWeighted` + `weightingPower` off
        every record (and raises a bare `Error: macroWeighted` if absent — this is
        what blocks a naive openpmd-api write). WarpX sets macroWeighted=0 on every
        record except `weighting` (which is the weight itself, macroWeighted=1), and
        weightingPower=1 on the additive/extensive records (momentum, charge, mass,
        weighting), 0 on position/positionOffset/id."""
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

    The handoff-IN reader for `linac_rest`: select `linac_sec1`'s EXIT dump (the
    largest-⟨z⟩ iteration by default — the captured coasting beam) and return it
    as a `pmd_beamphysics.ParticleGroup` (species "electron", t-coordinates) ready
    for `I.initial_particles`. The caller is responsible for `drift_to_t()` /
    zeroing z (Impact-T wants z==0 at injection) — kept out of here so the reader
    stays a pure format adapter.

    `iteration` : openPMD iteration to read; default = last (exit) dump.
    Returns a `pmd_beamphysics.ParticleGroup` with `.charge` == Σw·q_e of the dump.
    """
    from pmd_beamphysics import ParticleGroup
    from openpmd_viewer import OpenPMDTimeSeries

    ts = OpenPMDTimeSeries(particles_dir)
    if len(ts.iterations) == 0:
        raise RuntimeError(f"{particles_dir} has no iterations")
    it = ts.iterations[-1] if iteration is None else iteration

    x, y, z, ux, uy, uz, w = ts.get_particle(
        ["x", "y", "z", "ux", "uy", "uz", "w"], species=species, iteration=it)

    # ux/uy/uz are γβ (dimensionless). ParticleGroup wants px/py/pz in eV/c:
    # γβ · m_e c² [eV] = px [eV/c]. weight is a per-macro CHARGE [C] = w · q_e.
    data = dict(
        x=x, y=y, z=z,
        px=ux * MC2_EV, py=uy * MC2_EV, pz=uz * MC2_EV,
        t=np.zeros_like(x),
        weight=w * Q_E,
        status=np.ones_like(x, dtype=int),
        species="electron",
    )
    return ParticleGroup(data=data)
