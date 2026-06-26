"""Particle handoff IO shared across stages: read/write openPMD beams, build
`pmd_beamphysics.ParticleGroup`s, weighted downsampling, kinematics, the captured-core
cut between linac sections, the injector->linac iris scrape, and the Impact-T adapters.

Conventions: ux/uy/uz are dimensionless gamma*beta; ParticleGroup momentum px/py/pz is
eV/c (= gamma*beta*MC2_EV); `weight` is per-macroparticle CHARGE [C] (= count*q_e).
Stage-specific phase-space remaps (cathode 2D->RZ, gun->injector) live in their stage drivers.
"""

import json
import os

import numpy as np

from sim.helpers.tools import C_LIGHT, E_CHARGE as Q_E, M_E, MC2_EV


# ── openPMD particle series ──────────────────────────────────────────────────────
def open_particle_series(diag, stage_hint=None):
    """Open an openPMD particle series, raising if it has no iterations.
    `stage_hint` (e.g. "gun") tags the error with which upstream stage failed to write.
    """
    from openpmd_viewer import OpenPMDTimeSeries
    ts = OpenPMDTimeSeries(diag)
    if len(ts.iterations) == 0:
        hint = f" -- did the {stage_hint} stage run?" if stage_hint else ""
        raise RuntimeError(f"{diag} has no iterations{hint}")
    return ts


def anode_beam_mask(z, uz, gap_d, frac):
    """Forward-moving electrons in the top `frac` of the cathode gap -- the beam crossing the anode
    plane (the delivered flux). Excludes the dense near-cathode space-charge pileup and the reflected
    half of the over-injection that never exits. The cathode particle dumps are far sparser in time
    than a gap transit (~62 ps), so transits cannot be tracked as id-trajectory z-screens; the anode
    flux is taken as this crest-snapshot slab instead."""
    return (np.asarray(z) >= gap_d * (1.0 - frac)) & (np.asarray(uz) > 0.0)


def make_particle_group(x, y, z, ux, uy, uz, w):
    """ParticleGroup from RZ phase space: ux/uy/uz = gamma*beta, w = macroparticle count.

    Momentum px/py/pz = gamma*beta*MC2_EV [eV/c]; weight = count*q_e [C]; species "electron"
    (SINGULAR -- ParticleGroup's spelling; openPMD readers key the PLURAL "electrons").
    """
    from pmd_beamphysics import ParticleGroup
    n = x.size
    return ParticleGroup(data=dict(
        x=x, y=y, z=z,
        px=ux * MC2_EV, py=uy * MC2_EV, pz=uz * MC2_EV,
        t=np.zeros(n), weight=w * Q_E,
        status=np.ones(n, dtype=np.int64), species="electron"))


def downsample(arrays, w, max_part, rng):
    """Randomly thin to `max_part` macroparticles, rescaling weights to conserve charge.
    No-op if `max_part` is falsy or already <= the count. Returns (tuple of arrays, w).
    """
    n = w.size
    if not max_part or n <= max_part:
        return tuple(arrays), w
    sel = rng.choice(n, max_part, replace=False)
    return tuple(a[sel] for a in arrays), w[sel] * (n / max_part)


def resample(arrays, w, n_target, rng):
    """Resample to EXACTLY `n_target` macroparticles, conserving total charge.

    Downsamples without replacement when n > n_target; UPsamples WITH replacement when
    n < n_target (bootstrap — useful to refill a depleted beam to a fixed macroparticle
    count). Reweights the picks so sum(w) is preserved exactly (not just in expectation).
    No-op if `n_target` is falsy or already equal. Returns (tuple of arrays, w).

    Upsampling adds COINCIDENT duplicate macroparticles (no new phase-space information); it
    does NOT increase resolution. Safe only ahead of a stage that decorrelates them (the
    converter's stochastic showers) or with self-fields off — do not upsample into an SC-ON
    stage (cathode/gun/linac section 1), where duplicates would inject a spurious self-field.
    """
    n = w.size
    if not n_target or n == n_target:
        return tuple(arrays), w
    sel = rng.choice(n, n_target, replace=(n_target > n))
    w_sel = w[sel]
    return tuple(a[sel] for a in arrays), w_sel * (w.sum() / w_sel.sum())


def upsample_smeared(P, n_target, rng_seed=0, k_neighbors=8, smear=0.2,
                     smear_cols=("x", "y", "px", "py")):
    """Upsample a `ParticleGroup` to `n_target` macroparticles by KDE-style local smearing.

    Bootstrap-draws parents WITH replacement, then jitters each clone by `smear` x its distance
    to the `k_neighbors`-th nearest neighbour in std-whitened phase space, so the clones
    decorrelate and sample the LOCAL density. This is what `resample()` (plain bootstrap) cannot
    do: its coincident duplicates track identically under SC-off deterministic optics and add no
    statistics. Only `smear_cols` are jittered -- the default is the TRANSVERSE phase space
    (x, y, px, py), which governs aperture survival; pz/z/t are carried from the parent unchanged.
    Energy is thus APPROXIMATELY preserved (~0.5%) -- not exact, since total energy still depends on
    the smeared px/py -- so a clone of a low-pz / large-angle parent kept only by its transverse
    momentum CAN drop below an upstream KE floor; re-impose the cut downstream if it must hold.
    Total charge is preserved (uniform weight =
    q/n_target). Emittance inflation grows ~ smear^2, kept small by the per-particle (local)
    bandwidth -- dense regions get a small jitter, sparse halo a large one. No-op (returns P) if
    `n_target` is falsy or <= the current count. Returns a NEW group.

    SC-OFF ONLY: with self-fields on the synthetic clones would inject a spurious self-field.
    """
    from pmd_beamphysics import ParticleGroup
    from scipy.spatial import cKDTree
    n = P.n_particle
    if not n_target or n_target <= n:
        return P
    cols = ("x", "y", "z", "px", "py", "pz", "t")
    d = {c: np.asarray(getattr(P, c), dtype=float) for c in cols}
    feat = [c for c in smear_cols if d[c].std() > 0]           # smear only the requested live axes
    F = np.column_stack([d[c] for c in feat])
    mu, sd = F.mean(0), F.std(0)
    W = (F - mu) / sd                                           # per-axis std whitening
    dist, _ = cKDTree(W).query(W, k=min(k_neighbors + 1, n))
    h = dist[:, -1]                                             # whitened distance to the k-th neighbour
    rng = np.random.default_rng(rng_seed)
    par = rng.choice(n, n_target, replace=True)
    Fnew = (W[par] + smear * h[par][:, None] * rng.standard_normal((n_target, W.shape[1]))) * sd + mu
    out = {c: d[c][par].copy() for c in cols}                   # parents carry pz/z/t unchanged
    for j, c in enumerate(feat):
        out[c] = Fnew[:, j]
    q = float(P.charge)
    return ParticleGroup(data=dict(
        x=out["x"], y=out["y"], z=out["z"], px=out["px"], py=out["py"], pz=out["pz"],
        t=out["t"], weight=np.full(n_target, q / n_target),
        status=np.ones(n_target, dtype=int), species=P.species))


def beam_kinematics(ux, uy, uz, w):
    """(weighted mean v_z [m/s], weighted mean KE [keV]) from gamma*beta momenta."""
    gb = np.sqrt(1.0 + ux ** 2 + uy ** 2 + uz ** 2)        # gamma (ux/uy/uz are gamma*beta)
    v_beam = float(np.average(uz / gb, weights=w) * C_LIGHT)
    ke_mean_keV = float(np.average(gb - 1.0, weights=w) * MC2_EV / 1e3)
    return v_beam, ke_mean_keV


def load_warpx_exit_bunch(diag, label, max_part, rng_seed, z_inject, min_count=None,
                          core_ke_frac=0.5, resample_n=0):
    """Import an upstream WarpX section's EXIT beam (last well-populated dump) for the next
    linac section. Used by linac sections 2, 3 and 4 (no iris scrape -- that is the one-time
    injector->linac event at the section-1 entrance).

    Picks the last dump with >= `min_count` macroparticles (the captured beam coasting in the
    field-free exit drift, not a depleted boundary dump), keeps only the captured core
    (KE >= `core_ke_frac` * median KE), resizes the macroparticle count (reweighted), and shifts
    its tail to `z_inject`. `resample_n` > 0 forces EXACTLY that count (up- or down-sample, see
    `resample`); otherwise the core is downsampled to `max_part`. The core cut is essential: the
    section-exit dump trails a sparse slipping low-energy
    tail that lags the relativistic core by ~metres, not in the RF bucket -- genuinely lost
    between sections (same physics as the linac5-8 MIN_KE_MEV cut).

    Returns (bunch dict [gamma*beta momenta], v_beam [m/s], core <KE> [keV], info dict).
    info["exit_zmean_local_m"] is the read dump's <z> in the UPSTREAM local frame (for lab-z
    chaining via upstream_exit_lab_z); info["q_injected_C"] is the FULL exit charge.
    """
    ts = open_particle_series(diag, label)
    mc = min_count if min_count is not None else max(50, max_part // 50)
    it_exit = None
    for it in reversed([int(i) for i in ts.iterations]):
        z, = ts.get_particle(["z"], species="electrons", iteration=it)
        if len(z) >= mc:
            it_exit = it
            break
    if it_exit is None:
        raise RuntimeError(f"{diag}: no dump with >={mc} macroparticles -- did {label} run?")

    x, y, z, ux, uy, uz, w = ts.get_particle(
        ["x", "y", "z", "ux", "uy", "uz", "w"], species="electrons", iteration=it_exit)
    q_in = float(w.sum()) * Q_E                          # FULL exit charge (honest denominator)
    exit_zmean_local = float(np.average(z, weights=w))   # UPSTREAM local frame (lab-z chain)

    gb = np.sqrt(1.0 + ux ** 2 + uy ** 2 + uz ** 2)
    ke = (gb - 1.0) * MC2_EV                             # [eV]
    core = ke >= core_ke_frac * float(np.median(ke))
    if core.sum() < 50:
        raise RuntimeError(f"{diag}: only {int(core.sum())} core macroparticles above "
                           f"{core_ke_frac:g}x median KE -- upstream beam not accelerated?")
    x, y, z, ux, uy, uz, w = (a[core] for a in (x, y, z, ux, uy, uz, w))
    q_core = float(w.sum()) * Q_E

    rng = np.random.default_rng(rng_seed)
    if resample_n:                                       # fixed macroparticle count (up- or down-sample)
        (x, y, z, ux, uy, uz), w = resample((x, y, z, ux, uy, uz), w, resample_n, rng)
    else:
        (x, y, z, ux, uy, uz), w = downsample((x, y, z, ux, uy, uz), w, max_part, rng)
    z = z - z.min() + z_inject                           # core tail (smallest z) -> z_inject

    v_beam, ke_mean = beam_kinematics(ux, uy, uz, w)
    z_inject_mean = float(np.average(z, weights=w))
    sz = float(np.sqrt(np.average((z - z_inject_mean) ** 2, weights=w)))
    rmax = float(np.hypot(x, y).max())
    info = dict(it_exit=int(it_exit), n_injected=int(z.size), q_injected_C=q_in, q_core_C=q_core,
                core_frac=(q_core / q_in if q_in else 0.0),
                exit_zmean_local_m=exit_zmean_local, z_inject_mean_m=z_inject_mean,
                rmax_m=rmax, sigma_z_m=sz, ke_mean_keV=ke_mean)
    print(f"Injected {z.size} core macroparticles from {label} exit (iter {it_exit}; "
          f"{info['core_frac']*100:.0f}% of exit charge -- slipping tail dropped); "
          f"core <KE> {ke_mean/1e3:.2f} MeV, sigma_z {sz*1e3:.2f} mm, r_max {rmax*1e3:.2f} mm, "
          f"full-exit q {q_in*1e9:.4f} nC", flush=True)
    return dict(x=x, y=y, z=z, ux=ux, uy=uy, uz=uz, w=w), v_beam, ke_mean, info


def upstream_exit_lab_z(summary_path, exit_zmean_local, fallback=0.0):
    """Lab-frame z of an upstream linac section's exit, for chaining local frames.

    Reads the upstream section's injection_summary.json for its local->lab offset
    (z_inject_lab_m | z_handoff_m) - z_inject_mean_m, then adds the upstream exit dump's
    local <z>. Section 1 records z_handoff_m; sections 2/3 record z_inject_lab_m.
    """
    try:
        with open(summary_path) as fh:
            s = json.load(fh)
        z_inj_lab = s.get("z_inject_lab_m", s.get("z_handoff_m"))
        z0 = float(z_inj_lab) - float(s["z_inject_mean_m"])
        return z0 + exit_zmean_local
    except Exception:
        return exit_zmean_local + fallback


# ── injector -> linac iris collimation (multi-plane id scrape) ────────────────────
def pipe_violator_ids(ts, scan_iterations, collim_r, z_iris, species="electrons"):
    """Union of ids scraped by the pipe over `scan_iterations`.

    A particle is a violator if its own z >= `z_iris` and r = hypot(x, y) > `collim_r`
    in ANY scanned dump. Not a single cut -- the beam converges across the iris->handoff tail.
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
    """Boolean mask over `ids` (True = survives): id not in `violator_ids`."""
    ids = np.asarray(ids)
    if not violator_ids:
        return np.ones(ids.shape, dtype=bool)
    return ~np.isin(ids, np.fromiter(violator_ids, dtype=ids.dtype))


# ── Impact-T <-> WarpX-openPMD adapters (linac5-8) ───────────────────────────────
def read_warpx_dump(particles_dir, iteration=None, species="electrons"):
    """Read a WarpX-style openPMD particle dump into a `ParticleGroup` (handoff-IN reader
    for linac5-8; default last iteration = the upstream section's exit dump).
    """
    from pmd_beamphysics import ParticleGroup
    ts = open_particle_series(particles_dir)
    it = ts.iterations[-1] if iteration is None else iteration
    x, y, z, ux, uy, uz, w = ts.get_particle(
        ["x", "y", "z", "ux", "uy", "uz", "w"], species=species, iteration=it)
    return ParticleGroup(data=dict(
        x=x, y=y, z=z,
        px=ux * MC2_EV, py=uy * MC2_EV, pz=uz * MC2_EV,
        t=np.zeros_like(x), weight=w * Q_E,
        status=np.ones_like(x, dtype=int),
        species=species[:-1] if species.endswith("s") else species))


def write_openpmd_particles(pg, out_dir, iteration=0, time=0.0,
                            species="electrons", charge=-Q_E, mass=M_E):
    """Write a `ParticleGroup` to `out_dir` as a WarpX-style openPMD dump (handoff-OUT for
    linac5-8). Hand-rolled (not ParticleGroup.write, which emits openPMD 2.0 with a STRING
    extension openpmd-viewer rejects): replicate WarpX's byte-layout (openPMD 1.1.0, integer
    ED-PIC ext). `species` is the openPMD group key (PLURAL); `species`/`charge`/`mass` default to
    electrons, the converter passes positrons (`charge=+Q_E`). Records position [m], momentum
    [kg*m/s], weighting [count], charge, mass, id. Returns the written file path.
    """
    import openpmd_api as io
    os.makedirs(out_dir, exist_ok=True)
    n = pg.n_particle
    if n == 0:
        raise ValueError("write_openpmd_particles: ParticleGroup is empty")

    # px[eV/c]/MC2_EV = gamma*beta, then *m_e*c -> kg*m/s so viewer's mom/(mass*c) recovers u.
    mom_scale = M_E * C_LIGHT
    px = np.asarray(pg.px, dtype=np.float64) / MC2_EV * mom_scale
    py = np.asarray(pg.py, dtype=np.float64) / MC2_EV * mom_scale
    pz = np.asarray(pg.pz, dtype=np.float64) / MC2_EV * mom_scale
    x = np.asarray(pg.x, dtype=np.float64)
    y = np.asarray(pg.y, dtype=np.float64)
    z = np.asarray(pg.z, dtype=np.float64)
    w = np.asarray(pg.weight, dtype=np.float64) / Q_E    # PG.weight is CHARGE [C] -> count

    series = io.Series(os.path.join(out_dir, "openpmd_%06T.h5"), io.Access.create)
    series.set_openPMD("1.1.0")                           # viewer rejects 2.0 / STRING-ext
    series.set_openPMD_extension(1)                       # ED-PIC (integer)
    series.set_software("sim.linac5-8")
    series.set_particles_path("particles")

    it = series.iterations[int(iteration)]
    it.set_time(float(time)).set_dt(1.0).set_time_unit_SI(1.0)
    sp = it.particles[species]                           # PLURAL group key -- cross-stage contract

    dset_f = io.Dataset(np.dtype("float64"), [n])
    dset_i = io.Dataset(np.dtype("int64"), [n])

    def _tag(record, weighting_power, macro_weighted=0):
        # ED-PIC per-record attrs openpmd-viewer REQUIRES (else "Error: macroWeighted").
        record.set_attribute("macroWeighted", np.int32(macro_weighted))
        record.set_attribute("weightingPower", float(weighting_power))

    pos = sp["position"]
    pos.set_unit_dimension({io.Unit_Dimension.L: 1})
    _tag(pos, weighting_power=0)
    for comp, arr in (("x", x), ("y", y), ("z", z)):
        pos[comp].reset_dataset(dset_f)
        pos[comp].store_chunk(np.ascontiguousarray(arr))
        pos[comp].unit_SI = 1.0

    off = sp["positionOffset"]
    off.set_unit_dimension({io.Unit_Dimension.L: 1})
    _tag(off, weighting_power=0)
    zeros = np.zeros(n, dtype=np.float64)
    for comp in ("x", "y", "z"):
        off[comp].reset_dataset(dset_f)
        off[comp].store_chunk(np.ascontiguousarray(zeros))
        off[comp].unit_SI = 1.0

    mom = sp["momentum"]
    mom.set_unit_dimension({io.Unit_Dimension.M: 1, io.Unit_Dimension.L: 1,
                            io.Unit_Dimension.T: -1})
    _tag(mom, weighting_power=1)
    for comp, arr in (("x", px), ("y", py), ("z", pz)):
        mom[comp].reset_dataset(dset_f)
        mom[comp].store_chunk(np.ascontiguousarray(arr))
        mom[comp].unit_SI = 1.0

    wt = sp["weighting"][io.Mesh_Record_Component.SCALAR]
    sp["weighting"].set_unit_dimension({})
    _tag(sp["weighting"], weighting_power=1, macro_weighted=1)
    wt.reset_dataset(dset_f)
    wt.store_chunk(np.ascontiguousarray(w))
    wt.unit_SI = 1.0

    ch = sp["charge"][io.Mesh_Record_Component.SCALAR]
    sp["charge"].set_unit_dimension({io.Unit_Dimension.T: 1, io.Unit_Dimension.I: 1})
    _tag(sp["charge"], weighting_power=1)
    ch.reset_dataset(dset_f)
    ch.store_chunk(np.ascontiguousarray(np.full(n, charge, dtype=np.float64)))
    ch.unit_SI = 1.0

    ms = sp["mass"][io.Mesh_Record_Component.SCALAR]
    sp["mass"].set_unit_dimension({io.Unit_Dimension.M: 1})
    _tag(sp["mass"], weighting_power=1)
    ms.reset_dataset(dset_f)
    ms.store_chunk(np.ascontiguousarray(np.full(n, mass, dtype=np.float64)))
    ms.unit_SI = 1.0

    ids = np.asarray(pg["id"], dtype=np.int64) if "id" in pg else np.arange(1, n + 1, dtype=np.int64)
    sp["id"].set_unit_dimension({})
    _tag(sp["id"], weighting_power=0)
    idc = sp["id"][io.Mesh_Record_Component.SCALAR]
    idc.reset_dataset(dset_i)
    idc.store_chunk(np.ascontiguousarray(ids))
    idc.unit_SI = 1.0

    series.flush()
    del series                                           # close (flush on destruct)
    return os.path.join(out_dir, f"openpmd_{int(iteration):06d}.h5")
