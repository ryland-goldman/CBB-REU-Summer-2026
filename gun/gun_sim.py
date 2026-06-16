"""
CESR gun in WarpX (RZ): accelerate the cathode-emitted electrons through the
Poisson–Superfish gun field, with self-consistent space charge.

This is the second stage of the Cornell Linac chain modelled in WarpX. Stage 1
(`cathode/`) is the thermionic cathode operating at the Child–Langmuir
limit; here we take its emitted electrons and track them through the gun's
electrostatic accelerating field — the `CESR_gun.gdf` map scaled to 150 kV by
`build_gun_field.py` and applied as an external (electrode) field on the
particles, while WarpX's electromagnetostatic solver supplies the beam self-field
(electrostatic E plus the self magnetic B from the beam current, so the
relativistic magnetic pinch is included).

Geometry is RZ (cylindrical), matching the gun field map's native symmetry.

Run with (from the repo root, with `conda activate CBB`):
    python -c "import gun; gun.run()"               # build field map + sim + plots
    python -c "import gun; gun.run(plots=False)"    # build + sim only
    python -c "import gun; gun.plot()"              # plots only
Direct script invocation (`python gun/gun_sim.py`) does NOT work — this module
imports `pipeline._runner`, which is only on sys.path when launched from the
repo root (either via the facade above or `python -m gun.gun_sim`).

Beam source — see README: the cathode run is a continuous (DC) emitter, so the
weights in its last particle snapshot encode the steady-state population in
transit through the diode (~82 nC), not a bunch charge. We import the emitted
**phase-space distribution** (positions + momenta), remap the 2D (x, z) slab
into RZ by treating |x| as the radius r and smearing the particles uniformly in
azimuth — importance-resampling by r so the revolution supplies its 2πr Jacobian
(a uniform-in-x slab → a uniform-density disc, not a spurious 1/r on-axis cusp) —
and renormalize the total weight to a physical gun bunch charge
`BUNCH_CHARGE` (the CESR gun is pulse-grid gated). The full 82 nC injected as
one instantaneous bunch is unphysical — its radial space-charge field (~50 MV/m)
dwarfs the gun field and blows the beam apart before it accelerates.

Beam representation (`BEAM_RELEASE`, default "timed"): the CESR gun is gated by a 2 ns
grid pulse (cathode_master.in `twidth=2`), four gun-transit-times long, so the physical
beam is a long, low-density quasi-DC stream — the original GPT deck emits it that way via
`settdist(...,"t",...)`. Releasing the whole bunch at one instant ("snapshot") instead
over-concentrates the charge and over-states the space-charge force (the WarpX–GPT
150 kV-gun cross-code benchmark in CornellMisc/.../bench/writeup put the *controlled* beam-
representation term at ~28 % on emittance; here the over-dense snapshot ALSO blows the halo to
the wall — ≈81 % vs ≈100 % transmission — so the two modes' raw εn,x are computed on different
particle sets and their difference is confounded, see README). We release the imported
macroparticles over `PULSE_WIDTH` (uniform/flat-top — the real pulse has 30 V/ns edge ramps,
not modelled) via a per-step injection callback (the benchmark's warpx_tr.py technique), and
reconstruct the full exit beam for the injector by id-tracking it through the field-free pad
past the field map (`build_exit_handoff`). "snapshot" is kept for speed and back-compat. See
README → *Beam source* for the handoff and the injector-retuning caveat.
"""

import os
import shutil

import numpy as np
import pywarpx
from pywarpx import picmi, callbacks, particle_containers
from openpmd_viewer import OpenPMDTimeSeries

from pipeline._runner import run_step

c = picmi.constants.c
m_e = picmi.constants.m_e
q_e = picmi.constants.q_e
ep0 = picmi.constants.ep0

# ── Gun / field-map parameters (must match build_gun_field.py) ────────────────
GUN_FIELD = "gun/gun_field/gun_E.h5"
GUN_VOLTAGE = 150.0e3        # [V]
RMAX = 0.015                 # field-map R extent [m]
ZMAX_FIELD = 0.051765        # field-map Z extent [m] — the gun field ends here (Ez→0 at the
                             # map edge, verified) and is the physical exit plane.
# The RZ DOMAIN extends a field-free drift pad past the field map. WarpX zero-fills the
# applied field beyond the map (verified), so z ∈ (ZMAX_FIELD, ZMAX] is a field-free drift.
# This pad is REQUIRED by the timed-mode exit-beam reconstruction: build_exit_handoff samples
# each particle AFTER it clears the field, in field-free space, where the drift-to-a-common-
# time is emittance-preserving. Without the pad the domain ended exactly at the field map, so
# the reconstruction drifted still-in-field particles as if field-free and inflated εn,x ~8×
# (physics-review CRITICAL). ZPAD ≥ one inter-dump z-step (≈ v_exit·run_time/N_DIAGS ≈ 13 mm)
# so every exiting particle is caught in ≥1 field-free dump before it leaves the domain.
ZPAD = 0.020                 # field-free drift pad past the field map [m]
ZMAX = ZMAX_FIELD + ZPAD     # full RZ domain z-extent [m]

CATHODE_DIAG = "cathode/diags/particles"
BUNCH_CHARGE = 1.0e-9        # renormalized gun bunch charge [C] = 1 nC, matching the
                             # original LinacSim gpt_master.in total_charge = -1e-9;
                             # raw cathode snapshot is ~82 nC
RNG_SEED = 0

# ── Beam representation (snapshot vs time-release) ────────────────────────────
# The single largest unrealism the WarpX–GPT 150 kV-gun cross-code benchmark
# (CornellMisc/.../bench/writeup) identified: injecting the whole bunch at ONE instant
# ("snapshot") over-concentrates the charge and over-states the space-charge force vs
# releasing it over the real emission window ("time-release"), a ~28 % effect on εn,x in
# a comparable gun. The CESR gun is gated by a 2 ns grid pulse (cathode_master.in
# `twidth=2`), four gun-transit-times long, so the physical beam is a long, low-density,
# quasi-DC stream — exactly what the original GPT deck emits via `settdist(...,"t",...)`.
# We reproduce that by injecting the imported macroparticles over PULSE_WIDTH with a
# per-step `installbeforestep` callback (the same technique the benchmark's warpx_tr.py
# uses), instead of all at t=0.
BEAM_RELEASE = "timed"       # "timed" → release over PULSE_WIDTH (realistic, default);
                             # "snapshot" → all charge at t=0 (over-states space charge,
                             # kept for speed / back-compat with the old chain handoff).
PULSE_WIDTH = 2.0e-9         # grid-pulse emission window [s] (cathode_master.in twidth=2 ns).
                             # Flat-top model: emission times are drawn uniformly over
                             # [0, PULSE_WIDTH] (the real pulse has 30 V/ns edge ramps; a
                             # flat top is a documented first approximation — the dominant
                             # correction is the line-density drop, not the edge shape).
HANDOFF_DIR = "gun/diags/handoff"   # timed mode reconstructs the full released exit beam
                             # here (id-tracked exit-plane crossing across the volumetric
                             # dumps) for the injector; snapshot mode leaves it absent.

# ── Grid (RZ, single azimuthal mode — the gun field is m = 0) ─────────────────
# 128 (r) × 712 (z): the WarpX–GPT 150 kV-gun cross-code benchmark
# (CornellMisc/.../bench/writeup) found εn,x falls ~4 % and only converges once the
# grid resolves the near-cathode dynamics (their RZ study converged by nz≈720 over
# 55 mm); the old nz=384 over 51.77 mm (dz≈135 µm) sat on the unconverged side. nz=712
# over the padded 71.77 mm domain keeps dz≈101 µm (near-isotropic dz/dr≈0.86) across
# BOTH the field region and the field-free pad (uniform grid).
nr, nz = 128, 712            # divisible by the blocking factor (8)

# ── Diagnostics output directory ──────────────────────────────────────────────
DIAG_DIR = "gun/diags"

# ── Performance knobs (tunable via gun.config(...); see pipeline/run_pipeline.py) ─
# Defaults reproduce the original run exactly; lower them to trade accuracy for speed.
# Runtime ≈ nz² (per-step cost ∝ cells, and dz=ZMAX/nz ⇒ fewer derived steps as nz drops).
REQUIRED_PRECISION = 1e-5            # MLMG Poisson solve relative tolerance
SPACE_CHARGE = True                  # beam self-field (space charge) on/off. False →
                                     # warpx_do_not_deposit: the beam deposits no charge,
                                     # so only the applied gun field acts (no self-repulsion).
MAX_ITERS = None                     # MLMG iteration cap (None → PICMI default)
CFL = 0.4                            # dt = CFL · dz / v_exit
TRANSIT_MARGIN = 1.15                # run length = TRANSIT_MARGIN × gun-transit time
AVG_SPEED_FRAC = 0.6                 # bunch average speed as a fraction of v_exit
MAX_STEPS = 0                        # 0 → auto-derive from CFL/margins; >0 → fixed
N_DIAGS = 40                         # number of openPMD dumps over the run (≥20 keeps
                                     # space_charge.png's near-launch field snapshot)
MAX_PART = 0                         # 0/None → keep all cathode particles; >0 → cap


def load_cathode_bunch():
    """Import the last cathode snapshot and remap the (x, z) slab into RZ.

    Returns dict of x, y, z, ux, uy, uz, w arrays for ParticleListDistribution.
    """
    ts = OpenPMDTimeSeries(CATHODE_DIAG)
    if len(ts.iterations) == 0:
        raise RuntimeError(
            f"{CATHODE_DIAG} has no iterations — did the cathode stage run and "
            f"produce particles?")
    it = ts.iterations[-1]
    x, z, ux, uy, uz, w = ts.get_particle(
        ["x", "z", "ux", "uy", "uz", "w"], species="electrons", iteration=it,
    )

    rng = np.random.default_rng(RNG_SEED)
    r = np.abs(x)
    keep = r < RMAX
    if not keep.any():
        raise RuntimeError(
            f"no cathode particles with r < RMAX={RMAX} m; check RMAX or the "
            f"upstream cathode output")
    # Keep `xk` (the masked signed x) alongside the kept arrays so the radial-momentum
    # sign below survives the optional downsample (x[keep] would re-mask the full set).
    xk = x[keep]
    r, z, ux, uy, uz, w = (a[keep] for a in (r, z, ux, uy, uz, w))

    # Optionally downsample (reweighted to preserve total charge) to cap the cost.
    if MAX_PART and r.size > MAX_PART:
        sel = rng.choice(r.size, MAX_PART, replace=False)
        scale_w = r.size / MAX_PART
        xk, r, z, ux, uy, uz, w = (a[sel] for a in (xk, r, z, ux, uy, uz, w))
        w = w * scale_w

    # slab(x) → RZ disc: supply the 2πr revolution Jacobian that the naive r=|x|
    # map omits. A 2D Cartesian slab uniform in x has a flat dN/dr; revolving it
    # with r=|x| and unchanged weight yields areal density n(r) ∝ 1/r — a spurious
    # on-axis charge cusp that gives a radially-flat (nonlinear) self-field and
    # corrupts the σ_r, φ-well, and emittance this stage is meant to deliver.
    # Importance-resample (with replacement) with probability ∝ r·w (charge-correct;
    # ≡ ∝ r for the cathode's uniform weights), so dN/dr → r·dN/dr and the areal
    # density matches the cathode's true radial profile (a flat-top emitting strip
    # → a uniform-density disc). Drawing from the actual particles keeps weights
    # uniform (no weight-variance, so downstream RMS/emittance stay
    # unweighted-valid) and preserves the cathode-edge position–momentum correlations.
    if r.max() > 0.0:
        rw = r * w
        sel = rng.choice(r.size, r.size, replace=True, p=rw / rw.sum())
        xk, r, z, ux, uy, uz, w = (a[sel] for a in (xk, r, z, ux, uy, uz, w))

    theta = rng.uniform(0.0, 2.0 * np.pi, size=r.size)
    ct, st = np.cos(theta), np.sin(theta)

    # slab x -> radius; transverse momentum: radial = ux·sign(x), azimuthal = uy
    ur = ux * np.sign(np.where(xk == 0.0, 1.0, xk))
    ut = uy
    xpos = r * ct
    ypos = r * st
    uxn = ur * ct - ut * st
    uyn = ur * st + ut * ct
    zpos = np.clip(z, 0.0, ZMAX)

    # Renormalize weights so the imported distribution carries BUNCH_CHARGE.
    w = w * (BUNCH_CHARGE / (w.sum() * q_e))

    # Emission time per macroparticle. In "timed" mode the bunch is released over the
    # PULSE_WIDTH grid-pulse window (uniform → flat-top current); the macroparticles are
    # sorted by t so the per-step injection callback can walk them in one pass. In
    # "snapshot" mode every particle is emitted at t=0 (all charge present at once).
    if BEAM_RELEASE == "timed":
        t = rng.uniform(0.0, PULSE_WIDTH, size=r.size)
        order = np.argsort(t)
        t = t[order]
        xpos, ypos, zpos, uxn, uyn, uz, w = (
            a[order] for a in (xpos, ypos, zpos, uxn, uyn, uz, w))
    else:
        t = np.zeros(r.size)

    print(f"Imported {r.size} macroparticles from cathode (iter {it}); "
          f"renormalized to {BUNCH_CHARGE*1e9:.3f} nC, r ≤ {r.max()*1e3:.2f} mm; "
          f"release={BEAM_RELEASE}"
          + (f" over {PULSE_WIDTH*1e9:.1f} ns" if BEAM_RELEASE == "timed" else ""),
          flush=True)
    # openPMD ux/uy/uz are the dimensionless normalized momenta γβ; PICMI's
    # ParticleListDistribution wants proper velocity u = γβc in m/s, so ×c.
    # (Without this the beam is injected essentially at rest and the cathode's
    # thermal transverse momentum — hence its emittance — is lost; the energy
    # gain is insensitive because the cathode KE ≪ 150 keV gun voltage.)
    return dict(x=xpos, y=ypos, z=zpos, ux=uxn * c, uy=uyn * c, uz=uz * c, w=w, t=t)


def build_exit_handoff():
    """Reconstruct the full time-released exit beam and write it for the injector.

    Time-release makes the gun beam a ~2 ns quasi-DC stream whose ballistic z-extent
    (~v_exit·PULSE_WIDTH ≈ 0.4 m) is many times the gun domain, so no single volumetric
    snapshot can hold the whole released beam. We reconstruct it from the volumetric dumps
    by particle **id** (the `pipeline/collimator.py` idiom), sampling each particle in the
    FIELD-FREE pad past the field map so the reconstruction is emittance-correct:

      1. For each id, find its FIRST appearance with z ≥ ZMAX_FIELD — i.e. just after it
         clears the gun field into the field-free pad. Its (x, u) there IS the exit-plane
         phase space (the field is ≈0 beyond ZMAX_FIELD). Sampling in the pad — NOT at the
         particle's last in-field dump — is the fix for the εn,x inflation: a field-free
         ballistic drift preserves εn,x, but drifting a still-in-field particle as if
         field-free manufactures a spurious x–u correlation (physics-review CRITICAL).
      2. An id that NEVER reaches z ≥ ZMAX_FIELD did not exit: if it is still present in the
         final dump it is the un-flushed tail (run ended first); otherwise it was absorbed
         before the pad, i.e. scraped at the r = RMAX wall (the only other boundary) → a
         radial loss. Both are dropped and counted (this correctly identifies losses — the
         old r-at-last-dump test was a no-op because absorbed particles vanish at r<RMAX).
      3. Drift the kept (field-free) samples to a common reference time t_ref = max sample
         time — all motion is field-free so εn,x is preserved. early-emitted → exited
         early → drifted furthest → bunch HEAD (largest z); last-emitted → z≈ZMAX_FIELD →
         TAIL, matching the injector's `z − z.min() + Z_INJECT` (tail at the entrance).

    Writes an openPMD particle dump to HANDOFF_DIR via `pipeline.impact_io`. The injector's
    `load_gun_bunch` reads HANDOFF_DIR when present (else the volumetric snapshot).
    """
    from pmd_beamphysics import ParticleGroup
    from pipeline.impact_io import write_openpmd_particles

    MC2_EV = 510998.95069
    pdir = os.path.join(DIAG_DIR, "particles")
    ts = OpenPMDTimeSeries(pdir)
    if len(ts.iterations) == 0:
        print("  handoff: no volumetric dumps to reconstruct from — skipped", flush=True)
        return

    # Stack every (id, dump) row. Vectorized so the ~Npart×Ndumps rows never hit a loop.
    cols = {k: [] for k in ("id", "t", "x", "y", "z", "ux", "uy", "uz", "w")}
    final_it = ts.iterations[-1]
    for it, t_dump in zip(ts.iterations, ts.t):
        idp, x, y, z, ux, uy, uz, w = ts.get_particle(
            ["id", "x", "y", "z", "ux", "uy", "uz", "w"],
            species="electrons", iteration=it)
        if len(idp) == 0:
            continue
        cols["id"].append(idp.astype(np.int64)); cols["t"].append(np.full(len(idp), t_dump))
        cols["x"].append(x); cols["y"].append(y); cols["z"].append(z)
        cols["ux"].append(ux); cols["uy"].append(uy); cols["uz"].append(uz); cols["w"].append(w)
    if not cols["id"]:
        print("  handoff: no particles tracked — skipped", flush=True)
        return
    cat = {k: np.concatenate(v) for k, v in cols.items()}
    n_ids = np.unique(cat["id"]).size

    # First appearance in the FIELD-FREE pad (z ≥ ZMAX_FIELD) per id: restrict to pad rows,
    # sort by (id, t), take the earliest-t row of each id-group.
    pad = np.where(cat["z"] >= ZMAX_FIELD)[0]
    if pad.size == 0:
        print("  handoff: no particle reached the field-free pad (z ≥ "
              f"{ZMAX_FIELD*1e3:.1f} mm) — run longer? skipped", flush=True)
        return
    order = np.lexsort((cat["t"][pad], cat["id"][pad]))
    sid = cat["id"][pad][order]
    first = np.empty(sid.size, dtype=bool)
    first[0] = True
    first[1:] = sid[1:] != sid[:-1]                # True at each id-group start (earliest t)
    pick = pad[order[first]]
    t_s = cat["t"][pick]
    x, y, z = cat["x"][pick], cat["y"][pick], cat["z"][pick]
    ux, uy, uz, w = cat["ux"][pick], cat["uy"][pick], cat["uz"][pick], cat["w"][pick]
    n_exit = pick.size

    # Classify the ids that never reached the pad: still in the final dump → un-flushed tail;
    # else absorbed before the pad → radial (r=RMAX) loss.
    exited_ids = set(cat["id"][pick].tolist())
    final_ids = set(cat["id"][cat["t"] == ts.t[-1]].tolist()) if len(ts.t) else set()
    not_exited = [i for i in np.unique(cat["id"]).tolist() if i not in exited_ids]
    n_unflushed = sum(1 for i in not_exited if i in final_ids)
    n_radial = len(not_exited) - n_unflushed

    # Ballistic drift of the field-free samples to a common reference time t_ref = max(t_s).
    gamma = np.sqrt(1.0 + ux**2 + uy**2 + uz**2)
    vx, vy, vz = (c * ux / gamma, c * uy / gamma, c * uz / gamma)
    dtau = t_s.max() - t_s
    xh, yh, zh = x + vx * dtau, y + vy * dtau, z + vz * dtau

    pg = ParticleGroup(data=dict(
        x=xh, y=yh, z=zh,
        px=ux * MC2_EV, py=uy * MC2_EV, pz=uz * MC2_EV,   # γβ·m_ec² /c = γβ·(mc² in eV) [eV/c]
        t=np.zeros(xh.size),
        weight=w * q_e,                                   # macro count → charge [C]
        status=np.ones(xh.size, dtype=np.int64),
        species="electron",
    ))
    if os.path.isdir(HANDOFF_DIR):
        shutil.rmtree(HANDOFF_DIR)
    write_openpmd_particles(pg, HANDOFF_DIR, iteration=0, time=0.0)
    print(f"  handoff: {n_exit}/{n_ids} macroparticles cleared the field "
          f"({pg.charge*1e9:.3f} nC, transmission {100*n_exit/n_ids:.0f}%; "
          f"{n_radial} radial loss, {n_unflushed} un-flushed), "
          f"z-extent {(zh.max()-zh.min())*1e3:.0f} mm → {HANDOFF_DIR}", flush=True)


def main():
    grid = picmi.CylindricalGrid(
        number_of_cells=[nr, nz],
        n_azimuthal_modes=1,
        lower_bound=[0.0, 0.0],
        upper_bound=[RMAX, ZMAX],
        # r=0 must be "none" (axis); the electrode field is applied externally, so
        # the self-field solve just needs grounded z plates (dirichlet). The outer
        # radial wall is also dirichlet (grounded pipe at r=RMAX=15 mm, well outside
        # the r≲8 mm beam): the electromagnetostatic solver also does a vector-Poisson
        # solve for A, and the dominant A_z component (driven by the beam's j_z) would
        # have an all-Neumann, singular operator — and the MLMG bottom solve diverges —
        # if the outer wall were neumann; the dirichlet wall makes A_z well-posed.
        lower_boundary_conditions=["none", "dirichlet"],
        upper_boundary_conditions=["dirichlet", "dirichlet"],
        lower_boundary_conditions_particles=["none", "absorbing"],
        upper_boundary_conditions_particles=["absorbing", "absorbing"],
        warpx_blocking_factor=8,
    )

    # Electromagnetostatic solver for the beam self-field. In addition to the
    # electrostatic Poisson solve (∇²φ = -ρ/ε₀, E = -∇φ), warpx_magnetostatic=True
    # also solves the Coulomb-gauge vector potential from the beam current
    # (∇²A = -μ₀ j, B = ∇×A), so the self magnetic field is included. This supplies
    # the relativistic magnetic-pinch term qβ×B that partially cancels the radial
    # space-charge repulsion — the net transverse self-force is qE_r/γ² rather than the
    # pure-electrostatic qE_r (the ≈γ²=1.66× over-repulsion at the 146 keV exit the
    # plain labframe solver incurs). The magnetostatic MLMG solve is given the same
    # REQUIRED_PRECISION / MAX_ITERS knobs as the electrostatic solve via the explicit
    # warpx_magnetostatic_required_precision / warpx_magnetostatic_max_iters params.
    solver_kw = dict(grid=grid, method="Multigrid",
                     required_precision=REQUIRED_PRECISION,
                     warpx_magnetostatic=True,
                     warpx_magnetostatic_required_precision=REQUIRED_PRECISION,
                     warpx_self_fields_verbosity=0,           # silence ES MLMG per-iteration chatter
                     warpx_magnetostatic_verbosity=0)         # and the magnetostatic solve
    if MAX_ITERS:                                     # omit when None → PICMI default
        solver_kw["maximum_iterations"] = MAX_ITERS
        solver_kw["warpx_magnetostatic_max_iters"] = MAX_ITERS
    solver = picmi.ElectrostaticSolver(**solver_kw)

    # ── Applied gun field: the scaled CESR_gun.gdf map, read from file ────────
    # Applied directly to particles every step (the electrode field), on top of the
    # self-consistent space-charge field from the Poisson solve. PICMI has no class
    # for a tabulated particle-applied field, so set the raw WarpX inputs.
    pywarpx.particles.E_ext_particle_init_style = "read_from_file"
    pywarpx.particles.read_fields_from_path = GUN_FIELD
    pywarpx.particles.B_ext_particle_init_style = "none"

    bunch = load_cathode_bunch()
    timed = (BEAM_RELEASE == "timed")
    # snapshot: seed the species with the whole bunch (all at t=0). timed: seed with only
    # the earliest-emitted macroparticle (PICMI requires a non-empty initial distribution),
    # then inject the rest over the pulse via the per-step callback below. `sl` selects the
    # seed slice; the cathode arrays are already t-sorted in timed mode.
    sl = slice(0, 1) if timed else slice(None)
    electrons = picmi.Species(
        particle_type="electron",
        name="electrons",
        initial_distribution=picmi.ParticleListDistribution(
            x=bunch["x"][sl], y=bunch["y"][sl], z=bunch["z"][sl],
            ux=bunch["ux"][sl], uy=bunch["uy"][sl], uz=bunch["uz"][sl],
            weight=bunch["w"][sl],
        ),
        warpx_do_not_deposit=not SPACE_CHARGE,   # SPACE_CHARGE=False → no beam self-field
    )

    # ── Time step / duration ──────────────────────────────────────────────────
    # Exit kinetic energy ≈ 150 keV -> γ ≈ 1.29, β ≈ 0.63, v_exit ≈ 1.9e8 m/s.
    gamma = 1.0 + q_e * GUN_VOLTAGE / (m_e * c**2)
    v_exit = c * np.sqrt(1.0 - 1.0 / gamma**2)
    dz = ZMAX / nz
    dt = CFL * dz / v_exit
    # Steps for the bunch to just cross the full gun (average speed ~AVG_SPEED_FRAC·v_exit).
    # We stop as the beam reaches the exit: running longer empties the domain, and
    # the Multigrid self-field solve aborts when there is essentially no charge left.
    # Run length is sized on the time to clear the FIELD region (ZMAX_FIELD), NOT the full
    # padded domain: the run must STOP while the beam is still in the field-free pad, before
    # it drains out the padded ZMAX — an empty domain aborts the MLMG self-field solve
    # (`MLMG failed`). At PULSE_WIDTH + TRANSIT_MARGIN·transit_field the last-released particle
    # has cleared the field into the pad (caught in ≥1 field-free dump for the handoff) while
    # the bulk is still transiting, so the solve is never charge-starved. MAX_STEPS (module
    # constant, 0 = auto) overrides the derived value when set.
    transit_field = ZMAX_FIELD / (AVG_SPEED_FRAC * v_exit)
    run_time = (PULSE_WIDTH if timed else 0.0) + TRANSIT_MARGIN * transit_field
    max_steps = MAX_STEPS or int(run_time / dt)

    print(f"Gun: {GUN_VOLTAGE/1e3:.0f} kV  ->  γ={gamma:.3f}, β={v_exit/c:.3f}, "
          f"v_exit={v_exit:.2e} m/s", flush=True)
    print(f"dt = {dt:.3e} s, max_steps = {max_steps}"
          + (f" (release over {PULSE_WIDTH*1e9:.1f} ns + transit)" if timed else ""),
          flush=True)

    # ── Diagnostics (openPMD, HDF5) ───────────────────────────────────────────
    # Fresh diags: WarpX appends one openPMD file per dump, so re-running with a
    # different grid/step count would otherwise mix old and new iterations (whose
    # diag steps interleave) into one series — the plots then show a fan of
    # overlapping curves. diags are git-ignored and regenerated, so clearing is
    # safe. (Mirrors injector_sim.py / linac_sec1_sim.py.)
    if os.path.isdir(DIAG_DIR):
        shutil.rmtree(DIAG_DIR)

    period = max(1, max_steps // N_DIAGS)
    field_diag = picmi.FieldDiagnostic(
        name="fields",
        grid=grid,
        period=period,
        data_list=["phi", "rho", "E"],
        write_dir=DIAG_DIR,
        warpx_format="openpmd",
        warpx_openpmd_backend="h5",
    )
    part_diag = picmi.ParticleDiagnostic(
        name="particles",
        period=period,
        species=[electrons],
        data_list=["position", "momentum", "weighting"],
        write_dir=DIAG_DIR,
        warpx_format="openpmd",
        warpx_openpmd_backend="h5",
    )

    sim = picmi.Simulation(
        solver=solver,
        max_steps=max_steps,
        time_step_size=dt,
        verbose=0,                     # silence per-step "STEP N starts" — the tqdm bar is the progress display
        particle_shape="linear",
    )
    # ParticleListDistribution supplies the macroparticles explicitly, so this layout
    # (and its n_macroparticles_per_cell) is inert — the count is the imported list size.
    sim.add_species(
        electrons,
        layout=picmi.PseudoRandomLayout(n_macroparticles_per_cell=1, grid=grid),
    )
    sim.add_diagnostic(field_diag)
    sim.add_diagnostic(part_diag)

    # ── Time-release injection ────────────────────────────────────────────────
    # Build up the bunch over PULSE_WIDTH: each step, inject the macroparticles whose
    # emission time falls in the step window. Particle 0 is already seeded, so injection
    # starts at index 1. add_particles writes into the live container (the broken RZ
    # accessor is the READ path — `Component x does not exist` — not this write path,
    # verified by spike). ux/uy/uz are proper velocity γβc [m/s], matching the seed's
    # ParticleListDistribution units.
    if timed:
        state = {"next": 1, "step": 0}
        pc = [None]
        bt = bunch["t"]

        def _inject():
            if pc[0] is None:
                pc[0] = particle_containers.ParticleContainerWrapper("electrons")
            t_hi = (state["step"] + 1) * dt
            j = state["next"]
            k = j
            while k < bt.size and bt[k] < t_hi:
                k += 1
            if k > j:
                pc[0].add_particles(
                    x=bunch["x"][j:k], y=bunch["y"][j:k], z=bunch["z"][j:k],
                    ux=bunch["ux"][j:k], uy=bunch["uy"][j:k], uz=bunch["uz"][j:k],
                    w=bunch["w"][j:k])
                state["next"] = k
            state["step"] += 1

        callbacks.installbeforestep(_inject)

    print(f"\nRunning {max_steps} steps (diag every {period}) …")
    run_step(sim, max_steps, desc="gun")
    if timed:
        print(f"Injected {state['next']}/{bunch['t'].size} macroparticles over the pulse",
              flush=True)
        build_exit_handoff()
    print(f"\nDone. openPMD output → {DIAG_DIR}/{{fields,particles}}/")


if __name__ == "__main__":
    main()
