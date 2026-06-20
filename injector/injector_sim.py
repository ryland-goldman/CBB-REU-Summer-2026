"""
CESR injector in WarpX (RZ): the full LinacSim injector subsection in one
self-consistent space-charge run — two 214 MHz prebuncher cavities (Preb 2 reversed)
and three solenoid lenses (Lens 0A / Sol 0 / Lens 0E) — reading the gun exit beam and
handing a focused, velocity-bunched beam to linac_sec1 at z ≈ 2.03 m. Operating-point
constants below are config()-overridable.

See injector/README.md for physics, parameters, field maps, and gotchas.
"""

import os
import shutil
import time

import numpy as np
import pywarpx
import openpmd_api as io
from pywarpx import picmi
from openpmd_viewer import OpenPMDTimeSeries

from pipeline._runner import run_step
from pipeline.constants import MC2_EV


def _retry_io(fn, *args, tries=6, base=0.25, **kwargs):
    """Call an openPMD read, retrying a transient HDF5 "Inaccessible" open error.

    Backstop only — the production "Inaccessible" failure is fd exhaustion, fixed by
    raising RLIMIT_NOFILE (see README -> openPMD fd-leak gotcha); this retry does not
    rescue that. Re-raise after the last try.
    """
    for i in range(tries):
        try:
            return fn(*args, **kwargs)
        except io.Error:
            if i == tries - 1:
                raise
            time.sleep(base * 2 ** i)


# RF-drive constants live in the (pywarpx-free) build module as the single source of
# truth, so sim and plot_injector.py cannot drift apart.
from .build_injector_field import (
    Z_GAP_CENTER_1, Z_GAP_CENTER_2, V1J_KEV, F_RF, Q_L_1, Q_L_2,
    SOL_FILES, Z_HANDOFF,
)
from . import DEFAULT_OUTDIR

c = picmi.constants.c
m_e = picmi.constants.m_e
q_e = picmi.constants.q_e

# ── Field-map paths (must match build_injector_field.py) ───────────────────────
# Both cavities use the same FORWARD field; they differ only in lab-z placement
# (grid_global_offset baked per file, so Preb 2 needs its own file). The reversed
# install is a run-time phase, not a negated map. See README -> Field maps / Reversed install.
PREB1_FIELD = "injector/injector_field/preb1_EB.h5"   # forward field at Z_GAP_CENTER_1
PREB2_FIELD = "injector/injector_field/preb2_EB.h5"   # forward field at Z_GAP_CENTER_2

# Prefer the gun's reconstructed time-release exit beam (gun/diags/handoff) when present,
# else the legacy snapshot. See README -> Running (time-release handoff caveat).
GUN_DIAG = ("gun/diags/handoff" if os.path.isdir("gun/diags/handoff")
            else "gun/diags/particles")
Z_INJECT = 0.005                 # lab z where the bunch tail (smallest z) is placed [m]
MAX_PART = 50000                 # downsample the gun snapshot (reweighted) for speed
RNG_SEED = 0
CFL = 0.8                        # dt = CFL · Δz / v_beam

# ── Operating point — Prebuncher 1 (config()-overridable) ─────────────────────
PREB1_KW = 8.0                   # dissipated RF power [kW]
PREB1_Q = Q_L_1                  # loaded Q of prebuncher 1
PREB1_PHI_OFF = 0.0              # phase offset [deg] from the zc base (0 = centroid on zero-crossing)

# ── Operating point — Prebuncher 2 (reversed; config()-overridable) ───────────
PREB2_KW = 10.0                  # Prebuncher 2 design power [kW] (prebuncher2_input_power)
PREB2_Q = Q_L_2                  # loaded Q of prebuncher 2 (4300)
PREB2_PHI_OFF = 0.0              # phase offset [deg] from the zc base (0 = centroid on zero-crossing)
PREB2_REVERSED = True            # apply the reversed-install phase PREB2_REV_PHASE
# Reversed install ≡ +π in absolute drive phase (180° rotation flips Ez); in the zc+phi_off=0
# parametrization phi_off carries NO reversal info (unlike the old crest+GUI convention), so the
# +π is applied here explicitly. Both rev=0 and rev=π are energy-flat (centroid on a zero-crossing)
# — they differ only in the chirp SLOPE. Keep rev=π: it is the faithful reversed geometry AND lands
# the σ_z waist AT the 2.03 m handoff. rev=0 (a forward Preb-2) actually bunches HARDER — it
# over-compresses to an earlier waist (~1.64 m) that re-expands before the handoff. See README ->
# Reversed install.
PREB2_REV_PHASE = np.pi          # [rad] reversed install in the zc/centroid parametrization

# ── Solenoid lens currents [A] (config()-overridable; 0 disables a lens) ───────
# LinacSim GUI defaults (6/40/10 A); 1-A maps scale linearly with current.
# See README -> Solenoid lenses.
I_LENS0A = 6.0
I_SOL0 = 40.0
I_LENS0E = 10.0
# 0B/0C/0D default to the GUI 0 A (a 0-A lens is skipped below); config()-overridable.
I_LENS0B = 0.0
I_LENS0C = 0.0
I_LENS0D = 0.0

# ── Collimator (the faithful injector→linac iris/pipe) ────────────────────────
# 9.547 mm iris at z=1.922 m + pipe to 2.1 m. Applied POST-HOC (RZ build can't scrape r
# in-run) as a MULTI-PLANE id scrape, NOT a single 2.03 m cut: the beam CONVERGES across the
# 1.922→2.03 m tail, so a single cut would keep halo the real iris scrapes. The physical cut
# is the linac reader's scrape at injection; _report_collimated_handoff() below is a
# diagnostic print only. See README -> The 9.547 mm collimator.
COLLIM_R = 0.009547              # [m] iris/pipe radius (SLAC bore; gpt scatteriris)
COLLIM_Z = 1.922                # [m] iris start; the 9.547 mm pipe runs COLLIM_Z → ZMAX
COLLIMATE = True                # report the collimated handoff charge (set False to skip)

# Phase reference: "zc" (base=π/2) puts the bunch CENTROID on the RF zero-crossing, so the
# net mean-energy kick is zero and the cavity acts as a pure velocity buncher (the design
# goal). "crest" (base=π) is the legacy net-accelerating reference. See README -> RF drive.
PHASE = "zc"                     # default: zero-crossing (centroid-referenced) velocity bunching

# ── Performance knobs (tunable via injector.config(...); see run_pipeline.py) ──
# This stage dominates the pipeline. Do NOT coarsen NZ — convergence-bound, not cell-bound;
# speed it via CFL and MAX_ITERS/REQUIRED_PRECISION. See README -> Domain / grid.
REQUIRED_PRECISION = 1e-4        # MLMG relative tolerance (relaxed for the long-thin box)
MAX_ITERS = 500                  # MLMG iteration cap
SPACE_CHARGE = True              # beam self-field (space charge) on/off. False →
                                 # warpx_do_not_deposit: the beam deposits no charge, so
                                 # only the applied prebuncher/solenoid maps act.
MAX_STEPS = 0                    # 0 → auto-derive from transit; >0 → fixed
TRANSIT_MARGIN = 0.98            # stop just before the bunch centre reaches the exit (lands a dump ~2.03 m)
N_DIAGS = 60                     # number of openPMD dumps over the run

# ── Domain (RZ, single azimuthal mode — the cavity field is m = 0) ─────────────
# NZ=1664 gives dz=1.262 mm ⇒ 2.80:1 aspect (the ≈3:1 rule; below it the MLMG self-field
# solve diverges) and is ÷8 (blocking factor). Keep NR=80. See README -> Domain / grid.
RMAX = 0.036                     # covers the field-map bore (0–36.07 mm)
ZMAX = 2.10                      # full injector subsection (handoff at z≈2.03 m)
NR, NZ = 80, 1664                # dz=1.262 mm ⇒ 2.80:1 aspect; both ÷ blocking factor 8

OUTDIR = None                    # if None at main(), use DEFAULT_OUTDIR


def load_gun_bunch():
    """Import the gun's last snapshot (already RZ) and shift it to the entrance.

    Returns (dict for ParticleListDistribution, v_beam, mean KE [keV], z_centroid).
    ``z_centroid`` is the lab-z of the charge centroid after the entrance shift; the
    cavities are phased to put the CENTROID (not the tail) at the zero-crossing so the
    net mean-energy kick is zero. See README -> RF drive (centroid phase reference).
    """
    ts = OpenPMDTimeSeries(GUN_DIAG)
    if len(ts.iterations) == 0:
        raise RuntimeError(
            f"{GUN_DIAG} has no iterations — did the gun stage run and produce "
            f"particles?")
    it = ts.iterations[-1]
    x, y, z, ux, uy, uz, w = _retry_io(
        ts.get_particle,
        ["x", "y", "z", "ux", "uy", "uz", "w"], species="electrons", iteration=it,
    )
    # Downsample (reweighted to preserve total charge) to keep the run cheap.
    if z.size > MAX_PART:
        rng = np.random.default_rng(RNG_SEED)
        sel = rng.choice(z.size, MAX_PART, replace=False)
        scale_w = z.size / MAX_PART
        x, y, z, ux, uy, uz, w = (a[sel] for a in (x, y, z, ux, uy, uz, w))
        w = w * scale_w
    # Translate so the bunch *tail* (smallest z) sits at Z_INJECT (head is at larger z).
    z = z - z.min() + Z_INJECT

    # openPMD ux/uy/uz are the dimensionless normalized momenta γβ; PICMI's
    # ParticleListDistribution wants proper velocity u = γβc in m/s, so ×c.
    gb = np.sqrt(1.0 + ux**2 + uy**2 + uz**2)          # γ
    beta_z = uz / gb
    v_beam = float(np.average(beta_z, weights=w) * c)
    ke_mean = float(np.average(gb - 1.0, weights=w) * m_e * c**2 / q_e / 1e3)

    z_centroid = float(np.average(z, weights=w))
    print(f"Imported {z.size} macroparticles from gun (iter {it}); "
          f"z {z.min()*1e3:.1f}–{z.max()*1e3:.1f} mm, ⟨z⟩ {z_centroid*1e3:.1f} mm, "
          f"⟨KE⟩ {ke_mean:.1f} keV, v_beam {v_beam:.3e} m/s, "
          f"q {w.sum()*q_e*1e9:.3f} nC", flush=True)
    return (dict(x=x, y=y, z=z, ux=ux * c, uy=uy * c, uz=uz * c, w=w),
            v_beam, ke_mean, z_centroid)


def make_cavity(field_path, power, q_l, z_gap, v_at_gap, phi_off_deg, phase,
                omega, t_offset=0.0, rev_phase=0.0, z_ref=Z_INJECT):
    """Build one prebuncher cavity as a picmi.LoadAppliedField.

    Drives the 1-J map as a standing-wave TM mode: E ∝ scale·cos(ωt+φ),
    B ∝ scale·sin(ωt+φ). ``z_ref`` is the lab-z of the bunch CENTROID at the reference
    plane, so ``t_gap`` is the centroid's gap arrival and the zc base lands the centroid on
    the zero-crossing (net mean kick = 0). ``v_at_gap`` is the mean beam speed over
    ``z_ref → z_gap``; ``t_offset`` is the time already elapsed reaching ``z_ref`` (Preb 2
    uses a two-segment arrival to account for Preb-1's kick). ``rev_phase`` is added to φ.
    See README -> RF drive / Reversed install / Preb-2 timing caveat.

    Keep .10e precision on every term — ω·t truncation accumulates over the ~5 ns transit.
    """
    scale = float(np.sqrt(1e3 * q_l * power / (2.0 * np.pi * F_RF)))
    # Arrival time of the bunch centroid at this cavity's gap: time to z_ref + leg z_ref→z_gap.
    t_gap = t_offset + (z_gap - z_ref) / v_at_gap
    # Electron energy kick ΔW(t) ∝ -cos(ωt+φ) (on-axis Ez single-signed positive); zc base
    # = π/2 (centroid on zero-crossing -> net mean kick 0, pure buncher), crest base = π
    # (legacy net-accelerating). rev_phase carries the reversed install.
    base = np.pi / 2.0 if phase == "zc" else np.pi
    phi = -omega * t_gap + base + np.radians(phi_off_deg) + rev_phase
    e_time = f"{scale:.10e}*cos({omega:.10e}*t + ({phi:.10e}))"
    b_time = f"{scale:.10e}*sin({omega:.10e}*t + ({phi:.10e}))"
    fld = picmi.LoadAppliedField(
        read_fields_from_path=field_path, load_E=True, load_B=True,
        warpx_E_time_function=e_time, warpx_B_time_function=b_time,
    )
    return fld, scale, phi, t_gap


def _report_collimated_handoff(outdir):
    """Report the COLLIMATED handoff charge at the ~Z_HANDOFF plane for the sanity log.

    Applies the SAME multi-plane iris scrape the linac reader uses (pipeline.collimator):
    survivors / in-domain at the handoff. Diagnostic print only; the physical cut is the
    linac reader's at injection.
    """
    try:
        from openpmd_viewer import OpenPMDTimeSeries
        from pipeline.collimator import pipe_violator_ids, survivor_mask
        ts = OpenPMDTimeSeries(os.path.join(outdir, "particles"))
        recs = []
        for it in ts.iterations:
            z, w = _retry_io(ts.get_particle, ["z", "w"], species="electrons", iteration=it)
            if len(z) < 50:
                continue
            recs.append((it, float(np.average(z, weights=w))))
        if not recs:
            print("  collimated handoff: no populated snapshot near the plane", flush=True)
            return
        it_h, zm_h = min(recs, key=lambda t: abs(t[1] - Z_HANDOFF))
        idh, x, y, z, w = _retry_io(ts.get_particle, ["id", "x", "y", "z", "w"],
                                    species="electrons", iteration=it_h)
        q_dom = float(w.sum()) * q_e
        scan_iters = [it for it, zm in recs if (COLLIM_Z - 0.05) <= zm <= (Z_HANDOFF + 0.03)]
        violators = pipe_violator_ids(ts, scan_iters, COLLIM_R, COLLIM_Z)
        keep = survivor_mask(idh, violators)
        q_coll = float(w[keep].sum()) * q_e
        print(f"  COLLIMATED handoff (⟨z⟩={zm_h*1e3:.1f} mm, iris {COLLIM_R*1e3:.3f} mm, "
              f"multi-plane {len(scan_iters)} planes): {q_coll*1e9:.3f} nC survives the pipe / "
              f"{q_dom*1e9:.3f} nC in-domain = {100*q_coll/q_dom:.0f}% through the aperture",
              flush=True)
    except Exception as e:
        print(f"  collimated-handoff report unavailable: {e}", flush=True)


def main():
    outdir = OUTDIR or DEFAULT_OUTDIR

    # Fresh diags: WarpX appends per dump, so re-running would mix old/new iterations.
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)

    # Compute omega here (not at import) so a config(F_RF=...) override is honored.
    omega = 2.0 * np.pi * F_RF

    bunch, v_beam, ke_mean, z_centroid = load_gun_bunch()

    # ── Grid + electrostatic (self-field) solver ──────────────────────────────
    grid = picmi.CylindricalGrid(
        number_of_cells=[NR, NZ],
        n_azimuthal_modes=1,
        lower_bound=[0.0, 0.0],
        upper_bound=[RMAX, ZMAX],
        # r=0 is "none" (axis). Outer radial wall is dirichlet, NOT neumann: the EMS A_z
        # vector-Poisson (driven by j_z) is else all-Neumann SINGULAR and the MLMG bottom
        # solve diverges. At RMAX=36 mm self-field has decayed, so dynamics unaffected.
        lower_boundary_conditions=["none", "dirichlet"],
        upper_boundary_conditions=["dirichlet", "dirichlet"],
        lower_boundary_conditions_particles=["none", "absorbing"],
        upper_boundary_conditions_particles=["absorbing", "absorbing"],
        warpx_blocking_factor=8,
    )
    # Relativistic electromagnetostatic self-field solver (warpx_magnetostatic=True): adds the
    # Coulomb-gauge A solve so the qβ×B magnetic pinch gives net transverse self-force qE_r/γ²
    # (removes the ≈γ² lab-frame over-repulsion). See README -> Self-field solver.
    solver = picmi.ElectrostaticSolver(
        grid=grid, method="Multigrid", required_precision=REQUIRED_PRECISION,
        maximum_iterations=MAX_ITERS,
        warpx_magnetostatic=True,
        warpx_magnetostatic_required_precision=REQUIRED_PRECISION,
        warpx_magnetostatic_max_iters=MAX_ITERS,
        warpx_self_fields_verbosity=0,
        warpx_magnetostatic_verbosity=0,
    )

    # ── Applied fields ────────────────────────────────────────────────────────
    # B-only solenoid maps come BEFORE the RF cavities: picmi forces the global
    # E_ext_particle_init_style to "none" if the LAST-added field has load_E=False, so a
    # B-only map last would silently disable the accelerating E. 1-A maps scaled by current.
    applied = []
    for path, cur in [(SOL_FILES["LENS_0A"], I_LENS0A),
                      (SOL_FILES["LENS_0B"], I_LENS0B),
                      (SOL_FILES["LENS_0C"], I_LENS0C),
                      (SOL_FILES["LENS_0D"], I_LENS0D),
                      (SOL_FILES["SOL_0"],   I_SOL0),
                      (SOL_FILES["LENS_0E"], I_LENS0E)]:
        if cur != 0.0:
            applied.append(picmi.LoadAppliedField(
                read_fields_from_path=path, load_E=False, load_B=True,
                warpx_B_time_function=f"{cur:.8e}"))
            print(f"Solenoid {os.path.basename(path)}: I={cur:g} A", flush=True)

    # Prebuncher 1 (forward map) — centroid arrival uses v_beam over z_centroid→gap.
    fld1, scale1, phi1, t_gap1 = make_cavity(
        PREB1_FIELD, PREB1_KW, PREB1_Q, Z_GAP_CENTER_1, v_beam,
        PREB1_PHI_OFF, PHASE, omega, z_ref=z_centroid)
    if PREB1_KW > 0:
        applied.append(fld1)
    print(f"Preb 1: P={PREB1_KW:g} kW, Q={PREB1_Q}, scale={scale1:.3f}, "
          f"V_gap≈{scale1*V1J_KEV:.1f} kV, φ={phi1:.3f} rad, "
          f"t_gap={t_gap1*1e9:.3f} ns", flush=True)

    # Prebuncher 2 (reversed install, rev_phase=0). The Preb-2 time function is baked here
    # BEFORE WarpX integrates Preb 1, so we estimate the post-Preb-1 speed analytically (mean
    # kick) and time arrival in two segments. Valid only at the sub-threshold design point; a
    # hardened Preb-1 scan needs a two-pass run. See README -> Preb-2 timing caveat.
    if PREB2_KW > 0:
        base1 = np.pi / 2.0 if PHASE == "zc" else np.pi
        kick1 = -np.cos(base1 + np.radians(PREB1_PHI_OFF)) * scale1 * V1J_KEV
        ke_after1 = max(ke_mean + kick1, 1.0)
        gamma_a1 = 1.0 + ke_after1 / (m_e * c**2 / q_e / 1e3)
        v_after_preb1 = c * np.sqrt(1.0 - 1.0 / gamma_a1**2)
        rev_phase = PREB2_REV_PHASE if PREB2_REVERSED else 0.0
        fld2, scale2, phi2, t_gap2 = make_cavity(
            PREB2_FIELD, PREB2_KW, PREB2_Q, Z_GAP_CENTER_2, v_after_preb1,
            PREB2_PHI_OFF, PHASE, omega, t_offset=t_gap1, z_ref=Z_GAP_CENTER_1,
            rev_phase=rev_phase)
        applied.append(fld2)
        # Phase error vs. timing with bare injection β, for the sanity log.
        t_gap2_inj = (Z_GAP_CENTER_2 - z_centroid) / v_beam
        dphi_deg = np.degrees(omega * (t_gap2 - t_gap2_inj))
        print(f"Preb 2 (reversed): P={PREB2_KW:g} kW, Q={PREB2_Q}, scale={scale2:.3f}, "
              f"V_gap≈{scale2*V1J_KEV:.1f} kV, φ={phi2:.3f} rad, t_gap={t_gap2*1e9:.3f} ns "
              f"(two-segment: v_after_preb1={v_after_preb1:.3e} m/s from +{kick1:.1f} keV "
              f"Preb-1 kick; vs injection-β timing Δφ={dphi_deg:+.1f}°)", flush=True)

    # Ordering invariant: picmi sets the global E_ext_particle_init_style from the LAST-added
    # field, so when an RF cavity is present the last entry must load E. Skipped when there is
    # no RF E to protect (B-only baseline or empty drift).
    if any(getattr(f, "load_E", False) for f in applied):
        assert getattr(applied[-1], "load_E", False), (
            "last applied field must have load_E=True (an RF cavity), or the global E_ext "
            "style is forced to 'none' and the beam is unmodulated")

    electrons = picmi.Species(
        particle_type="electron", name="electrons",
        initial_distribution=picmi.ParticleListDistribution(
            x=bunch["x"], y=bunch["y"], z=bunch["z"],
            ux=bunch["ux"], uy=bunch["uy"], uz=bunch["uz"], weight=bunch["w"],
        ),
        warpx_do_not_deposit=not SPACE_CHARGE,   # SPACE_CHARGE=False → no beam self-field
    )

    # ── Time step / duration ──────────────────────────────────────────────────
    dz = ZMAX / NZ
    dt = CFL * dz / v_beam
    # Stop just before the bunch centre reaches the exit (margin < 1): once the beam clears
    # the absorbing boundary the domain empties and the Multigrid solve aborts. Size the
    # transit from the ACTUAL net kick (kick frac = -cos(base+phi_off), not the PHASE label)
    # via a 3-leg estimate with the real per-leg speed after EACH cavity (matters for the
    # net-accelerating crest path, where Preb-2's larger kick — if omitted — over-estimated
    # transit and over-ran the wall; at the energy-flat zc default both net kicks are ~0 so the
    # legs collapse to v_beam). The 0.98 margin stops short of the wall while landing a dump on
    # the 2.03 m plane.
    base1 = np.pi / 2.0 if PHASE == "zc" else np.pi
    MC2_KEV = MC2_EV / 1e3
    kick_frac1 = -np.cos(base1 + np.radians(PREB1_PHI_OFF))
    ke_after1 = max(ke_mean + (kick_frac1 * scale1 * V1J_KEV if PREB1_KW > 0 else 0.0), 1.0)
    v_after1 = c * np.sqrt(1.0 - 1.0 / (1.0 + ke_after1 / MC2_KEV) ** 2)
    if PREB2_KW > 0:
        kick_frac2 = -np.cos(base1 + np.radians(PREB2_PHI_OFF) + PREB2_REV_PHASE)
        ke_after2 = max(ke_after1 + kick_frac2 * scale2 * V1J_KEV, 1.0)
        v_after2 = c * np.sqrt(1.0 - 1.0 / (1.0 + ke_after2 / MC2_KEV) ** 2)
    else:
        ke_after2, v_after2 = ke_after1, v_after1
    v_after = v_after2          # speed AT the 2.03 m handoff (post-Preb-2); cadence-only below
    transit = ((Z_GAP_CENTER_1 - z_centroid) / v_beam
               + (Z_GAP_CENTER_2 - Z_GAP_CENTER_1) / v_after1
               + (ZMAX - Z_GAP_CENTER_2) / v_after2)
    n_steps = MAX_STEPS or int(TRANSIT_MARGIN * transit / dt)
    print(f"  Preb-1 net kick ≈ {kick_frac1*scale1*V1J_KEV:+.1f} keV (frac {kick_frac1:+.2f}); "
          f"⟨KE⟩ after Preb-1 ≈ {ke_after1:.1f} keV, after Preb-2 ≈ {ke_after2:.1f} keV", flush=True)
    print(f"dt = {dt:.3e} s, max_steps = {n_steps}, "
          f"RF period = {1/F_RF*1e9:.2f} ns ({1/F_RF/dt:.0f} steps/period)",
          flush=True)

    # ── Diagnostics (openPMD, HDF5) ───────────────────────────────────────────
    # Size `period` so dump spacing near the handoff is ≤ HANDOFF_DZ (the linac selector
    # picks the snapshot nearest ⟨z⟩=Z_HANDOFF), using v_after = post-Preb-2 speed. picmi
    # exposes only a uniform `period` (a z-station diagnostic isn't available — two same-name
    # diagnostics trip "Diagnostic attributes not consistent"). See README -> Domain / grid.
    HANDOFF_DZ = 0.008                              # [m] target dump spacing near 2.03 m
    period_handoff = max(1, int(HANDOFF_DZ / (v_after * dt)))
    period = min(max(1, n_steps // N_DIAGS), period_handoff)
    part_diag = picmi.ParticleDiagnostic(
        name="particles", period=period, species=[electrons],
        data_list=["position", "momentum", "weighting"],
        write_dir=outdir, warpx_format="openpmd", warpx_openpmd_backend="h5",
    )
    print(f"  diag period {period} steps (~{period*v_after*dt*1e3:.1f} mm near handoff; "
          f"≤{HANDOFF_DZ*1e3:.0f} mm so a dump lands near the {Z_HANDOFF*1e3:.0f} mm plane)",
          flush=True)

    sim = picmi.Simulation(
        solver=solver, max_steps=n_steps, time_step_size=dt,
        verbose=0,                     # silence per-step "STEP N starts" — tqdm is the display
        particle_shape="linear",
    )
    sim.add_species(
        electrons,
        layout=picmi.PseudoRandomLayout(n_macroparticles_per_cell=1, grid=grid),
    )
    for fld in applied:
        sim.add_applied_field(fld)
    sim.add_diagnostic(part_diag)

    print(f"\nRunning {n_steps} steps (diag every {period}) -> {outdir}/")
    run_step(sim, n_steps, desc="injector")
    print("\nDone.")

    # ── Collimator (9.547 mm iris at COLLIM_Z + pipe to ZMAX) ─────────────────
    # Diagnostic-only report of the multi-plane-collimated handoff charge (see the COLLIM_R
    # block above and README -> The 9.547 mm collimator). The physical cut is the linac
    # reader's scrape at injection.
    if COLLIMATE:
        _report_collimated_handoff(outdir)


if __name__ == "__main__":
    main()
