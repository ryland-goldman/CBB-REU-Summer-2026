"""
CESR injector in WarpX (RZ): the full LinacSim injector subsection in one
self-consistent space-charge drift — two 214 MHz prebuncher cavities (the second
installed reversed) and three static solenoid lenses (Lens 0A / Sol 0 / Lens 0E) —
handing a focused, velocity-bunched beam to `linac_sec1` at z ≈ 2.03 m.

Third stage of the Cornell Linac chain in WarpX:
    cathode (cathode/) -> gun (gun/) -> injector (this) -> linac_sec1.

The gun's exit beam (~146 keV, β≈0.63, ~1 nC, already RZ) is read from
`gun/diags/particles`, translated so it enters near z = 0, and tracked through the
injector cavities + lenses. Each cavity is the 1-J-normalised `prebuncher_25D`
field map (built by `build_injector_field.py`) driven as a standing-wave TM mode,
reproducing GPT's `Map25D_TM` convention:

    Er,Ez(t) = map · scale · cos(ω t + φ)
    Bφ(t)    = H   · scale · sin(ω t + φ)        (E and B 90° out of phase)

with f_RF = 18 × master RF = 214.18 MHz and scale = sqrt(1e3·Q·P / (2π f_RF))
from the loaded-Q / dissipated-power normalisation documented in
`reference/Linac Simulation Documentation/details.md`. Prebuncher 2 reuses the SAME
forward field (`preb2_EB.h5`, just placed at z=1.318 m); its *reversed* install is
encoded as a run-time phase, NOT a mirrored/negated map — and that phase is
PREB2_REV_PHASE=0, because crest-referencing the loaded field already absorbs the
geometric +π (see the long note at PREB2_REV_PHASE below). The operating point
constants below are config()-overridable.
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


def _retry_io(fn, *args, tries=6, base=0.25, **kwargs):
    """Call an openPMD read, retrying a transient HDF5 "Inaccessible" open error.

    NOTE: the production "OPEN_FILE failed ... Inaccessible" failure was fd
    exhaustion (openpmd-viewer leaks an fd per get_particle vs macOS's 256-fd
    default), now fixed by raising RLIMIT_NOFILE — see _runner._raise_fd_limit.
    This retry does NOT help that case (the fds stay spent). It is a backstop for
    a genuinely transient open: the in-sim handoff report opens a diag Series
    while WarpX's own diagnostic Series is still alive in this process (teardown
    is at process exit) and may be briefly releasing a just-flushed file. Re-raise
    after the last try so a genuinely missing file (or unfixed fd exhaustion)
    still surfaces.
    """
    for i in range(tries):
        try:
            return fn(*args, **kwargs)
        except io.Error:
            if i == tries - 1:
                raise
            time.sleep(base * 2 ** i)


# F_RF / Q_L / V1J_KEV / gap centres / phi-offsets live in the (pywarpx-free) build
# module as the single source of truth, so the sim and plot_injector.py cannot drift
# apart on the RF drive.
from .build_injector_field import (
    Z_GAP_CENTER_1, Z_GAP_CENTER_2, V1J_KEV, F_RF, Q_L_1, Q_L_2,
    PHI_OFF_1_DEG, PHI_OFF_2_DEG, SOL_FILES, Z_HANDOFF,
)
from . import DEFAULT_OUTDIR

c = picmi.constants.c
m_e = picmi.constants.m_e
q_e = picmi.constants.q_e

# ── Field-map paths (must match build_injector_field.py) ───────────────────────
# Both cavities use the FORWARD field values: prebuncher_25D has definite parity
# (Ez EVEN, Er ODD, Bφ EVEN about the gap; physics-measured corr ±0.9999), so the GPT
# `-1,0,0` reversed install of Preb 2 is a GLOBAL E/B sign flip ≡ +π in ABSOLUTE drive
# phase — but crest-referencing the loaded field absorbs it, so the APPLIED reversal
# phase is PREB2_REV_PHASE=0 (see the long note there). The reversal is in the time
# function, NOT a negated map. The asymmetric "z-reverse + negate Er/Bφ, keep Ez" build
# is a NO-OP for this parity (it reproduces the forward field) and is NOT used. The two
# files differ ONLY in lab-z placement (grid_global_offset, baked per openPMD file):
# preb1 at Z_GAP_CENTER_1, preb2 at Z_GAP_CENTER_2 — WarpX read_from_file takes a
# field's position from its file, so Preb 2 needs its own file despite identical values.
PREB1_FIELD = "injector/injector_field/preb1_EB.h5"   # forward field at Z_GAP_CENTER_1
PREB2_FIELD = "injector/injector_field/preb2_EB.h5"   # forward field at Z_GAP_CENTER_2; +π = reversal

GUN_DIAG = "gun/diags/particles"
Z_INJECT = 0.005                 # lab z where the bunch tail (smallest z) is placed [m]
MAX_PART = 50000                 # downsample the gun snapshot (reweighted) for speed
RNG_SEED = 0
CFL = 0.8                        # dt = CFL · Δz / v_beam

# ── Operating point — Prebuncher 1 (config()-overridable) ─────────────────────
# 8 kW is the faithful LinacSim default (prebuncher1_input_power) and is
# intentionally weak — single-cavity bunching is not the design; see README.md.
PREB1_KW = 8.0                   # dissipated RF power [kW]
PREB1_Q = Q_L_1                  # loaded Q of prebuncher 1
PREB1_PHI_OFF = PHI_OFF_1_DEG    # CREST-referenced GUI phase offset [deg] (-70)

# ── Operating point — Prebuncher 2 (reversed; config()-overridable) ───────────
PREB2_KW = 10.0                  # Prebuncher 2 design power [kW] (prebuncher2_input_power)
PREB2_Q = Q_L_2                  # loaded Q of prebuncher 2 (4300)
PREB2_PHI_OFF = PHI_OFF_2_DEG    # CREST-referenced GUI phase offset [deg] (-45)
PREB2_REVERSED = True            # apply the reversed-install phase PREB2_REV_PHASE
# Reversed-install phase added to Preb 2's time function. RESOLVED (physics-approved):
# PREB2_REV_PHASE = 0 — and geometry and the empirical run AGREE once you note where
# "crest" is defined.
#   - The geometric `-1,0,0` reversal IS a +π: a 180° rotation flips all three components
#     (Ez/Er/Bφ) given this map's Ez-even/Er-odd/Bφ-even parity = a global E,B sign flip =
#     a π phase shift, RELATIVE TO THE UN-REVERSED CAVITY AT THE SAME ABSOLUTE DRIVE PHASE.
#   - BUT we reference the drive phase to CREST (base = π = max(-cos) OF WHATEVER FIELD IS
#     LOADED). The field we load for Preb 2 is already the reversed (globally-flipped) one,
#     so base=π auto-lands on the REVERSED cavity's crest. Crest-referencing therefore
#     AUTO-ABSORBS the reversal; adding a separate +π moves 180° OFF that crest = double-
#     count → the debunching slope.
#   So in the (forward-map + crest-base + GUI φ_off) parametrization, rev_phase=0 IS the
#   geometric reversal (the +π and the crest-reference's built-in reversal cancel) — NOT its
#   absence. This rests on ONE assumption about the GUI's frame: that φ_off=−45° is referenced
#   to the AS-INSTALLED (reversed) cavity's crest. That assumption is not a standalone
#   geometric proof — it is ARBITRATED by the empirical kick-sign run below, which is the
#   decisive test. (rev_phase=+π would then double-count, debunching.)
# ARBITER — the Preb-2-only kick-sign run: REV_PHASE=0 bunches (compressive
# dchirp -0.33 keV/mm, tail gains); REV_PHASE=+π decelerates (-67 keV, no bunching).
# DO NOT "fix" this back to +π — that re-introduces the double-count. (Knob retained for a
# future map whose loaded drive phase is NOT the as-installed crest.)
PREB2_REV_PHASE = 0.0            # [rad] faithful reversed install in this crest-referenced parametrization

# ── Solenoid lens currents [A] (config()-overridable; 0 disables a lens) ───────
# LinacSim GUI defaults: Lens 0A 6 A, Sol 0 40 A, Lens 0E 10 A. These provide the
# transverse focusing that keeps the beam inside the bore through the 2.03 m handoff
# (the physical fix for the 68% radial scrape that the old linac-solenoid hack stood
# in for). The 1-A-normalised maps scale linearly with current.
I_LENS0A = 6.0
I_SOL0 = 40.0
I_LENS0E = 10.0

# ── Collimator (the faithful injector→linac iris/pipe) ────────────────────────
# LinacSim gpt_master.in: a scatteriris of radius 9.547 mm at z=1.922 m, followed by
# a 9.547 mm beam pipe to 2.1 m. So past 1.922 m the aperture is the SLAC ~9.55 mm
# bore — particles outside it are scraped. MECHANISM (the plan's accepted-minimal
# option): this is applied as a RADIAL CUT, not an in-run particle scrape — this
# pywarpx RZ build's particle-position SoA accessors raise "Component x does not exist",
# so an afterstep weight-zeroing callback is not available here. Two pieces:
#   1. _report_collimated_handoff() below prints the collimated handoff charge
#      (r ≤ COLLIM_R at the ~2.03 m dump) for the sanity log — a DIAGNOSTIC only.
#      The injector run itself is NOT collimated in-flight (its space charge over the
#      0.1 m COLLIM_Z→handoff tail includes the soon-to-be-scraped halo — a small,
#      late, β≈0.7 approximation).
#   2. the PHYSICAL cut is the linac reader's r ≤ RMAX=9.547 mm at injection (BORE_R),
#      which is equivalent to the continuous pipe because the envelope grows
#      monotonically over the 0.1 m tail (a particle inside 9.547 mm at 2.03 m was
#      inside the whole pipe; one outside hit the wall before 2.03 m).
# Do NOT widen the linac RMAX to contain a re-expanded envelope — that accepts charge
# the real iris scrapes and inflates capture.
COLLIM_R = 0.009547              # [m] iris/pipe radius (SLAC bore; gpt scatteriris)
COLLIM_Z = 1.922                # [m] iris start; the 9.547 mm pipe runs COLLIM_Z → ZMAX
COLLIMATE = True                # report the collimated handoff charge (set False to skip)

# Phase reference. "crest" (base = π) is the FAITHFUL LinacSim reference: the GUI
# phi_off values (-70 Preb-1, -45 Preb-2) are crest-referenced, so the operating
# point is base=π + phi_off. "zc" (base = π/2) is the bare zero-crossing reference
# kept ONLY for the exploratory power/phase scan (use with phi_off=0).
PHASE = "crest"                  # faithful default: crest base + GUI phi_off

# ── Performance knobs (tunable via injector.config(...); see run_pipeline.py) ──
# This stage dominates the pipeline. Do NOT coarsen NZ to go faster: this long-thin
# box is convergence-bound, not cell-bound, so the MLMG solve is slower per step at
# low NZ (and under-resolves the ~1 mm bunch) — see the README. Speed it via CFL
# (fewer steps) and MAX_ITERS/REQUIRED_PRECISION (cheaper solve).
REQUIRED_PRECISION = 1e-4        # MLMG relative tolerance (relaxed for the long-thin box)
MAX_ITERS = 500                  # MLMG iteration cap
MAX_STEPS = 0                    # 0 → auto-derive from transit; >0 → fixed
TRANSIT_MARGIN = 0.97            # stop just before the bunch centre reaches the exit
N_DIAGS = 60                     # number of openPMD dumps over the run

# ── Domain (RZ, single azimuthal mode — the cavity field is m = 0) ─────────────
# ZMAX = 2.10 m is the LinacSim prebuncher-subsection ZSTOP; the z≈2.03 m handoff
# plane (Z_acc_1) sits just inside, with a field-free exit drift past it so the
# handoff beam coasts. NR=80 (dr=0.45 mm) keeps the RF map's 36 mm bore resolved —
# do NOT copy the linac's NR=16. NZ=1664 gives dz=1.262 mm ⇒ dz/dr = 2.80:1 (the
# ≈3:1 cell-aspect rule for this long-thin box; below it the MLMG self-field solve
# diverges) and 1664 is divisible by the blocking factor 8. Do NOT coarsen NZ to go
# faster — this box is convergence-bound, not cell-bound (see README / CLAUDE.md);
# speed it via CFL and MAX_ITERS/REQUIRED_PRECISION instead.
RMAX = 0.036                     # covers the field-map bore (0–36.07 mm)
ZMAX = 2.10                      # full injector subsection (handoff at z≈2.03 m)
NR, NZ = 80, 1664                # dz=1.262 mm ⇒ 2.80:1 aspect; both ÷ blocking factor 8

OUTDIR = None                    # if None at main(), use DEFAULT_OUTDIR


def load_gun_bunch():
    """Import the gun's last snapshot (already RZ) and shift it to the entrance.

    Returns (dict for ParticleListDistribution, v_beam, mean KE [keV]).
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

    print(f"Imported {z.size} macroparticles from gun (iter {it}); "
          f"z {z.min()*1e3:.1f}–{z.max()*1e3:.1f} mm, "
          f"⟨KE⟩ {ke_mean:.1f} keV, v_beam {v_beam:.3e} m/s, "
          f"q {w.sum()*q_e*1e9:.3f} nC", flush=True)
    return (dict(x=x, y=y, z=z, ux=ux * c, uy=uy * c, uz=uz * c, w=w),
            v_beam, ke_mean)


def make_cavity(field_path, power, q_l, z_gap, v_at_gap, phi_off_deg, phase,
                omega, t_offset=0.0, rev_phase=0.0, z_ref=Z_INJECT):
    """Build one prebuncher cavity as a picmi.LoadAppliedField.

    The cavity drives the (raw 1-J) map at ``field_path`` as a standing-wave TM
    mode: E ∝ scale·cos(ωt+φ), B ∝ scale·sin(ωt+φ). ``scale`` = sqrt(1e3·Q·P/(2πf))
    sets the amplitude from dissipated power; ``φ`` phases the bunch at the cavity
    gap. The GUI ``phi_off_deg`` is **crest-referenced**, so the faithful operating
    point uses a crest base (``phase="crest"`` → base = π) plus ``phi_off_deg``.
    ``phase="zc"`` (base = π/2) is the bare zero-crossing reference kept only for the
    exploratory scan.

    ``v_at_gap`` is the mean beam speed over the leg ``z_ref → z_gap``; ``t_offset`` is
    the time already elapsed reaching ``z_ref``. Prebuncher 1 uses z_ref=Z_INJECT,
    t_offset=0, v_at_gap=v_beam. Prebuncher 2 uses z_ref=Z_GAP_CENTER_1, t_offset=t_gap1,
    v_at_gap=v_after_preb1 — a TWO-SEGMENT arrival (v_beam to Z1, then the post-Preb-1 β
    over Z1→Z2), so the +~15 keV Preb-1 kick that speeds the beam over the inter-cavity
    drift is accounted for (cuts the constant-injection-β phase error from ~10° to ~few°).

    REVERSED INSTALL (``rev_phase``, Preb 2): GPT's `-1,0,0` is a 180° rotation. For a
    standing-wave TM map it is a sign flip of the rotation-odd field components, which
    can map onto a time-phase shift on the cos/sin drive — but the exact value (±π or
    already folded into the GUI -45° offset) is uncertain and is resolved EMPIRICALLY
    by the Preb-2-only kick-sign diagnostic (caller passes PREB2_REV_PHASE). ``rev_phase``
    is added to φ. Keep .10e precision on every term — ω·t truncation accumulates over
    the ~5 ns transit at 214 MHz.
    """
    scale = float(np.sqrt(1e3 * q_l * power / (2.0 * np.pi * F_RF)))
    # Arrival time of the bunch tail at this cavity's gap: time to z_ref + leg z_ref→z_gap.
    t_gap = t_offset + (z_gap - z_ref) / v_at_gap
    # The energy kick of an electron is ΔW(t) ∝ -cos(ω t + φ) (on-axis Ez is
    # single-signed positive). The GUI phi_off is CREST-referenced, so:
    #   crest: base = π    (faithful path; phi_off measured from crest)
    #   zc:    base = π/2  (bare zero-crossing; exploratory scan only)
    # The reversed install adds +π (global E/B sign flip); base/phi_off are NOT
    # otherwise touched by the reversal.
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

    Reads the snapshot nearest ⟨z⟩ = Z_HANDOFF and applies the 9.547 mm iris/pipe cut
    (r ≤ COLLIM_R), printing the charge that survives the aperture vs the in-domain
    charge — the real transmission through the injector→linac iris. This is a diagnostic
    print only; the physical cut at injection is applied by the linac reader (BORE_R).
    """
    try:
        from openpmd_viewer import OpenPMDTimeSeries
        ts = OpenPMDTimeSeries(os.path.join(outdir, "particles"))
        best, bd = None, 9e9
        for it in ts.iterations:
            x, y, z, w = _retry_io(ts.get_particle, ["x", "y", "z", "w"],
                                   species="electrons", iteration=it)
            if len(z) < 50:
                continue
            zm = float(np.average(z, weights=w))
            if abs(zm - Z_HANDOFF) < bd:
                bd, best = abs(zm - Z_HANDOFF), (zm, x, y, z, w)
        if best is None:
            print("  collimated handoff: no populated snapshot near the plane", flush=True)
            return
        zm, x, y, z, w = best
        r = np.hypot(x, y)
        q_dom = float(w.sum()) * q_e
        q_coll = float(w[r <= COLLIM_R].sum()) * q_e
        print(f"  COLLIMATED handoff (⟨z⟩={zm*1e3:.1f} mm, iris {COLLIM_R*1e3:.3f} mm): "
              f"{q_coll*1e9:.3f} nC within the iris / {q_dom*1e9:.3f} nC in-domain "
              f"= {100*q_coll/q_dom:.0f}% through the aperture", flush=True)
    except Exception as e:
        print(f"  collimated-handoff report unavailable: {e}", flush=True)


def main():
    outdir = OUTDIR or DEFAULT_OUTDIR

    # Fresh diags: WarpX appends one openPMD file per dump, so re-running the same
    # case would otherwise mix old and new iterations into one series and corrupt the
    # focus/σ_z analysis. diags are git-ignored and regenerated, so clearing is safe.
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)

    # Compute omega here, not at import, so a config(F_RF=...) override is honored
    # (an import-time module constant would be frozen before the override lands).
    omega = 2.0 * np.pi * F_RF

    bunch, v_beam, ke_mean = load_gun_bunch()

    # ── Grid + electrostatic (self-field) solver ──────────────────────────────
    grid = picmi.CylindricalGrid(
        number_of_cells=[NR, NZ],
        n_azimuthal_modes=1,
        lower_bound=[0.0, 0.0],
        upper_bound=[RMAX, ZMAX],
        lower_boundary_conditions=["none", "dirichlet"],
        upper_boundary_conditions=["neumann", "dirichlet"],
        lower_boundary_conditions_particles=["none", "absorbing"],
        upper_boundary_conditions_particles=["absorbing", "absorbing"],
        warpx_blocking_factor=8,
    )
    solver = picmi.ElectrostaticSolver(
        grid=grid, method="Multigrid", required_precision=REQUIRED_PRECISION,
        maximum_iterations=MAX_ITERS, warpx_self_fields_verbosity=0,
    )

    # ── Applied fields ────────────────────────────────────────────────────────
    # Solenoid B-only maps come BEFORE the RF cavities: picmi forces the global
    # E_ext_particle_init_style to "none" if the LAST-added LoadAppliedField has
    # load_E=False, so a B-only map last would silently disable the accelerating E.
    # The 1-A-normalised maps are scaled by a constant time function = the current.
    applied = []
    for path, cur in [(SOL_FILES["LENS_0A"], I_LENS0A),
                      (SOL_FILES["SOL_0"],   I_SOL0),
                      (SOL_FILES["LENS_0E"], I_LENS0E)]:
        if cur != 0.0:
            applied.append(picmi.LoadAppliedField(
                read_fields_from_path=path, load_E=False, load_B=True,
                warpx_B_time_function=f"{cur:.8e}"))
            print(f"Solenoid {os.path.basename(path)}: I={cur:g} A", flush=True)

    # Prebuncher 1 (forward map) — arrival uses v_beam.
    fld1, scale1, phi1, t_gap1 = make_cavity(
        PREB1_FIELD, PREB1_KW, PREB1_Q, Z_GAP_CENTER_1, v_beam,
        PREB1_PHI_OFF, PHASE, omega)
    if PREB1_KW > 0:
        applied.append(fld1)
    print(f"Preb 1: P={PREB1_KW:g} kW, Q={PREB1_Q}, scale={scale1:.3f}, "
          f"V_gap≈{scale1*V1J_KEV:.1f} kV, φ={phi1:.3f} rad, "
          f"t_gap={t_gap1*1e9:.3f} ns", flush=True)

    # Prebuncher 2 (reversed install, rev_phase=0 — see the PREB2_REV_PHASE note above for
    # why 0, not +π, is the faithful reversal in this crest-referenced parametrization).
    # ARRIVAL TIMING (two-segment): t_gap2 = t_gap1 + (Z2−Z1)/v_after_preb1. The Preb-2 time
    # function is baked here, BEFORE WarpX integrates Preb 1, so the true post-Preb-1 β is not
    # yet known. The faithful crest-base Preb-1 imparts a mean +~15 keV that SPEEDS the beam
    # over the 0.534→1.318 m inter-cavity drift; timing Preb-2 with the bare injection β would
    # mis-time arrival by ≈ −13° at 214 MHz. So we estimate the post-Preb-1 speed ANALYTICALLY
    # from the same mean-kick fraction the transit uses (−cos(base+φ_off)·scale1·V1J) and time
    # Preb-2 in two segments (v_beam to Z1, then v_after_preb1 over Z1→Z2), cutting the residual
    # to ~few°. This is an analytic estimate of the MEAN kick, not the true post-cavity β
    # distribution: valid ONLY at the sub-threshold design point; a hardened Preb-1 power scan
    # desyncs the Preb-2 reference and needs a two-pass run (read the diagnostic β, rebuild
    # Preb-2). See injector/README.md.
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
        t_gap2_inj = (Z_GAP_CENTER_2 - Z_INJECT) / v_beam
        dphi_deg = np.degrees(omega * (t_gap2 - t_gap2_inj))
        print(f"Preb 2 (reversed): P={PREB2_KW:g} kW, Q={PREB2_Q}, scale={scale2:.3f}, "
              f"V_gap≈{scale2*V1J_KEV:.1f} kV, φ={phi2:.3f} rad, t_gap={t_gap2*1e9:.3f} ns "
              f"(two-segment: v_after_preb1={v_after_preb1:.3e} m/s from +{kick1:.1f} keV "
              f"Preb-1 kick; vs injection-β timing Δφ={dphi_deg:+.1f}°)", flush=True)

    # Enforce the ordering invariant: picmi sets the *global* E_ext_particle_init_style
    # from the LAST-added field, so when an RF cavity (load_E) IS present the LAST entry
    # must load E — else a trailing B-only solenoid silently disables the accelerating E.
    # The guard only matters when there's an RF E field to protect: a baseline with only
    # B-only solenoids (no RF) has no E to disable, and a pure drift (empty list) has no
    # field at all — both legitimately skip it.
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
    )

    # ── Time step / duration ──────────────────────────────────────────────────
    dz = ZMAX / NZ
    dt = CFL * dz / v_beam
    # Stop just before the bunch centre reaches the exit (margin < 1): once the beam
    # clears the absorbing boundary the domain empties and the Multigrid solve aborts.
    # Size the transit from the ACTUAL net energy kick at the gap, not the PHASE label:
    # the faithful operating point is PHASE="crest" base + phi_off=-70 (a partial kick),
    # NOT full on-crest acceleration, so keying off PHASE=="crest" would assume the full
    # scale·V1J gain and stop far too early. The on-axis energy kick of an electron is
    # ΔW ∝ -cos(ω t_gap + φ); at the gap ω t_gap + φ = base + phi_off (the -ω t_gap
    # cancels), so the net kick fraction is -cos(base + radians(phi_off)). The kick can
    # decelerate (negative), which SLOWS the beam, so the transit must be the LONGER of
    # the kicked and unkicked estimates to guarantee the run still spans the box.
    base1 = np.pi / 2.0 if PHASE == "zc" else np.pi
    kick_frac1 = -np.cos(base1 + np.radians(PREB1_PHI_OFF))
    ke_after = ke_mean + kick_frac1 * scale1 * V1J_KEV
    ke_after = max(ke_after, 1.0)                       # guard: keep γ real if over-decel
    gamma_a = 1.0 + ke_after / (m_e * c**2 / q_e / 1e3)
    v_after = c * np.sqrt(1.0 - 1.0 / gamma_a**2)
    transit_kicked = ((Z_GAP_CENTER_1 - Z_INJECT) / v_beam
                      + (ZMAX - Z_GAP_CENTER_1) / v_after)
    transit_coast = (ZMAX - Z_INJECT) / v_beam
    transit = max(transit_kicked, transit_coast)       # never stop short of the box
    n_steps = MAX_STEPS or int(TRANSIT_MARGIN * transit / dt)
    print(f"  Preb-1 net kick ≈ {kick_frac1*scale1*V1J_KEV:+.1f} keV "
          f"(frac {kick_frac1:+.2f}); ⟨KE⟩ after ≈ {ke_after:.1f} keV", flush=True)
    print(f"dt = {dt:.3e} s, max_steps = {n_steps}, "
          f"RF period = {1/F_RF*1e9:.2f} ns ({1/F_RF/dt:.0f} steps/period)",
          flush=True)

    # ── Diagnostics (openPMD, HDF5) ───────────────────────────────────────────
    # The dump cadence (period) must be fine enough that one snapshot lands within a
    # few mm of the z≈2.03 m handoff plane — the linac selector picks the snapshot
    # nearest ⟨z⟩=Z_HANDOFF, so a coarse cadence that straddles the plane by tens of mm
    # would hand off an off-plane beam. The bunch advances v_after·dt ≈ 1.2 mm/step near
    # the handoff, so we size `period` to keep the handoff-region spacing ≤ HANDOFF_DZ
    # (≈8 mm), then take the finer of that and the N_DIAGS cadence. (picmi exposes only a
    # single uniform `period`; a true z-station / multi-interval diagnostic isn't
    # available through this picmi build — two same-name ParticleDiagnostics trip the
    # "Diagnostic attributes not consistent" assertion, and `warpx_intervals` is rejected
    # by picmistandard — so a uniformly fine cadence is the portable mechanism.)
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
    # The faithful injector→linac aperture (gpt scatteriris 9.547 mm at 1.922 m + a
    # 9.547 mm pipe to 2.1 m). It is applied as a RADIAL CUT on the openPMD diagnostic
    # snapshots (r > COLLIM_R for z ≥ COLLIM_Z → scraped), NOT as an in-run particle
    # scrape: this pywarpx RZ build's particle-position SoA accessors raise "Component x
    # does not exist" (the radial position is the AMReX particle position, not a named
    # real comp), so an afterstep weight-zeroing callback is not reliably available here.
    # Because the pipe HOLDS 9.547 mm continuously from COLLIM_Z to ZMAX and the envelope
    # grows monotonically over that 0.1 m tail, the radial cut at the 2.03 m handoff
    # plane is EQUIVALENT to the continuous pipe (any particle outside 9.547 mm at 2.03 m
    # hit the pipe wall before reaching it; any inside was inside all along). The only
    # approximation is the self-field of the scraped halo over the ~0.1 m COLLIM_Z→handoff
    # tail (a small, late-stage, near-relativistic correction). The linac reader applies
    # the same 9.547 mm cut at injection (BORE_R), so the handoff beam IS collimated; do
    # NOT widen the linac RMAX to contain an uncollimated 36 mm envelope. The collimated
    # handoff charge is reported by collimated_handoff_charge() below for the sanity log.
    if COLLIMATE:
        _report_collimated_handoff(outdir)


if __name__ == "__main__":
    main()
