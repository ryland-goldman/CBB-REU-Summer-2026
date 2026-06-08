"""
Per-section field-scale calibration + validation gates for the ``linac_rest`` Impact-T stage
(Tasks 5 and the §5 validation gates of the plan).

The seven travelling-wave sections (2–8) reuse one S-band field *shape* (``rfdata4–7``); ALL
per-section physics lives in the calibrated ``rf_field_scale`` (the solrf element has no scalar
gradient input, so the scale is fit, not computed analytically — plan §1/§4).

This module is import-only (no module-level side effects) so ``linac_rest_sim.py`` can call:

    scales = calibrate_sections(I, P_in, power_mw=POWER_MW)   # Task 5
    gates  = validate_run(I, P_in, power_mw=POWER_MW)         # §5 gates, after I.run()

Conventions established while writing this (verified against the installed lume-impact 0.x):

* ``autophase_and_scale``'s ``metric="mean_energy"`` is **total** energy in eV (rest mass
  included), NOT kinetic energy. The calibration ``target`` is therefore
  ``(mean total energy entering the section) + ΔE_target``, computed per section by tracking the
  beam through the already-calibrated upstream sections.
* The 4 solrf sub-elements of one section are linked with a ``ControlGroup`` (``add_group``) with
  ``absolute=True`` and ``factors = [1, 1/sin(β₀d), 1/sin(β₀d), 1]`` so setting the group's
  ``rf_field_scale`` value S sets entrance/exit = S and the two body cells = S/sin(β₀d). That
  **preserves the template's 1/sin body ratio across the group** for free (ControlGroup
  ``set_absolute`` multiplies each element's attribute by its factor).
* Phase is **pinned on-crest (θ₀ = 0)** — at β > 0.999 the synchronous-phase reference is trivial
  and per-section re-phasing risks drifting off-crest to a lower scale that still hits the
  (sub-maximal) target. We calibrate **scale only** (a 1-D ``brentq`` on the group value), reusing
  the ``autophase_and_scale`` scale-step semantics. This honours the plan's "fix θ₀ = 0, calibrate
  scale only" while still using the lume-impact machinery (``add_group`` + the same root-find as
  ``autophase_and_scale_brent2``'s scale step).
"""

import math

from scipy.optimize import brentq

from . import build_linac_rest_lattice as L

ELECTRON_REST_MEV = 0.51099895069  # electron rest energy [MeV]


# ── element-name scheme ───────────────────────────────────────────────────────
# T3 (build_linac_rest_lattice) OWNS the element naming and exposes it as
# ``section_group_names(index)`` → ["sec{N}_entrance","sec{N}_body_1","sec{N}_body_2",
# "sec{N}_exit"] (N=index+2). We delegate to that single source of truth so the calibration
# can never drift from the actual deck names (the body cells are body_1/body_2, with underscores).
def section_ele_names(index):
    """(entrance, body_1, body_2, exit) solrf element names for section ``index`` (0-based)."""
    return tuple(L.section_group_names(index))


def section_group_name(index):
    return f"sec{index + 2}_scale"


def _ensure_section_group(I, index):
    """Create (idempotently) the rf_field_scale ControlGroup over a section's 4 solrf cells.

    factors = [1, 1/sin(β₀d), 1/sin(β₀d), 1] with absolute=True ⇒ setting the group value S sets
    entrance/exit scale = S and body scale = S/sin(β₀d), preserving the template body ratio.
    Returns the group name.
    """
    gname = section_group_name(index)
    if gname in getattr(I, "group", {}):
        return gname
    inv_sin = 1.0 / L.SIN_BETA0_D
    I.add_group(
        gname,
        ele_names=list(section_ele_names(index)),
        var_name="rf_field_scale",
        factors=[1.0, inv_sin, inv_sin, 1.0],
        absolute=True,
    )
    return gname


def _set_group_scale(I, gname, value):
    """Set a section's ControlGroup field-scale value.

    lume-impact's ``Impact.__setitem__`` requires the ``"name:attribute"`` form; for a group the
    attribute IS the group's ``var_name`` (``rf_field_scale`` here). A bare ``I[gname] = value``
    raises (it tries to split on ':'), so always go through this helper.
    """
    I[f"{gname}:rf_field_scale"] = float(value)


def _set_section_phase(I, index, phase_deg):
    """Pin every solrf sub-element of a section to the on-crest driven phase (with the template's
    fixed inter-line offsets entrance, +30°, +90°, exit added on top of phase_deg)."""
    entrance, body1, body2, exit_ = section_ele_names(index)
    I[entrance]["theta0_deg"] = phase_deg + 0.0
    I[body1]["theta0_deg"] = phase_deg + 30.0
    I[body2]["theta0_deg"] = phase_deg + 90.0
    I[exit_]["theta0_deg"] = phase_deg + 0.0



def _parabolic_peak(phases, gains, step):
    """Parabolic-refine the max of (phases, gains) sampled on a uniform `step` grid.

    `phases` need not be sorted/wrapped; we bracket the sampled max and unwrap the ±step
    neighbours onto a monotone axis before fitting. Returns the refined phase (deg, [0,360)).
    """
    k = max(range(len(phases)), key=lambda j: gains[j])
    p1, g1 = phases[k], gains[k]
    p0 = phases[(k - 1) % len(phases)]
    p2 = phases[(k + 1) % len(phases)]
    g0 = gains[(k - 1) % len(phases)]
    g2 = gains[(k + 1) % len(phases)]
    if p0 > p1:
        p0 -= 360.0
    if p2 < p1:
        p2 += 360.0
    denom = (g0 - 2 * g1 + g2)
    phase_star = p1 - 0.5 * step * (g2 - g0) / denom if denom != 0 else p1
    return phase_star % 360.0


def _find_crest_phase(Ic, index, gname, P_entrance, s1, probe_scale,
                      seed_phase=None, coarse_step=15.0, fine_half_window=18.0,
                      fine_step=9.0):
    """Find the base phase [deg] that puts section ``index`` on crest (max ΔE) for the bunch's
    actual arrival time.

    Impact-T ``theta0_deg`` is ABSOLUTE (referenced to t=0), so each downstream section's crest is
    at a different base phase (the later-arriving bunch sees a shifted RF phase). Two modes:

    * ``seed_phase=None`` (the first section): full 0–360° coarse scan (24 pts at 15°) + parabolic
      refine — no prior to anchor on.
    * ``seed_phase`` given (downstream sections): the crest moves with arrival time by
      −2πf·Δt_arrival, so the caller seeds with (prev crest + analytic arrival shift); we only do a
      tight LOCAL scan (±``fine_half_window`` at ``fine_step`` ≈ 5 pts) + parabolic refine around
      that seed. This cuts the per-section crest cost from ~24 tracks to ~5 (the dominant runtime).

    The crest phase is ~scale-independent, so the scale fit afterward keeps the section on crest.
    Returns the crest base phase (deg, [0,360)).
    """
    _set_group_scale(Ic, gname, probe_scale)
    e_in = P_entrance["mean_energy"]

    def de_at(phase):
        _set_section_phase(Ic, index, float(phase))
        P = Ic.track(P_entrance, s=s1)
        return ((P["mean_energy"] - e_in) if P else -e_in)

    if seed_phase is None:
        phases = [coarse_step * k for k in range(int(round(360.0 / coarse_step)))]
        gains = [de_at(ph) for ph in phases]
        return _parabolic_peak(phases, gains, coarse_step)

    # Local refine around the analytic seed.
    n_side = int(round(fine_half_window / fine_step))
    offsets = [fine_step * j for j in range(-n_side, n_side + 1)]
    phases = [(seed_phase + d) % 360.0 for d in offsets]
    gains = [de_at(ph) for ph in phases]
    # If the seed missed (peak at a window edge), fall back to the full coarse scan once.
    k = max(range(len(gains)), key=lambda j: gains[j])
    if k == 0 or k == len(gains) - 1:
        phases = [coarse_step * m for m in range(int(round(360.0 / coarse_step)))]
        gains = [de_at(ph) for ph in phases]
        return _parabolic_peak(phases, gains, coarse_step)
    return _parabolic_peak(phases, gains, fine_step)


def calibrate_sections(I, P_in, power_mw=None, scale_range=(5e6, 90e6),
                       rtol=2e-3, probe_scale=2.0e7, np_calib=400, verbose=True,
                       bar=None):
    """Calibrate each section to its energy-gain target ΔE_target(P_op), ON CREST.

    Per section i (lab order), sequential. To be IMMUNE to cross-section state bleed (a reused
    Impact copy re-links its ControlGroups to rebuilt ele dicts after each track(), so a group
    added early goes stale and later scale writes silently no-op), each section gets a FRESH copy
    with all UPSTREAM fitted (scale, crest-phase) re-applied, and the beam is tracked from z=0:
      1. fresh Ic = I.copy(); apply upstream fitted scales+phases; zero this section and all
         downstream sections; track P_in from 0 to this section's entrance s0 → entrance energy;
      2. **find the crest base phase** — the absolute θ₀ maximizing ΔE for the bunch's actual
         arrival (Impact-T θ₀ is absolute ⇒ crest ≠ 0° downstream; coarse scan + parabolic refine);
      3. at crest, brentq the section's rf_field_scale ControlGroup to hit
         target_total = entrance energy + ΔE_target,i(P_op);
      4. store the fitted (scale, crest-phase); apply onto the live deck I.

    O(N²) tracks (re-track upstream each section) but fully correct — SC-off tracks are cheap.

    PERF: the calibration is decimated and seeded so it doesn't dominate runtime:
    * ``np_calib`` — calibrate on a DECIMATED bunch (~400 macroparticles): ⟨KE⟩/crest are mean
      quantities, insensitive to particle count, so the fit is unchanged but every track is
      cheaper. The final ``I.run()`` (caller) uses the full ``Np``.
    * crest SEEDING — section 0 does a full coarse phase scan; sections ≥1 seed the crest from the
      previous section's crest + the analytic arrival-time shift (−2πf·Δt_arrival, since at β>0.999
      ΔE(φ)=A·cos(φ−φc) with φc moving as the bunch's RF arrival phase) and only LOCAL-refine
      (~5 pts) — cutting the dominant 24-pt/section scan to ~5 pts.

    Returns a list of dicts (index, name, scale, crest_phase_deg, entrance_ke_mev, target_de_mev,
    achieved_de_mev, err_frac) and leaves ``I`` with every section's calibrated scale AND crest
    phase applied (ready for the full ``I.run()``).

    ``bar`` (optional tqdm) ticks once per section so the calibration phase shows progress like
    the WarpX stages; per-section summary lines route through ``bar.write`` so they don't scramble
    the bar. When ``bar`` is None the lines just ``print`` (standalone / direct call).
    """
    def _emit(msg):
        if bar is not None:
            bar.write(msg)
        else:
            print(msg)
    p = L.POWER_MW if power_mw is None else power_mw
    from impact.lattice import ele_bounds
    results = []
    fitted = []   # (scale, crest_phase_deg) per calibrated upstream section

    # Decimate the calibration bunch (mean ⟨KE⟩/crest are count-insensitive; final run uses full Np).
    if np_calib and P_in.n_particle > np_calib:
        stride = P_in.n_particle // np_calib
        P_cal = P_in[::stride]
    else:
        P_cal = P_in

    omega = 2.0 * math.pi * L.RF_FREQ_HZ
    prev_crest = None        # previous section's crest base phase [deg]
    prev_t_arr = None        # previous section's bunch arrival mean_t [s]

    for i in range(L.N_SECTIONS):
        # 1) Fresh copy so the ControlGroups link to THIS copy's current ele dicts (no stale link).
        Ic = I.copy()
        Ic.verbose = False
        Ic.configure()
        gnames = [_ensure_section_group(Ic, j) for j in range(L.N_SECTIONS)]
        for j in range(L.N_SECTIONS):
            if j < i:                                   # upstream: re-apply fitted scale + phase
                sc, ph = fitted[j]
                _set_section_phase(Ic, j, ph)
                _set_group_scale(Ic, gnames[j], sc)
            else:                                       # this section + downstream: OFF
                _set_section_phase(Ic, j, L.PHASE_DEG)
                _set_group_scale(Ic, gnames[j], 0.0)

        names = set(section_ele_names(i))
        s0, s1 = ele_bounds([e for e in Ic.lattice if e.get("name") in names])
        gname = gnames[i]

        # Track from z=0 through the calibrated upstream sections to this section's entrance.
        # GUARD: never call track with s ≤ 0 — Impact-T treats stop=0.0 as "no early stop" and runs
        # the ENTIRE lattice, so for section 2 (s0=0) with a ~0/slightly-negative input mean_z (the
        # z-zeroing in load leaves mean_z ≈ -1e-15) the naive `mean_z < s0` would track the whole
        # 36 m line and return the line-exit beam. For the first section the beam is already at the
        # entrance, so skip the track; otherwise track only if genuinely upstream of s0.
        if s0 > 1e-9 and P_cal["mean_z"] < s0 - 1e-9:
            P_entrance = Ic.track(P_cal, s=s0)
        else:
            P_entrance = P_cal
        e_entrance = P_entrance["mean_energy"]              # total energy [eV]
        ke_entrance_mev = (e_entrance / 1e6) - ELECTRON_REST_MEV

        de_target_mev = L.section_de_target(i, p)
        target_total = e_entrance + de_target_mev * 1e6     # total energy target [eV]

        # 2) Find the crest base phase for this section at the bunch's arrival (θ₀ is absolute).
        # Seed from the previous section's crest + the analytic arrival-time shift (−ω·Δt_arr,
        # deg) so sections ≥1 only LOCAL-refine; section 0 does the full coarse scan (seed=None).
        t_arr = P_entrance["mean_t"]
        seed = None
        if prev_crest is not None:
            shift_deg = -math.degrees(omega * (t_arr - prev_t_arr))
            seed = (prev_crest + shift_deg) % 360.0
        crest_phase = _find_crest_phase(Ic, i, gname, P_entrance, s1, probe_scale,
                                        seed_phase=seed)
        _set_section_phase(Ic, i, crest_phase)
        prev_crest, prev_t_arr = crest_phase, t_arr

        # 3) At crest, calibrate the field scale to hit ΔE_target. Guard the bracket: if the scaled
        # energy doesn't change across [lo, hi] the group link is broken — fail loudly with context.
        def gain_minus_target(S):
            _set_group_scale(Ic, gname, S)
            P = Ic.track(P_entrance, s=s1)
            en = P["mean_energy"] if P else 0.0
            return en / target_total - 1.0

        f_lo, f_hi = gain_minus_target(scale_range[0]), gain_minus_target(scale_range[1])
        if f_lo == f_hi:
            raise RuntimeError(
                f"section {i + 2}: tracked energy does not vary with rf_field_scale "
                f"(f(lo)=f(hi)={f_lo:.4g}) — ControlGroup '{gname}' link is stale or the deck "
                f"isn't picking up the scale; cannot bracket the scale fit.")
        if f_lo * f_hi > 0:
            raise RuntimeError(
                f"section {i + 2}: ΔE target {de_target_mev:.1f} MeV unreachable in scale range "
                f"{scale_range} at crest {crest_phase:.1f}° (f(lo)={f_lo:.3g}, f(hi)={f_hi:.3g}).")
        S_fit = brentq(gain_minus_target, scale_range[0], scale_range[1],
                       maxiter=40, rtol=rtol)

        # Achieved gain at the fitted scale.
        _set_group_scale(Ic, gname, S_fit)
        P_exit = Ic.track(P_entrance, s=s1)
        achieved_de_mev = (P_exit["mean_energy"] - e_entrance) / 1e6
        err_frac = achieved_de_mev / de_target_mev - 1.0

        fitted.append((S_fit, crest_phase))
        _ensure_section_group(I, i)
        _set_section_phase(I, i, crest_phase)   # write crest phase onto the live deck
        _set_group_scale(I, gname, S_fit)       # write the calibrated scale onto the live deck

        rec = {
            "index": i, "name": L.SECTIONS[i]["name"], "scale": S_fit,
            "crest_phase_deg": crest_phase,
            "entrance_ke_mev": ke_entrance_mev, "target_de_mev": de_target_mev,
            "achieved_de_mev": achieved_de_mev, "err_frac": err_frac,
        }
        results.append(rec)
        if verbose:
            _emit(f"  sec {i + 2} {rec['name']:<6}  KE_in={ke_entrance_mev:7.2f} MeV  "
                  f"crest={crest_phase:6.1f}°  scale={S_fit:.4e}  ΔE target={de_target_mev:6.2f}  "
                  f"achieved={achieved_de_mev:6.2f} MeV  ({err_frac * 100:+.2f}%)")
        if bar is not None:
            bar.set_postfix_str(f"sec {i + 2}, {ke_entrance_mev + achieved_de_mev:.0f} MeV")
            bar.update(1)

    if verbose:
        worst = max(abs(r["err_frac"]) for r in results)
        _emit(f"  calibration worst |error| = {worst * 100:.2f}%  (gate: ±3%)")
    return results


# ── §5 validation gates ──────────────────────────────────────────────────────
# (The MIN_KE_MEV model-validity cut lives in linac_rest_sim.load_sec1_core — the single source
# of MIN_KE_MEV — not here; calibration only consumes the already-cut captured-core beam.)
def _beta_from_ke_mev(ke_mev):
    g = 1.0 + ke_mev / ELECTRON_REST_MEV
    return math.sqrt(max(0.0, 1.0 - 1.0 / (g * g)))


def validate_run(I, P_in, power_mw=None, calib=None, require_gates=False):
    """Compute + print the plan §5 validation gates after a full ``I.run()``.

    Gates:
      1. per-section ΔE vs ΔE_target — ±3% (from ``calib`` if supplied).
      2. cumulative exit ⟨KE⟩ ≈ measured ⟨KE⟩_in + Σ ΔE_target (predicted ≈307 MeV @11 MW;
         achieved ≈308 MeV — 307.97 survivors through the real bore / 309.2 full-beam).
      3. σ_KE absolute (NOT conserved — grows ~3.9× from second-order crest curvature) and
         relative (still shrinks, ⟨KE⟩ grows faster). Diagnostic only.
      4. normalized emittance εn,x / εn,y in vs out — diagnostic. NOTE: quads-OFF εn is NOT
         conserved; the recorded vs-z εn sawtooths ~2.4× — a fort.10N norm_emit artifact at
         bore/section crossings (σ_x stays smooth across the jumps ⇒ not physical growth).
      5. beam reached final zedge: I.stat("mean_z")[-1] ≈ Σ lattice length (catches Ntstep
         truncation, which falsely reports finished=True).
      6. min captured KE ⇒ β > 0.999 (justifies the rigid-crest no-slip assumption).
      7. transmission against the real bore (count-based n_out/n_in) — a no-focusing LOWER BOUND
         (~78.5%), NOT a prediction; the real FODO contains the beam. Diagnostic.

    Returns a dict of computed values + a boolean per gate. Prints a clear PASS/FAIL line each.

    If ``require_gates`` is True, raises ``AssertionError`` on any FAILED **hard** gate (the four
    the team locked as must-pass: per-section ΔE ±3%, cumulative exit ⟨KE⟩ ±3%, beam reached the
    final zedge, min-captured-KE β > 0.999). The **soft/diagnostic** gates (σ_KE evolution, εn
    growth, transmission) are always print-only — they characterize the beam but are not pass/fail
    criteria (εn/transmission are only meaningful with quad optics / aperture scraping, which are
    OFF for the headline). T6 calls this with ``require_gates=True`` so a bad run fails loudly.
    """
    p = L.POWER_MW if power_mw is None else power_mw
    P_out = I.particles["final_particles"]

    ke_in_mev = (P_in["mean_energy"] / 1e6) - ELECTRON_REST_MEV
    ke_out_mev = (P_out["mean_energy"] / 1e6) - ELECTRON_REST_MEV
    sum_de = sum(L.section_de_target(i, p) for i in range(L.N_SECTIONS))
    expected_out = ke_in_mev + sum_de

    # σ_KE (absolute, MeV). ParticleGroup energy is total; σ of KE == σ of total energy.
    sig_ke_in = P_in["sigma_energy"] / 1e6
    sig_ke_out = P_out["sigma_energy"] / 1e6

    enx_in, eny_in = P_in["norm_emit_x"], P_in["norm_emit_y"]
    enx_out, eny_out = P_out["norm_emit_x"], P_out["norm_emit_y"]

    mean_z_reached = I.stat("mean_z")[-1]
    z_expected = L.total_lattice_length_m()

    # min captured KE (over surviving particles) → β.
    ke_min_mev = (P_out["energy"].min() / 1e6) - ELECTRON_REST_MEV
    beta_min = _beta_from_ke_mev(ke_min_mev)

    # Envelope-in-bore soft gate (NEW, T6): does the transverse RMS envelope stay inside the
    # real tapered bore along the whole line? σ_x/σ_y are RMS, not the beam edge, so we test a
    # BEAM_EDGE_SIGMA=3σ multiple (a meaningful edge-vs-bore proxy) against the NARROWEST bore
    # radius (entrance + exit of every section). This is the meaningful check that the FODO
    # actually CONTAINS the beam (vs. a few core particles squeaking through). It is print-only /
    # soft — K1 is the guessed-strength placeholder FODO (A→T undocumented), so it must never
    # gate the energy headline. NOT a "FAIL-when-quads-OFF" liveness test: the quads-OFF beam
    # transmits ~78.5% with adiabatic damping, so its RMS σ may stay under the bore even while
    # distribution tails scrape — this gate can legitimately PASS quads-OFF (liveness is shown by
    # the figure: bounded oscillation quads-ON vs. monotonic RMS rise quads-OFF).
    BEAM_EDGE_SIGMA = 3.0
    sx = I.stat("sigma_x")
    sy = I.stat("sigma_y")
    max_env_m = BEAM_EDGE_SIGMA * float(max(sx.max(), sy.max())) if len(sx) else float("nan")
    min_bore_m = min(min(L.section_bore_radii(i)) for i in range(L.N_SECTIONS))

    # Transmission from the macroparticle COUNT, NOT charge. Count-based is the authoritative
    # measure (it's what sim.main() records); a charge ratio would only equal it because main()
    # re-imposes q_out = q_core·(n_out/n_in) before calling this — a hidden ordering dependency.
    # Computing from counts here keeps the gate correct regardless of when charge is re-imposed.
    n_in = P_in.n_particle
    n_out = P_out.n_particle
    transmission = (n_out / n_in) if n_in else 0.0

    gates = {}
    gates["calib_within_3pct"] = (
        all(abs(r["err_frac"]) <= 0.03 for r in calib) if calib else None
    )
    gates["exit_ke_mev"] = ke_out_mev
    gates["expected_exit_ke_mev"] = expected_out
    gates["exit_ke_within_3pct"] = abs(ke_out_mev / expected_out - 1.0) <= 0.03
    gates["sigma_ke_in_mev"] = sig_ke_in
    gates["sigma_ke_out_mev"] = sig_ke_out
    gates["rel_spread_in"] = sig_ke_in / ke_in_mev if ke_in_mev else float("nan")
    gates["rel_spread_out"] = sig_ke_out / ke_out_mev if ke_out_mev else float("nan")
    gates["emit_x_growth"] = enx_out / enx_in - 1.0 if enx_in else float("nan")
    gates["emit_y_growth"] = eny_out / eny_in - 1.0 if eny_in else float("nan")
    gates["mean_z_reached_m"] = mean_z_reached
    gates["mean_z_expected_m"] = z_expected
    gates["mean_z_ok"] = abs(mean_z_reached - z_expected) <= 0.05 * z_expected
    gates["min_ke_mev"] = ke_min_mev
    gates["beta_min"] = beta_min
    gates["beta_min_ok"] = beta_min > 0.999
    gates["transmission"] = transmission
    gates["max_envelope_m"] = max_env_m
    gates["min_bore_m"] = min_bore_m
    gates["envelope_in_bore"] = (max_env_m < min_bore_m) if max_env_m == max_env_m else None

    def mark(ok):
        return "PASS" if ok else "FAIL" if ok is not None else "n/a "

    print("\n── linac_rest validation gates (§5) ──")
    print(f"[{mark(gates['calib_within_3pct'])}] 1. per-section ΔE within ±3% of target")
    print(f"[{mark(gates['exit_ke_within_3pct'])}] 2. exit ⟨KE⟩ = {ke_out_mev:.1f} MeV "
          f"(expected {expected_out:.1f} = {ke_in_mev:.1f} in + {sum_de:.1f} Σ ΔE_target @ {p} MW)")
    print(f"[ -- ] 3. σ_KE  in {sig_ke_in:.2f} → out {sig_ke_out:.2f} MeV; "
          f"rel spread {gates['rel_spread_in']*100:.1f}% → {gates['rel_spread_out']*100:.1f}%")
    print(f"[ -- ] 4. εn growth  x {gates['emit_x_growth']*100:+.1f}%  "
          f"y {gates['emit_y_growth']*100:+.1f}%  (diagnostic; the quads-OFF ~2.4× is a fort.10N "
          f"εn artifact at bore/section crossings, σ_x smooth ⇒ not physical)")
    print(f"[{mark(gates['mean_z_ok'])}] 5. beam reached z = {mean_z_reached:.2f} m "
          f"(Σ lattice {z_expected:.2f} m — catches Ntstep truncation)")
    print(f"[{mark(gates['beta_min_ok'])}] 6. min captured KE {ke_min_mev:.1f} MeV ⇒ "
          f"β_min = {beta_min:.5f} (>0.999 ⇒ no-slip OK)")
    print(f"[ -- ] 7. transmission {transmission*100:.1f}% "
          f"(meaningful only with aperture scraping; else ~100% tautology)")
    print(f"[{mark(gates['envelope_in_bore'])}] 8. envelope-in-bore: 3σ_max = "
          f"{max_env_m*1e3:.2f} mm vs min bore {min_bore_m*1e3:.2f} mm "
          f"(soft — does the FODO contain the RMS envelope?)")

    # Hard gates the team locked as must-pass (concrete asserts for T6 to verify the run).
    if require_gates:
        if calib is not None:
            worst = max((abs(r["err_frac"]) for r in calib), default=0.0)
            assert gates["calib_within_3pct"], (
                f"gate 1: per-section ΔE off by up to {worst*100:.2f}% (>3%)")
        assert gates["exit_ke_within_3pct"], (
            f"gate 2: exit ⟨KE⟩ {ke_out_mev:.1f} MeV vs expected {expected_out:.1f} MeV "
            f"({(ke_out_mev/expected_out-1)*100:+.2f}%, >3%)")
        assert gates["mean_z_ok"], (
            f"gate 5: beam reached z={mean_z_reached:.2f} m, expected ~{z_expected:.2f} m "
            f"(Ntstep truncation? off by {(mean_z_reached/z_expected-1)*100:+.1f}%)")
        assert gates["beta_min_ok"], (
            f"gate 6: min captured KE {ke_min_mev:.1f} MeV ⇒ β_min={beta_min:.5f} (≤0.999, "
            f"no-slip assumption violated)")
    return gates
