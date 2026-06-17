"""Per-section field-scale calibration + §5 validation gates for the ``linac_rest`` Impact-T stage.

``calibrate_sections`` fits each TW section's ``rf_field_scale`` to its ΔE target on crest;
``validate_run`` computes/prints the §5 gates after ``I.run()``. Import-only (no side effects).

mean_energy metric is TOTAL energy in eV (rest mass included), NOT kinetic.
See linac_rest/README.md for physics, the energy budget, conventions, and gotchas.
"""

import math

from scipy.optimize import brentq

from . import build_linac_rest_lattice as L

ELECTRON_REST_MEV = 0.51099895069  # electron rest energy [MeV]


# Delegate naming to build_linac_rest_lattice (single source of truth for deck element names).
def section_ele_names(index):
    """(entrance, body_1, body_2, exit) solrf element names for section ``index`` (0-based)."""
    return tuple(L.section_group_names(index))


def section_group_name(index):
    return f"sec{index + 2}_scale"


def _ensure_section_group(I, index):
    """Create (idempotently) the rf_field_scale ControlGroup over a section's 4 solrf cells.

    factors=[1, 1/sin(β₀d), 1/sin(β₀d), 1], absolute=True ⇒ group value S sets entrance/exit=S,
    body=S/sin(β₀d), preserving the template body ratio. Returns the group name.
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

    Impact.__setitem__ needs the "name:attribute" form (attribute = the group var_name); a bare
    I[gname]=value raises (splits on ':'), so always go through this helper.
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
    """Find the base phase [deg] that puts section ``index`` on crest (max ΔE) at the bunch arrival.

    Impact-T theta0_deg is ABSOLUTE (t=0 reference) ⇒ each downstream section crests at a different
    base phase. seed_phase=None ⇒ full 0–360° coarse scan; seed given ⇒ tight local scan around it.
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

    Per section i (lab order), sequential: fresh Ic=I.copy() with upstream fitted scales+phases,
    this + downstream sections zeroed; track from z=0 to entrance; find crest base phase; brentq
    the rf_field_scale group to entrance energy + ΔE_target; store and apply onto the live deck I.

    A FRESH copy per section is required: a reused Impact copy re-links its ControlGroups to
    rebuilt ele dicts after each track(), so an early-added group goes stale and later scale
    writes silently no-op. ``np_calib`` decimates the bunch (mean ⟨KE⟩/crest are count-insensitive)
    and crest seeding (sections ≥1 from prev crest + arrival-time shift) keep cost off the runtime.

    Returns per-section dicts and leaves ``I`` fully calibrated (ready for ``I.run()``). Optional
    tqdm ``bar`` ticks once per section; summary lines route through ``bar.write`` (else ``print``).
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
        # Fresh copy so the ControlGroups link to THIS copy's current ele dicts (no stale link).
        Ic = I.copy()
        Ic.verbose = False
        Ic.configure()
        gnames = [_ensure_section_group(Ic, j) for j in range(L.N_SECTIONS)]
        for j in range(L.N_SECTIONS):
            if j < i:                                   # upstream: re-apply fitted scale + phase
                sc, ph = fitted[j]
                _set_section_phase(Ic, j, ph)
                _set_group_scale(Ic, gnames[j], sc)
            else:                                       # this section + downstream: off
                _set_section_phase(Ic, j, L.PHASE_DEG)
                _set_group_scale(Ic, gnames[j], 0.0)

        names = set(section_ele_names(i))
        s0, s1 = ele_bounds([e for e in Ic.lattice if e.get("name") in names])
        gname = gnames[i]

        # GUARD: never track with s ≤ 0 — Impact-T treats stop=0.0 as "no early stop" and runs the
        # ENTIRE lattice (section 2 has s0=0, and load z-zeroing leaves mean_z ≈ -1e-15).
        if s0 > 1e-9 and P_cal["mean_z"] < s0 - 1e-9:
            P_entrance = Ic.track(P_cal, s=s0)
        else:
            P_entrance = P_cal
        e_entrance = P_entrance["mean_energy"]              # total energy [eV]
        ke_entrance_mev = (e_entrance / 1e6) - ELECTRON_REST_MEV

        de_target_mev = L.section_de_target(i, p)
        target_total = e_entrance + de_target_mev * 1e6     # total energy target [eV]

        # Crest base phase at the bunch arrival (θ₀ is absolute). Seed sections ≥1 from prev crest
        # + arrival-time shift (−ω·Δt_arr) for a local refine; section 0 does the full coarse scan.
        t_arr = P_entrance["mean_t"]
        seed = None
        if prev_crest is not None:
            shift_deg = -math.degrees(omega * (t_arr - prev_t_arr))
            seed = (prev_crest + shift_deg) % 360.0
        crest_phase = _find_crest_phase(Ic, i, gname, P_entrance, s1, probe_scale,
                                        seed_phase=seed)
        _set_section_phase(Ic, i, crest_phase)
        prev_crest, prev_t_arr = crest_phase, t_arr

        # At crest, calibrate the field scale to hit ΔE_target. Guard the bracket: no energy change
        # across [lo, hi] means the group link is broken — fail loudly with context.
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


# §5 validation gates. (The MIN_KE_MEV model-validity cut lives in linac_rest_sim.load_sec1_core,
# not here; calibration only consumes the already-cut captured-core beam.)
def _beta_from_ke_mev(ke_mev):
    g = 1.0 + ke_mev / ELECTRON_REST_MEV
    return math.sqrt(max(0.0, 1.0 - 1.0 / (g * g)))


def validate_run(I, P_in, power_mw=None, calib=None, require_gates=False):
    """Compute + print the §5 validation gates after a full ``I.run()``. See README -> "Validation
    gates (§5)" for the gate list and rationale.

    Returns a dict of computed values + a boolean per gate, printing a PASS/FAIL line each.
    ``require_gates=True`` raises AssertionError on any failed HARD gate (per-section ΔE ±3%,
    exit ⟨KE⟩ ±3%, beam reached final zedge, min-captured-KE β > 0.999); soft/diagnostic gates
    (σ_KE, εn, transmission, envelope) stay print-only.
    """
    p = L.POWER_MW if power_mw is None else power_mw
    P_out = I.particles["final_particles"]

    ke_in_mev = (P_in["mean_energy"] / 1e6) - ELECTRON_REST_MEV
    ke_out_mev = (P_out["mean_energy"] / 1e6) - ELECTRON_REST_MEV
    sum_de = sum(L.section_de_target(i, p) for i in range(L.N_SECTIONS))
    expected_out = ke_in_mev + sum_de

    # σ_KE (absolute, MeV); ParticleGroup energy is total, so σ of KE == σ of total energy.
    sig_ke_in = P_in["sigma_energy"] / 1e6
    sig_ke_out = P_out["sigma_energy"] / 1e6

    enx_in, eny_in = P_in["norm_emit_x"], P_in["norm_emit_y"]
    enx_out, eny_out = P_out["norm_emit_x"], P_out["norm_emit_y"]

    mean_z_reached = I.stat("mean_z")[-1]
    z_expected = L.total_lattice_length_m()

    # min captured KE (over surviving particles) → β.
    ke_min_mev = (P_out["energy"].min() / 1e6) - ELECTRON_REST_MEV
    beta_min = _beta_from_ke_mev(ke_min_mev)

    # Envelope-in-bore soft gate: 3σ_max envelope vs the narrowest section bore. Soft/print-only
    # (K1 is the guessed placeholder FODO ⇒ must never gate the energy headline).
    BEAM_EDGE_SIGMA = 3.0
    sx = I.stat("sigma_x")
    sy = I.stat("sigma_y")
    max_env_m = BEAM_EDGE_SIGMA * float(max(sx.max(), sy.max())) if len(sx) else float("nan")
    min_bore_m = min(min(L.section_bore_radii(i)) for i in range(L.N_SECTIONS))

    # Transmission from the macroparticle COUNT, not charge — correct regardless of when main()
    # re-imposes q_out = q_core·(n_out/n_in) (a charge ratio would carry that ordering dependency).
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

    # Hard gates (must-pass) asserted only when require_gates.
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
