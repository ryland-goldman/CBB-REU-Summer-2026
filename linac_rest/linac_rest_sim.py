"""
linac_rest stage (Cornell Linac sections 2–8) main():  handoff IN (captured core from
linac_sec1) → build/calibrate the in-memory Impact-T deck → I.run() → §5 gates →
openPMD handoff OUT + injection_summary.json.

See linac_rest/README.md for physics, operating point, calibration, and gotchas.
"""

import inspect
import json
import os
import shutil
import threading

import numpy as np

from pipeline.impact_io import read_warpx_dump, write_openpmd_particles
from pipeline._impact_runner import terminal_progress
from . import build_linac_rest_lattice as L
from . import calibration as cal
from . import DEFAULT_OUTDIR

MC2_MEV = 0.51099895069          # electron rest energy [MeV]

# Read from the helper's signature default so the recorded μ can't drift from what
# fodo_quad_gradients() actually uses (the sim calls it with the default).
_FODO_PHASE_ADV_DEG = inspect.signature(L.fodo_quad_gradients).parameters["phase_adv_deg"].default

# ── Upstream input ────────────────────────────────────────────────────────────
LINAC_SEC1_DIAG = "linac_sec1/diags/main/particles"  # last/exit dump (the captured ~25 MeV beam)
LINAC_SEC1_SUMMARY = "linac_sec1/diags/main/injection_summary.json"  # for the true-injected denom
Z_INJECT_LAB_M = None            # lab-z of sec-1 exit (None ⇒ measure from the dump at run time)

MIN_KE_MEV = 12.0                # captured-core energy cut [MeV]; β≈0.99917 ⇒ rigid-crest no-slip holds

# ── Operating point (tunable via linac_rest.config(...)) ──────────────────────
POWER_MW = 11.0                  # RF input power [MW] per section — ONE convention (sec-1 point)
PHASE_DEG = 0.0                  # on-crest synchronous phase (β > 0.999 ⇒ no per-section rephasing)

# ── Performance / deck knobs ──────────────────────────────────────────────────
Np = 4000                        # macroparticles tracked in the final I.run() (downsample core to this)
Np_calib = 400                   # decimated bunch for per-section calibration (0 ⇒ full Np)
Ntstep = 200000                  # Impact-T step cap (sized for ~36 m at Dt≈2e-12; mean_z asserted)
Dt = 2.0e-12                     # time step [s]
Nxyz = 16                        # SC mesh per axis (power of 2; used only when SPACE_CHARGE)
SPACE_CHARGE = False             # beam self-field. False (headline) ⇒ Bcurr=0; True ⇒ Bcurr=q·Bfreq
                                 # (exploratory/unvalidated; calibration stays SC-free, SC on run only)
DRIFT_M = None                   # inter-section drift override [m] (None ⇒ build default 0.4)
QUADS_ON = False                 # headline: quads OFF (K1 = 0). True ⇒ exploratory FODO.
QUAD_K = None                    # per-section quad b1_gradient [T/m] (exploratory; None ⇒ zeros)
RNG_SEED = 0
REQUIRE_GATES = True             # assert the hard §5 gates so a bad run fails loudly

OUTDIR = None                    # if None at main(), use DEFAULT_OUTDIR


def _sec1_lab_z0():
    """linac_sec1's local→lab z offset (z_handoff_m − z_inject_mean_m), mirroring
    plot_chain._apply_linac_z0. Falls back to 2.03 m if the sec-1 summary lacks the fields.
    """
    try:
        with open(LINAC_SEC1_SUMMARY) as fh:
            s = json.load(fh)
        return float(s["z_handoff_m"]) - float(s["z_inject_mean_m"])
    except Exception:
        return 2.03


def load_sec1_core():
    """Read linac_sec1's exit beam, keep the captured core, return (ParticleGroup, info).

    Drops the low-energy un-captured tail (KE < MIN_KE_MEV), downsamples to Np (reweighted to
    preserve the surviving charge), ``drift_to_t()`` + zeroes z for Impact-T injection. The
    ParticleGroup carries the captured-core charge (no renormalisation).
    """
    P = read_warpx_dump(LINAC_SEC1_DIAG)             # species "electron" (singular), t-coords, last dump
    n_all = P.n_particle
    # sec-1 dump ⟨z⟩ is in sec-1's LOCAL frame; add its local→lab offset so z_inject_lab is the
    # LAB z where linac_rest abuts sec-1 (else plot_chain places the segment ~1.9 m too early).
    z_local = float(P["mean_z"])
    z_inject_lab = (z_local + _sec1_lab_z0()) if Z_INJECT_LAB_M is None else Z_INJECT_LAB_M
    q_exit = float(P["charge"])                      # all sec-1 exit charge (pre-core-cut, honest denom)

    ke_mev = (P.energy - MC2_MEV * 1e6) / 1e6
    core = ke_mev >= MIN_KE_MEV
    if core.sum() < 50:
        raise RuntimeError(
            f"only {int(core.sum())} sec-1 particles above MIN_KE_MEV={MIN_KE_MEV} MeV — "
            f"capture cut too aggressive or upstream beam not accelerated")
    Pc = P[core]
    q_core = float(Pc.charge)

    # Downsample the core to Np (reweighted to preserve the core charge).
    if Pc.n_particle > Np:
        rng = np.random.default_rng(RNG_SEED)
        sel = rng.choice(Pc.n_particle, Np, replace=False)
        Pc = Pc[sel]
        Pc.weight = Pc.weight * (q_core / float(Pc.charge))   # restore total core charge

    # Impact-T injects at a common time with z == 0: drift to mean t, then translate z to 0.
    Pc.drift_to_t(Pc["mean_t"])
    Pc.z = Pc.z - Pc["mean_z"]

    ke_in = float(Pc["mean_energy"] / 1e6 - MC2_MEV)
    ke_min = float(Pc.energy.min() / 1e6 - MC2_MEV)
    beta_min = cal._beta_from_ke_mev(ke_min)
    info = dict(
        n_sec1_exit=int(n_all), n_core=int(Pc.n_particle),
        q_sec1_exit_C=q_exit, q_core_C=q_core,
        q_after_cut_C=q_core,                        # survivors of the KE≥MIN_KE_MEV cut
        q_dropped_lowKE_C=(q_exit - q_core),         # the dropped low-energy (β<0.999) tail (counted as loss)
        core_charge_frac=(q_core / q_exit if q_exit else 0.0),
        min_ke_mev_cut=float(MIN_KE_MEV), ke_in_mev=ke_in,
        ke_min_core_mev=ke_min, beta_min_core=beta_min,
        z_inject_lab_m=z_inject_lab,
    )
    print(f"sec-1 exit: {n_all} parts, {q_exit*1e12:.1f} pC; captured core (KE≥{MIN_KE_MEV} MeV): "
          f"{Pc.n_particle} parts, {q_core*1e12:.1f} pC ({info['core_charge_frac']*100:.1f}% of "
          f"exit charge). ⟨KE⟩_in {ke_in:.2f} MeV, min-core KE {ke_min:.2f} MeV "
          f"(β_min={beta_min:.5f}), inject lab-z {z_inject_lab:.3f} m", flush=True)
    if beta_min <= 0.999:
        print(f"  WARNING: min-core β {beta_min:.5f} ≤ 0.999 — raise MIN_KE_MEV to keep the "
              f"rigid-crest no-slip assumption (plan §5 gate 6).", flush=True)
    return Pc, info


def _stat_vs_z(I, n=200):
    """Thin Impact-T's I.stat(...) z-arrays to ~`n` samples for the vs-z plots (write_beam
    dumps are off). σ_KE uses σ_gamma·mc2 (sigma_energy is not a stat key).
    """
    zc = I.stat("mean_z")
    if len(zc) == 0:
        return {}
    idx = np.unique(np.linspace(0, len(zc) - 1, min(n, len(zc))).astype(int))
    out = {"z_m": zc[idx].tolist(),
           "ke_mev": (I.stat("mean_kinetic_energy")[idx] / 1e6).tolist(),
           "sigma_ke_mev": (I.stat("sigma_gamma")[idx] * MC2_MEV).tolist(),
           "sigma_x_m": I.stat("sigma_x")[idx].tolist(),
           "sigma_y_m": I.stat("sigma_y")[idx].tolist(),
           "norm_emit_x": I.stat("norm_emit_x")[idx].tolist(),
           "norm_emit_y": I.stat("norm_emit_y")[idx].tolist()}
    return out


def _watch_fort18(I, bar, stop_evt, poll=0.3):
    """Drive the run progress bar from Impact-T's ``fort.18`` column 2 (reference ``mean_z`` [m]),
    polled live from the run workdir. Best-effort: a partial/missing line just skips a tick.
    """
    f18 = None
    while not stop_evt.is_set():
        try:
            if f18 is None and I.path:
                cand = os.path.join(I.path, "fort.18")
                if os.path.exists(cand):
                    f18 = cand
            if f18:
                with open(f18, "rb") as fh:
                    tail = fh.read().rsplit(b"\n", 2)
                line = next((ln for ln in reversed(tail) if ln.strip()), b"")
                z = float(line.split()[1])
                if bar.total and z > bar.n:
                    bar.update(min(z, bar.total) - bar.n)
        except Exception:
            pass
        stop_evt.wait(poll)


def _run_with_progress(I, total_len):
    """Run the full Impact-T deck with a fort.18-driven z-progress bar (0 → total_len [m])."""
    with terminal_progress(total=round(float(total_len), 2),
                           desc="linac_rest: track", unit="m") as rbar:
        stop_evt = threading.Event()
        watcher = threading.Thread(target=_watch_fort18, args=(I, rbar, stop_evt),
                                   daemon=True)
        watcher.start()
        try:
            I.run()
        finally:
            stop_evt.set()
            watcher.join(timeout=1.0)


def _write_outputs(I, outdir, inj):
    """Write the surviving ParticleGroups as WarpX-layout openPMD slices (sorted by ⟨z⟩) plus
    injection_summary.json. Group charges were already re-imposed in main().
    """
    part_dir = os.path.join(outdir, "particles")
    os.makedirs(part_dir, exist_ok=True)

    slices = []
    for name, pg in I.particles.items():
        if pg is None or pg.n_particle < 50:
            continue
        slices.append((float(pg["mean_z"]), pg))
    slices.sort(key=lambda t: t[0])
    if not slices:
        raise RuntimeError("Impact-T produced no usable particle groups")

    for it, (_zc, pg) in enumerate(slices):
        write_openpmd_particles(pg, part_dir, iteration=it, time=float(pg["mean_t"]))

    with open(os.path.join(outdir, "injection_summary.json"), "w") as fh:
        json.dump(inj, fh, indent=2)


def main():
    outdir = OUTDIR or DEFAULT_OUTDIR
    # Fresh diags (regenerated, git-ignored): clear so a rerun doesn't mix old iterations.
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir, exist_ok=True)

    # ── Handoff IN: captured core from linac_sec1 ─────────────────────────────
    P_in, core_info = load_sec1_core()

    # ── Build the CALIBRATION deck — ALWAYS quads-OFF, SC-free (bcurr=0) ──────
    # The fit's mean_energy is transverse-independent only on-axis; a quads-OFF, zero-aperture,
    # SC-free deck keeps the fit clean and identical to the headline. `bcurr` (Bcurr=|q|·Bfreq ⇒
    # Q=|q| on the SC mesh) applies only to the FINAL run deck below.
    bcurr = abs(P_in.charge) * L.RF_FREQ_HZ if SPACE_CHARGE else 0.0
    I_cal, total_len = L.build_impact(
        power_mw=POWER_MW, phase_deg=PHASE_DEG, drift_m=DRIFT_M,
        np_particles=P_in.n_particle, dt=Dt, ntstep=Ntstep, nxyz=Nxyz,
        quads_on=False, quad_k=None, bcurr=0.0)
    I_cal.initial_particles = P_in
    I_cal.configure()
    print(f"Deck: {L.N_SECTIONS} TW sections, Σ {total_len:.2f} m, P={POWER_MW:g} MW, "
          f"on-crest θ₀={PHASE_DEG:g}°, SC {'on (Bcurr=%.4g A)' % bcurr if SPACE_CHARGE else 'off'}, "
          f"quads {'ON' if QUADS_ON else 'OFF (K1=0)'} → {outdir}/", flush=True)

    # ── Per-section scale calibration — on the quads-OFF deck ─────────────────
    print("Calibrating per-section field scale to ΔE_target (on-crest, scale-only)…",
          flush=True)
    with terminal_progress(total=L.N_SECTIONS, desc="linac_rest: calibrate",
                           unit="sec") as cbar:
        calib = cal.calibrate_sections(I_cal, P_in, power_mw=POWER_MW, np_calib=Np_calib,
                                       bar=cbar)
    # `calib` carries per-section {scale, crest_phase_deg} — both re-applied to the run deck.
    calibrated_scales = [r["scale"] for r in calib]

    # ── Assemble the RUN deck ─────────────────────────────────────────────────
    applied_quad_k = [0.0] * L.N_SECTIONS    # placed b1_gradient [T/m] (all-zero quads-OFF)
    if QUADS_ON:
        # Build a FRESH quads-ON deck (calibrated scales baked in, FODO gradients + quad/drift
        # bore radius active); the measured ke_in energy-scales the FODO.
        quad_k = (QUAD_K if QUAD_K is not None
                  else L.fodo_quad_gradients(ke_in_mev=core_info["ke_in_mev"]))
        applied_quad_k = list(quad_k)        # the per-quad T/m actually placed on the run deck
        I, total_len = L.build_impact(
            power_mw=POWER_MW, phase_deg=PHASE_DEG, drift_m=DRIFT_M,
            np_particles=P_in.n_particle, dt=Dt, ntstep=Ntstep, nxyz=Nxyz,
            quads_on=True, quad_k=quad_k, scales=calibrated_scales, bcurr=bcurr)
        # Re-apply calibration via the ControlGroup, not just the build-time element scales: the
        # rf_field_scale group is absolute=True defaulting 0, so adding it + configure() would
        # overwrite the baked-in scales with 0 (silent no-acceleration). Set the group scale AND
        # the absolute crest phase per section (θ₀ is ABSOLUTE), then configure once.
        for r in calib:
            gname = cal._ensure_section_group(I, r["index"])
            cal._set_section_phase(I, r["index"], r["crest_phase_deg"])
            cal._set_group_scale(I, gname, r["scale"])
        I.initial_particles = P_in
        I.configure()
        # Guard the FODO was actually placed (catches a silent quads-off build); quad3a is the
        # lead half of the first H/V doublet.
        assert I.ele["quad3a"]["b1_gradient"] != 0.0, (
            "QUADS_ON run deck has zero b1_gradient on quad3a — FODO gradients did not apply")
        print(f"Run deck: quads ON, H/V-doublet lead-pole b1_gradient [T/m] = "
              f"{[round(quad_k[i], 4) for i in range(L.N_SECTIONS - 1)]} "
              f"(placeholder optics — guessed K1, A→T undocumented, H/V doublet (±g halves), "
              f"nominal μ={_FODO_PHASE_ADV_DEG:g}°)", flush=True)
    else:
        # Headline: the run deck IS the quads-OFF calibrated deck. SC-ON re-imposes the run Bcurr
        # on the SC-free-calibrated deck and re-configures (calibration stayed SC-free).
        I = I_cal
        if bcurr:
            I.header["Bcurr"] = bcurr
            I.configure()

    # ── Full run ──────────────────────────────────────────────────────────────
    print(f"Running Impact-T ({L.N_SECTIONS} sections, Ntstep={Ntstep})…", flush=True)
    _run_with_progress(I, total_len)
    if not I.finished or I.error:
        raise RuntimeError(f"Impact-T did not finish cleanly (finished={I.finished}, "
                           f"error={I.error})")

    # ── Transmission from MACRO COUNT, measured BEFORE re-imposing charge ──────
    # n_out/n_in on the macro count (uniform per-macro weight) is the only honest transmission;
    # computing it from charge AFTER the re-impose below would force 1.0 and mask aperture loss.
    P_out = I.particles["final_particles"]
    n_in = int(P_in.n_particle)
    n_out = int(P_out.n_particle)
    transmission = (n_out / n_in) if n_in else 0.0
    q_core = float(P_in["charge"])
    q_out = q_core * transmission                       # physically transmitted core charge

    # ── Re-impose the physical charge for the openPMD `weighting` (SC-OFF loses it) ──
    # Impact-T returns a default 1 C normalisation; rescale each group to q_core × (group n /
    # n_in). Output-only — transmission was already measured from counts above.
    for _name, _pg in I.particles.items():
        if _pg is not None and _pg.n_particle > 0:
            _pg.charge = q_core * (_pg.n_particle / n_in)

    # ── Validation gates (§5) — assert the hard gates ─────────────────────────
    gates = cal.validate_run(I, P_in, power_mw=POWER_MW, calib=calib,
                             require_gates=REQUIRE_GATES)

    # ── Handoff OUT: openPMD + summary ────────────────────────────────────────
    # Record sec-1's q_injected_C (the injector→linac handoff charge) for the chain capture narrative.
    sec1_true_injected = None
    if os.path.exists(LINAC_SEC1_SUMMARY):
        try:
            sec1_true_injected = json.load(open(LINAC_SEC1_SUMMARY)).get("q_injected_C")
        except Exception:
            sec1_true_injected = None

    inj = dict(
        # HONEST capture denominator: the FULL sec-1 captured charge at the handoff, NOT the
        # post-cut core — so _beam_summary's q_out/q_injected counts the dropped tail + in-run loss.
        q_injected_C=core_info["q_sec1_exit_C"],     # full sec-1 exit charge (honest denominator)
        q_core_injected_C=core_info["q_core_C"],     # of that, the captured core actually tracked
        z_inject_lab_m=core_info["z_inject_lab_m"],  # lab-z of injection (Impact-T z is local)
        z_inject_local_m=0.0,                        # Impact-T local z at injection (beam zeroed)
        total_lattice_length_m=float(total_len),
        power_mw=float(POWER_MW), phase_deg=float(PHASE_DEG),
        quads_on=bool(QUADS_ON),
        # Per-quad b1_gradient [T/m] ACTUALLY PLACED (all-zero quads-OFF); length-N_SECTIONS but
        # only [0..N_SECTIONS-2] (Q2–Q7) placed — the Q8 trailing entry is never installed.
        quad_k=[float(k) for k in applied_quad_k],
        quad_phase_adv_deg=float(_FODO_PHASE_ADV_DEG),
        # Which aperture the transmission was measured against; xyrad_m is the containment-box
        # half-width (kept just above the bore so the bore is the binding aperture).
        bore_aperture_on=bool(L.BORE_APERTURE_ON or QUADS_ON),
        xyrad_m=float(L.XYRAD_M),
        ke_in_mev=core_info["ke_in_mev"],
        ke_out_mev=float(gates["exit_ke_mev"]),
        expected_ke_out_mev=float(gates["expected_exit_ke_mev"]),
        mean_z_reached_m=float(gates["mean_z_reached_m"]),
        beta_min_core=core_info["beta_min_core"],
        # Transmission from MACRO COUNT (n_out/n_in), measured before re-imposing charge.
        n_core_in=n_in, n_out=n_out,
        transmission_core=transmission,              # n_out/n_in — honest (1.0 only if no scrape)
        # Soft envelope-in-bore gate: 3σ_max RMS envelope vs narrowest bore [m] + PASS/FAIL.
        # Never gated (guessed K1); can legitimately FAIL quads-OFF (no-focusing beam diverges).
        max_envelope_m=float(gates["max_envelope_m"]),
        min_bore_m=float(gates["min_bore_m"]),
        envelope_in_bore=gates["envelope_in_bore"],
        q_out_C=q_out,                               # q_core × transmission (physically transmitted)
        core_charge_frac_of_sec1_exit=core_info["core_charge_frac"],
        q_after_cut_C=core_info["q_after_cut_C"],
        q_dropped_lowKE_C=core_info["q_dropped_lowKE_C"],
        sec1_true_injected_C=sec1_true_injected,
        n_sec1_exit=core_info["n_sec1_exit"], n_core=core_info["n_core"],
        min_ke_mev_cut=core_info["min_ke_mev_cut"],
        calibration=[{k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                      for k, v in r.items()} for r in calib],
        stat_vs_z=_stat_vs_z(I),
    )
    _write_outputs(I, outdir, inj)
    print(f"\nDone. Exit ⟨KE⟩ {gates['exit_ke_mev']:.1f} MeV "
          f"(expected {gates['expected_exit_ke_mev']:.1f}); beam reached "
          f"{gates['mean_z_reached_m']:.2f}/{total_len:.2f} m. → {outdir}/", flush=True)


if __name__ == "__main__":
    main()
