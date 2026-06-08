"""
Cornell Linac **Sections 2–8** in Impact-T (the ``linac_rest`` stage): accelerate the
captured ~25 MeV beam from ``linac_sec1`` through seven S-band travelling-wave sections to
≈308 MeV at the default 11 MW klystron point (307.97 survivors through the real bore).

Fifth and final stage of the Cornell Linac chain:
    cathode → gun → injector → linac_sec1 → linac_rest (this).

Unlike the four WarpX stages, this is an **Impact-T** run driven through lume-impact
(external serial ``ImpactTexe``), so it runs IN-PROCESS via ``ImpactStage`` (no pywarpx
global-geometry binding). The deck is assembled in-memory by ``build_linac_rest_lattice``
(seven TW sections reusing the ``rfdata4–7`` field shape, SC off, quads K1=0 for the
headline), each section's field scale calibrated to its ΔE_target by
``calibration.calibrate_sections``.

Handoff IN (``linac_sec1`` → ParticleGroup):
  * read ``linac_sec1/diags/main/particles`` LAST (exit) dump — the captured coasting
    beam (``pipeline.impact_io.read_warpx_dump``);
  * keep only the **captured core** (KE ≥ ``MIN_KE_MEV``): the sec-1 exit dump carries a
    low-energy un-captured tail (~16% of charge below 10 MeV) whose β < 0.999 would break
    the rigid-crest no-slip assumption across the 36 m line and is not part of the beam the
    downstream optics transport. The discarded tail is recorded in the summary;
  * ``drift_to_t()`` then zero z (Impact-T injects at z == 0, t-coordinates);
  * set ``I.initial_particles`` (carries the surviving charge — no renormalisation).

Handoff OUT (Impact-T → openPMD for plot_chain + summary):
  * ``P_out = I.particles["final_particles"]`` → ``write_openpmd_particles`` writes the exact
    WarpX openPMD layout (species ``"electrons"``, ``ux=γβ``, ``w``=count) into
    ``diags/main/particles/`` so ``plot_chain``/``_beam_summary``/``build_moment_table`` read
    it with the same ``get_particle([...], species="electrons")`` call as every WarpX stage;
  * Impact-T's only particle dumps (with write_beam off) are the initial and final beams,
    so a 2-iteration series (injected core + ≈308 MeV exit) is written for plot_chain; the
    per-section vs-z evolution lives in ``injection_summary.json`` ``stat_vs_z`` (I.stat), not
    in particle slices;
  * ``injection_summary.json`` records ``q_injected_C`` (the captured-core charge read from
    sec-1, so ``_beam_summary`` reports ~100% within-stage and the chain capture narrative
    stays coherent) and ``z_inject_lab_m`` (the lab-z the beam was injected at — Impact-T
    output z is local-frame; the cross-stage z0 shift lives in ``plot_chain``).
"""

import json
import os
import shutil

import numpy as np

from pipeline.impact_io import read_warpx_dump, write_openpmd_particles
from . import build_linac_rest_lattice as L
from . import calibration as cal
from . import DEFAULT_OUTDIR

MC2_MEV = 0.51099895069          # electron rest energy [MeV]

# ── Upstream input ────────────────────────────────────────────────────────────
# linac_sec1's exit dump (its captured ~25 MeV coasting beam). The reader picks the
# last/exit iteration by default — the captured beam, NOT a handoff-plane nearest-⟨z⟩ pick
# (mirrors linac_sec1's _exit_row rows[-1] for the linac exit).
LINAC_SEC1_DIAG = "linac_sec1/diags/main/particles"
LINAC_SEC1_SUMMARY = "linac_sec1/diags/main/injection_summary.json"  # for the true-injected denom
# z the sec-1 captured beam exits at (its exit-dump ⟨z⟩ ≈ 3.31 m lab). Recorded into the
# linac_rest summary as z_inject_lab_m so plot_chain's _apply_linac_rest_z0 can place the
# local-frame Impact-T z in the lab without a hard-coded offset.
Z_INJECT_LAB_M = None            # measured from the read-in dump at run time (None ⇒ measure)

# Captured-core energy cut [MeV]: drop the low-energy un-captured tail so the injected beam
# is the captured core (β > 0.999 ⇒ rigid-crest no-slip holds). β > 0.999 needs KE > 10.92 MeV
# (γ > 22.37); 12 MeV (β ≈ 0.99917) keeps a margin and retains ~88% of the sec-1 exit charge.
# The discarded tail is the un-captured halo, not part of the transported beam. Configurable.
MIN_KE_MEV = 12.0

# ── Operating point (tunable via linac_rest.config(...)) ──────────────────────
POWER_MW = 11.0                  # RF input power [MW] per section — ONE convention (sec-1 point)
PHASE_DEG = 0.0                  # on-crest synchronous phase (β > 0.999 ⇒ no per-section rephasing)

# ── Performance / deck knobs ──────────────────────────────────────────────────
Np = 4000                        # macroparticles tracked in the final I.run() (downsample core to this)
Np_calib = 400                   # decimated bunch for per-section calibration (⟨KE⟩/crest are mean
                                 # quantities, so a coarse bunch fits the scale+phase fine and the
                                 # ~240 calibration tracks run ~10× faster); 0 ⇒ calibrate on full Np
Ntstep = 200000                  # Impact-T step cap (sized for ~36 m at Dt≈2e-12; mean_z asserted)
Dt = 2.0e-12                     # time step [s]
Nxyz = 16                        # SC mesh per axis (unused — SC off; power of 2)
DRIFT_M = None                   # inter-section drift override [m] (None ⇒ build default 0.4)
QUADS_ON = False                 # headline: quads OFF (K1 = 0). True ⇒ exploratory FODO.
QUAD_K = None                    # per-section quad b1_gradient [T/m] (exploratory; None ⇒ zeros)
RNG_SEED = 0
REQUIRE_GATES = True             # T6: assert the hard §5 gates so a bad run fails loudly

OUTDIR = None                    # if None at main(), use DEFAULT_OUTDIR


def _sec1_lab_z0():
    """linac_sec1's local→lab z offset, mirroring plot_chain._apply_linac_z0.

    sec-1's diagnostics are in its own local frame; plot_chain shifts them to lab by
    z0 = z_handoff_m − z_inject_mean_m (recorded in sec-1's injection_summary.json). We add the
    same offset to the sec-1 dump ⟨z⟩ so linac_rest's recorded z_inject_lab_m is in the lab
    frame and the segments abut without overlap. Falls back to Z_HANDOFF (2.03 m) if the
    sec-1 summary lacks the fields (older runs), matching plot_chain's own fallback.
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
    P = read_warpx_dump(LINAC_SEC1_DIAG)             # species "electron", t-coords, last dump
    n_all = P.n_particle
    # LAB-frame injection z: the sec-1 dump's ⟨z⟩ is in sec-1's LOCAL frame; plot_chain shifts
    # sec-1 to lab by z0 = z_handoff_m − z_inject_mean_m (from sec-1's injection_summary). Add
    # that same offset so z_inject_lab_m is the LAB z where the sec-1 beam exits — i.e. exactly
    # where linac_rest abuts sec-1 (else plot_chain's _apply_linac_rest_z0 would place the
    # segment ~1.9 m too early, overlapping sec-1's exit-drift dumps).
    z_local = float(P["mean_z"])
    z_inject_lab = (z_local + _sec1_lab_z0()) if Z_INJECT_LAB_M is None else Z_INJECT_LAB_M
    q_exit = float(P["charge"])                      # all sec-1 exit charge (pre-core-cut)

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

    # Impact-T injects at a common time with z == 0: drift every particle to the mean t,
    # then translate z so the bunch sits at z == 0.
    Pc.drift_to_t(Pc["mean_t"])
    Pc.z = Pc.z - Pc["mean_z"]

    ke_in = float(Pc["mean_energy"] / 1e6 - MC2_MEV)
    ke_min = float(Pc.energy.min() / 1e6 - MC2_MEV)
    beta_min = cal._beta_from_ke_mev(ke_min)
    info = dict(
        n_sec1_exit=int(n_all), n_core=int(Pc.n_particle),
        q_sec1_exit_C=q_exit, q_core_C=q_core,
        # Explicit model-validity-cut bookkeeping (per physicist's policy): the FULL sec-1
        # charge is the honest denominator; the dropped sub-MIN_KE_MEV tail is counted as loss.
        q_after_cut_C=q_core,                        # survivors of the KE≥MIN_KE_MEV cut
        q_dropped_lowKE_C=(q_exit - q_core),         # the dropped low-energy (β<0.999) tail
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
    """Downsample Impact-T's continuous z-stats into a small along-z table for the plots.

    The per-section evolution (⟨KE⟩, σ_KE, σ_x, ε_n) comes from Impact-T's I.stat(...) arrays
    (thousands of points over the line — write_beam dumps are off, see build_impact), thinned
    to ~`n` evenly-spaced samples and stored in injection_summary.json so plot_linac_rest can
    render vs-z without re-running. σ_KE uses σ_gamma·mc2 (sigma_energy is not a stat key).
    """
    zc = I.stat("mean_z")
    if len(zc) == 0:
        return {}
    idx = np.unique(np.linspace(0, len(zc) - 1, min(n, len(zc))).astype(int))
    out = {"z_m": zc[idx].tolist(),
           "ke_mev": (I.stat("mean_kinetic_energy")[idx] / 1e6).tolist(),
           "sigma_ke_mev": (I.stat("sigma_gamma")[idx] * MC2_MEV).tolist(),
           "sigma_x_m": I.stat("sigma_x")[idx].tolist(),
           "norm_emit_x": I.stat("norm_emit_x")[idx].tolist(),
           "norm_emit_y": I.stat("norm_emit_y")[idx].tolist()}
    return out


def _write_outputs(I, outdir, inj):
    """Write the Impact-T result as WarpX-layout openPMD + the injection summary JSON.

    Writes the exit beam (final_particles) — and any other surviving ParticleGroups — as
    openPMD slices for the cross-stage handoff (plot_chain sorts by ⟨z⟩, tolerates a single
    exit dump). The per-section vs-z evolution lives in inj["stat_vs_z"] (Impact-T I.stat),
    not in particle slices.

    The output groups' charges were already re-imposed in main() — each scaled to
    q_core × (group n / n_in) so the openPMD `weighting` sums to the physically transmitted
    charge (Impact-T with Bcurr=0 returns a default 1 C normalisation; transmission was
    measured from macro count BEFORE that rescale, so the rescale can't mask aperture loss).
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

    # ── Build the chained 7-section deck ──────────────────────────────────────
    I, total_len = L.build_impact(
        power_mw=POWER_MW, phase_deg=PHASE_DEG, drift_m=DRIFT_M,
        np_particles=P_in.n_particle, dt=Dt, ntstep=Ntstep, nxyz=Nxyz,
        quads_on=QUADS_ON, quad_k=QUAD_K)
    I.initial_particles = P_in
    I.configure()
    print(f"Deck: {L.N_SECTIONS} TW sections, Σ {total_len:.2f} m, P={POWER_MW:g} MW, "
          f"on-crest θ₀={PHASE_DEG:g}°, SC off, quads {'ON' if QUADS_ON else 'OFF (K1=0)'} "
          f"→ {outdir}/", flush=True)

    # ── Per-section scale calibration (Task 5) ────────────────────────────────
    print("Calibrating per-section field scale to ΔE_target (on-crest, scale-only)…",
          flush=True)
    calib = cal.calibrate_sections(I, P_in, power_mw=POWER_MW, np_calib=Np_calib)

    # ── Full run ──────────────────────────────────────────────────────────────
    print(f"Running Impact-T ({L.N_SECTIONS} sections, Ntstep={Ntstep})…", flush=True)
    I.run()
    if not I.finished or I.error:
        raise RuntimeError(f"Impact-T did not finish cleanly (finished={I.finished}, "
                           f"error={I.error})")

    # ── Transmission from MACRO COUNT (BEFORE re-imposing charge) ─────────────
    # Measure transmission as n_out/n_in on the macroparticle COUNT, while the output still
    # carries Impact-T's own weights — this is the only honest transmission. (Computing it
    # from charge AFTER re-imposing q_core below would force 1.0, masking any aperture scrape:
    # the reviewer's bug.) With Bcurr=0 Impact-T carries no charge, but the per-macro WEIGHT
    # is uniform, so n_out/n_in is the true surviving fraction whether or not the bore aperture
    # scrapes. q_out is then q_core × that fraction — the physically transmitted charge.
    P_out = I.particles["final_particles"]
    n_in = int(P_in.n_particle)
    n_out = int(P_out.n_particle)
    transmission = (n_out / n_in) if n_in else 0.0
    q_core = float(P_in["charge"])
    q_out = q_core * transmission                       # physically transmitted core charge

    # ── Re-impose the physical beam charge for the openPMD output (SC OFF loses it) ──
    # Impact-T returns the output with a default 1 C normalisation; rescale each output group's
    # per-macro weights so the group sums to its physically transmitted charge (q_core × the
    # group's own surviving fraction vs n_in). This is for the openPMD `weighting` record only;
    # the transmission number above was already measured from counts, so this can't mask loss.
    for _name, _pg in I.particles.items():
        if _pg is not None and _pg.n_particle > 0:
            _pg.charge = q_core * (_pg.n_particle / n_in)

    # ── Validation gates (§5) — assert the hard gates (Task 6) ────────────────
    gates = cal.validate_run(I, P_in, power_mw=POWER_MW, calib=calib,
                             require_gates=REQUIRE_GATES)

    # ── Handoff OUT: openPMD + summary ────────────────────────────────────────
    # The true end-to-end denominator (cathode-injected) lives in the upstream summaries;
    # record sec-1's captured charge for the chain capture narrative.
    sec1_true_injected = None
    if os.path.exists(LINAC_SEC1_SUMMARY):
        try:
            sec1_true_injected = json.load(open(LINAC_SEC1_SUMMARY)).get("q_injected_C")
        except Exception:
            sec1_true_injected = None

    inj = dict(
        # HONEST capture denominator (team ruling): the FULL linac_sec1 captured charge that
        # arrived at the handoff, NOT the post-cut core. So _beam_summary's q_out/q_injected
        # counts the dropped low-energy tail (the MIN_KE_MEV model-validity cut) AND the in-run
        # aperture loss as within-stage loss — normalizing to the core charge would overstate it.
        q_injected_C=core_info["q_sec1_exit_C"],     # full sec-1 exit charge (honest denominator)
        q_core_injected_C=core_info["q_core_C"],     # of that, the captured core actually tracked
        z_inject_lab_m=core_info["z_inject_lab_m"],  # lab-z of injection (Impact-T z is local)
        z_inject_local_m=0.0,                        # Impact-T local z at injection (beam zeroed)
        total_lattice_length_m=float(total_len),
        power_mw=float(POWER_MW), phase_deg=float(PHASE_DEG),
        quads_on=bool(QUADS_ON),
        # Aperture provenance: which aperture the transmission was measured against. The real
        # tapered bore is ON for the headline (bore_aperture_on True or quads_on); xyrad_m is the
        # containment-box half-width (kept just above the bore so the bore is the binding aperture).
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
