"""
SLAC Linac Section 1 in WarpX (RZ): capture the injector beam in a 3 m
2π/3 traveling-wave accelerating structure (on-crest gain ~35 MeV at the default
11 MW), with self-consistent space charge. Transverse focusing is NOT here — it is
applied upstream in the injector stage by the real lenses at their true lab z, which
hand the linac a beam already focused and collimated to the 9.547 mm iris.

Fourth stage of the Cornell Linac chain in WarpX:
    cathode → gun → injector → linac_sec1 (this).

The injector's focused, velocity-bunched beam is read from `injector/diags/main/particles`
at the **z ≈ 2.03 m handoff plane** (the dump whose ⟨z⟩ is nearest Z_HANDOFF — see
load_injector_bunch; NOT the old min-σ_z / max-in-bore selector), translated to the
injection point, and tracked through:

  * a short injection drift (collimation is the multi-plane iris scrape applied in
    load_injector_bunch — only the particles that pass the 9.547 mm iris/pipe from
    z=1.922 m onward are injected; the halo removed is what the real aperture scrapes), and
  * the SLAC Section-1 traveling wave.

The traveling wave is synthesised exactly as the reference GPT model does — the sum
of two standing-wave maps 90° apart (`build_linac_sec1_field.py`):

    E(t) = scale·[ map1·cos(ωt+φ) + map2·cos(ωt+φ+π/2) ]
    B(t) = scale·[ map1·sin(ωt+φ) + map2·sin(ωt+φ+π/2) ]   (= forward traveling wave)

with f_RF = 2856 MHz (S-band) and scale = sqrt(POWER_MW / 1e-3) from the 1-kW field
normalisation. Each map is one `picmi.LoadAppliedField`; WarpX sums them on the
particles. The synchronous phase is undocumented, so `PHASE_DEG` is an offset scanned
to find the acceptance/crest (see run_pipeline.py).
"""

import os
import json
import shutil
import numpy as np
from pywarpx import picmi
from openpmd_viewer import OpenPMDTimeSeries

from pipeline._runner import run_step
from pipeline.collimator import pipe_violator_ids, survivor_mask   # multi-plane iris scrape
from .build_linac_sec1_field import Z_STRUCT, RMAX, BORE_R, V1KW_KEV   # keep field/phasing/domain in sync
from . import DEFAULT_OUTDIR

c = picmi.constants.c
m_e = picmi.constants.m_e
q_e = picmi.constants.q_e
MC2_KEV = m_e * c**2 / q_e / 1e3        # electron rest energy [keV] ≈ 511

# ── Field maps (built by build_linac_sec1_field.py) ───────────────────────────
# Only the two SLAC quadrature RF maps — the in-linac solenoid was removed in the
# injector upgrade (focusing now lives in the injector at the lenses' true lab z).
RF1_FIELD = "linac_sec1/linac_sec1_field/linac_rf1.h5"
RF2_FIELD = "linac_sec1/linac_sec1_field/linac_rf2.h5"

F_RF = 2856.0e6                  # SLAC S-band [Hz] (Linac_RF in details.md)
RF_NORM_MW = 0.001               # field-map power normalisation (1 kW)

# ── Upstream input ────────────────────────────────────────────────────────────
# The injector hands off a focused, velocity-bunched beam already collimated to the
# 9.547 mm iris at the z ≈ 2.03 m handoff plane (Z_acc_1). The selector below picks the
# snapshot whose ⟨z⟩ is nearest that plane (NOT min-σ_z / max-in-bore — those would land
# on an off-plane snapshot now that a real longitudinal waist forms near the handoff).
INJECTOR_DIAG = "injector/diags/main/particles"      # focused, bunched injector beam
Z_HANDOFF = 2.03                 # [m] injector→linac handoff plane (Z_acc_1)
COLLIM_Z = 1.922                 # [m] iris start (LinacSim scatteriris); the 9.547 mm pipe runs
                                 #     COLLIM_Z→wall. RMAX (imported) is the iris radius. The
                                 #     injector→linac collimation is a multi-plane id scrape from this
                                 #     plane (pipeline.collimator), NOT a single cut at the 2.03 m
                                 #     handoff — the beam CONVERGES through the tail, so a 2.03 m cut
                                 #     would pass halo the real 1.922 m iris scrapes.
Z_INJECT = 0.005                 # lab z where the bunch head is placed [m]
MAX_PART = 50000                 # downsample the injected snapshot (reweighted)
RNG_SEED = 0

# ── Operating point (tunable via linac_sec1.config(...)) ──────────────────────
# Matched to the original LinacSim gpt_master.in section-1 GUI defaults:
# sec1_input_power = 11.0 MW, relative phase phi_sec1_off = 0. (Transverse focusing is
# upstream in the injector now; the linac carries no solenoid.)
POWER_MW = 11.0                  # RF input power [MW] (sec1_input_power)
PHASE_DEG = 0.0                  # injection RF phase offset [deg] (= phi_sec1_off; crest found empirically)

# ── Performance / domain knobs ────────────────────────────────────────────────
CFL = 0.5                        # dt = CFL · Δz / v_inject
REQUIRED_PRECISION = 1e-4        # MLMG relative tolerance (space charge is a small perturbation)
SPACE_CHARGE = True              # beam self-field (space charge) on/off. False →
                                 # warpx_do_not_deposit: only the applied RF maps act. Keep TRUE:
                                 # the self-field is LARGEST at the ~220 keV injection handoff (the
                                 # beam is captured to relativistic energy over the stage, exiting
                                 # ~25 MeV), where the lab-frame ES solver even overstates it ~γ², so
                                 # SC-off is a real change at the low-energy front — diagnostic-only.
MAX_ITERS = 200                  # MLMG iteration cap
MAX_STEPS = 0                    # 0 → auto-derive from transit; >0 → fixed
TRANSIT_MARGIN = 1.0             # transit already targets Z_END short of the wall
N_DIAGS = 60                     # number of openPMD dumps over the run

# ── Domain (RZ, single azimuthal mode — the maps are m = 0) ───────────────────
# 0.1 m injection drift + 3.0 m structure (entrance at Z_STRUCT=0.1 m, exit ≈3.12 m)
# + a field-free exit drift so the bunch coasts (not absorbed) at the last dump.
ZMAX = 3.50                      # [m]
# NR=16 over RMAX (9.547 mm) → dr=0.597 mm; with NZ=1664 (dz≈2.1 mm) the cells are
# ≈3.5:1. The long-thin box needs cells within the ≈3:1 rule for MLMG to converge;
# 3.5:1 still converges here (the self field is a small perturbation on the captured
# beam). If it diverges, raise NR (keeping it ÷ blocking factor) rather than coarsening
# NZ. (RMAX dropped from 12 mm to the 9.547 mm bore/iris in the injector upgrade.)
NR, NZ = 16, 1664

OUTDIR = None                    # if None at main(), use DEFAULT_OUTDIR


def load_injector_bunch():
    """Import the injector beam at the z ≈ 2.03 m handoff plane and shift it to entry.

    SELECTOR (changed in the injector upgrade): pick the diagnostic dump whose bunch
    ⟨z⟩ is closest to Z_HANDOFF = 2.03 m (Z_acc_1) — NOT min-σ_z and NOT max-in-bore
    charge. With the two-cavity + solenoid injector the bunch now forms a real
    longitudinal waist near the handoff, so the old max-q_bore selector would land on an
    early, debunched snapshot and silently discard the bunching; and min-σ_z would pick
    the upstream waist (~1.45 m), not the handoff. The injector places a dump within
    ~1 mm of 2.03 m (fine-cadence handoff diagnostic), so nearest-⟨z⟩ lands on the plane.

    Returns (dict for ParticleListDistribution, v_beam, mean KE [keV], inj summary).
    """
    ts = OpenPMDTimeSeries(INJECTOR_DIAG)
    if len(ts.iterations) == 0:
        raise RuntimeError(
            f"{INJECTOR_DIAG} has no iterations — did the injector stage run?")

    # One light pass to find the well-populated dump nearest the handoff plane. The
    # min-count floor scales with the requested statistics (drops only near-empty
    # boundary dumps); the relative n ≥ 0.8·nmax gate keeps the selection on a
    # representative snapshot, not a depleted late dump that happens to sit near 2.03 m.
    min_count = max(50, MAX_PART // 50)
    recs = []
    for it in ts.iterations:
        z, w = ts.get_particle(["z", "w"], species="electrons", iteration=it)
        if len(z) < min_count:
            continue
        recs.append((it, len(z), float(np.average(z, weights=w))))
    if not recs:
        raise RuntimeError(
            f"{INJECTOR_DIAG}: no snapshot with ≥{min_count} macroparticles")
    nmax = max(n for _, n, _ in recs)
    cands = [(it, zm) for it, n, zm in recs if n >= 0.8 * nmax]
    # Nearest ⟨z⟩ to the handoff plane.
    it_handoff, zm_handoff = min(cands, key=lambda t: abs(t[1] - Z_HANDOFF))
    if abs(zm_handoff - Z_HANDOFF) > 0.02:
        print(f"  WARNING: nearest injector dump to the {Z_HANDOFF*1e3:.0f} mm handoff is at "
              f"⟨z⟩={zm_handoff*1e3:.0f} mm ({abs(zm_handoff-Z_HANDOFF)*1e3:.0f} mm off) — "
              f"the handoff diagnostic may be too coarse; the injected beam is off-plane.",
              flush=True)

    idh, x, y, z, ux, uy, uz, w = ts.get_particle(
        ["id", "x", "y", "z", "ux", "uy", "uz", "w"], species="electrons", iteration=it_handoff)
    # TRUE injected charge = the full beam delivered to the handoff (all r) — the honest
    # capture denominator. Recorded BEFORE the iris scrape below.
    q_inj = float(w.sum()) * q_e

    # ── Iris/pipe collimation: multi-plane id scrape (the physical injector→linac cut) ──
    # The real 9.547 mm scatteriris is at COLLIM_Z=1.922 m with a 9.547 mm pipe to the wall.
    # The Sol 0 / Lens 0E telescope focuses the beam HARD across the 1.922→2.03 m tail (it
    # CONVERGES — measured in-iris ~38 % @1.92 m → ~93 % @2.03 m), so a single radial cut at
    # the 2.03 m handoff (what this used to do) would pass halo the real iris scrapes and
    # overstate transmission. Instead scrape by tracking particle ids across every dump that
    # samples the pipe: a particle outside the iris at ANY plane z ≥ COLLIM_Z hit the wall.
    # See pipeline.collimator. Only the survivors are injected into the linac.
    scan_iters = [it for it, _n, zm in recs if (COLLIM_Z - 0.05) <= zm <= (Z_HANDOFF + 0.03)]
    violators = pipe_violator_ids(ts, scan_iters, RMAX, COLLIM_Z)
    keep = survivor_mask(idh, violators)
    idh, x, y, z, ux, uy, uz, w = (a[keep] for a in (idh, x, y, z, ux, uy, uz, w))
    r = np.hypot(x, y)
    q_dom = float(w.sum()) * q_e                         # iris survivors carried to the handoff (in-iris)
    q_bore = float(w[r <= BORE_R].sum()) * q_e           # of those, within the 9.55 mm RF bore

    # Downsample the SURVIVORS (reweighted to preserve their charge) for the run.
    if z.size > MAX_PART:
        rng = np.random.default_rng(RNG_SEED)
        sel = rng.choice(z.size, MAX_PART, replace=False)
        scale_w = z.size / MAX_PART        # reweight to preserve total charge
        x, y, z, ux, uy, uz, w = (a[sel] for a in (x, y, z, ux, uy, uz, w))
        w = w * scale_w
    # Translate so the bunch tail (smallest z; the beam travels +z, so the head leads
    # at larger z) sits at Z_INJECT.
    z = z - z.min() + Z_INJECT

    gb = np.sqrt(1.0 + ux**2 + uy**2 + uz**2)          # γ
    v_beam = float(np.average(uz / gb, weights=w) * c)
    ke_mean = float(np.average(gb - 1.0, weights=w) * MC2_KEV)
    sz = float(np.sqrt(np.average((z - np.average(z, weights=w)) ** 2, weights=w)))
    rmax = float(np.hypot(x, y).max())
    inj = dict(it_handoff=int(it_handoff), z_handoff_m=float(zm_handoff),
               n_injected=int(z.size), q_injected_C=q_inj, q_in_domain_C=q_dom,
               q_in_bore_C=q_bore, z_inject_mean_m=float(np.average(z, weights=w)),
               rmax_m=rmax, sigma_z_m=sz, ke_mean_keV=ke_mean)
    print(f"Injected {z.size} macroparticles (iris survivors) from injector handoff "
          f"(iter {it_handoff}, ⟨z⟩={zm_handoff*1e3:.0f} mm); ⟨KE⟩ {ke_mean:.1f} keV, "
          f"σ_z {sz*1e3:.2f} mm, r_max {rmax*1e3:.2f} mm, v_beam {v_beam:.3e} m/s, "
          f"true-injected q {q_inj*1e9:.4f} nC", flush=True)
    print(f"  multi-plane iris scrape at the {RMAX*1e3:.3f} mm aperture (z≥{COLLIM_Z*1e3:.0f} mm, "
          f"{len(scan_iters)} planes): {q_dom/q_inj*100:.1f}% of the handoff charge survives the "
          f"pipe ({q_dom*1e12:.1f} pC into the bore), {q_bore/q_inj*100:.1f}% within the "
          f"{BORE_R*1e3:.2f} mm RF bore — the rest hit the iris/pipe wall (the real aperture).",
          flush=True)
    # openPMD ux/uy/uz are γβ; PICMI wants proper velocity γβc in m/s → ×c.
    return (dict(x=x, y=y, z=z, ux=ux * c, uy=uy * c, uz=uz * c, w=w),
            v_beam, ke_mean, inj)


def main():
    outdir = OUTDIR or DEFAULT_OUTDIR
    # Fresh diags: WarpX appends one openPMD file per dump, so a rerun of the same
    # case (e.g. a scan point) would otherwise mix old and new iterations. The diags
    # are git-ignored and regenerated, so clearing the case dir is safe.
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    # Compute omega in main() (not at module import) so a config(F_RF=...) override,
    # which lands after import, is honoured rather than frozen at the import-time value.
    omega = 2.0 * np.pi * F_RF

    bunch, v_beam, ke_mean, inj = load_injector_bunch()
    z_center = float(np.average(bunch["z"], weights=bunch["w"]))
    # Persist the true injected charge (+ in-domain/in-bore breakdown) so the plotter reports
    # capture against what was actually injected, not the post-scrape first-dump charge. WarpX
    # creates `outdir` when it writes the first diagnostic; create it up front for this sidecar.
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "injection_summary.json"), "w") as fh:
        json.dump(inj, fh, indent=2)

    # ── RF amplitude + phase ──────────────────────────────────────────────────
    scale = float(np.sqrt(POWER_MW / RF_NORM_MW))
    # Phase referenced to the bunch-centre arrival at the structure entrance; the
    # absolute crest is undocumented (capture from β≈0.63 into a β_phase=1 wave), so
    # PHASE_DEG is scanned and the headline run uses the scan maximum.
    t_in = (Z_STRUCT - z_center) / v_beam
    phi = -omega * t_in + np.deg2rad(PHASE_DEG)
    phi2 = phi + np.pi / 2.0          # field2 quadrature offset
    print(f"Case: P={POWER_MW:g} MW (scale={scale:.1f}), phase_off={PHASE_DEG:g}°, "
          f"f_RF={F_RF/1e6:.0f} MHz → {outdir}/  (focusing is upstream in the injector)",
          flush=True)

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
        maximum_iterations=MAX_ITERS, warpx_self_fields_verbosity=0)

    # ── Applied fields: the two quadrature RF maps, summed ────────────────────
    # Constants are baked into the AMReX parser strings (LoadAppliedField forwards
    # extra kwargs to picmistandard, which rejects them). The picmi layer registers
    # each as a uniquely named ext_fieldN and WarpX sums them on the particles. The
    # in-linac solenoid was removed (focusing is upstream in the injector), so the
    # linac's only applied fields are the two RF maps — both load_E=True.
    applied = [
        picmi.LoadAppliedField(
            read_fields_from_path=RF1_FIELD, load_E=True, load_B=True,
            warpx_E_time_function=f"{scale:.8e}*cos({omega:.10e}*t + ({phi:.8e}))",
            warpx_B_time_function=f"{scale:.8e}*sin({omega:.10e}*t + ({phi:.8e}))"),
        picmi.LoadAppliedField(
            read_fields_from_path=RF2_FIELD, load_E=True, load_B=True,
            warpx_E_time_function=f"{scale:.8e}*cos({omega:.10e}*t + ({phi2:.8e}))",
            warpx_B_time_function=f"{scale:.8e}*sin({omega:.10e}*t + ({phi2:.8e}))"),
    ]
    # Guard the E-init-style invariant: the LAST applied field must load E (an RF map),
    # else the global E_ext_particle_init_style is forced to "none" and the beam coasts
    # unaccelerated. Unconditional — there is always an RF map last (no B-only field now).
    assert getattr(applied[-1], "load_E", False), (
        "last applied field must have load_E=True (an RF map), or the global E_ext "
        "style is forced to 'none' and the beam coasts unaccelerated")

    electrons = picmi.Species(
        particle_type="electron", name="electrons",
        initial_distribution=picmi.ParticleListDistribution(
            x=bunch["x"], y=bunch["y"], z=bunch["z"],
            ux=bunch["ux"], uy=bunch["uy"], uz=bunch["uz"], weight=bunch["w"]),
        warpx_do_not_deposit=not SPACE_CHARGE)   # SPACE_CHARGE=False → no beam self-field

    # ── Time step / duration ──────────────────────────────────────────────────
    # dt is sized at the (slowest) injection velocity. The beam accelerates to β→1,
    # so the transit is estimated in segments: injection drift + a short capture
    # length (β: in→≈c) + the relativistic remainder. We target a stop plane Z_END a
    # little short of ZMAX so the bunch finishes in the field-free exit drift, NOT at
    # the absorbing z-wall — once the domain empties the MLMG self-field solve aborts.
    # All scan phases share this length (gain_keV uses the on-crest max), so the
    # fastest (crest) case is the binding one for the no-abort margin.
    dz = ZMAX / NZ
    dt = CFL * dz / v_beam
    beta_in = v_beam / c
    gain_keV = scale * V1KW_KEV                     # ≈ on-crest gain [keV] (1-kW V1kW from the maps)
    gamma_hi = 1.0 + (ke_mean + gain_keV) / MC2_KEV
    beta_hi = float(np.sqrt(1.0 - 1.0 / gamma_hi**2))
    z_end = ZMAX - 0.20                            # stop in the drift, clear of the wall
    L_cap = 0.40                                   # capture length over which β: in→hi
    beta_cap = 0.5 * (beta_in + beta_hi)
    transit = ((Z_STRUCT - z_center) / v_beam
               + min(L_cap, max(0.0, z_end - Z_STRUCT)) / (beta_cap * c)
               + max(0.0, z_end - Z_STRUCT - L_cap) / (beta_hi * c))
    n_steps = MAX_STEPS or int(TRANSIT_MARGIN * transit / dt)
    print(f"dt = {dt:.3e} s, max_steps = {n_steps}, RF period {1/F_RF*1e9:.3f} ns "
          f"({1/F_RF/dt:.0f} steps/period), β_in={beta_in:.3f}, "
          f"est. on-crest gain {gain_keV/1e3:.1f} MeV", flush=True)

    # ── Diagnostics (openPMD, HDF5) ───────────────────────────────────────────
    period = max(1, n_steps // N_DIAGS)
    part_diag = picmi.ParticleDiagnostic(
        name="particles", period=period, species=[electrons],
        data_list=["position", "momentum", "weighting"],
        write_dir=outdir, warpx_format="openpmd", warpx_openpmd_backend="h5")

    sim = picmi.Simulation(
        solver=solver, max_steps=n_steps, time_step_size=dt,
        verbose=0, particle_shape="linear")
    sim.add_species(
        electrons,
        layout=picmi.PseudoRandomLayout(n_macroparticles_per_cell=1, grid=grid))
    for fld in applied:
        sim.add_applied_field(fld)
    sim.add_diagnostic(part_diag)

    print(f"\nRunning {n_steps} steps (diag every {period}) → {outdir}/")
    run_step(sim, n_steps, desc="linac_sec1")
    print("\nDone.")


if __name__ == "__main__":
    main()
