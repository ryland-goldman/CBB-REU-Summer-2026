"""
SLAC Linac Section 1 in WarpX (RZ): capture the prebuncher beam in a 3 m
2π/3 traveling-wave accelerating structure (from ~137 keV; on-crest gain ~35 MeV at the
default 11 MW), with solenoid focusing and self-consistent space charge. At the
faithful-to-LinacSim 40 A / 11 MW default the weak focusing captures only ~2 %.

Fourth stage of the Cornell Linac chain in WarpX:
    cathode → gun → prebuncher → linac_sec1 (this).

The prebuncher's bunched beam is read from `prebuncher/diags/<case>/particles` at
its **tightest longitudinal focus** (min σ_z snapshot — the only point where the
beam both is bunched and still fits the 9.5 mm structure bore), translated to the
injection point, and tracked through:

  * a short injection drift,
  * a solenoid focusing channel (the `SOL_0`/`LENS_*` map × current `I_SOL`), and
  * the SLAC Section-1 traveling wave.

The traveling wave is synthesised exactly as the reference GPT model does — the sum
of two standing-wave maps 90° apart (`build_linac_sec1_field.py`):

    E(t) = scale·[ map1·cos(ωt+φ) + map2·cos(ωt+φ+π/2) ]
    B(t) = scale·[ map1·sin(ωt+φ) + map2·sin(ωt+φ+π/2) ]   (= forward traveling wave)

with f_RF = 2856 MHz (S-band) and scale = sqrt(POWER_MW / 1e-3) from the 1-kW field
normalisation. Each map is one `picmi.LoadAppliedField`; WarpX sums them (and the
static solenoid field) on the particles. The synchronous phase is undocumented, so
`PHASE_DEG` is an offset scanned to find the acceptance/crest (see run_pipeline.py).
"""

import os
import json
import shutil
import numpy as np
from pywarpx import picmi
from openpmd_viewer import OpenPMDTimeSeries

from pipeline._runner import run_step
from .build_linac_sec1_field import Z_STRUCT, RMAX, BORE_R, SOL_Z, V1KW_KEV   # keep field/phasing/domain in sync
from . import DEFAULT_OUTDIR

c = picmi.constants.c
m_e = picmi.constants.m_e
q_e = picmi.constants.q_e
MC2_KEV = m_e * c**2 / q_e / 1e3        # electron rest energy [keV] ≈ 511

# ── Field maps (built by build_linac_sec1_field.py) ───────────────────────────
RF1_FIELD = "linac_sec1/linac_sec1_field/linac_rf1.h5"
RF2_FIELD = "linac_sec1/linac_sec1_field/linac_rf2.h5"
SOL_FIELD = "linac_sec1/linac_sec1_field/linac_sol.h5"

F_RF = 2856.0e6                  # SLAC S-band [Hz] (Linac_RF in details.md)
RF_NORM_MW = 0.001               # field-map power normalisation (1 kW)

# ── Upstream input ────────────────────────────────────────────────────────────
PREBUNCH_DIAG = "prebuncher/diags/P8_zc/particles"   # bunched gun→prebuncher beam
Z_INJECT = 0.005                 # lab z where the bunch head is placed [m]
Z_FOCUS_MIN = 0.30               # only seek the bunch focus past the cavity (drift), so
                                 # the pre-modulation injection snapshots near z=0 (which
                                 # are tighter but un-bunched) are not mistaken for it [m]
MAX_PART = 50000                 # downsample the injected snapshot (reweighted)
RNG_SEED = 0

# ── Operating point (tunable via linac_sec1.config(...)) ──────────────────────
# Matched to the original LinacSim gpt_master.in section-1 GUI defaults:
# sec1_input_power = 11.0 MW, Sol 0 current = 40 A, relative phase phi_sec1_off = 0.
POWER_MW = 11.0                  # RF input power [MW] (sec1_input_power)
PHASE_DEG = 0.0                  # injection RF phase offset [deg] (= phi_sec1_off; crest found empirically)
# Solenoid current [A] (0 → focusing off), from the original Sol 0 current = 40 A.
# NOTE: at 40 A the focusing is far weaker than the I≈1000 A the rebuild previously
# used to reach ~97% capture into the 9.5 mm bore (capture was ~4% unfocused), so
# expect low capture at this faithful-to-original current. Set I_SOL via config to scan.
I_SOL = 40.0

# ── Performance / domain knobs ────────────────────────────────────────────────
CFL = 0.5                        # dt = CFL · Δz / v_inject
REQUIRED_PRECISION = 1e-4        # MLMG relative tolerance (space charge is a small perturbation)
MAX_ITERS = 200                  # MLMG iteration cap
MAX_STEPS = 0                    # 0 → auto-derive from transit; >0 → fixed
TRANSIT_MARGIN = 1.0             # transit already targets Z_END short of the wall
N_DIAGS = 60                     # number of openPMD dumps over the run

# ── Domain (RZ, single azimuthal mode — the maps are m = 0) ───────────────────
# 0.1 m injection drift + 3.0 m structure (entrance at Z_STRUCT=0.1 m, exit ≈3.12 m)
# + a field-free exit drift so the bunch coasts (not absorbed) at the last dump.
ZMAX = 3.50                      # [m]
# NR=16 over RMAX (12 mm) → dr=0.75 mm; with NZ=1664 (dz≈2.1 mm) the cells are ≈2.8:1,
# matching the prebuncher's long-thin box where MLMG converges. Coarser/anisotropic
# cells (e.g. NR=32) make the Poisson solve diverge ("MLMG failed").
NR, NZ = 16, 1664

OUTDIR = None                    # if None at main(), use DEFAULT_OUTDIR


def load_prebuncher_bunch():
    """Import the prebuncher beam at its min-σ_z focus and shift it to the entrance.

    Returns (dict for ParticleListDistribution, v_beam, mean KE [keV]).
    """
    ts = OpenPMDTimeSeries(PREBUNCH_DIAG)
    if len(ts.iterations) == 0:
        raise RuntimeError(
            f"{PREBUNCH_DIAG} has no iterations — did the prebuncher stage run?")

    # One light pass (z, w only) to locate the bunch focus: the smallest σ_z among
    # well-populated snapshots (a snapshot that has lost particles can have a
    # spuriously small σ_z, so require ≥ 80% of the peak count).
    # Floor scales with the requested statistics rather than an absolute 1000, so a
    # low-MAX_PART upstream run isn't silently rejected. At the default MAX_PART=50000
    # this is 1000 — identical to the original behaviour. (The "well-populated" gate is
    # the relative n ≥ 0.8·nmax below; this floor only drops near-empty dumps.)
    min_count = max(50, MAX_PART // 50)
    recs = []
    for it in ts.iterations:
        x, y, z, w = ts.get_particle(["x", "y", "z", "w"], species="electrons", iteration=it)
        if len(z) < min_count:
            continue
        zm = np.average(z, weights=w)
        sz = np.sqrt(np.average((z - zm) ** 2, weights=w))
        # Charge that fits the structure bore — the only charge the linac can capture; the
        # rest scrapes the iris/domain wall at injection no matter how strong the solenoid.
        # This is the raw-snapshot weight sum, used only to RANK snapshots for selection; the
        # sidecar's q_in_bore_C below is the post-downsample value in Coulombs for the chosen one.
        q_bore = float(w[np.hypot(x, y) <= BORE_R].sum())
        recs.append((it, len(z), zm, sz, q_bore))
    if not recs:
        raise RuntimeError(
            f"{PREBUNCH_DIAG}: no snapshot with ≥{min_count} macroparticles")
    nmax = max(n for _, n, _, _, _ in recs)
    cands = [(it, sz, qb) for it, n, zm, sz, qb in recs if n >= 0.8 * nmax and zm > Z_FOCUS_MIN]
    if not cands:                      # fallback: no snapshot past the gate
        cands = [(it, sz, qb) for it, n, zm, sz, qb in recs if n >= 0.8 * nmax]
    # Bore-aware focus: maximize the in-bore (capturable) charge, breaking ties toward the
    # tightest longitudinal focus (−σ_z). q_bore is a continuous float, so the −σ_z tie-break
    # only engages on a bit-identical in-bore charge; in practice the pick is just "max in-bore
    # charge". For the present low-power (8 kW) prebuncher beam — which never bunches and whose
    # σ_r grows monotonically over the drift — this lands on the earliest/least-expanded post-gate
    # snapshot, i.e. the most charge the linac can actually capture.
    #
    # CAVEAT: this optimizes the TRANSVERSE bore fit, not the LONGITUDINAL bunch. If the
    # prebuncher is driven into its bunching regime (≳160 kW), the min-σ_z waist forms LATE in
    # the drift (~1.26 m) where the beam is also radially over-expanded, so max-q_bore would
    # prefer an EARLY, debunched snapshot and silently discard the velocity-bunching. The shipped
    # 8 kW default has no such waist (the two criteria agree); the guard below warns if a future
    # operating point makes them disagree so the dropped bunch focus is not silent.
    it_focus = max(cands, key=lambda t: (t[2], -t[1]))[0]
    it_szmin = min(cands, key=lambda t: t[1])[0]        # the old pure min-σ_z bunching focus
    if it_szmin != it_focus:
        sz_focus = next(sz for it, sz, _ in cands if it == it_focus)
        sz_min = next(sz for it, sz, _ in cands if it == it_szmin)
        if sz_focus > 1.2 * sz_min:
            print(f"  WARNING: bore-aware focus (iter {it_focus}, σ_z {sz_focus*1e3:.2f} mm) is "
                  f"NOT the longitudinal bunch focus (iter {it_szmin}, σ_z {sz_min*1e3:.2f} mm) — "
                  f"max-in-bore-charge selection is discarding a {sz_focus/sz_min:.1f}× tighter "
                  f"bunch. Expected only if the prebuncher is in its bunching regime; the injected "
                  f"beam will be radially contained but longitudinally debunched.", flush=True)

    x, y, z, ux, uy, uz, w = ts.get_particle(
        ["x", "y", "z", "ux", "uy", "uz", "w"], species="electrons", iteration=it_focus)
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
    r = np.hypot(x, y)
    rmax = float(r.max())
    # Injection-loss bookkeeping. WarpX silently drops particles initialised at r > RMAX
    # before the first diagnostic dump, so the first dump's charge is already the post-scrape
    # (in-domain) charge. Capture must therefore be reported against the TRUE injected charge
    # recorded here — otherwise the dominant loss (radial scraping at t=0) is invisible.
    q_inj = float(w.sum()) * q_e                         # total injected [C]
    q_dom = float(w[r <= RMAX].sum()) * q_e              # within the 12 mm domain (survives step 0)
    q_bore = float(w[r <= BORE_R].sum()) * q_e           # within the 9.55 mm RF bore (capturable)
    inj = dict(it_focus=int(it_focus), n_injected=int(z.size),
               q_injected_C=q_inj, q_in_domain_C=q_dom, q_in_bore_C=q_bore,
               z_inject_mean_m=float(np.average(z, weights=w)),
               rmax_m=rmax, sigma_z_m=sz, ke_mean_keV=ke_mean)
    print(f"Injected {z.size} macroparticles from prebuncher focus (iter {it_focus}); "
          f"⟨KE⟩ {ke_mean:.1f} keV, σ_z {sz*1e3:.2f} mm, r_max {rmax*1e3:.2f} mm, "
          f"v_beam {v_beam:.3e} m/s, q {q_inj*1e9:.4f} nC", flush=True)
    print(f"  injection radial fit: {q_dom/q_inj*100:.1f}% of charge within RMAX={RMAX*1e3:.0f} mm "
          f"({q_dom*1e12:.1f} pC enters the domain), {q_bore/q_inj*100:.1f}% within the "
          f"{BORE_R*1e3:.2f} mm bore — the rest scrapes the wall at injection.", flush=True)
    if q_dom < 0.9 * q_inj:
        print(f"  WARNING: {(1-q_dom/q_inj)*100:.0f}% of the injected beam starts outside the "
              f"radial domain and is lost at step 0; reported capture is vs the {q_inj*1e9:.4f} nC "
              f"injected, not the {q_dom*1e9:.4f} nC that enters.", flush=True)
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

    bunch, v_beam, ke_mean, inj = load_prebuncher_bunch()
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
          f"I_sol={I_SOL:g} A, f_RF={F_RF/1e6:.0f} MHz → {outdir}/", flush=True)

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

    # ── Applied fields: two quadrature RF maps (+ static solenoid), summed ─────
    # Constants are baked into the AMReX parser strings (LoadAppliedField forwards
    # extra kwargs to picmistandard, which rejects them). The picmi layer registers
    # each as a uniquely named ext_fieldN and WarpX sums them on the particles.
    applied = []
    # IMPORTANT ordering: the picmi LoadAppliedField sets the *global*
    # `E_ext_particle_init_style="none"` for any field with load_E=False, and the
    # field initialised LAST wins. So the solenoid (B only, load_E=False) must be
    # added BEFORE the RF maps — otherwise it disables the accelerating E field for
    # the whole run. With an RF map (load_E=True) last, E stays "read_from_file".
    if I_SOL != 0.0:
        applied.append(picmi.LoadAppliedField(
            read_fields_from_path=SOL_FIELD, load_E=False, load_B=True,
            warpx_B_time_function=f"{I_SOL:.8e}"))     # static: B = per-A map × I_SOL
    applied += [
        picmi.LoadAppliedField(
            read_fields_from_path=RF1_FIELD, load_E=True, load_B=True,
            warpx_E_time_function=f"{scale:.8e}*cos({omega:.10e}*t + ({phi:.8e}))",
            warpx_B_time_function=f"{scale:.8e}*sin({omega:.10e}*t + ({phi:.8e}))"),
        picmi.LoadAppliedField(
            read_fields_from_path=RF2_FIELD, load_E=True, load_B=True,
            warpx_E_time_function=f"{scale:.8e}*cos({omega:.10e}*t + ({phi2:.8e}))",
            warpx_B_time_function=f"{scale:.8e}*sin({omega:.10e}*t + ({phi2:.8e}))"),
    ]
    # Enforce the ordering invariant the comment above relies on: picmi sets the *global*
    # E_ext_particle_init_style from the LAST-added field, so the last entry MUST load E
    # (an RF map) — else a B-only field (the solenoid) silently disables the accelerating
    # E and the beam just coasts at ~137 keV with no error. Guard it so a future reorder
    # fails loudly instead of producing a wrong, silent run.
    assert getattr(applied[-1], "load_E", False), (
        "last applied field must have load_E=True (an RF map), or the global E_ext "
        "style is forced to 'none' and the beam coasts unaccelerated")

    electrons = picmi.Species(
        particle_type="electron", name="electrons",
        initial_distribution=picmi.ParticleListDistribution(
            x=bunch["x"], y=bunch["y"], z=bunch["z"],
            ux=bunch["ux"], uy=bunch["uy"], uz=bunch["uz"], weight=bunch["w"]))

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
