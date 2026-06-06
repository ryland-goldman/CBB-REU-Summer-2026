"""
SLAC Linac Section 1 in WarpX (RZ): capture the bunched prebuncher beam in a 3 m
2π/3 traveling-wave accelerating structure and take it from ~148 keV to ~37 MeV,
with solenoid focusing and self-consistent space charge.

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
import shutil
import numpy as np
from pywarpx import picmi
from openpmd_viewer import OpenPMDTimeSeries

from pipeline._runner import run_step
from .build_linac_sec1_field import Z_STRUCT, RMAX, SOL_Z   # keep field/phasing/domain in sync
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
OMEGA = 2.0 * np.pi * F_RF
RF_NORM_MW = 0.001               # field-map power normalisation (1 kW)

# ── Upstream input ────────────────────────────────────────────────────────────
PREBUNCH_DIAG = "prebuncher/diags/P800_zc/particles"   # bunched gun→prebuncher beam
Z_INJECT = 0.005                 # lab z where the bunch head is placed [m]
Z_FOCUS_MIN = 0.30               # only seek the bunch focus past the cavity (drift), so
                                 # the pre-modulation injection snapshots near z=0 (which
                                 # are tighter but un-bunched) are not mistaken for it [m]
MAX_PART = 50000                 # downsample the injected snapshot (reweighted)
RNG_SEED = 0

# ── Operating point (tunable via linac_sec1.config(...)) ──────────────────────
POWER_MW = 15.0                  # RF input power [MW]  (~37 MeV on crest)
PHASE_DEG = 0.0                  # injection RF phase offset [deg] (scanned for crest)
# Solenoid current [A] (0 → focusing off). The prebuncher beam arrives diverging
# (it lacked focusing over its 1.3 m drift), so strong focusing is needed to keep it
# in the bore: capture rises from ~4% (I=0) to ~95% (I≈1000 A → 0.15 T peak), so the
# default is the strongly focused, on-crest operating point. Set I_SOL=0 (via config)
# to run the unfocused case.
I_SOL = 1000.0

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
    recs = []
    for it in ts.iterations:
        z, w = ts.get_particle(["z", "w"], species="electrons", iteration=it)
        if len(z) < 1000:
            continue
        zm = np.average(z, weights=w)
        sz = np.sqrt(np.average((z - zm) ** 2, weights=w))
        recs.append((it, len(z), zm, sz))
    if not recs:
        raise RuntimeError(f"{PREBUNCH_DIAG}: no snapshot with >1000 macroparticles")
    nmax = max(n for _, n, _, _ in recs)
    cands = [(it, sz) for it, n, zm, sz in recs if n >= 0.8 * nmax and zm > Z_FOCUS_MIN]
    if not cands:                      # fallback: no snapshot past the gate
        cands = [(it, sz) for it, n, zm, sz in recs if n >= 0.8 * nmax]
    it_focus, _ = min(cands, key=lambda t: t[1])

    x, y, z, ux, uy, uz, w = ts.get_particle(
        ["x", "y", "z", "ux", "uy", "uz", "w"], species="electrons", iteration=it_focus)
    if z.size > MAX_PART:
        rng = np.random.default_rng(RNG_SEED)
        sel = rng.choice(z.size, MAX_PART, replace=False)
        scale_w = z.size / MAX_PART        # reweight to preserve total charge
        x, y, z, ux, uy, uz, w = (a[sel] for a in (x, y, z, ux, uy, uz, w))
        w = w * scale_w
    # Translate so the bunch head (smallest z) sits at Z_INJECT.
    z = z - z.min() + Z_INJECT

    gb = np.sqrt(1.0 + ux**2 + uy**2 + uz**2)          # γ
    v_beam = float(np.average(uz / gb, weights=w) * c)
    ke_mean = float(np.average(gb - 1.0, weights=w) * MC2_KEV)
    sz = float(np.sqrt(np.average((z - np.average(z, weights=w)) ** 2, weights=w)))
    rmax = float(np.sqrt(x**2 + y**2).max())
    print(f"Injected {z.size} macroparticles from prebuncher focus (iter {it_focus}); "
          f"⟨KE⟩ {ke_mean:.1f} keV, σ_z {sz*1e3:.2f} mm, r_max {rmax*1e3:.2f} mm, "
          f"v_beam {v_beam:.3e} m/s, q {w.sum()*q_e*1e9:.4f} nC", flush=True)
    # openPMD ux/uy/uz are γβ; PICMI wants proper velocity γβc in m/s → ×c.
    return (dict(x=x, y=y, z=z, ux=ux * c, uy=uy * c, uz=uz * c, w=w),
            v_beam, ke_mean)


def main():
    outdir = OUTDIR or DEFAULT_OUTDIR
    # Fresh diags: WarpX appends one openPMD file per dump, so a rerun of the same
    # case (e.g. a scan point) would otherwise mix old and new iterations. The diags
    # are git-ignored and regenerated, so clearing the case dir is safe.
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    # Recompute OMEGA here (not just at import) so a config(F_RF=...) override stays
    # consistent — the module-level OMEGA is frozen at import, before the override lands.
    omega = 2.0 * np.pi * F_RF

    bunch, v_beam, ke_mean = load_prebuncher_bunch()
    z_center = float(np.average(bunch["z"], weights=bunch["w"]))

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
    gain_keV = scale * 331.0                       # ≈ on-crest gain [keV] (1-kW V1kW≈331 kV)
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
    field_diag = picmi.FieldDiagnostic(
        name="fields", grid=grid, period=period,
        data_list=["phi", "rho", "E"],
        write_dir=outdir, warpx_format="openpmd", warpx_openpmd_backend="h5")
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
    sim.add_diagnostic(field_diag)
    sim.add_diagnostic(part_diag)

    print(f"\nRunning {n_steps} steps (diag every {period}) → {outdir}/")
    run_step(sim, n_steps, desc="linac_sec1")
    print("\nDone.")


if __name__ == "__main__":
    main()
