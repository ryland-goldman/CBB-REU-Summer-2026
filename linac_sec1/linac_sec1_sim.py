"""
SLAC Linac Section 1 in WarpX (RZ): capture the injector beam in a 3 m 2π/3
traveling-wave accelerating structure with self-consistent space charge. Reads the
injector handoff beam, applies the iris scrape, tracks through the two quadrature RF
maps, and writes openPMD diagnostics.

See linac_sec1/README.md for physics, parameters, and gotchas.
"""

import os
import json
import shutil
import numpy as np
from pywarpx import picmi
from openpmd_viewer import OpenPMDTimeSeries

from pipeline._runner import run_step
from pipeline.collimator import pipe_violator_ids, survivor_mask
from .build_linac_sec1_field import Z_STRUCT, RMAX, BORE_R, V1KW_KEV
from . import DEFAULT_OUTDIR

c = picmi.constants.c
m_e = picmi.constants.m_e
q_e = picmi.constants.q_e
MC2_KEV = m_e * c**2 / q_e / 1e3        # electron rest energy [keV] ≈ 511

# ── Field maps (built by build_linac_sec1_field.py): two SLAC quadrature RF maps ──
RF1_FIELD = "linac_sec1/linac_sec1_field/linac_rf1.h5"
RF2_FIELD = "linac_sec1/linac_sec1_field/linac_rf2.h5"

F_RF = 2856.0e6                  # SLAC S-band [Hz] (Linac_RF in details.md)
RF_NORM_MW = 0.001               # field-map power normalisation (1 kW)

# ── Upstream input ────────────────────────────────────────────────────────────
INJECTOR_DIAG = "injector/diags/main/particles"
Z_HANDOFF = 2.03                 # [m] injector→linac handoff plane (Z_acc_1)
COLLIM_Z = 1.922                 # [m] iris start; iris scrape is multi-plane from here, NOT a
                                 #     single cut at the 2.03 m handoff (beam converges through tail)
Z_INJECT = 0.005                 # lab z where the bunch head is placed [m]
MAX_PART = 50000                 # downsample the injected snapshot (reweighted)
RNG_SEED = 0

# ── Operating point (tunable via linac_sec1.config(...)) ──────────────────────
POWER_MW = 11.0                  # RF input power [MW] (sec1_input_power)
PHASE_DEG = 0.0                  # injection RF phase offset [deg]; crest found empirically

# ── Performance / domain knobs ────────────────────────────────────────────────
CFL = 0.5                        # dt = CFL · Δz / v_inject
REQUIRED_PRECISION = 1e-4        # MLMG relative tolerance (space charge is a small perturbation)
SPACE_CHARGE = True              # beam self-field; False → warpx_do_not_deposit (RF maps only)
MAX_ITERS = 200                  # MLMG iteration cap
MAX_STEPS = 0                    # 0 → auto-derive from transit; >0 → fixed
TRANSIT_MARGIN = 1.0             # transit already targets Z_END short of the wall
N_DIAGS = 60                     # number of openPMD dumps over the run

# ── Domain (RZ, single azimuthal mode — the maps are m = 0) ───────────────────
# ZMAX includes a field-free exit drift so the bunch coasts (not absorbed) at the last dump.
ZMAX = 3.50                      # [m]
# Keep cells near the ≈3:1 rule or the MLMG self-field solve aborts; if it diverges raise
# NR (÷ blocking factor) rather than coarsening NZ.
NR, NZ = 16, 1664

OUTDIR = None                    # if None at main(), use DEFAULT_OUTDIR


def load_injector_bunch():
    """Import the injector beam at the z ≈ 2.03 m handoff plane and shift it to entry.

    Selects the dump whose bunch ⟨z⟩ is nearest Z_HANDOFF (not min-σ_z / max-in-bore).
    Returns (dict for ParticleListDistribution, v_beam, mean KE [keV], inj summary).
    """
    ts = OpenPMDTimeSeries(INJECTOR_DIAG)
    if len(ts.iterations) == 0:
        raise RuntimeError(
            f"{INJECTOR_DIAG} has no iterations — did the injector stage run?")

    # One light pass to find the well-populated dump nearest the handoff plane; the
    # n ≥ 0.8·nmax gate avoids a depleted late dump that happens to sit near 2.03 m.
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
    it_handoff, zm_handoff = min(cands, key=lambda t: abs(t[1] - Z_HANDOFF))
    if abs(zm_handoff - Z_HANDOFF) > 0.02:
        print(f"  WARNING: nearest injector dump to the {Z_HANDOFF*1e3:.0f} mm handoff is at "
              f"⟨z⟩={zm_handoff*1e3:.0f} mm ({abs(zm_handoff-Z_HANDOFF)*1e3:.0f} mm off) — "
              f"the handoff diagnostic may be too coarse; the injected beam is off-plane.",
              flush=True)

    idh, x, y, z, ux, uy, uz, w = ts.get_particle(
        ["id", "x", "y", "z", "ux", "uy", "uz", "w"], species="electrons", iteration=it_handoff)
    # TRUE injected charge = full beam at the handoff (all r), the capture denominator;
    # recorded BEFORE the iris scrape below.
    q_inj = float(w.sum()) * q_e

    # Iris/pipe collimation: multi-plane id scrape (a particle outside the iris at ANY
    # plane z ≥ COLLIM_Z hit the wall), NOT a single 2.03 m cut. See pipeline.collimator.
    scan_iters = [it for it, _n, zm in recs if (COLLIM_Z - 0.05) <= zm <= (Z_HANDOFF + 0.03)]
    violators = pipe_violator_ids(ts, scan_iters, RMAX, COLLIM_Z)
    keep = survivor_mask(idh, violators)
    idh, x, y, z, ux, uy, uz, w = (a[keep] for a in (idh, x, y, z, ux, uy, uz, w))
    r = np.hypot(x, y)
    q_dom = float(w.sum()) * q_e                         # in-iris survivors
    q_bore = float(w[r <= BORE_R].sum()) * q_e           # of those, within the RF bore

    # Downsample the survivors, reweighted to preserve total charge.
    if z.size > MAX_PART:
        rng = np.random.default_rng(RNG_SEED)
        sel = rng.choice(z.size, MAX_PART, replace=False)
        scale_w = z.size / MAX_PART
        x, y, z, ux, uy, uz, w = (a[sel] for a in (x, y, z, ux, uy, uz, w))
        w = w * scale_w
    # Translate the bunch tail (smallest z; beam travels +z) to Z_INJECT.
    z = z - z.min() + Z_INJECT

    gb = np.sqrt(1.0 + ux**2 + uy**2 + uz**2)          # γ (ux/uy/uz are γβ)
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
    # openPMD ux/uy/uz are γβ; PICMI wants proper velocity γβc [m/s] → ×c.
    return (dict(x=x, y=y, z=z, ux=ux * c, uy=uy * c, uz=uz * c, w=w),
            v_beam, ke_mean, inj)


def main():
    outdir = OUTDIR or DEFAULT_OUTDIR
    # Fresh diags: WarpX appends one file per dump, so a rerun would mix iterations.
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    # omega in main() (not at import) so a config(F_RF=...) override is honoured.
    omega = 2.0 * np.pi * F_RF

    bunch, v_beam, ke_mean, inj = load_injector_bunch()
    z_center = float(np.average(bunch["z"], weights=bunch["w"]))
    # Persist the true injected charge so the plotter reports capture against it, not the
    # post-scrape first dump. WarpX makes outdir at first diag; make it now for this sidecar.
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "injection_summary.json"), "w") as fh:
        json.dump(inj, fh, indent=2)

    # ── RF amplitude + phase ──────────────────────────────────────────────────
    scale = float(np.sqrt(POWER_MW / RF_NORM_MW))
    # Phase referenced to bunch-centre arrival at the structure entrance (absolute crest
    # undocumented; PHASE_DEG scanned).
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

    # ── Applied fields: the two quadrature RF maps, summed by WarpX ───────────
    # Constants are baked into the AMReX parser strings (LoadAppliedField rejects extra kwargs).
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
    # PICMI guard: the LAST applied field must load_E, else the global E_ext init style
    # is forced "none" and the beam coasts unaccelerated.
    assert getattr(applied[-1], "load_E", False), (
        "last applied field must have load_E=True (an RF map), or the global E_ext "
        "style is forced to 'none' and the beam coasts unaccelerated")

    electrons = picmi.Species(
        particle_type="electron", name="electrons",
        initial_distribution=picmi.ParticleListDistribution(
            x=bunch["x"], y=bunch["y"], z=bunch["z"],
            ux=bunch["ux"], uy=bunch["uy"], uz=bunch["uz"], weight=bunch["w"]),
        warpx_do_not_deposit=not SPACE_CHARGE)

    # ── Time step / duration ──────────────────────────────────────────────────
    # dt sized at the slowest (injection) velocity. Transit estimated in segments
    # (drift + capture + relativistic remainder) to a stop plane short of ZMAX so the
    # bunch finishes in the field-free exit drift, NOT at the absorbing wall (an empty
    # domain aborts the MLMG solve). gain_keV uses the on-crest max — the binding case.
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
