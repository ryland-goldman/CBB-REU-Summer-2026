"""
CESR prebuncher in WarpX (RZ): velocity-modulate the gun-exit beam with the RF
TM cavity field and let it bunch in the downstream drift, with self-consistent
space charge.

Third stage of the Cornell Linac chain in WarpX:
    cathode (cathode/) -> gun (gun/) -> prebuncher (this).

The gun's exit beam (~148 keV, β≈0.63, 0.1 nC, already RZ) is read from
`gun/diags/particles`, translated so it enters near z = 0, and tracked
through the prebuncher cavity. The cavity is the 1-J-normalised `prebuncher_25D`
field map (built by `build_prebuncher_field.py`) driven as a standing-wave TM
mode, reproducing GPT's `Map25D_TM` convention:

    Er,Ez(t) = map · scale · cos(ω t + φ)
    Bφ(t)    = H   · scale · sin(ω t + φ)        (E and B 90° out of phase)

with f_RF = 18 × master RF = 214.18 MHz and scale = sqrt(1e3·Q·P / (2π f_RF))
from the loaded-Q (Q = 3000) / dissipated-power (P) normalisation documented in
`reference/Linac Simulation Documentation/details.md`. The operating point
(POWER_KW, PHASE, OUTDIR) is set at module level below. main() also accepts the
same names as keyword overrides for direct/script use (e.g.
`prebuncher_sim.main(power=400, phase="zc", outdir=...)`); the pipeline instead
applies config() values by setattr on the module before calling bare main().
"""

import os
import shutil

import numpy as np
import pywarpx
from pywarpx import picmi
from openpmd_viewer import OpenPMDTimeSeries

from pipeline._runner import run_step
# F_RF/Q_L/V1J_KEV live in the (pywarpx-free) build module as the single source of
# truth, so the sim and plot_prebuncher.py cannot drift apart on the RF drive.
from .build_prebuncher_field import Z_GAP_CENTER, V1J_KEV, F_RF, Q_L

c = picmi.constants.c
m_e = picmi.constants.m_e
q_e = picmi.constants.q_e

# ── Prebuncher / field-map parameters (must match build_prebuncher_field.py) ───
PREBUNCH_FIELD = "prebuncher/prebuncher_field/prebuncher_EB.h5"
# F_RF (214.18 MHz), Q_L (3000), and V1J_KEV are imported from build_prebuncher_field.

GUN_DIAG = "gun/diags/particles"
Z_INJECT = 0.005                 # lab z where the bunch tail (smallest z) is placed [m]
MAX_PART = 50000                 # downsample the gun snapshot (reweighted) for speed
RNG_SEED = 0
CFL = 0.8                        # dt = CFL · Δz / v_beam

# ── Performance knobs (tunable via prebuncher.config(...); see run_pipeline.py) ─
# This stage dominates the pipeline (≈75% of runtime). Unlike the gun, do NOT coarsen
# NZ to go faster: this long-thin box is convergence-bound, not cell-bound, so the MLMG
# solve is ~1.37× slower per step at NZ=512 than NZ=1024 (and dz=2.5 mm then under-
# resolves the ~1 mm bunch) — see the README. Speed it via CFL (fewer steps) and
# MAX_ITERS/REQUIRED_PRECISION (cheaper solve). Defaults reproduce the original run.
REQUIRED_PRECISION = 1e-4        # MLMG relative tolerance (relaxed for the long-thin box)
MAX_ITERS = 500                  # MLMG iteration cap
MAX_STEPS = 0                    # 0 → auto-derive from transit; >0 → fixed
TRANSIT_MARGIN = 0.97            # stop just before the bunch centre reaches the exit
N_DIAGS = 60                     # number of openPMD dumps over the run

# ── Domain (RZ, single azimuthal mode — the cavity field is m = 0) ─────────────
RMAX = 0.036                     # covers the field-map bore (0–36.07 mm)
ZMAX = 1.30                      # entrance drift + 305 mm cavity + bunching drift
NR, NZ = 80, 1024                # divisible by the blocking factor (8)

# ── Operating point (the two undocumented prebuncher inputs; see details.md) ───
# POWER_KW = 0 runs the drift-only baseline (no cavity). OUTDIR defaults to None
# and is derived from POWER_KW/PHASE in main() — see prebuncher.resolve_outdir()
# so the parent process can compute the same path without importing pywarpx.
from . import DEFAULT_POWER_KW, DEFAULT_PHASE
POWER_KW = DEFAULT_POWER_KW      # dissipated RF power [kW]  (V_gap ≈ scale·439 kV)
PHASE = DEFAULT_PHASE            # "zc" = zero-crossing (bunching), "crest" = max gain
OUTDIR = None                    # if None at main(), derived from POWER_KW/PHASE


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
    x, y, z, ux, uy, uz, w = ts.get_particle(
        ["x", "y", "z", "ux", "uy", "uz", "w"], species="electrons", iteration=it,
    )
    # Downsample (reweighted to preserve total charge) to keep the scan cheap.
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


def main(power=None, phase=None, outdir=None, nz=None, zmax=None, max_steps=0):
    """Run one prebuncher case.

    All arguments default to the module-level constants (POWER_KW, PHASE, OUTDIR,
    NZ, ZMAX); pass keyword overrides to run a different operating point without
    editing the file.
    """
    if power is None:   power = POWER_KW
    if phase is None:   phase = PHASE
    if outdir is None:  outdir = OUTDIR
    if nz is None:      nz = NZ
    if zmax is None:    zmax = ZMAX
    if outdir is None:
        from . import _derive_outdir
        outdir = _derive_outdir(power, phase)

    # Fresh diags: WarpX appends one openPMD file per dump, so re-running the same case
    # (e.g. a power/phase scan point) would otherwise mix old and new iterations into one
    # series and corrupt the focus/σ_z analysis. diags are git-ignored and regenerated,
    # so clearing the case dir is safe. (Mirrors linac_sec1_sim.py.)
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)

    # Compute omega here, not at import, so a config(F_RF=...) override is honored
    # (an import-time module constant would be frozen before the override lands).
    omega = 2.0 * np.pi * F_RF

    bunch, v_beam, ke_mean = load_gun_bunch()

    # ── RF amplitude and phase ────────────────────────────────────────────────
    # scale = sqrt(stored_energy / 1 J),  stored_energy = 1e3·Q·P/(2π f_RF).
    # Power is in kW (GPT convention), hence the 1e3 kW->W factor.
    scale = float(np.sqrt(1e3 * Q_L * power / (2.0 * np.pi * F_RF)))
    # Time at which the bunch *tail* (placed at Z_INJECT) reaches the gap — the phase
    # reference. (Tail, not centre: the centre is ~σ_z/2 ≈ 0.5 mm farther, a <1° RF-phase
    # shift at 214 MHz, negligible; plot_prebuncher.py uses the same Z_INJECT reference so
    # the drawn waveform matches the run.) The energy kick of an electron is
    # ΔW(t) ∝ -cos(ω t + φ) (on-axis Ez is single-signed positive), so bunching needs the
    # *tail* (later arrival) to gain energy:
    #   d/dt[-cos] = ω sin(ω t_gap + φ) > 0  ->  sin(ω t_gap + φ) = +1.
    # zc:    φ = -ω t_gap + π/2  (zero net kick, max +slope -> velocity bunching)
    # crest: φ = -ω t_gap + π    (-cos = +1 -> maximum energy gain, little bunching)
    t_gap = (Z_GAP_CENTER - Z_INJECT) / v_beam
    if phase == "zc":
        phi = -omega * t_gap + np.pi / 2.0
    else:
        phi = -omega * t_gap + np.pi
    print(f"Case: P={power:g} kW, phase={phase}, scale={scale:.3f}, "
          f"f_RF={F_RF/1e6:.2f} MHz, φ={phi:.3f} rad, t_gap={t_gap*1e9:.3f} ns",
          flush=True)

    # ── Grid + electrostatic (self-field) solver ──────────────────────────────
    grid = picmi.CylindricalGrid(
        number_of_cells=[NR, nz],
        n_azimuthal_modes=1,
        lower_bound=[0.0, 0.0],
        upper_bound=[RMAX, zmax],
        lower_boundary_conditions=["none", "dirichlet"],
        upper_boundary_conditions=["neumann", "dirichlet"],
        lower_boundary_conditions_particles=["none", "absorbing"],
        upper_boundary_conditions_particles=["absorbing", "absorbing"],
        warpx_blocking_factor=8,
    )
    # MLMG converges slowly in this long, thin (1.3 m × 36 mm) box. The 0.1 nC
    # space-charge field is a small perturbation on the 148 keV beam, so a
    # relative precision of 1e-4 is ample and far cheaper than 1e-5.
    solver = picmi.ElectrostaticSolver(
        grid=grid, method="Multigrid", required_precision=REQUIRED_PRECISION,
        maximum_iterations=MAX_ITERS, warpx_self_fields_verbosity=0,
    )

    # ── Cavity field: 1-J map × scale × cos/sin(ω t + φ), read from file ───────
    # E and B share one openPMD path; E uses cos, B uses sin (TM standing wave).
    # Constants are baked into the parser strings: the LoadAppliedField wrapper
    # forwards extra kwargs to the picmistandard base class, which rejects them.
    # P = 0 is the drift-only baseline (no cavity) used to isolate bunching.
    prebuncher = None
    if power > 0:
        e_time = f"{scale:.10e}*cos({omega:.10e}*t + ({phi:.10e}))"
        b_time = f"{scale:.10e}*sin({omega:.10e}*t + ({phi:.10e}))"
        prebuncher = picmi.LoadAppliedField(
            read_fields_from_path=PREBUNCH_FIELD,
            load_E=True, load_B=True,
            warpx_E_time_function=e_time,
            warpx_B_time_function=b_time,
        )

    electrons = picmi.Species(
        particle_type="electron", name="electrons",
        initial_distribution=picmi.ParticleListDistribution(
            x=bunch["x"], y=bunch["y"], z=bunch["z"],
            ux=bunch["ux"], uy=bunch["uy"], uz=bunch["uz"], weight=bunch["w"],
        ),
    )

    # ── Time step / duration ──────────────────────────────────────────────────
    dz = zmax / nz
    dt = CFL * dz / v_beam
    # Stop just before the bunch *centre* reaches the exit (margin < 1): once the
    # beam clears the absorbing boundary the domain empties and the Multigrid solve
    # aborts. margin = 0.97 keeps the bunch in-domain while still capturing the full
    # cavity + bunching drift (all foci sit at ≤ 1.16 m < 0.97·1.30 m). At the crest
    # the cavity accelerates the beam (+V_gap ≈ scale·V1J), so the post-gap drift is
    # crossed faster — size that segment with the accelerated speed.
    if phase == "crest":
        ke_final = ke_mean + scale * V1J_KEV
        gamma_f = 1.0 + ke_final / (m_e * c**2 / q_e / 1e3)
        v_final = c * np.sqrt(1.0 - 1.0 / gamma_f**2)
        transit = ((Z_GAP_CENTER - Z_INJECT) / v_beam
                   + (zmax - Z_GAP_CENTER) / v_final)
    else:
        transit = (zmax - Z_INJECT) / v_beam
    n_steps = max_steps or MAX_STEPS or int(TRANSIT_MARGIN * transit / dt)
    print(f"dt = {dt:.3e} s, max_steps = {n_steps}, "
          f"RF period = {1/F_RF*1e9:.2f} ns ({1/F_RF/dt:.0f} steps/period)",
          flush=True)

    # ── Diagnostics (openPMD, HDF5) ───────────────────────────────────────────
    period = max(1, n_steps // N_DIAGS)
    part_diag = picmi.ParticleDiagnostic(
        name="particles", period=period, species=[electrons],
        data_list=["position", "momentum", "weighting"],
        write_dir=outdir, warpx_format="openpmd", warpx_openpmd_backend="h5",
    )

    sim = picmi.Simulation(
        solver=solver, max_steps=n_steps, time_step_size=dt,
        verbose=0,                     # silence per-step "STEP N starts" — the tqdm bar is the progress display
        particle_shape="linear",
    )
    # ParticleListDistribution supplies the macroparticles explicitly, so this layout
    # (and its n_macroparticles_per_cell) is inert — the count is the imported (downsampled) list size.
    sim.add_species(
        electrons,
        layout=picmi.PseudoRandomLayout(n_macroparticles_per_cell=1, grid=grid),
    )
    if prebuncher is not None:
        sim.add_applied_field(prebuncher)
    sim.add_diagnostic(part_diag)

    print(f"\nRunning {n_steps} steps (diag every {period}) -> {outdir}/")
    run_step(sim, n_steps, desc="prebuncher")
    print("\nDone.")


if __name__ == "__main__":
    main()
