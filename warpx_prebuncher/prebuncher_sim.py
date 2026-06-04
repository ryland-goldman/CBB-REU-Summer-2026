"""
CESR prebuncher in WarpX (RZ): velocity-modulate the gun-exit beam with the RF
TM cavity field and let it bunch in the downstream drift, with self-consistent
space charge.

Third stage of the Cornell Linac chain in WarpX:
    cathode (warpx_cathode/) -> gun (warpx_gun/) -> prebuncher (this).

The gun's exit beam (~148 keV, β≈0.63, 0.1 nC, already RZ) is read from
`warpx_gun/diags/particles`, translated so it enters near z = 0, and tracked
through the prebuncher cavity. The cavity is the 1-J-normalised `prebuncher_25D`
field map (built by `build_prebuncher_field.py`) driven as a standing-wave TM
mode, reproducing GPT's `Map25D_TM` convention:

    Er,Ez(t) = map · scale · cos(ω t + φ)
    Bφ(t)    = H   · scale · sin(ω t + φ)        (E and B 90° out of phase)

with f_RF = 18 × master RF = 214.18 MHz and scale = sqrt(1e3·Q·P / (2π f_RF))
from the loaded-Q (Q = 3000) / dissipated-power (P) normalisation documented in
`reference/Linac Simulation Documentation/details.md`. P and the RF phase are the
two undocumented operating-point inputs; this script takes them on the CLI so a
shell loop over it (one --outdir per power) can sweep power and compare phases.

WarpX's electrostatic Multigrid solver supplies the beam self-field; the cavity
field is applied to the particles via `picmi.LoadAppliedField` (read_from_file +
time-dependent E/B scaling functions).

Run one case with:
    conda run -n CBB python warpx_prebuncher/prebuncher_sim.py \
        --power 40 --phase zc --outdir warpx_prebuncher/diags/P40_zc
"""

import argparse
import numpy as np
import pywarpx
from pywarpx import picmi
from openpmd_viewer import OpenPMDTimeSeries

from build_prebuncher_field import Z_GAP_CENTER   # keep field/phasing in sync

c = picmi.constants.c
m_e = picmi.constants.m_e
q_e = picmi.constants.q_e

# ── Prebuncher / field-map parameters (must match build_prebuncher_field.py) ───
PREBUNCH_FIELD = "warpx_prebuncher/prebuncher_field/prebuncher_EB.h5"
F_RF = 499.7645e6 / 42 * 18      # 18 × master RF = 214.18 MHz (details.md)
OMEGA = 2.0 * np.pi * F_RF
Q_L = 3000                       # loaded Q of prebuncher 1 (details.md)

GUN_DIAG = "warpx_gun/diags/particles"
Z_INJECT = 0.005                 # lab z where the bunch head is placed [m]
MAX_PART = 50000                 # downsample the gun snapshot (reweighted) for speed
RNG_SEED = 0
CFL = 0.8                        # dt = CFL · Δz / v_beam

# ── Domain (RZ, single azimuthal mode — the cavity field is m = 0) ─────────────
RMAX = 0.036                     # covers the field-map bore (0–36.07 mm)
ZMAX = 1.30                      # entrance drift + 305 mm cavity + bunching drift
nr, nz = 80, 1024                # divisible by the blocking factor (8)

# ── Single-case operating point (set here, or override on the CLI) ─────────────
# P and the RF phase are the two undocumented prebuncher inputs (details.md).
# P = 0 runs the drift-only baseline (no cavity).
POWER_W = 800.0                  # dissipated RF power [W]  (V_gap ≈ scale·430 kV)
PHASE = "zc"                     # "zc" = zero-crossing (bunching), "crest" = max gain
OUTDIR = "warpx_prebuncher/diags/P800_zc"


def parse_args():
    p = argparse.ArgumentParser(description="WarpX CESR prebuncher (one case).")
    p.add_argument("--power", type=float, default=POWER_W,
                   help="dissipated RF power P [W] -> field scale (0 = drift baseline)")
    p.add_argument("--phase", choices=["zc", "crest"], default=PHASE,
                   help="zc = zero-crossing (velocity bunching); "
                        "crest = max energy gain")
    p.add_argument("--outdir", default=OUTDIR,
                   help="openPMD diagnostics output directory")
    p.add_argument("--nz", type=int, default=nz, help="cells in z")
    p.add_argument("--zmax", type=float, default=ZMAX, help="domain length [m]")
    p.add_argument("--max-steps", type=int, default=0,
                   help="override step count (0 = auto from transit time)")
    return p.parse_args()


def load_gun_bunch():
    """Import the gun's last snapshot (already RZ) and shift it to the entrance.

    Returns (dict for ParticleListDistribution, v_beam, mean KE [keV]).
    """
    ts = OpenPMDTimeSeries(GUN_DIAG)
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
    # Translate so the bunch *head* (smallest z) sits at Z_INJECT.
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


def main():
    args = parse_args()
    bunch, v_beam, ke_mean = load_gun_bunch()

    # ── RF amplitude and phase ────────────────────────────────────────────────
    # scale = sqrt(stored_energy / 1 J),  stored_energy = 1e3·Q·P/(2π f_RF).
    scale = float(np.sqrt(1e3 * Q_L * args.power / (2.0 * np.pi * F_RF)))
    # Time at which the bunch centre reaches the gap. The energy kick of an
    # electron is ΔW(t) ∝ -cos(ω t + φ) (on-axis Ez is single-signed positive),
    # so bunching needs the *tail* (later arrival) to gain energy:
    #   d/dt[-cos] = ω sin(ω t_gap + φ) > 0  ->  sin(ω t_gap + φ) = +1.
    # zc:    φ = -ω t_gap + π/2  (zero net kick, max +slope -> velocity bunching)
    # crest: φ = -ω t_gap + π    (-cos = +1 -> maximum energy gain, little bunching)
    t_gap = (Z_GAP_CENTER - Z_INJECT) / v_beam
    if args.phase == "zc":
        phi = -OMEGA * t_gap + np.pi / 2.0
    else:
        phi = -OMEGA * t_gap + np.pi
    print(f"Case: P={args.power:g} W, phase={args.phase}, scale={scale:.3f}, "
          f"f_RF={F_RF/1e6:.2f} MHz, φ={phi:.3f} rad, t_gap={t_gap*1e9:.3f} ns",
          flush=True)

    # ── Grid + electrostatic (self-field) solver ──────────────────────────────
    grid = picmi.CylindricalGrid(
        number_of_cells=[nr, args.nz],
        n_azimuthal_modes=1,
        lower_bound=[0.0, 0.0],
        upper_bound=[RMAX, args.zmax],
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
        grid=grid, method="Multigrid", required_precision=1e-4,
        maximum_iterations=500, warpx_self_fields_verbosity=0,
    )

    # ── Cavity field: 1-J map × scale × cos/sin(ω t + φ), read from file ───────
    # E and B share one openPMD path; E uses cos, B uses sin (TM standing wave).
    # Constants are baked into the parser strings: the LoadAppliedField wrapper
    # forwards extra kwargs to the picmistandard base class, which rejects them.
    # P = 0 is the drift-only baseline (no cavity) used to isolate bunching.
    prebuncher = None
    if args.power > 0:
        e_time = f"{scale:.10e}*cos({OMEGA:.10e}*t + ({phi:.10e}))"
        b_time = f"{scale:.10e}*sin({OMEGA:.10e}*t + ({phi:.10e}))"
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
    dz = args.zmax / args.nz
    dt = CFL * dz / v_beam
    # Stop just before the bunch *centre* reaches the exit (margin < 1): once the
    # beam clears the absorbing boundary the domain empties and the Multigrid solve
    # aborts. margin = 0.97 keeps the bunch in-domain while still capturing the full
    # cavity + bunching drift (all foci sit at ≤ 1.16 m < 0.97·1.30 m). At the crest
    # the cavity accelerates the beam (+V_gap ≈ scale·V1J), so the post-gap drift is
    # crossed faster — size that segment with the accelerated speed.
    V1J_KEV = 430.2                                  # 1-J effective gap voltage
    if args.phase == "crest":
        ke_final = ke_mean + scale * V1J_KEV
        gamma_f = 1.0 + ke_final / (m_e * c**2 / q_e / 1e3)
        v_final = c * np.sqrt(1.0 - 1.0 / gamma_f**2)
        transit = ((Z_GAP_CENTER - Z_INJECT) / v_beam
                   + (args.zmax - Z_GAP_CENTER) / v_final)
    else:
        transit = (args.zmax - Z_INJECT) / v_beam
    max_steps = args.max_steps or int(0.97 * transit / dt)
    print(f"dt = {dt:.3e} s, max_steps = {max_steps}, "
          f"RF period = {1/F_RF*1e9:.2f} ns ({1/F_RF/dt:.0f} steps/period)",
          flush=True)

    # ── Diagnostics (openPMD, HDF5) ───────────────────────────────────────────
    period = max(1, max_steps // 60)
    field_diag = picmi.FieldDiagnostic(
        name="fields", grid=grid, period=period,
        data_list=["phi", "rho", "E"],
        write_dir=args.outdir, warpx_format="openpmd", warpx_openpmd_backend="h5",
    )
    part_diag = picmi.ParticleDiagnostic(
        name="particles", period=period, species=[electrons],
        data_list=["position", "momentum", "weighting"],
        write_dir=args.outdir, warpx_format="openpmd", warpx_openpmd_backend="h5",
    )

    sim = picmi.Simulation(
        solver=solver, max_steps=max_steps, time_step_size=dt,
        verbose=1, particle_shape="linear",
    )
    sim.add_species(
        electrons,
        layout=picmi.PseudoRandomLayout(n_macroparticles_per_cell=1, grid=grid),
    )
    if prebuncher is not None:
        sim.add_applied_field(prebuncher)
    sim.add_diagnostic(field_diag)
    sim.add_diagnostic(part_diag)

    print(f"\nRunning {max_steps} steps (diag every {period}) -> {args.outdir}/")
    sim.step(max_steps)
    print("\nDone.")


if __name__ == "__main__":
    main()
