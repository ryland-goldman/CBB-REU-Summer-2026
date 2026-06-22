"""Shared beam-handoff helpers for the WarpX stages (gun, injector, linac_sec1).

The single source of truth for the idioms every stage's `load_*_bunch` / handoff
path repeated: building a `pmd_beamphysics.ParticleGroup` from RZ phase space,
weighted downsampling, mean velocity/energy from γβ momenta, and the RF
`warpx_E/B_time_function` strings. See pipeline/README.md → Shared modules.

Conventions (CLAUDE.md): ux/uy/uz are dimensionless γβ; momentum is eV/c
(= γβ·MC2_EV); ParticleGroup `weight` is per-macroparticle CHARGE [C] (= count·q_e).
"""

import numpy as np

from pipeline.constants import C_LIGHT, E_CHARGE as Q_E, MC2_EV


def open_particle_series(diag, stage_hint=None):
    """Open an openPMD particle series, raising if it has no iterations.

    `stage_hint` (e.g. "gun") tags the error with which upstream stage failed to write.
    """
    from openpmd_viewer import OpenPMDTimeSeries
    ts = OpenPMDTimeSeries(diag)
    if len(ts.iterations) == 0:
        hint = f" — did the {stage_hint} stage run?" if stage_hint else ""
        raise RuntimeError(f"{diag} has no iterations{hint}")
    return ts


def make_particle_group(x, y, z, ux, uy, uz, w):
    """ParticleGroup from RZ phase space: ux/uy/uz = γβ, w = macroparticle count.

    Momentum px/py/pz = γβ·MC2_EV [eV/c]; weight = count·q_e [C]; species "electron"
    (SINGULAR — ParticleGroup's spelling; openPMD readers key the PLURAL "electrons").
    """
    from pmd_beamphysics import ParticleGroup
    n = x.size
    return ParticleGroup(data=dict(
        x=x, y=y, z=z,
        px=ux * MC2_EV, py=uy * MC2_EV, pz=uz * MC2_EV,
        t=np.zeros(n), weight=w * Q_E,
        status=np.ones(n, dtype=np.int64), species="electron"))


def downsample(arrays, w, max_part, rng):
    """Randomly thin to `max_part` macroparticles, rescaling weights to conserve charge.

    No-op (returns the inputs) if `max_part` is falsy or already ≤ the count. Draws one
    `rng.choice` from the passed generator, so the caller controls the RNG stream/seed.
    Returns (tuple of downsampled arrays, rescaled w).
    """
    n = w.size
    if not max_part or n <= max_part:
        return tuple(arrays), w
    sel = rng.choice(n, max_part, replace=False)
    return tuple(a[sel] for a in arrays), w[sel] * (n / max_part)


def beam_kinematics(ux, uy, uz, w):
    """(weighted mean v_z [m/s], weighted mean KE [keV]) from γβ momenta."""
    gb = np.sqrt(1.0 + ux**2 + uy**2 + uz**2)              # γ (ux/uy/uz are γβ)
    v_beam = float(np.average(uz / gb, weights=w) * C_LIGHT)
    ke_mean_keV = float(np.average(gb - 1.0, weights=w) * MC2_EV / 1e3)
    return v_beam, ke_mean_keV


def rf_time_functions(scale, omega, phi, amp_prec=10, phase_prec=10):
    """(E, B) `warpx_*_time_function` strings for a standing-wave TM cavity drive.

    E ∝ scale·cos(ωt+φ), B ∝ scale·sin(ωt+φ). ω keeps .10e (its truncation accumulates
    over the ~5 ns transit); `amp_prec`/`phase_prec` set the scale/φ precision per caller.
    """
    e = f"{scale:.{amp_prec}e}*cos({omega:.10e}*t + ({phi:.{phase_prec}e}))"
    b = f"{scale:.{amp_prec}e}*sin({omega:.10e}*t + ({phi:.{phase_prec}e}))"
    return e, b
