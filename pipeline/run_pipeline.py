"""End-to-end Cornell Linac beam simulation in WarpX.

Orchestrates the three stages in order from one driver process:

    cathode -> gun -> prebuncher

Each stage's build/plot run in-process, but its WarpX simulation runs in a fresh
Python subprocess (`pipeline._launch_sim`) to sidestep pywarpx's per-process
geometry binding (cathode 2D -> gun/prebuncher RZ would otherwise trip a diagnostic
state assertion). Each stage exposes `config(**kwargs)` (override module-level
constants) and `run()` (build field map if any, simulate, plot). Stage execution is
wrapped by the shared runner in `pipeline/_runner.py`, which drives the tqdm bar from
an afterstep callback and redirects WarpX's per-step stdout to the pipeline log file,
writing a full DEBUG-level transcript to `pipeline/logs/pipeline_<ts>.log`.

Run with:
    conda activate CBB
    python pipeline/run_pipeline.py
"""

import os
import sys

# Set OMP_NUM_THREADS BEFORE any pywarpx import (read by OpenMP at load time).
# Default 1: these stages run fastest single-threaded — the grids are small and
# the MLMG Poisson solve is memory-bandwidth bound, so OpenMP threads contend for
# the same memory bus and add fork/join + barrier overhead without speeding up the
# solve (measured: full chain ~1.1 min at OMP=1; OMP=14 showed only ~450% CPU and
# no gain). Keep this single-threaded; only raise OMP_THREADS for the much larger
# original-config grids, where per-thread work outgrows the overhead.
os.environ.setdefault("OMP_NUM_THREADS", os.environ.get("OMP_THREADS", "1"))

# Run from the repo root so each stage's hard-coded relative paths resolve.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(_ROOT)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import cathode
import gun
import prebuncher
import linac_sec1

from pipeline._runner import setup_logging, _cl, _BOLD, _RESET

# ── Operating-point overrides (physics; defaults live in the stage modules) ──
cathode.config(V_anode=50.0)
gun.config(GUN_VOLTAGE=150e3, BUNCH_CHARGE=0.1e-9)
prebuncher.config(POWER_W=800, PHASE="zc")
linac_sec1.config(POWER_MW=15.0)            # SLAC Section 1 input power (~37 MeV on crest)

# ── Performance knobs (accuracy ↔ speed). Full knob list, runtime split, and the
#    reason NZ must stay at 1024: see pipeline/README.md § Configuration. ──
# Balanced profile: ACTIVE (~1.7×, ~5 min). Comment these 3 lines for the baseline.
cathode.config(PPC=6, REQUIRED_PRECISION=3e-5)
gun.config(nz=256, MAX_PART=50000, REQUIRED_PRECISION=1e-4)
prebuncher.config(CFL=0.95, MAX_ITERS=150, REQUIRED_PRECISION=1e-3)
# Conservative (~1.3×, near-identical):
# gun.config(MAX_PART=80000, REQUIRED_PRECISION=1e-4)
# prebuncher.config(REQUIRED_PRECISION=2e-4, MAX_ITERS=400)
# Aggressive (~2.2×, looser space-charge solve):
# cathode.config(nz=48, PPC=4, REQUIRED_PRECISION=5e-5, MAX_STEPS=1200)
# gun.config(nz=192, MAX_PART=40000, REQUIRED_PRECISION=3e-4, N_DIAGS=20)
# prebuncher.config(CFL=0.97, MAX_ITERS=80, REQUIRED_PRECISION=3e-3, N_DIAGS=20)
# linac_sec1: the demo() driver (called below) sets its own grid per case — a coarse
# RF-phase acceptance scan (scan_nz) + a full-resolution headline/focus-off (full_nz).
# Tune cost via linac_sec1.demo(scan_nz=..., full_nz=..., phases=...).


def _beam_summary(diag, label, unit="keV"):
    """Report the final bunch from the last snapshot of `diag` (console + log).

    `unit` selects the kinetic-energy scale ("keV" for the prebuncher exit, "MeV"
    for the linac exit). Also reports the captured-charge fraction (last/first
    snapshot) when both are available.
    """
    try:
        import numpy as np
        from openpmd_viewer import OpenPMDTimeSeries
        ts = OpenPMDTimeSeries(os.path.join(diag, "particles"))
        its = list(ts.iterations)
        q0 = None
        if its:
            _, _, _, _, w0 = ts.get_particle(
                ["z", "ux", "uy", "uz", "w"], species="electrons", iteration=its[0])
            q0 = w0.sum()
        z = None
        for it in reversed(its):
            z, ux, uy, uz, w = ts.get_particle(
                ["z", "ux", "uy", "uz", "w"], species="electrons", iteration=it)
            if len(z) > 50:
                break
        if z is None or len(z) <= 50:
            _cl(f"\n(final-beam summary [{label}]: no snapshot with >50 macroparticles — "
                "the beam may have cleared the domain)")
            return
        fac = 0.51099895e3 if unit == "keV" else 0.51099895
        ke = (np.sqrt(1 + ux**2 + uy**2 + uz**2) - 1) * fac
        zm = np.average(z, weights=w)
        sz = np.sqrt(np.average((z - zm) ** 2, weights=w))
        km = np.average(ke, weights=w)
        dk = np.sqrt(np.average((ke - km) ** 2, weights=w))
        cap = f"   captured {w.sum()/q0*100:.0f}%" if q0 else ""
        _cl(f"\n{_BOLD}Final beam{_RESET} ({label}, {len(z)} macroparticles):")
        _cl(f"      ⟨z⟩ = {zm*1e3:.0f} mm   σ_z = {sz*1e3:.3f} mm   "
            f"⟨KE⟩ = {km:.1f} {unit}   σ_KE = {dk:.2f} {unit}   "
            f"q = {w.sum()*1.602176634e-19*1e9:.3f} nC{cap}")
    except Exception as e:
        import logging
        _cl(f"    (final-beam summary [{label}] unavailable: {e})", level=logging.WARNING)


def main():
    log_path = setup_logging()
    import time
    t0 = time.time()
    _cl("=" * 72)
    _cl(" Cornell Linac WarpX pipeline:  cathode -> gun -> prebuncher -> linac_sec1")
    _cl(f" OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS', '?')}")
    _cl("=" * 72)
    print(f" log: {log_path}")

    cathode.run()
    gun.run()
    prebuncher.run()
    linac_sec1.demo()           # RF-phase acceptance scan + headline + focus-off + plots

    _beam_summary(prebuncher.resolve_outdir(), "prebuncher exit", "keV")
    _beam_summary(linac_sec1.resolve_outdir(), "linac_sec1 exit", "MeV")

    _cl("\n" + "=" * 72)
    _cl(f" Pipeline complete in {(time.time()-t0)/60:.1f} min.")
    _cl(" Figures: cathode/results/, gun/results/, prebuncher/results/, linac_sec1/results/")
    _cl("=" * 72)
    print(f" log: {log_path}")


if __name__ == "__main__":
    main()
