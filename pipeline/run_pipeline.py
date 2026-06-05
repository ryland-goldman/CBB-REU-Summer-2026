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

from pipeline._runner import setup_logging, _cl, _BOLD, _RESET

# ── Operating-point overrides (physics; defaults live in the stage modules) ──
cathode.config(V_anode=50.0)
gun.config(GUN_VOLTAGE=150e3, BUNCH_CHARGE=0.1e-9)
prebuncher.config(POWER_W=800, PHASE="zc")

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


def _final_beam_summary(diag):
    """Report the final bunch from the last prebuncher snapshot (console + log)."""
    try:
        import numpy as np
        from openpmd_viewer import OpenPMDTimeSeries
        ts = OpenPMDTimeSeries(os.path.join(diag, "particles"))
        z = None
        for it in reversed(ts.iterations):
            z, ux, uy, uz, w = ts.get_particle(
                ["z", "ux", "uy", "uz", "w"], species="electrons", iteration=it)
            if len(z) > 50:
                break
        if z is None or len(z) <= 50:
            _cl("\n(final-beam summary: no snapshot with >50 macroparticles — "
                "the beam may have cleared the domain)")
            return
        ke = (np.sqrt(1 + ux**2 + uy**2 + uz**2) - 1) * 0.51099895e3
        zm = np.average(z, weights=w)
        sz = np.sqrt(np.average((z - zm) ** 2, weights=w))
        km = np.average(ke, weights=w)
        dk = np.sqrt(np.average((ke - km) ** 2, weights=w))
        _cl(f"\n{_BOLD}Final beam{_RESET} (prebuncher exit, {len(z)} macroparticles):")
        _cl(f"      ⟨z⟩ = {zm*1e3:.0f} mm   σ_z = {sz*1e3:.3f} mm   "
            f"⟨KE⟩ = {km:.1f} keV   σ_KE = {dk:.2f} keV   "
            f"q = {w.sum()*1.602176634e-19*1e9:.3f} nC")
    except Exception as e:
        import logging
        _cl(f"    (final-beam summary unavailable: {e})", level=logging.WARNING)


def main():
    log_path = setup_logging()
    import time
    t0 = time.time()
    _cl("=" * 72)
    _cl(" Cornell Linac WarpX pipeline:  cathode -> gun -> prebuncher")
    _cl(f" OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS', '?')}")
    _cl("=" * 72)
    print(f" log: {log_path}")

    cathode.run()
    gun.run()
    prebuncher.run()

    _final_beam_summary(prebuncher.resolve_outdir())

    _cl("\n" + "=" * 72)
    _cl(f" Pipeline complete in {(time.time()-t0)/60:.1f} min.")
    _cl(" Figures: cathode/results/, gun/results/, prebuncher/results/")
    _cl("=" * 72)
    print(f" log: {log_path}")


if __name__ == "__main__":
    main()
