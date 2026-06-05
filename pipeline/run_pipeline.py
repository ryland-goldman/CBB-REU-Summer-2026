"""End-to-end Cornell Linac beam simulation in WarpX.

Runs the three stages in order in a single Python process:

    cathode -> gun -> prebuncher

Each stage exposes `config(**kwargs)` (override module-level constants) and `run()`
(build field map if any, simulate, plot). Stage execution is wrapped by the shared
runner in `pipeline/_runner.py`, which captures WarpX's stdout for tqdm progress
bars and writes a full DEBUG-level transcript to `pipeline/logs/pipeline_<ts>.log`.

Run with:
    conda activate CBB
    python pipeline/run_pipeline.py
"""

import os
import sys

# Set OMP_NUM_THREADS BEFORE any pywarpx import (read by OpenMP at load time).
os.environ.setdefault("OMP_NUM_THREADS", os.environ.get("OMP_THREADS", "6"))

# Run from the repo root so each stage's hard-coded relative paths resolve.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(_ROOT)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import cathode
import gun
import prebuncher

from pipeline._runner import setup_logging, _cl, _BOLD, _RESET

# ── Operating-point overrides (optional; defaults live in the stage modules) ──
cathode.config(V_anode=50.0)
gun.config(GUN_VOLTAGE=150e3, BUNCH_CHARGE=0.1e-9)
prebuncher.config(POWER_W=800, PHASE="zc", OUTDIR="prebuncher/diags")


def _final_beam_summary(diag):
    """Report the final bunch from the last prebuncher snapshot (console + log)."""
    try:
        import numpy as np
        from openpmd_viewer import OpenPMDTimeSeries
        ts = OpenPMDTimeSeries(os.path.join(diag, "particles"))
        for it in reversed(ts.iterations):
            z, ux, uy, uz, w = ts.get_particle(
                ["z", "ux", "uy", "uz", "w"], species="electrons", iteration=it)
            if len(z) > 50:
                break
        ke = (np.sqrt(1 + ux**2 + uy**2 + uz**2) - 1) * 0.51099895e3
        zm = np.average(z, weights=w)
        sz = np.sqrt(np.average((z - zm) ** 2, weights=w))
        km = np.average(ke, weights=w)
        dk = np.sqrt(np.average((ke - km) ** 2, weights=w))
        _cl(f"\n{_BOLD}Final beam{_RESET} (prebuncher exit, {len(z)} macroparticles):")
        _cl(f"      ⟨z⟩ = {zm*1e3:.0f} mm   σ_z = {sz*1e3:.3f} mm   "
            f"⟨KE⟩ = {km:.1f} keV   σ_KE = {dk:.2f} keV   "
            f"q = {w.sum()*1.602e-19*1e9:.3f} nC")
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

    from prebuncher import prebuncher_sim
    _final_beam_summary(prebuncher_sim.OUTDIR)

    _cl("\n" + "=" * 72)
    _cl(f" Pipeline complete in {(time.time()-t0)/60:.1f} min.")
    _cl(" Figures: cathode/results/, gun/results/, prebuncher/results/")
    _cl("=" * 72)
    print(f" log: {log_path}")


if __name__ == "__main__":
    main()
