"""End-to-end Cornell Linac beam simulation in WarpX.

Orchestrates the four stages in order from one driver process:

    cathode -> gun -> injector -> linac_sec1

Each stage's build/plot run in-process, but its WarpX simulation runs in a fresh
Python subprocess (`pipeline._launch_sim`) to sidestep pywarpx's per-process
geometry binding (cathode 2D -> gun/injector/linac_sec1 RZ would otherwise trip a
diagnostic state assertion). Each stage exposes `config(**kwargs)` (override module-level
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
# An explicitly-set OMP_THREADS always wins, even if OMP_NUM_THREADS was already
# exported (otherwise setdefault would silently ignore OMP_THREADS); else default 1.
if "OMP_THREADS" in os.environ:
    os.environ["OMP_NUM_THREADS"] = os.environ["OMP_THREADS"]
else:
    os.environ.setdefault("OMP_NUM_THREADS", "1")

# Disable HDF5 file locking BEFORE any openPMD/h5py import. HDF5 reads this env
# var once at library init (first HDF5 call) — the stage imports below trigger
# that, so setting it later (e.g. inside _runner._prepare_environment) is too
# late for this parent process and the in-process plot step still fails. Without
# it, the post-run handoff report and the plot step hit
# "IO Task OPEN_FILE failed ... Inaccessible" opening openPMD .h5 files WarpX has
# only just flushed/closed — the files are intact; macOS HDF5's default locking
# just refuses them in that window. Reads never need the lock. _launch_sim.py
# sets the same flag so the spawned sim subprocess honors it from its own start.
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

# Run from the repo root so each stage's hard-coded relative paths resolve.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(_ROOT)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import cathode
import gun
import injector
import linac_sec1

from pipeline._runner import setup_logging, _cl, _BOLD, _RESET

# ── Operating-point overrides (physics; defaults live in the stage modules) ──
# Matched to the original LinacSim input files (reference/Linac Simulation
# Documentation/input_files/): cathode_master.in + gpt_master.in (gun, injector,
# section-1) GUI defaults.
cathode.config(V_anode=60.0)                          # Vpulse = 60 V
gun.config(GUN_VOLTAGE=150e3, BUNCH_CHARGE=1.0e-9)    # total_charge = -1e-9 (1 nC)
injector.config(PREB1_KW=8, PREB2_KW=10, PHASE="crest")  # preb1=8 kW, preb2=10 kW; crest base + GUI phi_off
# Sol 0 (40 A) and Lens 0A/0E live on the injector now (I_SOL0/I_LENS0A/I_LENS0E, faithful
# defaults in injector_sim) — the linac no longer has a solenoid, only RF power + phase.
linac_sec1.config(POWER_MW=11.0)                      # sec1_input_power = 11 MW (PHASE_DEG=0 default)

# ── Performance knobs (accuracy ↔ speed). Full knob list, runtime split, and the
#    reason the injector NZ must stay at 1664: see pipeline/README.md § Configuration. ──
# Balanced profile: ACTIVE (~1.7×, ~5 min). Comment these 3 lines for the baseline.
cathode.config(PPC=6, REQUIRED_PRECISION=3e-5)
gun.config(nz=256, MAX_PART=50000, REQUIRED_PRECISION=1e-4)
injector.config(CFL=0.95, MAX_ITERS=150, REQUIRED_PRECISION=1e-3)
# Conservative (~1.3×, near-identical):
# gun.config(MAX_PART=80000, REQUIRED_PRECISION=1e-4)
# injector.config(REQUIRED_PRECISION=2e-4, MAX_ITERS=400)
# Aggressive (~2.2×, looser space-charge solve):
# cathode.config(nz=48, PPC=4, REQUIRED_PRECISION=5e-5, MAX_STEPS=1200)
# gun.config(nz=192, MAX_PART=40000, REQUIRED_PRECISION=3e-4, N_DIAGS=20)
# injector.config(CFL=0.97, MAX_ITERS=80, REQUIRED_PRECISION=3e-3, N_DIAGS=20)
# linac_sec1.config(NZ=1024, CFL=0.6)   # coarser/faster linac run (default NZ=1664, ~40 s)


def _beam_summary(diag, label, unit="keV"):
    """Report the final bunch from the last snapshot of `diag` (console + log).

    `unit` selects the kinetic-energy scale ("keV" for the injector exit, "MeV"
    for the linac exit). Also reports the captured-charge fraction (last/first
    snapshot) when both are available.
    """
    try:
        import numpy as np
        from openpmd_viewer import OpenPMDTimeSeries
        ts = OpenPMDTimeSeries(os.path.join(diag, "particles"))
        its = list(ts.iterations)
        # Capture denominator: prefer the TRUE injected charge the sim records in
        # injection_summary.json (the linac drops r>RMAX particles before the first dump, so
        # the first-dump charge already hides the injection loss). Fall back to the first dump
        # for stages without a sidecar (e.g. the injector exit).
        q0 = None
        summ_path = os.path.join(diag, "injection_summary.json")
        if os.path.isfile(summ_path):
            import json
            with open(summ_path) as fh:
                q0 = json.load(fh)["q_injected_C"] / 1.602176634e-19
        elif its:
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
    _cl(" Cornell Linac WarpX pipeline:  cathode -> gun -> injector -> linac_sec1")
    _cl(f" OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS', '?')}")
    _cl("=" * 72)
    print(f" log: {log_path}")

    cathode.run()
    gun.run()
    injector.run()
    linac_sec1.run()            # SLAC Section 1: ~26 MeV captured at 11 MW (~18% of true injected; γ² lower bound)

    _beam_summary(injector.resolve_outdir(), "injector exit", "keV")
    _beam_summary(linac_sec1.resolve_outdir(), "linac_sec1 exit", "MeV")

    # Cross-stage figures (in-process, no pywarpx): one moment table per stage →
    # chain_evolution / emittance_budget / transmission_waterfall / scorecard in results/.
    try:
        import pipeline.plot_chain as plot_chain   # submodule, not the pipeline.plot_chain() fn
        plot_chain.main()
    except Exception as e:
        import logging
        _cl(f"    (cross-stage figures unavailable: {e})", level=logging.WARNING)

    _cl("\n" + "=" * 72)
    _cl(f" Pipeline complete in {(time.time()-t0)/60:.1f} min.")
    _cl(" Figures: cathode/results/, gun/results/, injector/results/, linac_sec1/results/, "
        "results/ (cross-stage)")
    _cl("=" * 72)
    print(f" log: {log_path}")


if __name__ == "__main__":
    main()
