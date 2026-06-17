"""End-to-end Cornell Linac pipeline: cathode -> gun -> injector -> linac_sec1 -> linac_rest.

Imports each stage facade and calls config()/run() in order, then renders the
cross-stage figures. See pipeline/README.md for the run command, configuration,
performance knobs, and the subprocess-isolation rationale.
"""

import os
import sys

# Set OMP_NUM_THREADS BEFORE any pywarpx import (OpenMP reads it at load time).
# Explicit OMP_THREADS wins over an already-exported OMP_NUM_THREADS; else default 1.
if "OMP_THREADS" in os.environ:
    os.environ["OMP_NUM_THREADS"] = os.environ["OMP_THREADS"]
else:
    os.environ.setdefault("OMP_NUM_THREADS", "1")

# Set before any openPMD/h5py import (HDF5 latches locking at library init).
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
import linac_rest

from pipeline._runner import setup_logging, _cl, _BOLD, _RESET

# ── Operating-point overrides (physics; defaults live in the stage modules) ──
cathode.config(V_anode=60.0)
gun.config(GUN_VOLTAGE=150e3, BUNCH_CHARGE=1.0e-9)
# Crest-referenced phase convention (base=π, kick ∝ −cos(base+φ_off)): null = φ_off −90°, crest = φ_off 0°.
injector.config(PREB1_KW=8, PREB2_KW=10, PREB1_PHI_OFF=-90, PREB2_PHI_OFF=0, PHASE="crest")
linac_sec1.config(POWER_MW=11.0)
linac_rest.config(POWER_MW=11.0)

# ── Space charge (beam self-field) per stage; see pipeline/README.md § Configuration.
#    Uncomment to override the per-stage defaults:
# cathode.config(SPACE_CHARGE=True)
# gun.config(SPACE_CHARGE=True)
# injector.config(SPACE_CHARGE=True)
# linac_sec1.config(SPACE_CHARGE=True)
# linac_rest.config(SPACE_CHARGE=True)   # turn ON Impact-T space charge (headline is OFF)

# ── Performance knobs (accuracy ↔ speed); see pipeline/README.md § Configuration. ──
# Balanced profile: ACTIVE. Comment these 3 lines for the baseline.
cathode.config(PPC=6, REQUIRED_PRECISION=3e-5)
gun.config(nz=384, MAX_PART=50000, REQUIRED_PRECISION=1e-4)
injector.config(CFL=0.95, MAX_ITERS=150, REQUIRED_PRECISION=1e-3)
# Conservative (~1.3×, near-identical):
# gun.config(MAX_PART=80000, REQUIRED_PRECISION=1e-4)
# injector.config(REQUIRED_PRECISION=2e-4, MAX_ITERS=400)
# Aggressive (~2.2×, looser space-charge solve):
# cathode.config(nz=48, PPC=4, REQUIRED_PRECISION=5e-5, MAX_STEPS=1200)
# gun.config(nz=192, MAX_PART=40000, REQUIRED_PRECISION=3e-4, N_DIAGS=20)
# injector.config(CFL=0.97, MAX_ITERS=80, REQUIRED_PRECISION=3e-3, N_DIAGS=20)
# linac_sec1.config(NZ=1024, CFL=0.6)   # coarser/faster linac run (default NZ=1664, ~40 s)
# Np = tracked macroparticle count; Ntstep sized for the ~36 m line.
linac_rest.config(Np=4000, Ntstep=200000)
# Exploratory FODO (headline stays quads OFF); see linac_rest/README.md. Leave commented.
# linac_rest.config(QUADS_ON=True)                     # exploratory FODO


def _beam_summary(diag, label, unit="keV"):
    """Report the final bunch from the last snapshot of `diag` (console + log).

    `unit` is the KE scale ("keV"/"MeV"). Charge fraction is "captured" vs the
    sidecar's true injected charge, else "transmitted" vs the first snapshot.
    """
    try:
        import numpy as np
        from openpmd_viewer import OpenPMDTimeSeries
        ts = OpenPMDTimeSeries(os.path.join(diag, "particles"))
        its = list(ts.iterations)
        # Prefer the TRUE injected charge from injection_summary.json: the linac drops
        # r>RMAX particles before the first dump, so first-dump charge hides injection loss.
        q0 = None
        q0_label = "captured"
        summ_path = os.path.join(diag, "injection_summary.json")
        if os.path.isfile(summ_path):
            import json
            with open(summ_path) as fh:
                q0 = json.load(fh)["q_injected_C"] / 1.602176634e-19
        elif its:
            # First-dump fallback: a within-stage transmission, not a capture fraction.
            q0_label = "transmitted"
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
        cap = f"   {q0_label} {w.sum()/q0*100:.0f}%" if q0 else ""
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
    _cl(" Cornell Linac pipeline:  cathode -> gun -> injector -> linac_sec1 -> linac_rest")
    _cl(f" OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS', '?')}")
    _cl("=" * 72)
    print(f" log: {log_path}")

    cathode.run()
    gun.run()
    injector.run()
    linac_sec1.run()
    linac_rest.run()

    _beam_summary(injector.resolve_outdir(), "injector exit", "keV")
    _beam_summary(linac_sec1.resolve_outdir(), "linac_sec1 exit", "MeV")
    _beam_summary(linac_rest.resolve_outdir(), "linac exit (8 sections)", "MeV")

    try:
        import pipeline.plot_chain as plot_chain   # submodule, not the pipeline.plot_chain() fn
        plot_chain.main()
    except Exception as e:
        import logging
        _cl(f"    (cross-stage figures unavailable: {e})", level=logging.WARNING)

    _cl("\n" + "=" * 72)
    _cl(f" Pipeline complete in {(time.time()-t0)/60:.1f} min.")
    _cl(" Figures: cathode/results/, gun/results/, injector/results/, linac_sec1/results/, "
        "linac_rest/results/, results/ (cross-stage)")
    _cl("=" * 72)
    print(f" log: {log_path}")


if __name__ == "__main__":
    main()
