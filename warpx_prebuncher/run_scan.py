"""
Driver for the WarpX prebuncher study: sweep RF power and compare phases.

The two undocumented prebuncher operating-point inputs are the dissipated power P
(which sets the field scale via scale = sqrt(1e3·Q·P / (2π f_RF)), Q = 3000) and
the RF phase. Per the plan we scan several powers and run BOTH the zero-crossing
(velocity-bunching) and on-crest phases for each, so `plot_prebuncher.py` can
show bunching factor / min bunch length vs. power for both phases.

Each (power, phase) case is a separate WarpX run — AMReX/WarpX cannot be cleanly
re-initialised within one process. The cases are run CONCURRENTLY: WarpX's MLMG
Poisson solve (the per-step bottleneck here) scales poorly past a few OpenMP
threads, so running several cases at once with fewer threads each uses the 14
cores far better than one case on all threads. `WORKERS` cases run at a time, each
pinned to `THREADS_PER` OpenMP threads (via OMP_NUM_THREADS); each case streams to
its own `diags/P{power}_{phase}.log`.

Run with:
    conda run -n CBB python warpx_prebuncher/run_scan.py
Override concurrency:
    WORKERS=5 THREADS_PER=2 conda run -n CBB python warpx_prebuncher/run_scan.py
"""

import os
import subprocess
import sys
import time
import numpy as np

# Map-derived constants (see build_prebuncher_field.py output).
F_RF = 499.7645e6 / 42 * 18      # 214.18 MHz
Q_L = 3000
V1J = 430.2e3                    # 1-J transit-weighted effective gap voltage [V]

# The gun beam carries an intrinsic +1.40 keV/mm (debunching) energy chirp; the
# zero-crossing cavity adds −3.05 keV/mm per unit field scale, so the net chirp is
# 1.40 − 3.05·scale and bunching needs scale ≳ 0.46 (≈95 W, V_gap ≈ 205 kV). These
# powers span weak→strong bunching with foci ≈ 1.16, 0.60, 0.44, 0.36 m from the
# gap (z = 0.20 m), all inside the 1.30 m domain.
POWERS = [160.0, 300.0, 500.0, 800.0]   # [W]
PHASES = ["zc", "crest"]

DIAG_ROOT = "warpx_prebuncher/diags"

# Concurrency: WORKERS cases at once, THREADS_PER OpenMP threads each.
try:
    NCORE = len(os.sched_getaffinity(0))       # Linux
except AttributeError:
    NCORE = os.cpu_count() or 8
WORKERS = int(os.environ.get("WORKERS", min(4, max(1, NCORE // 3))))
THREADS_PER = int(os.environ.get("THREADS_PER", max(1, NCORE // WORKERS)))


def scale_of(power):
    return float(np.sqrt(1e3 * Q_L * power / (2.0 * np.pi * F_RF)))


def main():
    print(f"f_RF = {F_RF/1e6:.2f} MHz, Q = {Q_L}, "
          f"1-J effective gap voltage = {V1J/1e3:.0f} kV")
    print(f"concurrency: {WORKERS} cases × {THREADS_PER} OpenMP threads "
          f"({NCORE} cores)\n")
    print(f"{'P [W]':>6}  {'scale':>6}  {'V_gap [kV]':>10}")
    for power in POWERS:
        s = scale_of(power)
        print(f"{power:6.0f}  {s:6.3f}  {s*V1J/1e3:10.1f}")
    print()

    # P = 0 is the drift-only baseline (no cavity) used to isolate bunching.
    cases = [(0.0, "zc")] + [(p, ph) for p in POWERS for ph in PHASES]
    env = dict(os.environ, OMP_NUM_THREADS=str(THREADS_PER))
    running = {}        # Popen -> (tag, logfile handle)
    pending = list(cases)
    failed = []

    def launch(power, phase):
        tag = "P0_drift" if power == 0 else f"P{int(power)}_{phase}"
        outdir = os.path.join(DIAG_ROOT, tag)
        log = open(os.path.join(DIAG_ROOT, f"{tag}.log"), "w")
        os.makedirs(DIAG_ROOT, exist_ok=True)
        cmd = [sys.executable, "warpx_prebuncher/prebuncher_sim.py",
               "--power", str(power), "--phase", phase, "--outdir", outdir]
        print(f"  launch {tag}  (P={power:g} W, {phase})", flush=True)
        p = subprocess.Popen(cmd, env=env, stdout=log, stderr=subprocess.STDOUT)
        running[p] = (tag, log)

    while pending or running:
        while pending and len(running) < WORKERS:
            launch(*pending.pop(0))
        time.sleep(2.0)
        for p in list(running):
            if p.poll() is not None:
                tag, log = running.pop(p)
                log.close()
                status = "ok" if p.returncode == 0 else f"FAILED ({p.returncode})"
                if p.returncode != 0:
                    failed.append(tag)
                print(f"  done   {tag}  [{status}]", flush=True)

    print(f"\nAll cases finished. {'failures: ' + ', '.join(failed) if failed else 'no failures.'}")
    print("Plot with:\n    conda run -n CBB python warpx_prebuncher/plot_prebuncher.py")


if __name__ == "__main__":
    main()
