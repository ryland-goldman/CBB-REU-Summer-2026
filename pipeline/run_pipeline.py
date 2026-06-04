"""
End-to-end Cornell Linac beam simulation in WarpX.

Runs the full chain in order, each stage as a subprocess, with a live progress bar
and the key physics output for every component:

    1. Cathode diode         warpx_cathode/cathode_diode.py      (SCL emission)
    2. Gun field map          warpx_gun/build_gun_field.py        (CESR_gun.gdf -> openPMD)
    3. Gun acceleration       warpx_gun/gun_sim.py                (RZ, -> ~148 keV)
    4. Prebuncher field map   warpx_prebuncher/build_prebuncher_field.py
    5. Prebuncher bunching    warpx_prebuncher/prebuncher_sim.py  (RZ RF cavity)
    6. Plots (optional)       each stage's plot_*.py

Each downstream stage reads the previous stage's openPMD output, so the order is
fixed. Simulation parameters are set in the CONFIG block below; the prebuncher
operating point (power, phase) is passed through to `prebuncher_sim.py`.

Output:
  * Terminal — a live tqdm progress bar per simulation stage plus the key prints.
  * Log file — pipeline/logs/pipeline_<timestamp>.log: the same information with
    the progress bar replaced by periodic step/rate/ETA lines, per-stage durations
    and return codes, and extra DEBUG detail (full commands, environment, and the
    WarpX output normally hidden from the console). MLMG iteration spam is dropped.

Run with:
    conda activate CBB
    python pipeline/run_pipeline.py
"""

import logging
import os
import platform
import re
import subprocess
import sys
import time

from tqdm import tqdm

# Run from the repository root so the stage scripts' relative paths resolve.
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

# ════════════════════════════════════════════════════════════════════════════
#  CONFIG — simulation parameters
# ════════════════════════════════════════════════════════════════════════════
RUN_CATHODE = True            # stage 1
RUN_GUN = True                # stages 2–3
RUN_PREBUNCHER = True         # stages 4–5
MAKE_PLOTS = True             # stage 6

# Prebuncher operating point (the two undocumented inputs; see warpx_prebuncher).
PREBUNCHER_POWER_W = 800.0    # dissipated RF power [W]  (0 = drift-only baseline)
PREBUNCHER_PHASE = "zc"       # "zc" = zero-crossing bunching, "crest" = max energy gain

# OpenMP threads per WarpX run. The MLMG Poisson solve is memory-bandwidth bound,
# so all 14 cores is actually slower than a moderate count; 6 is a good default.
OMP_THREADS = int(os.environ.get("OMP_THREADS", "6"))
# ════════════════════════════════════════════════════════════════════════════

# ANSI styling, but only on an interactive terminal — otherwise (e.g. redirected
# to a log file) escape codes show up as literal "[1m"/"[0m" noise.
_TTY = sys.stdout.isatty()
BOLD = "\033[1m" if _TTY else ""
GREEN = "\033[32m" if _TTY else ""
YELLOW = "\033[33m" if _TTY else ""
RESET = "\033[0m" if _TTY else ""
ANSI = re.compile(r"\033\[[0-9;]*m")

# Lines worth surfacing on the console from the (otherwise noisy) WarpX stdout.
INFO = re.compile(
    r"Imported|Case:|Diode|Child|v_th|v_final|dt =|max_steps|Running|On-axis|"
    r"Gun:|Scaled|Wrote|Done|Prebuncher map|Peak|Cavity gap|renormalized|"
    r"γ|β|nC|keV|MV/m|kV")
STEP = re.compile(r"STEP\s+(\d+)\s+starts")
TOTAL = re.compile(r"(?:max_steps\s*=\s*|Running\s+)(\d+)")
DROP = re.compile(r"MLMG: Iteration")     # convergence spam — keep out of the log too

log = logging.getLogger("pipeline")       # configured in setup_logging()


def setup_logging():
    """Create a timestamped DEBUG log file under pipeline/logs/ and return its path."""
    logdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "logs")
    os.makedirs(logdir, exist_ok=True)
    path = os.path.join(logdir, time.strftime("pipeline_%Y%m%d_%H%M%S.log"))
    log.setLevel(logging.DEBUG)
    log.propagate = False
    log.handlers.clear()
    fh = logging.FileHandler(path, mode="w", encoding="utf-8")
    fh.setFormatter(logging.Formatter("%(asctime)s  %(levelname)-5s  %(message)s", "%H:%M:%S"))
    log.addHandler(fh)
    return path


def cl(msg="", level=logging.INFO):
    """Print to the console (with ANSI) and log the plain text to the file."""
    print(msg)
    log.log(level, ANSI.sub("", msg).rstrip() or " ")


def run_stage(idx, n, title, script, args, is_sim):
    """Run one stage as a subprocess; tqdm bar on the console, step lines in the log."""
    tag = f"[{idx}/{n}]"
    cl(f"\n{BOLD}{tag} {title}{RESET}")
    cmd = [sys.executable, script, *args]
    cl(f"    $ {os.path.basename(script)} {' '.join(args)}".rstrip())
    env = dict(os.environ, OMP_NUM_THREADS=str(OMP_THREADS))
    log.debug(f"    cmd: {' '.join(cmd)}")
    log.debug(f"    OMP_NUM_THREADS={OMP_THREADS}  cwd={os.getcwd()}")
    t0 = time.time()

    proc = subprocess.Popen(cmd, env=env, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, bufsize=1)
    bar, total = None, None
    log_every = next_log = 0
    for line in proc.stdout:
        line = line.rstrip()
        if not line or DROP.search(line):
            continue
        if is_sim:
            m = TOTAL.search(line)
            if m and bar is None:
                total = int(m.group(1))
                log_every = max(1, total // 20)        # ~20 progress lines in the log
                next_log = log_every
                bar = tqdm(total=total, unit="step", desc=f"    {title}",
                           ncols=88, leave=True, dynamic_ncols=False, disable=not _TTY)
            s = STEP.search(line)
            if s and bar is not None:
                step = min(int(s.group(1)), total)
                bar.n = step
                bar.refresh()
                if step >= next_log:                    # periodic step/rate/ETA -> log only
                    el = time.time() - t0
                    rate = step / el if el else 0.0
                    eta = (total - step) / rate if rate else 0.0
                    log.info(f"    progress: step {step}/{total} ({100*step/total:3.0f}%)  "
                             f"elapsed {el:5.0f}s  {rate:5.1f} step/s  eta {eta:4.0f}s")
                    next_log += log_every
                continue
        # Informative physics lines -> console + log; other WarpX lines -> log DEBUG only.
        if (not is_sim) or INFO.search(line):
            (bar.write if bar else print)(f"      {line}")
            log.info(f"      {ANSI.sub('', line)}")
        else:
            log.debug(f"      {line}")
    proc.wait()
    if bar is not None:
        bar.n = total
        bar.refresh()
        bar.close()
    dt = time.time() - t0
    ok = proc.returncode == 0
    flag = f"{GREEN}✓{RESET}" if ok else f"{YELLOW}⚠ exit {proc.returncode}{RESET}"
    cl(f"    {flag}  {title}  ({dt:5.1f} s)")
    log.info(f"    stage {idx}/{n} return code = {proc.returncode}, duration = {dt:.1f} s")
    return ok


def final_beam_summary(diag):
    """Report the final bunch from the last prebuncher snapshot (console + log)."""
    try:
        import numpy as np
        from openpmd_viewer import OpenPMDTimeSeries
        ts = OpenPMDTimeSeries(os.path.join(diag, "particles"))
        for it in reversed(ts.iterations):           # last dump with a real beam
            z, ux, uy, uz, w = ts.get_particle(
                ["z", "ux", "uy", "uz", "w"], species="electrons", iteration=it)
            if len(z) > 50:
                break
        ke = (np.sqrt(1 + ux**2 + uy**2 + uz**2) - 1) * 0.51099895e3
        zm = np.average(z, weights=w)
        sz = np.sqrt(np.average((z - zm) ** 2, weights=w))
        km = np.average(ke, weights=w)
        dk = np.sqrt(np.average((ke - km) ** 2, weights=w))
        cl(f"\n{BOLD}Final beam{RESET} (prebuncher exit, {len(z)} macroparticles):")
        cl(f"      ⟨z⟩ = {zm*1e3:.0f} mm   σ_z = {sz*1e3:.3f} mm   "
           f"⟨KE⟩ = {km:.1f} keV   σ_KE = {dk:.2f} keV   q = {w.sum()*1.602e-19*1e9:.3f} nC")
    except Exception as e:
        cl(f"    (final-beam summary unavailable: {e})", level=logging.WARNING)


def main():
    log_path = setup_logging()

    stages = []
    if RUN_CATHODE:
        stages.append(("Cathode diode (space-charge-limited emission)",
                       "warpx_cathode/cathode_diode.py", [], True))
    if RUN_GUN:
        stages.append(("Gun field map (CESR_gun.gdf -> openPMD)",
                       "warpx_gun/build_gun_field.py", [], False))
        stages.append(("Gun acceleration (RZ, -> ~148 keV)",
                       "warpx_gun/gun_sim.py", [], True))
    if RUN_PREBUNCHER:
        stages.append(("Prebuncher field map (prebuncher_25D.gdf -> openPMD)",
                       "warpx_prebuncher/build_prebuncher_field.py", [], False))
        stages.append((f"Prebuncher RF bunching (P={PREBUNCHER_POWER_W:g} W, "
                       f"phase={PREBUNCHER_PHASE})",
                       "warpx_prebuncher/prebuncher_sim.py",
                       ["--power", str(PREBUNCHER_POWER_W),
                        "--phase", PREBUNCHER_PHASE], True))
    if MAKE_PLOTS:
        if RUN_CATHODE:
            stages.append(("Cathode plots", "warpx_cathode/plot_cathode.py", [], False))
        if RUN_GUN:
            stages.append(("Gun plots", "warpx_gun/plot_gun.py", [], False))
        if RUN_PREBUNCHER:
            stages.append(("Prebuncher plots", "warpx_prebuncher/plot_prebuncher.py", [], False))

    n = len(stages)
    cl("=" * 72)
    cl(" Cornell Linac WarpX pipeline:  cathode -> gun -> prebuncher")
    cl(f" {n} stages, OMP_NUM_THREADS={OMP_THREADS}, "
       f"prebuncher P={PREBUNCHER_POWER_W:g} W / {PREBUNCHER_PHASE}")
    cl("=" * 72)
    print(f" log: {log_path}")
    # Environment / config detail for the log only.
    log.debug(f"python {sys.version.split()[0]} on {platform.platform()}")
    log.debug(f"executable {sys.executable}")
    log.debug(f"conda env {os.environ.get('CONDA_DEFAULT_ENV', '?')}, cwd {os.getcwd()}")
    log.debug(f"config: cathode={RUN_CATHODE} gun={RUN_GUN} prebuncher={RUN_PREBUNCHER} "
              f"plots={MAKE_PLOTS}")

    t0 = time.time()
    failures = []
    for i, (title, script, args, is_sim) in enumerate(stages, 1):
        if not run_stage(i, n, title, script, args, is_sim):
            failures.append(title)

    if RUN_PREBUNCHER:
        tag = ("P0_drift" if PREBUNCHER_POWER_W == 0
               else f"P{int(PREBUNCHER_POWER_W)}_{PREBUNCHER_PHASE}")
        final_beam_summary(f"warpx_prebuncher/diags/{tag}")

    cl("\n" + "=" * 72)
    mins = (time.time() - t0) / 60
    if failures:
        cl(f" Pipeline finished in {mins:.1f} min with warnings: {', '.join(failures)}",
           level=logging.WARNING)
    else:
        cl(f" Pipeline complete in {mins:.1f} min — all {n} stages OK.")
    cl(" Figures: warpx_{cathode,gun,prebuncher}/results/")
    cl("=" * 72)
    print(f" log: {log_path}")


if __name__ == "__main__":
    main()
