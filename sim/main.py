"""End-to-end pipeline driver: cathode -> gun -> injector -> linac1/2/3 -> linac4-8,
then the cross-stage figures.

Each stage runs as a fresh subprocess (pywarpx binds one geometry per interpreter, so the
WarpX stages MUST be isolated; Impact-T runs the same way for uniformity). A stage's sim and
its plotter are separate subprocess calls. The Impact-T stage and the plotters pipe stdout to
logs/pipeline/log_<date>.log with the bar on stderr; the WarpX stages keep stdout on the terminal
(lume-warpx puts its bar on fd1 and routes the engine output to the log via PIPELINE_LOG_PATH).

All tuning lives in the per-stage config/*.yaml -- there is no config() override layer.
Run from the repo root:  python sim/main.py
"""

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import subprocess
import time

from sim.helpers.tools import REPO_ROOT, MC2_EV, E_CHARGE

# (label, sim script, plot script, extra args, exit-diag dir, KE unit)
STAGES = [
    ("cathode",  "sim/cathode.py",  "sim/plot/cathode.py",  [],    "logs/diags/cathode",            "keV"),
    ("gun",      "sim/gun.py",      "sim/plot/gun.py",      [],    "logs/diags/gun",                "keV"),
    ("injector", "sim/injector.py", "sim/plot/injector.py", [],    "logs/diags/injector/main",      "keV"),
    ("linac1",   "sim/linac1-3.py", "sim/plot/linac1-3.py", ["1"], "logs/diags/linac1-3/sec1/main", "MeV"),
    ("linac2",   "sim/linac1-3.py", "sim/plot/linac1-3.py", ["2"], "logs/diags/linac1-3/sec2/main", "MeV"),
    ("linac3",   "sim/linac1-3.py", "sim/plot/linac1-3.py", ["3"], "logs/diags/linac1-3/sec3/main", "MeV"),
    ("linac4-8", "sim/linac4-8.py", "sim/plot/linac4-8.py", [],    "logs/diags/linac4-8/main",      "MeV"),
]

_lf = None


def say(msg=""):
    """Print to the terminal and append to the run log."""
    print(msg, flush=True)
    if _lf is not None:
        _lf.write(msg + "\n")
        _lf.flush()


def run_subprocess(argv, title, fatal=True, warpx=False):
    """Run one stage step as a subprocess: stdout -> log, stderr (progress bar) -> terminal.

    `fatal=True` (the simulations) aborts the pipeline on a non-zero exit; `fatal=False` (the
    plotters) only warns -- a figure bug must not discard the completed physics, which is on disk.

    `warpx=True`: lume-warpx puts its bar on a dup of fd1 (disabled if that fd is not a tty) and
    redirects the engine's stdout/stderr to PIPELINE_LOG_PATH itself during the step -- so leave the
    child's stdout on the terminal (fd1 must stay a tty for the bar) and hand it the log path instead
    of piping stdout to the log here (which would make fd1 a file and silently disable the bar).
    """
    say(f"\n> {title}")
    if _lf is not None:
        _lf.flush()
    env = dict(os.environ)
    env.setdefault("OMP_NUM_THREADS", env.get("OMP_THREADS", "1"))
    env["HDF5_USE_FILE_LOCKING"] = "FALSE"
    env["PYTHONPATH"] = REPO_ROOT + os.pathsep + env.get("PYTHONPATH", "")
    out = None if (warpx and _lf is not None) else (_lf or None)
    if warpx and _lf is not None:
        env["PIPELINE_LOG_PATH"] = os.path.abspath(_lf.name)
    t0 = time.time()
    rc = subprocess.run([sys.executable, *argv], cwd=REPO_ROOT, env=env,
                        stdout=out, stderr=None).returncode
    dt = time.time() - t0
    flag = "ok" if rc == 0 else f"FAILED (exit {rc})"
    say(f"  {flag}  {title}  ({dt:.1f} s)")
    if rc != 0 and fatal:
        raise RuntimeError(f"{title} exited with code {rc}")
    return rc == 0


def beam_summary(diag, label, unit="keV"):
    """Report the final bunch from `diag` (prefers the true injected charge sidecar for capture %)."""
    try:
        import json
        import numpy as np
        from openpmd_viewer import OpenPMDTimeSeries
        ts = OpenPMDTimeSeries(os.path.join(diag, "particles"))
        its = list(ts.iterations)
        # "end-to-end" denominator = full upstream injected charge, so the % folds BOTH the
        # captured-core cut and in-transit loss (the per-stage sidecars split them via core_frac).
        q0, q0_label = None, "end-to-end capture"
        summ = os.path.join(diag, "injection_summary.json")
        sdata = json.load(open(summ)) if os.path.isfile(summ) else {}
        if sdata.get("q_injected_C") is not None:    # linac sidecars carry it; injector's does not
            q0 = sdata["q_injected_C"] / E_CHARGE
        elif its:
            q0_label = "transmitted"
            _, _, _, _, w0 = ts.get_particle(["z", "ux", "uy", "uz", "w"],
                                             species="electrons", iteration=its[0])
            q0 = w0.sum()
        z = None
        for it in reversed(its):
            z, ux, uy, uz, w = ts.get_particle(["z", "ux", "uy", "uz", "w"],
                                               species="electrons", iteration=it)
            if len(z) > 50:
                break
        if z is None or len(z) <= 50:
            say(f"  ({label}: no snapshot with >50 macroparticles)")
            return
        fac = MC2_EV / 1e3 if unit == "keV" else MC2_EV / 1e6
        ke = (np.sqrt(1 + ux**2 + uy**2 + uz**2) - 1) * fac
        zm = np.average(z, weights=w)
        sz = np.sqrt(np.average((z - zm) ** 2, weights=w))
        km = np.average(ke, weights=w)
        dk = np.sqrt(np.average((ke - km) ** 2, weights=w))
        cap = f"   {q0_label} {w.sum()/q0*100:.0f}%" if q0 else ""
        say(f"  {label}: <z> {zm*1e3:.0f} mm   sigma_z {sz*1e3:.3f} mm   "
            f"<KE> {km:.1f} {unit}   sigma_KE {dk:.2f} {unit}   "
            f"q {w.sum()*E_CHARGE*1e9:.3f} nC{cap}")
    except Exception as e:
        say(f"  ({label} summary unavailable: {e})")


def main():
    global _lf
    os.chdir(REPO_ROOT)
    os.makedirs("logs/pipeline", exist_ok=True)
    log_path = os.path.join("logs", "pipeline", time.strftime("log_%Y%m%d_%H%M%S.log"))
    _lf = open(log_path, "a", buffering=1, encoding="utf-8")

    t0 = time.time()
    say("=" * 72)
    say(" Cornell Linac pipeline:  cathode -> gun -> injector -> linac1/2/3 -> linac4-8")
    say(f" log: {log_path}")
    say("=" * 72)

    for label, sim, plot, args, _diag, _unit in STAGES:
        run_subprocess([sim, *args], f"{label}: simulation", warpx=(label != "linac4-8"))
        run_subprocess([plot, *args], f"{label}: plots", fatal=False)

    say("\n" + "-" * 72)
    for label, _sim, _plot, _args, diag, unit in STAGES:
        if label in ("cathode", "gun"):
            continue                       # source/low-energy: capture % not meaningful
        beam_summary(diag, f"{label} exit", unit)

    say("\n" + "=" * 72)
    say(f" Pipeline complete in {(time.time()-t0)/60:.1f} min.")
    say(" Figures: logs/plots/<stage>/")
    say(f" log: {log_path}")
    say("=" * 72)
    _lf.close()


if __name__ == "__main__":
    main()
