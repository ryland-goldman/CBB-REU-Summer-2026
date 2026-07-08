"""End-to-end pipeline driver: cathode -> gun -> injector -> linac1/2/3/4 -> converter -> linac5-8,
then the cross-stage figures. See README.md.

Run from the repo root: python sim/main.py
"""

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import subprocess
import time

from sim.helpers.tools import REPO_ROOT, MC2_EV, E_CHARGE, out_root
from sim.helpers.sandbox import make_out_dir

# (label, sim script, plot script, extra args, autophase argv (or None), exit-diag dir, KE unit)
STAGES = [
    ("cathode",  "sim/cathode.py",  "sim/plot/cathode.py",  [],    None,                          "logs/diags/cathode",            "keV"),
    ("gun",      "sim/gun.py",      "sim/plot/gun.py",      [],    None,                          "logs/diags/gun",                "keV"),
    ("injector", "sim/injector.py", "sim/plot/injector.py", [],    None,                          "logs/diags/injector/main",      "keV"),
    ("linac1",   "sim/linac1-4.py", "sim/plot/linac1-4.py", ["1"], ["sim/autophase.py", "1"],     "logs/diags/linac1-4/sec1/main", "MeV"),
    ("linac2",   "sim/linac1-4.py", "sim/plot/linac1-4.py", ["2"], ["sim/autophase.py", "2"],     "logs/diags/linac1-4/sec2/main", "MeV"),
    ("linac3",   "sim/linac1-4.py", "sim/plot/linac1-4.py", ["3"], ["sim/autophase.py", "3"],     "logs/diags/linac1-4/sec3/main", "MeV"),
    ("linac4",   "sim/linac1-4.py", "sim/plot/linac1-4.py", ["4"], ["sim/autophase.py", "4"],     "logs/diags/linac1-4/sec4/main", "MeV"),
    ("converter","sim/converter.py","sim/plot/converter.py", [],    None,                          "logs/diags/converter/main",     "MeV"),
    ("linac5-8", "sim/linac5-8.py", "sim/plot/linac5-8.py", [],    ["sim/autophase_impact.py"],   "logs/diags/linac5-8/main",      "MeV"),
]

_lf = None


def say(msg=""):
    print(msg, flush=True)
    if _lf is not None:
        _lf.write(msg + "\n")
        _lf.flush()


def run_subprocess(argv, title, fatal=True, warpx=False):
    """`fatal=False` (plotters) only warns -- a figure bug must not discard completed physics on disk.

    `warpx=True`: lume-warpx puts its bar on a dup of fd1 (disabled if not a tty) and redirects the
    engine's stdout/stderr to PIPELINE_LOG_PATH itself -- leave the child's stdout on the terminal
    (fd1 must stay a tty for the bar) instead of piping it to the log here (that would disable the bar).
    """
    say(f"\n> {title}")
    if _lf is not None:
        _lf.flush()
    out_dir = out_root()
    env = dict(os.environ)
    env.setdefault("OMP_NUM_THREADS", env.get("OMP_THREADS", "1"))
    env["HDF5_USE_FILE_LOCKING"] = "FALSE"
    env["PYTHONPATH"] = REPO_ROOT + os.pathsep + env.get("PYTHONPATH", "")
    if out_dir != REPO_ROOT:                    # only a real sandbox; a plain run leaves it unset so
        env["LINACSIM_OUT_DIR"] = out_dir       # the fieldmap skip-guard stays off (truncate-rebuild)
    out = None if (warpx and _lf is not None) else (_lf or None)
    if warpx and _lf is not None:
        env["PIPELINE_LOG_PATH"] = os.path.abspath(_lf.name)
    t0 = time.time()
    # argv[0] is a repo-relative script path; resolve it against REPO_ROOT (the code lives there),
    # but run with cwd=out_dir so config/logs I/O lands in the sandbox.
    cmd = [sys.executable, os.path.join(REPO_ROOT, argv[0]), *argv[1:]]
    rc = subprocess.run(cmd, cwd=out_dir, env=env,
                        stdout=out, stderr=None).returncode
    dt = time.time() - t0
    flag = "ok" if rc == 0 else f"FAILED (exit {rc})"
    say(f"  {flag}  {title}  ({dt:.1f} s)")
    if rc != 0 and fatal:
        raise RuntimeError(f"{title} exited with code {rc}")
    return rc == 0


def _fmt_charge(q):
    """Charge with an auto-scaled SI prefix (nC/pC/fC) so sub-nC bunches keep their figures."""
    aq = abs(q)
    if aq >= 1e-9:
        return f"{q*1e9:.3f} nC"
    if aq >= 1e-12:
        return f"{q*1e12:.3f} pC"
    return f"{q*1e15:.3f} fC"


def beam_summary(diag, label, unit="keV"):
    """Report the final bunch from `diag` (prefers the true injected charge sidecar for capture %)."""
    try:
        import json
        import numpy as np
        from openpmd_viewer import OpenPMDTimeSeries
        ts = OpenPMDTimeSeries(os.path.join(diag, "particles"))
        sp = ts.avail_species[0] if ts.avail_species else "electrons"   # converter/linac5-8 write "positrons"
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
                                             species=sp, iteration=its[0])
            q0 = w0.sum()
        z = None
        for it in reversed(its):
            z, ux, uy, uz, w = ts.get_particle(["z", "ux", "uy", "uz", "w"],
                                               species=sp, iteration=it)
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
            f"q {_fmt_charge(w.sum()*E_CHARGE)}{cap}")
    except Exception as e:
        say(f"  ({label} summary unavailable: {e})")


def _select_stages(argv):
    """Slice STAGES by optional `--from <label>` (start) and `--to <label>` (inclusive end) -- a
    partial chain that skips frozen upstream stages (and their autophase) and/or stops early. No
    flags => the full chain (unchanged)."""
    labels = [s[0] for s in STAGES]
    lo, hi = 0, len(STAGES)
    if "--from" in argv:
        start = argv[argv.index("--from") + 1]
        if start not in labels:
            sys.exit(f"--from: unknown stage {start!r} (choices: {', '.join(labels)})")
        lo = labels.index(start)
    if "--to" in argv:
        end = argv[argv.index("--to") + 1]
        if end not in labels:
            sys.exit(f"--to: unknown stage {end!r} (choices: {', '.join(labels)})")
        hi = labels.index(end) + 1
    return STAGES[lo:hi]


def main():
    global _lf
    stages = _select_stages(sys.argv)
    no_plots = "--no-plots" in sys.argv          # optimizer runs read the diags, never the figures
    out_dir = out_root()
    # Build a not-yet-populated standalone sandbox; defer to a caller (Xopt) that already wrote
    # overrides into <out_dir>/config (an unconditional copytree would wipe them).
    if out_dir != REPO_ROOT and not os.path.isdir(f"{out_dir}/config"):
        make_out_dir(out_dir)
    os.chdir(out_dir)
    os.makedirs("logs/pipeline", exist_ok=True)
    log_path = os.path.join("logs", "pipeline", time.strftime("log_%Y%m%d_%H%M%S.log"))
    _lf = open(log_path, "a", buffering=1, encoding="utf-8")

    t0 = time.time()
    say("=" * 72)
    say(" Cornell Linac pipeline:  cathode -> gun -> injector -> linac1/2/3/4 -> converter -> linac5-8")
    say(f" log: {log_path}")
    say("=" * 72)

    for label, sim, plot, args, autophase, _diag, _unit in stages:
        if autophase:
            # Fatal: a stale/garbage crest would silently invalidate the stage.
            run_subprocess(autophase, f"{label}: autophase")
        run_subprocess([sim, *args], f"{label}: simulation", warpx=(label not in ("linac5-8", "converter")))
        if not no_plots:
            run_subprocess([plot, *args], f"{label}: plots", fatal=False)

    say("\n" + "-" * 72)
    for label, _sim, _plot, _args, _ap, diag, unit in stages:
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
