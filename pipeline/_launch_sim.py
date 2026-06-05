"""Subprocess entry point for a single stage simulation.

Invoked by `Stage._run_sim_subprocess` as

    python -m pipeline._launch_sim <sim_module_dotted_path>

with the optional `STAGE_CONFIG_JSON` env var carrying a JSON dict of parameter
overrides to apply by setattr on the sim module before its `main()` runs.

Each stage runs in its own fresh interpreter because pywarpx binds globally to
one geometry (2D/RZ/3D) at first .so load and caches diagnostic state per
diagnostic name — running cathode (2D) then gun (RZ) in one process trips
`AssertionError: Diagnostic attributes not consistent`. A subprocess per stage
sidesteps both.
"""

import importlib
import json
import os
import sys


def _silence_finalize():
    """Redirect fd 1/2 to the pipeline log file (or /dev/null).

    Called synchronously right after the sim's `main()` returns — once the
    child has nothing else user-visible to print, AMReX's TinyProfiler /
    memory tables (dumped during pywarpx teardown) end up in the log
    instead of the parent's terminal. atexit doesn't work here because
    pywarpx registered its finalize hook first; LIFO order means AMReX
    would dump before our atexit got a chance to redirect.
    """
    target = os.environ.get("PIPELINE_LOG_PATH") or os.devnull
    try:
        fd = os.open(target, os.O_WRONLY | os.O_APPEND)
        try:
            os.dup2(fd, 1)
            os.dup2(fd, 2)
        finally:
            os.close(fd)
    except Exception:
        pass


def main():
    if len(sys.argv) < 2:
        print("usage: python -m pipeline._launch_sim <dotted.sim.module>",
              file=sys.stderr)
        sys.exit(2)

    cwd = os.getcwd()
    if cwd not in sys.path:
        sys.path.insert(0, cwd)

    sim_module_path = sys.argv[1]
    params = json.loads(os.environ.get("STAGE_CONFIG_JSON", "{}"))

    sim = importlib.import_module(sim_module_path)
    for key, value in params.items():
        if hasattr(sim, key):
            setattr(sim, key, value)

    # Run sim.main() guarded so AMReX/pywarpx teardown chatter (TinyProfiler,
    # memory tables) is redirected to the log on the failure path too — Python
    # interpreter shutdown will trigger AMReX finalize next either way, and we
    # don't want it on the parent's terminal between stages.
    try:
        sim.main()
    finally:
        try: sys.stdout.flush()
        except Exception: pass
        try: sys.stderr.flush()
        except Exception: pass
        _silence_finalize()


if __name__ == "__main__":
    main()
