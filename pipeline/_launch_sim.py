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
import resource
import sys

# Raise the open-file-descriptor soft limit so the post-run handoff report can
# loop over the full diag series without exhausting macOS's default 256 — see
# _runner._raise_fd_limit for the full rationale (openpmd-viewer leaks an fd per
# get_particle, so the open fails with "OPEN_FILE failed ... Inaccessible" at a
# fixed dump count). A stage launched via the parent inherits the parent's
# raised limit; set it here too for a directly-invoked `python -m
# pipeline._launch_sim`. The disabled HDF5 lock is a minor extra mitigation.
try:
    _soft, _hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    _want = min(_hard, max(_soft, 16384))
    if _want > _soft:
        resource.setrlimit(resource.RLIMIT_NOFILE, (_want, _hard))
except (ValueError, OSError):
    pass
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")


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

    # Attach this child's `pipeline` logger to the parent's log file (append) so
    # run_step's `progress:` lines and any other child log.* calls land in the
    # shared log. Without this the child logger has no handler and those records
    # are silently dropped (run_step only ever executes here, in the child).
    log_path = os.environ.get("PIPELINE_LOG_PATH")
    if log_path:
        from pipeline._runner import setup_logging
        setup_logging(log_path)

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
