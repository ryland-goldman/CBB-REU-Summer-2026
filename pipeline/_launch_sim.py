"""Fresh-interpreter entry point for one stage sim (`python -m pipeline._launch_sim <module>`).

Applies `STAGE_CONFIG_JSON` overrides by setattr, then runs the module's `main()`.
See pipeline/README.md for the subprocess-isolation rationale and the run flow.
"""

import importlib
import json
import os
import resource
import sys

# Raise the fd soft limit (openpmd-viewer leaks an fd per get_particle); set here too
# for a directly-invoked launch. See _runner._raise_fd_limit.
try:
    _soft, _hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    _want = min(_hard, max(_soft, 16384))
    if _want > _soft:
        resource.setrlimit(resource.RLIMIT_NOFILE, (_want, _hard))
except (ValueError, OSError):
    pass
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")


def _silence_finalize():
    """Redirect fd 1/2 to the pipeline log (or /dev/null) so AMReX/pywarpx teardown
    chatter lands in the log, not the parent's terminal.

    Called synchronously after main() returns — atexit can't be used, pywarpx
    registered its finalize hook first and LIFO order would dump it before ours ran.
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

    # Attach this child's `pipeline` logger to the parent's log file; without it the
    # child logger has no handler and run_step's records are silently dropped.
    log_path = os.environ.get("PIPELINE_LOG_PATH")
    if log_path:
        from pipeline._runner import setup_logging
        setup_logging(log_path)

    sim = importlib.import_module(sim_module_path)
    for key, value in params.items():
        if hasattr(sim, key):
            setattr(sim, key, value)

    # finally-guard the redirect so teardown chatter is silenced on the failure path too.
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
