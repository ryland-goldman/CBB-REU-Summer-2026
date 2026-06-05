"""Stage shim + progress/log helpers shared by the cathode/gun/prebuncher facades.

Each top-level package (cathode/, gun/, prebuncher/) instantiates a `Stage` in
its __init__.py and re-exports `config`, `run`, `plot`. The Stage object:

  * Sets `OMP_NUM_THREADS` once (read by OpenMP at WarpX library load) and
    chdirs to the repo root so each stage's hard-coded relative paths resolve.
  * Applies `config(**kwargs)` overrides by setattr on the underlying build/sim/plot
    modules. Soft-apply: keys that don't match any module attribute are silently
    ignored (the sim module lives in a subprocess, so strict cross-phase
    validation isn't possible — cross-check `<stage>/*.py` if a config call
    seems to have no effect).
  * Runs **builds and plots in-process** — they don't touch pywarpx.
  * Runs the **simulation in a fresh subprocess** via `pipeline._launch_sim`.
    pywarpx binds globally to one geometry (2D/RZ/3D) at first .so load and
    caches diagnostic state by name; chaining cathode (2D) → gun (RZ) →
    prebuncher (RZ) in a single interpreter trips
    `AssertionError: Diagnostic attributes not consistent for "fields"`.
    A child interpreter per stage sidesteps both.
  * Writes a structured pipeline log (banner, timing, exceptions) to
    `pipeline/logs/pipeline_<ts>.log`. The child inherits this path via the
    `PIPELINE_LOG_PATH` env var and `run_step()` redirects WarpX's per-step
    output into it, so the tqdm bar (writing to the child's inherited stdout
    = parent terminal) updates on a clean line.

Logging is initialised once per process and reused across stages.
"""

import importlib
import json
import logging
import os
import re
import subprocess
import sys
import time

_PIPELINE_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_PIPELINE_DIR)
_LOG_DIR = os.path.join(_PIPELINE_DIR, "logs")

_TTY = sys.stdout.isatty()
_BOLD = "\033[1m" if _TTY else ""
_GREEN = "\033[32m" if _TTY else ""
_YELLOW = "\033[33m" if _TTY else ""
_RESET = "\033[0m" if _TTY else ""
_ANSI = re.compile(r"\033\[[0-9;]*m")

log = logging.getLogger("pipeline")
_log_path = None


def setup_logging():
    """Initialise the per-process pipeline log file. Idempotent."""
    global _log_path
    if _log_path is not None:
        return _log_path
    os.makedirs(_LOG_DIR, exist_ok=True)
    _log_path = os.path.join(_LOG_DIR, time.strftime("pipeline_%Y%m%d_%H%M%S.log"))
    log.setLevel(logging.DEBUG)
    log.propagate = False
    log.handlers.clear()
    fh = logging.FileHandler(_log_path, mode="w", encoding="utf-8")
    fh.setFormatter(logging.Formatter("%(asctime)s  %(levelname)-5s  %(message)s",
                                      "%H:%M:%S"))
    log.addHandler(fh)
    log.debug(f"python {sys.version.split()[0]} on {sys.platform}")
    log.debug(f"executable {sys.executable}")
    log.debug(f"conda env {os.environ.get('CONDA_DEFAULT_ENV', '?')}, cwd {os.getcwd()}")
    return _log_path


def _cl(msg="", level=logging.INFO):
    """Print to the console (with ANSI) AND log the plain text."""
    print(msg)
    log.log(level, _ANSI.sub("", msg).rstrip() or " ")


def run_step(sim, nsteps, desc):
    """Run `sim.step(nsteps)` with a tqdm progress bar driven by WarpX's
    `afterstep` callback. WarpX's own stdout/stderr is redirected to the
    pipeline log file for the duration of the step, so the bar updates on a
    clean terminal line (rather than being scrolled off by WarpX's init banner
    and post-step warnings). Stage sim files call this in place of bare
    `sim.step(...)`.
    """
    from pywarpx.callbacks import installcallback, uninstallcallback
    from tqdm import tqdm as _tqdm

    # The bar writes to a saved duplicate of fd 1, so it still hits the real
    # terminal even after we redirect fd 1/2.
    bar_fd = os.dup(1)
    bar_file = os.fdopen(bar_fd, "w", buffering=1, closefd=False)
    bar = _tqdm(total=nsteps, unit="step", desc=desc,
                ncols=88, leave=True, file=bar_file, disable=not _TTY)
    tick = lambda: bar.update(1)
    installcallback("afterstep", tick)

    # In a subprocess child, PIPELINE_LOG_PATH points at the parent's log.
    target = os.environ.get("PIPELINE_LOG_PATH") or _log_path or os.devnull
    redir_fd = os.open(target, os.O_WRONLY | os.O_APPEND)
    saved_out = os.dup(1)
    saved_err = os.dup(2)
    try:
        os.dup2(redir_fd, 1)
        os.dup2(redir_fd, 2)
    finally:
        os.close(redir_fd)

    try:
        sim.step(nsteps)
    finally:
        try: sys.stdout.flush()
        except Exception: pass
        try: sys.stderr.flush()
        except Exception: pass
        os.dup2(saved_out, 1)
        os.dup2(saved_err, 2)
        os.close(saved_out)
        os.close(saved_err)
        try:
            uninstallcallback("afterstep", tick)
        except Exception:
            pass
        if bar.total and bar.n < bar.total:
            bar.n = bar.total
            bar.refresh()
        bar.close()
        try: bar_file.flush()
        except Exception: pass
        os.close(bar_fd)


def _prepare_environment():
    """Set OMP_NUM_THREADS (before any pywarpx import) and chdir to the repo root."""
    os.environ.setdefault("OMP_NUM_THREADS",
                          os.environ.get("OMP_THREADS", "6"))
    if os.getcwd() != _REPO_ROOT:
        os.chdir(_REPO_ROOT)
    if _REPO_ROOT not in sys.path:
        sys.path.insert(0, _REPO_ROOT)


class Stage:
    """Facade for one accelerator stage. See cathode/gun/prebuncher __init__.py."""

    def __init__(self, name, sim_module, plot_module, build_module=None):
        self.name = name
        self._sim_path = sim_module
        self._plot_path = plot_module
        self._build_path = build_module
        self._params = {}

    # ── Public API ───────────────────────────────────────────────────────────
    def config(self, **kwargs):
        """Stage parameter overrides applied at the next run()/plot()."""
        self._params.update(kwargs)

    def run(self, plots=True):
        """Build field map (if any), simulate, then plot (unless plots=False).

        The build and plot phases run in-process (no pywarpx involved). The
        simulation phase is run in a fresh Python subprocess because pywarpx
        binds globally to one geometry at first .so load and caches
        diagnostic state by name — sharing one interpreter across stages
        trips `AssertionError: Diagnostic attributes not consistent for ...`.
        """
        _prepare_environment()
        setup_logging()
        build = self._load(self._build_path)
        if build is not None:
            self._apply_params(build)
            self._run_step(f"{self.name}: field map", build.main, is_sim=False)
        self._run_sim_subprocess()
        if plots:
            self.plot()

    def _run_sim_subprocess(self):
        title = f"{self.name}: simulation"
        _cl(f"\n{_BOLD}▶ {title}{_RESET}")
        log.info(f"    subprocess: {self._sim_path}.main()  cwd={os.getcwd()}")
        env = os.environ.copy()
        if _log_path:
            env["PIPELINE_LOG_PATH"] = _log_path
        if self._params:
            env["STAGE_CONFIG_JSON"] = json.dumps(self._params)
        cmd = [sys.executable, "-m", "pipeline._launch_sim", self._sim_path]
        t0 = time.time()
        try:
            result = subprocess.run(cmd, env=env, cwd=_REPO_ROOT)
            ok = result.returncode == 0
            err = None if ok else RuntimeError(
                f"{self.name} subprocess exited with code {result.returncode}")
        except Exception as e:
            ok, err = False, e
            log.exception(f"{title} subprocess raised {type(e).__name__}: {e}")
        dt = time.time() - t0
        flag = (f"{_GREEN}✓{_RESET}" if ok
                else f"{_YELLOW}⚠ exit {(err and getattr(err, 'args', [''])[0]) or 'fail'}{_RESET}")
        _cl(f"    {flag}  {title}  ({dt:5.1f} s)")
        log.info(f"    {title}: ok={ok}, duration = {dt:.1f} s")
        if not ok:
            raise err

    def plot(self):
        """Generate figures from the stage's existing diagnostics."""
        _prepare_environment()
        setup_logging()
        plot_mod = self._load(self._plot_path)
        self._apply_params(plot_mod)
        self._run_step(f"{self.name}: plots", plot_mod.main, is_sim=False)

    # ── Internals ────────────────────────────────────────────────────────────
    def _load(self, dotted):
        return importlib.import_module(dotted) if dotted else None

    def _apply_params(self, *modules):
        """Soft-apply config kwargs onto the given modules (no-op when a key
        is absent). Strict validation is impossible across phases because the
        sim module lives in a separate subprocess — typos in `config()` keys
        will be silently ignored, so cross-check against the parameter names
        at the top of `<stage>/*.py` if a config call doesn't seem to take."""
        for mod in modules:
            if mod is None:
                continue
            for key, value in self._params.items():
                if hasattr(mod, key):
                    setattr(mod, key, value)

    def _run_step(self, title, func, is_sim):
        _cl(f"\n{_BOLD}▶ {title}{_RESET}")
        log.info(f"    {func.__module__}.main()  cwd={os.getcwd()}")
        t0 = time.time()
        ok, err = True, None
        try:
            func()
        except Exception as e:
            ok, err = False, e
            log.exception(f"{title} raised {type(e).__name__}: {e}")
        dt = time.time() - t0
        flag = (f"{_GREEN}✓{_RESET}" if ok
                else f"{_YELLOW}⚠ {type(err).__name__}{_RESET}")
        _cl(f"    {flag}  {title}  ({dt:5.1f} s)")
        log.info(f"    {title}: ok={ok}, duration = {dt:.1f} s")
        if not ok:
            raise err
