"""Stage shim + progress/log helpers shared by each <stage>/ facade
(cathode, gun, prebuncher, linac_sec1, …).

Each top-level package (cathode/, gun/, prebuncher/, linac_sec1/, …) instantiates
a `Stage` in its __init__.py and re-exports `config`, `run`, `plot`. The Stage object:

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
import importlib.util
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


def setup_logging(path=None):
    """Initialise the per-process pipeline log file. Idempotent.

    With `path` given, attach to that existing file in APPEND mode instead of
    creating a fresh timestamped one. This is the child-subprocess case:
    `_launch_sim` passes PIPELINE_LOG_PATH so the child's `progress:` lines and
    any other `log.*` calls (emitted by run_step) join the parent's log — without
    it the child's `pipeline` logger has no handler and those records are silently
    dropped (and a "w" reopen would truncate the parent's file).
    """
    global _log_path
    if _log_path is not None:
        return _log_path
    if path:
        _log_path, mode = path, "a"
    else:
        os.makedirs(_LOG_DIR, exist_ok=True)
        _log_path = os.path.join(_LOG_DIR, time.strftime("pipeline_%Y%m%d_%H%M%S.log"))
        mode = "w"
    log.setLevel(logging.DEBUG)
    log.propagate = False
    log.handlers.clear()
    fh = logging.FileHandler(_log_path, mode=mode, encoding="utf-8")
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
    and post-step warnings). The same callback also emits periodic
    `progress: step N/total (%) — elapsed / rate / ETA` lines to the log (~20
    over the run), so non-TTY runs (CI, nohup, redirected output) — where tqdm
    is disabled — still get progress, and the log retains it regardless of TTY.
    Stage sim files call this in place of bare `sim.step(...)`.
    """
    from pywarpx.callbacks import installcallback, uninstallcallback
    from tqdm import tqdm as _tqdm

    # The bar writes to a saved duplicate of fd 1, so it still hits the real
    # terminal even after we redirect fd 1/2.
    bar_fd = os.dup(1)
    bar_file = os.fdopen(bar_fd, "w", buffering=1, closefd=False)
    bar = _tqdm(total=nsteps, unit="step", desc=desc,
                ncols=88, leave=True, file=bar_file, disable=not _TTY)
    t0 = time.time()
    log_every = max(1, nsteps // 20)        # ~20 progress lines in the log
    state = {"step": 0, "next_log": log_every}

    def tick():
        bar.update(1)
        state["step"] += 1
        step = state["step"]
        if step >= state["next_log"] or step >= nsteps:
            el = time.time() - t0
            rate = step / el if el else 0.0
            eta = (nsteps - step) / rate if rate else 0.0
            log.info(f"    progress: step {step}/{nsteps} ({100*step/nsteps:3.0f}%)  "
                     f"elapsed {el:5.0f}s  {rate:5.1f} step/s  eta {eta:4.0f}s")
            state["next_log"] += log_every

    # From here on every acquired fd / installed callback is released in the
    # finally, so a failure mid-setup (e.g. os.open on a bad PIPELINE_LOG_PATH)
    # can't leak the bar fd or leave the afterstep hook installed.
    saved_out = saved_err = None
    try:
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

        sim.step(nsteps)
    finally:
        try: sys.stdout.flush()
        except Exception: pass
        try: sys.stderr.flush()
        except Exception: pass
        if saved_out is not None:
            os.dup2(saved_out, 1)
            os.close(saved_out)
        if saved_err is not None:
            os.dup2(saved_err, 2)
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


def _module_top_level_names(dotted):
    """Return the set of top-level assignment targets in a module's source,
    without importing it. Used to validate config() keys against the sim
    module (which we cannot import in-parent without breaking subprocess
    isolation). Best-effort: on any failure returns an empty set, so unknown
    keys still get flagged — safer to over-warn than to miss a typo."""
    try:
        import ast
        spec = importlib.util.find_spec(dotted)
        if spec is None or spec.origin is None:
            return set()
        with open(spec.origin, "r", encoding="utf-8") as f:
            tree = ast.parse(f.read(), filename=spec.origin)
        names = set()

        def _add_target(t):
            # Recurse into tuple/list targets so unpacked module constants
            # (e.g. `nr, nz = 96, 384`) are recorded, not just bare names.
            if isinstance(t, ast.Name):
                names.add(t.id)
            elif isinstance(t, (ast.Tuple, ast.List)):
                for elt in t.elts:
                    _add_target(elt)

        for node in tree.body:
            if isinstance(node, ast.Assign):
                for t in node.targets:
                    _add_target(t)
            elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
                names.add(node.target.id)
        return names
    except Exception:
        return set()


def _prepare_environment():
    """Set OMP_NUM_THREADS (before any pywarpx import) and chdir to the repo root.

    Default 1: these stages are fastest single-threaded — small grids + a
    memory-bandwidth-bound MLMG solve mean OpenMP threads only contend for the
    memory bus and add barrier overhead. Keep single-threaded; see the OMP note
    in run_pipeline.py / CLAUDE.md.
    """
    os.environ.setdefault("OMP_NUM_THREADS",
                          os.environ.get("OMP_THREADS", "1"))
    if os.getcwd() != _REPO_ROOT:
        os.chdir(_REPO_ROOT)
    if _REPO_ROOT not in sys.path:
        sys.path.insert(0, _REPO_ROOT)


class Stage:
    """Facade for one accelerator stage. See each <stage>/__init__.py
    (cathode, gun, prebuncher, linac_sec1, …)."""

    def __init__(self, name, sim_module, plot_module, build_module=None):
        self.name = name
        self._sim_path = sim_module
        self._plot_path = plot_module
        self._build_path = build_module
        self._params = {}

    # ── Public API ───────────────────────────────────────────────────────────
    def config(self, **kwargs):
        """Stage parameter overrides applied at the next run()/plot().

        **Cumulative:** keys accumulate across calls (`dict.update`) and are
        never auto-cleared, so a key set once persists into every later
        run() until overwritten. A scan loop that varies, say, POWER_W but
        also writes per-point OUTDIRs must set OUTDIR every iteration (as the
        documented prebuncher scan does) — set it once and the stale value
        leaks into subsequent runs, silently overwriting the same directory.
        """
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
        plot_mod = self._load(self._plot_path)
        recognized = self._apply_params(build, plot_mod)
        self._warn_unknown_params(recognized)
        if build is not None:
            self._run_step(f"{self.name}: field map", build.main)
        self._run_sim_subprocess()
        if plots:
            self._run_step(f"{self.name}: plots", plot_mod.main)

    def _run_sim_subprocess(self):
        title = f"{self.name}: simulation"
        _cl(f"\n{_BOLD}▶ {title}{_RESET}")
        log.info(f"    subprocess: {self._sim_path}.main()  cwd={os.getcwd()}")
        env = os.environ.copy()
        if _log_path:
            env["PIPELINE_LOG_PATH"] = _log_path
        cmd = [sys.executable, "-m", "pipeline._launch_sim", self._sim_path]
        t0 = time.time()
        rc = None
        try:
            # json.dumps is inside the try so a non-serializable config value
            # is reported through the same ok/err path as a sim failure.
            if self._params:
                env["STAGE_CONFIG_JSON"] = json.dumps(self._params)
            result = subprocess.run(cmd, env=env, cwd=_REPO_ROOT)
            rc = result.returncode
            ok = rc == 0
            err = None if ok else RuntimeError(
                f"{self.name} subprocess exited with code {rc}")
        except Exception as e:
            ok, err = False, e
            log.exception(f"{title} subprocess raised {type(e).__name__}: {e}")
        dt = time.time() - t0
        flag = (f"{_GREEN}✓{_RESET}" if ok
                else f"{_YELLOW}⚠ exit {rc if rc is not None else type(err).__name__}{_RESET}")
        _cl(f"    {flag}  {title}  ({dt:5.1f} s)")
        log.info(f"    {title}: ok={ok}, duration = {dt:.1f} s")
        if not ok:
            raise err

    def plot(self):
        """Generate figures from the stage's existing diagnostics."""
        _prepare_environment()
        setup_logging()
        plot_mod = self._load(self._plot_path)
        recognized = self._apply_params(plot_mod)
        self._warn_unknown_params(recognized)
        self._run_step(f"{self.name}: plots", plot_mod.main)

    # ── Internals ────────────────────────────────────────────────────────────
    def _load(self, dotted):
        return importlib.import_module(dotted) if dotted else None

    def _apply_params(self, *modules):
        """Soft-apply config kwargs onto the given modules (no-op when a key
        is absent). The sim module lives in a subprocess, so the child checks
        its own keys (see pipeline._launch_sim) and the parent only knows
        about build/plot — `recognized` here is the union across all modules
        the parent sees, used by `_warn_unknown_params` to flag typos."""
        recognized = set()
        for mod in modules:
            if mod is None:
                continue
            for key, value in self._params.items():
                if hasattr(mod, key):
                    setattr(mod, key, value)
                    recognized.add(key)
        return recognized

    def _warn_unknown_params(self, recognized):
        """Warn about config() keys that matched no attribute on any module
        the parent loaded (build/plot) AND no attribute on the sim module.
        The sim is in a subprocess, so we introspect its source for top-level
        bindings rather than importing it here (which would defeat the
        per-stage subprocess isolation)."""
        unknown = set(self._params) - set(recognized)
        if not unknown:
            return
        sim_names = _module_top_level_names(self._sim_path)
        unknown -= sim_names
        if unknown:
            msg = (f"{self.name}: config() keys ignored (no matching attribute "
                   f"on build/sim/plot): {sorted(unknown)}")
            log.warning(msg)
            _cl(f"    {_YELLOW}⚠ {msg}{_RESET}", level=logging.WARNING)

    def _run_step(self, title, func):
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
