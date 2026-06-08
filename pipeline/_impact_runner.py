"""In-process stage shim for the Impact-T `linac_rest/` stage.

`ImpactStage` mirrors the public surface of `pipeline._runner.Stage`
(`config`/`run`/`plot`/`_params`/`_warn_unknown_params`) but runs **entirely
in-process** — `build.main()`, `sim.main()`, `plot.main()` are plain function
calls, NOT a `_launch_sim` subprocess with a `run_step`/`afterstep` tqdm bar.
Impact-T has no per-step Python callback to drive a bar from, so `sim.main()`
uses the `terminal_progress` helper here (a calibrate bar ticked per section, a
final-run bar driven by a fort.18 watcher); `_run_step` exposes the saved
terminal fd as `_TERMINAL_FD` so those bars survive the sim-phase stdout
redirect, mirroring `run_step`'s saved-fd trick.

Why no subprocess (unlike every WarpX stage): the WarpX stages each spawn a
fresh interpreter because pywarpx binds globally to one geometry (2D/RZ/3D) at
first `.so` load. Impact-T has no such global — it is an external serial exe
(`ImpactTexe`) driven through lume-impact, which runs its OWN child process and
captures the exe's stdout into `Impact.log`. So `linac_rest` can run in the
parent interpreter without tripping the diagnostic-consistency assertion.

What it STILL reuses from `_runner` (the integration contract — do not drop):
  * `_prepare_environment()` — chdir to repo root (the stage's relative paths
    `linac_rest/diags`, `linac_rest/rfdata` resolve from there) AND raise
    RLIMIT_NOFILE to 16384. The fd raise matters here too: `impact_io` loops
    openPMD dumps on the handoff IN/OUT and `plot_chain` already hit the macOS
    256-fd wall (openpmd-viewer leaks an fd per `get_particle`).
  * `setup_logging()` — the shared pipeline log.
  * an fd-level stdout/stderr redirect around `sim.main()` so any lume-impact /
    `ImpactTexe` chatter (lume prints the run script + an interactive progress
    line when `verbose`; the exe banner) lands in the pipeline log rather than
    the parent terminal mid-pipeline. Our own `_cl` status lines bypass it via a
    saved duplicate fd, exactly like `run_step`'s tqdm bar.

`_warn_unknown_params` matches `Stage` exactly: live-check the build + plot
modules (imported here) AND AST-introspect the sim module's top-level names
(without importing it — keeps a `config()` key that targets a build/sim-module
section-table constant from spuriously warning).
"""

import contextlib
import importlib
import logging
import os
import sys
import time

from pipeline._runner import (
    _prepare_environment, setup_logging, _module_top_level_names,
    _cl, log, _BOLD, _GREEN, _YELLOW, _RESET, _TTY,
)

# A dup of the real terminal's fd 1, set by `_run_step` for the duration of the
# redirected sim phase (where fd 1/2 point at the capture file). `terminal_progress`
# reads it so a tqdm bar reaches the terminal even while the exe's stdout is captured
# — the same saved-fd trick `pipeline._runner.run_step` uses for the WarpX bar. None
# outside a redirected step (e.g. a direct `python -m linac_rest.linac_rest_sim`), in
# which case the bar falls back to fd 1.
_TERMINAL_FD = None


@contextlib.contextmanager
def terminal_progress(total=None, desc="", unit="it"):
    """A tqdm bar that reaches the real terminal even while the sim phase has fd 1/2
    redirected to the capture file (see `_run_step(redirect=True)` / `_TERMINAL_FD`).

    Used by the Impact-T `linac_rest` stage to show progress (calibration sections,
    final-run z) the way the WarpX stages' `run_step` bar does — they install an
    `afterstep` callback, but Impact-T is an opaque external exe, so the bars here are
    driven by the calibration loop (per section) and a fort.18 watcher (final run).
    Disabled on a non-TTY (matches `run_step`). Closes its own fd on exit.
    """
    from tqdm import tqdm as _tqdm
    fd = _TERMINAL_FD if _TERMINAL_FD is not None else 1
    bar_fd = os.dup(fd)
    bar_file = os.fdopen(bar_fd, "w", buffering=1, closefd=False)
    bar = _tqdm(total=total, desc=desc, unit=unit, ncols=88, leave=True,
                file=bar_file, disable=not _TTY)
    try:
        yield bar
    finally:
        if bar.total and bar.n < bar.total:
            bar.n = bar.total
            bar.refresh()
        bar.close()
        try: bar_file.flush()
        except Exception: pass
        os.close(bar_fd)


class ImpactStage:
    """In-process facade for the Impact-T `linac_rest/` stage.

    Same public surface as `pipeline._runner.Stage`; see this module's docstring
    for why it does NOT use the subprocess launcher.
    """

    def __init__(self, name, build_module, sim_module, plot_module):
        self.name = name
        self._build_path = build_module
        self._sim_path = sim_module
        self._plot_path = plot_module
        self._params = {}

    # ── Public API (mirrors Stage) ─────────────────────────────────────────────
    def config(self, **kwargs):
        """Cumulative parameter overrides applied at the next run()/plot().

        Same semantics as `Stage.config`: keys accumulate (`dict.update`) and
        persist into later runs until overwritten — a scan loop must set OUTDIR
        every iteration or the stale value leaks.
        """
        self._params.update(kwargs)

    def run(self, plots=True):
        """Build the deck, run Impact-T, then plot (unless plots=False).

        All in-process. build/sim/plot modules are imported here (no subprocess
        isolation needed — no pywarpx global geometry). Config overrides are
        applied to all three; unknown keys are warned about (build/plot live-
        checked, sim AST-checked).
        """
        _prepare_environment()
        setup_logging()
        build = self._load(self._build_path)
        sim = self._load(self._sim_path)
        plot_mod = self._load(self._plot_path)
        recognized = self._apply_params(build, sim, plot_mod)
        self._warn_unknown_params(recognized)
        if build is not None:
            self._run_step(f"{self.name}: lattice + rfdata", build.main)
        self._run_step(f"{self.name}: simulation", sim.main, redirect=True)
        if plots:
            self._run_step(f"{self.name}: plots", plot_mod.main)

    def plot(self):
        """Generate figures from the stage's existing diagnostics."""
        _prepare_environment()
        setup_logging()
        plot_mod = self._load(self._plot_path)
        recognized = self._apply_params(plot_mod)
        self._warn_unknown_params(recognized)
        self._run_step(f"{self.name}: plots", plot_mod.main)

    # ── Internals ──────────────────────────────────────────────────────────────
    def _load(self, dotted):
        return importlib.import_module(dotted) if dotted else None

    def _apply_params(self, *modules):
        """Soft-apply config kwargs onto the given modules (no-op when a key is
        absent on a module). Unlike `Stage`, the sim module IS imported here, so
        sim-targeting keys are applied directly (not just AST-recognised) — but
        we still AST-check in `_warn_unknown_params` so the warning logic matches
        `Stage` and tolerates a sim that fails to import in `plot()`-only mode."""
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
        """Warn about config() keys that matched no attribute on build/plot (and,
        when imported, sim) AND no top-level name in the sim module's source.

        AST-checking the sim source (in addition to the live `hasattr` in
        `_apply_params`) mirrors `Stage` exactly and covers `plot()`-only calls
        where the sim module isn't imported."""
        unknown = set(self._params) - set(recognized)
        if not unknown:
            return
        unknown -= _module_top_level_names(self._sim_path)
        if unknown:
            msg = (f"{self.name}: config() keys ignored (no matching attribute "
                   f"on build/sim/plot): {sorted(unknown)}")
            log.warning(msg)
            _cl(f"    {_YELLOW}⚠ {msg}{_RESET}", level=logging.WARNING)

    def _run_step(self, title, func, redirect=False):
        """Run `func()` with timing + the shared ok/✓ / ⚠ console+log line.

        With `redirect=True` (the sim phase), fd 1/2 are captured for the
        duration of `func()` so lume-impact / `ImpactTexe` output doesn't scroll
        the parent terminal; the captured text is then replayed into the pipeline
        log THROUGH the logger after `func()` returns. The status line is written
        to a saved duplicate of fd 1 so it still reaches the real terminal — same
        trick `run_step` uses for the tqdm bar. Mirrors `Stage._run_step`'s
        raise-on-failure contract.

        Why capture into a SEPARATE temp file rather than redirecting fd 1/2 at
        the pipeline log directly: the `setup_logging` FileHandler holds that log
        file open with its OWN buffered write offset, so a raw `os.write` (the
        exe is a real subprocess writing to its inherited fd) appended at EOF gets
        clobbered the moment the handler next flushes from its stale offset (the
        exe banner silently vanishes). Capturing to a throwaway file and then
        emitting it via `log.info` keeps a single writer (the handler) on the log.
        """
        import tempfile
        # Status banner first (before any redirect), so it lands on the terminal.
        _cl(f"\n{_BOLD}▶ {title}{_RESET}")
        log.info(f"    {func.__module__}.main()  cwd={os.getcwd()}")
        t0 = time.time()
        ok, err = True, None

        global _TERMINAL_FD
        saved_out = saved_err = cap_fd = None
        cap_path = None
        if redirect:
            try:
                cap_fd, cap_path = tempfile.mkstemp(prefix="linac_rest_sim_", suffix=".log")
                saved_out, saved_err = os.dup(1), os.dup(2)
                os.dup2(cap_fd, 1)
                os.dup2(cap_fd, 2)
                # Expose the saved terminal fd so `terminal_progress` bars (calibration,
                # final-run watcher) inside func() still reach the real terminal.
                _TERMINAL_FD = saved_out
            except Exception:
                # If the capture can't be set up, fall back to running un-redirected
                # rather than failing the whole stage. Restore fd 1/2 from any saved
                # duplicates FIRST — a partial dup2 (fd 1 redirected, fd 2 raised) would
                # otherwise leave fd 1 dangling on the about-to-be-closed temp file.
                if saved_out is not None:
                    try: os.dup2(saved_out, 1)
                    except Exception: pass
                if saved_err is not None:
                    try: os.dup2(saved_err, 2)
                    except Exception: pass
                for fd in (cap_fd, saved_out, saved_err):
                    if fd is not None:
                        try: os.close(fd)
                        except Exception: pass
                if cap_path is not None:
                    try: os.unlink(cap_path)
                    except Exception: pass
                cap_fd = saved_out = saved_err = cap_path = None
                _TERMINAL_FD = None
        try:
            func()
        except Exception as e:
            ok, err = False, e
            log.exception(f"{title} raised {type(e).__name__}: {e}")
        finally:
            if redirect:
                _TERMINAL_FD = None     # the bars are done; saved_out is about to close
                try: sys.stdout.flush()
                except Exception: pass
                try: sys.stderr.flush()
                except Exception: pass
                if saved_out is not None:
                    os.dup2(saved_out, 1); os.close(saved_out)
                if saved_err is not None:
                    os.dup2(saved_err, 2); os.close(saved_err)
                if cap_fd is not None:
                    try: os.close(cap_fd)
                    except Exception: pass
                # Replay the captured exe/sim output INTO the log via the logger,
                # so the FileHandler stays the single writer (no offset clobber).
                if cap_path is not None:
                    try:
                        with open(cap_path, "r", errors="replace") as fh:
                            captured = fh.read().rstrip()
                        if captured:
                            log.info(f"    ── {title} output ──\n{captured}")
                    except Exception:
                        pass
                    finally:
                        try: os.unlink(cap_path)
                        except Exception: pass
        dt = time.time() - t0
        flag = (f"{_GREEN}✓{_RESET}" if ok
                else f"{_YELLOW}⚠ {type(err).__name__}{_RESET}")
        _cl(f"    {flag}  {title}  ({dt:5.1f} s)")
        log.info(f"    {title}: ok={ok}, duration = {dt:.1f} s")
        if not ok:
            raise err
