"""In-process stage shim for the Impact-T `linac_rest/` stage.

`ImpactStage` mirrors `pipeline._runner.Stage`'s public surface but runs
build/sim/plot in-process (no `_launch_sim` subprocess): Impact-T is an external
exe with no pywarpx global-geometry binding to isolate. It reuses `_runner`'s
`_prepare_environment()` (repo-root chdir + RLIMIT_NOFILE raise) and
`setup_logging()`, and redirects the sim phase's stdout into the pipeline log.

See linac_rest/README.md and pipeline/README.md for the design rationale.
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

# Dup of the real terminal's fd 1 during a redirected sim phase, so a tqdm bar can
# reach the terminal while fd 1/2 point at the capture file. None outside a redirect.
_TERMINAL_FD = None


@contextlib.contextmanager
def terminal_progress(total=None, desc="", unit="it"):
    """A tqdm bar that reaches the real terminal even while the sim phase has fd 1/2
    redirected to the capture file (see `_run_step(redirect=True)` / `_TERMINAL_FD`).
    Disabled on a non-TTY. Closes its own fd on exit.
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
    """In-process facade for the Impact-T `linac_rest/` stage (mirrors
    `pipeline._runner.Stage`'s public surface; see this module's docstring)."""

    def __init__(self, name, build_module, sim_module, plot_module):
        self.name = name
        self._build_path = build_module
        self._sim_path = sim_module
        self._plot_path = plot_module
        self._params = {}

    # ── Public API (mirrors Stage) ─────────────────────────────────────────────
    def config(self, **kwargs):
        """Cumulative parameter overrides applied at the next run()/plot().

        Keys accumulate (`dict.update`) and persist until overwritten — a scan
        loop must set OUTDIR every iteration or the stale value leaks.
        """
        self._params.update(kwargs)

    def run(self, plots=True):
        """Build the deck, run Impact-T, then plot (unless plots=False), all
        in-process. Config overrides apply to build/sim/plot; unknown keys warn."""
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
        absent). Unlike `Stage`, the sim module IS imported here, so sim keys are
        applied directly — but still AST-checked in `_warn_unknown_params` to
        match `Stage` and tolerate a sim that won't import in `plot()`-only mode."""
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
        """Warn about config() keys matching no attribute on build/sim/plot AND no
        top-level name in the sim source. The AST check (vs the live `hasattr` in
        `_apply_params`) matches `Stage` and covers `plot()`-only calls."""
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
        """Run `func()` with timing + the shared ok/✓ / ⚠ console+log line;
        raise-on-failure like `Stage._run_step`.

        With `redirect=True` (sim phase), capture fd 1/2 to a SEPARATE temp file,
        then replay it into the log via `log.info` after `func()` returns. The
        temp file is required: redirecting fd 1/2 straight at the pipeline log
        would race the `setup_logging` FileHandler's buffered write offset (a raw
        `os.write` at EOF gets clobbered on the handler's next flush), so a single
        writer (the handler) must own the log. The status line goes to a saved dup
        of fd 1 so it still reaches the terminal.
        """
        import tempfile
        # Status banner before any redirect, so it lands on the terminal.
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
                # Expose the saved terminal fd so `terminal_progress` bars in func() reach it.
                _TERMINAL_FD = saved_out
            except Exception:
                # Capture setup failed: fall back to un-redirected. Restore fd 1/2 FIRST
                # — a partial dup2 would leave fd 1 dangling on the about-to-close temp file.
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
                # Replay captured output via the logger (keeps the FileHandler sole writer).
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
