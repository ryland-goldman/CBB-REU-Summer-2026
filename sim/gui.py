"""Tk control panel for the Cornell Linac chain: a per-stage pipeline runner, a
static-figure browser, and an interactive beam-properties explorer over the existing
openPMD dumps (nothing re-simulated).

Three surfaces in one window:
  - left  : one card per stage (cathode … linac5-8) with [edit config] [run] [plot]
            and, for the linac stages, [autophase] — each shells out to the same
            sim/plot/autophase scripts sim/main.py drives, streaming output to the console.
  - right : a notebook — "Beam Explorer" (Trends / 1D / 2D over a stage's ⟨z⟩-ordered
            screens) and "Plots" (the PNGs under logs/plots/<stage>/).
  - bottom: a console mirroring every subprocess's stdout/stderr.

Run from the repo root in the CBB env:  python sim/gui.py
See docs/ for the per-stage physics; the stage-local-z / σ_z-vs-σ_t conventions match sim/main.py.
"""

import os
import re
import sys
import pty
import glob
import queue
import base64
import threading
import subprocess
import warnings
from io import BytesIO

warnings.filterwarnings("ignore")

# Run from the repo root so the stage-relative diagnostic/script paths resolve.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(_ROOT)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

# openpmd-viewer leaks an fd per get_particle; raise RLIMIT_NOFILE so a full-stage
# browse doesn't hit the fd wall. prepare_env() also pins OMP=1 / HDF5 locking. Best-effort.
try:
    from sim.helpers.tools import prepare_env
    prepare_env()
except Exception:
    pass

import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from openpmd_viewer import OpenPMDTimeSeries
from pmd_beamphysics import ParticleGroup

from sim.helpers.loadparticles import make_particle_group

# ── Stages, in chain order ───────────────────────────────────────────────────
# Dumps store positions [m] and momenta u = γβ. The cathode is 2D (x–z, no y); the
# rest are RZ. `sim`/`plot`/`autophase` are argv lists passed to the current interpreter;
# `config` is the YAML edited by [edit config]; `plots` is the figure glob for the Plots tab.
STAGES = [
    {"name": "Cathode",          "config": "config/cathode.yaml",   "geom": "2d",
     "diag": "logs/diags/cathode/particles",
     "sim": ["sim/cathode.py"],          "plot": ["sim/plot/cathode.py"],
     "autophase": None,                  "plots": "logs/plots/cathode/*.png"},
    {"name": "Gun",              "config": "config/gun.yaml",       "geom": "rz",
     "diag": "logs/diags/gun/particles",
     "sim": ["sim/gun.py"],              "plot": ["sim/plot/gun.py"],
     "autophase": None,                  "plots": "logs/plots/gun/*.png"},
    {"name": "Injector",         "config": "config/injector.yaml",  "geom": "rz",
     "diag": "logs/diags/injector/main/particles",
     "sim": ["sim/injector.py"],         "plot": ["sim/plot/injector.py"],
     "autophase": None,                  "plots": "logs/plots/injector/*.png"},
    {"name": "Linac Section 1",  "config": "config/linac1.yaml",    "geom": "rz",
     "diag": "logs/diags/linac1-4/sec1/main/particles",
     "sim": ["sim/linac1-4.py", "1"],    "plot": ["sim/plot/linac1-4.py", "1"],
     "autophase": ["sim/autophase.py", "1"], "plots": "logs/plots/linac1-4/sec1_*.png"},
    {"name": "Linac Section 2",  "config": "config/linac2.yaml",    "geom": "rz",
     "diag": "logs/diags/linac1-4/sec2/main/particles",
     "sim": ["sim/linac1-4.py", "2"],    "plot": ["sim/plot/linac1-4.py", "2"],
     "autophase": ["sim/autophase.py", "2"], "plots": "logs/plots/linac1-4/sec2_*.png"},
    {"name": "Linac Section 3",  "config": "config/linac3.yaml",    "geom": "rz",
     "diag": "logs/diags/linac1-4/sec3/main/particles",
     "sim": ["sim/linac1-4.py", "3"],    "plot": ["sim/plot/linac1-4.py", "3"],
     "autophase": ["sim/autophase.py", "3"], "plots": "logs/plots/linac1-4/sec3_*.png"},
    {"name": "Linac Section 4",  "config": "config/linac4.yaml",    "geom": "rz",
     "diag": "logs/diags/linac1-4/sec4/main/particles",
     "sim": ["sim/linac1-4.py", "4"],    "plot": ["sim/plot/linac1-4.py", "4"],
     "autophase": ["sim/autophase.py", "4"], "plots": "logs/plots/linac1-4/sec4_*.png"},
    {"name": "Converter",        "config": "config/converter.yaml", "geom": "rz",
     "diag": "logs/diags/converter/main/particles",
     "sim": ["sim/converter.py"],        "plot": ["sim/plot/converter.py"],
     "autophase": None,                  "plots": "logs/plots/converter/*.png"},
    {"name": "Linac 5–8",        "config": "config/linac5-8.yaml",  "geom": "rz",
     "diag": "logs/diags/linac5-8/main/particles",
     "sim": ["sim/linac5-8.py"],         "plot": ["sim/plot/linac5-8.py"],
     "autophase": ["sim/autophase_impact.py"], "plots": "logs/plots/linac5-8/*.png"},
]

# ── Per-particle variables: ParticleGroup key → (label, SI→display scale) ────
VARS = {
    "x":              ("x [mm]",            1e3),
    "y":              ("y [mm]",            1e3),
    "z":              ("z [mm]",            1e3),
    "r":              ("r [mm]",            1e3),
    "px":             ("px [keV/c]",        1e-3),
    "py":             ("py [keV/c]",        1e-3),
    "pz":             ("pz [keV/c]",        1e-3),
    "pr":             ("pr [keV/c]",        1e-3),
    "xp":             ("x' [mrad]",         1e3),
    "yp":             ("y' [mrad]",         1e3),
    "energy":         ("energy [MeV]",      1e-6),
    "kinetic_energy": ("KE [MeV]",          1e-6),
    "gamma":          ("gamma",             1.0),
}
VARS_2D_ONLY = {"y", "py", "yp"}   # hidden when the active stage is the 2D cathode

# ── Trend Y options: label → (stat key(s), axis label, scale) ────────────────
# Bunch length is σ_z (NOT σ_t — the dumps are time snapshots).
TRENDS = {
    "Beam size σ_x, σ_y":   (["sigma_x", "sigma_y"],          "σ [mm]",         1e3),
    "Bunch length σ_z":     (["sigma_z"],                     "σ_z [mm]",       1e3),
    "Norm. emittance x, y": (["norm_emit_x", "norm_emit_y"],  "ε_n [mm·mrad]",  1e6),
    "Mean kinetic energy":  (["mean_kinetic_energy"],         "⟨KE⟩ [MeV]",     1e-6),
    "Energy spread σ_E":    (["sigma_energy"],                "σ_E [keV]",      1e-3),
    "Charge":               (["charge"],                      "q [nC]",         1e9),
    "Trajectory ⟨x⟩, ⟨y⟩":  (["mean_x", "mean_y"],            "⟨pos⟩ [mm]",     1e3),
}


# ═════════════════════════════════════════════════════════════════════════════
# Data layer: lazy per-stage loader with caching (one StageData per diag path).
# ═════════════════════════════════════════════════════════════════════════════
class StageData:
    """One stage's openPMD series: cached ParticleGroups and a ⟨z⟩-ordered screen list."""

    def __init__(self, stage):
        self.name = stage["name"]
        self.path = stage["diag"]
        self.geom = stage["geom"]
        self.ts = OpenPMDTimeSeries(self.path)
        self.species = self.ts.avail_species[0] if self.ts.avail_species else "electrons"
        self.iterations = list(self.ts.iterations)
        self.screens = None          # filled by build_screen_list(); list of (it, mean_z)
        self._pg_cache = {}          # iteration -> ParticleGroup
        self._trend_cache = {}       # trend-label -> (z[N], {key: vals[N]})
        self._range_cache = {}       # var key -> (lo, hi) raw-unit data range over ALL screens

    def build_screen_list(self, progress=None):
        """Populate `self.screens` = [(iteration, mean_z), …] sorted by ⟨z⟩.

        Reads only z/w per dump. Dumps with <2 particles are skipped (boundary/empty).
        """
        if self.screens is not None:
            return self.screens
        out = []
        n = len(self.iterations)
        for i, it in enumerate(self.iterations):
            try:
                z, w = self.ts.get_particle(["z", "w"], species=self.species, iteration=it)
            except Exception:
                continue
            if len(z) >= 2 and w.sum() > 0:
                out.append((it, float(np.average(z, weights=w))))
            if progress:
                progress(i + 1, n)
        out.sort(key=lambda r: r[1])
        self.screens = out
        return out

    def particle_group(self, iteration):
        if iteration in self._pg_cache:
            return self._pg_cache[iteration]
        if self.geom == "rz":
            x, y, z, ux, uy, uz, w = self.ts.get_particle(
                ["x", "y", "z", "ux", "uy", "uz", "w"], species=self.species, iteration=iteration)
        else:   # 2D cathode: no y / uy. Fill zeros so ParticleGroup is well-formed.
            x, z, ux, uz, w = self.ts.get_particle(
                ["x", "z", "ux", "uz", "w"], species=self.species, iteration=iteration)
            y = np.zeros_like(x)
            uy = np.zeros_like(x)
        P = make_particle_group(x, y, z, ux, uy, uz, w)       # γβ → eV/c, count → charge [C]
        # Bounded LRU: keep 16 most-recent dumps so a full-stage sweep doesn't pin RAM.
        if len(self._pg_cache) > 16:
            self._pg_cache.pop(next(iter(self._pg_cache)))
        self._pg_cache[iteration] = P
        return P

    def trend(self, label, progress=None):
        """Return (z[N], {stat_key: values[N]}) for a TRENDS entry, cached per label."""
        if label in self._trend_cache:
            return self._trend_cache[label]
        keys, _, _ = TRENDS[label]
        self.build_screen_list(progress)
        zs, series = [], {k: [] for k in keys}
        n = len(self.screens)
        for i, (it, _zmean) in enumerate(self.screens):
            P = self.particle_group(it)
            zs.append(P["mean_z"])
            for k in keys:
                try:
                    series[k].append(P[k])
                except Exception:
                    series[k].append(np.nan)
            if progress:
                progress(i + 1, n)
        result = (np.array(zs), {k: np.array(v) for k, v in series.items()})
        self._trend_cache[label] = result
        return result

    def var_range(self, key, progress=None):
        """Return (lo, hi) of `key` in raw units across ALL screens (locks fixed axes).

        Zero/negative-weight macroparticles are ignored.
        """
        if key in self._range_cache:
            return self._range_cache[key]
        self.build_screen_list(progress)
        lo, hi = np.inf, -np.inf
        n = len(self.screens)
        for i, (it, _zmean) in enumerate(self.screens):
            P = self.particle_group(it)
            try:
                v = P[key]
            except Exception:
                v = None
            if v is not None and len(v):
                m = P.weight > 0
                v = v[m] if m.any() else v
                if len(v):
                    lo = min(lo, float(np.min(v)))
                    hi = max(hi, float(np.max(v)))
            if progress:
                progress(i + 1, n)
        if not np.isfinite(lo):
            lo, hi = 0.0, 1.0
        self._range_cache[key] = (lo, hi)
        return self._range_cache[key]

    def cached_range(self, key):
        """The already-computed (lo, hi) for `key`, or None if var_range hasn't run yet."""
        return self._range_cache.get(key)


def has_dumps(diag):
    """True if `diag` holds at least one openPMD file (the dir alone may be empty)."""
    exts = ("h5", "bp", "bp4", "bp5", "sst", "json", "toml")
    return any(glob.glob(os.path.join(diag, f"*.{e}")) for e in exts)


class BouncyScroll:
    """Inertial, rubber-banding vertical scroll for a frame inside a Canvas, with a live
    scrollbar. Tk has no native overscroll, so we drive the inner frame's y-coordinate
    ourselves: wheel ticks add velocity, friction decays it, and an edge spring lets the
    content overshoot the top/bottom and settle back (the macOS rubber-band feel).
    """
    FRICTION = 0.80       # per-frame velocity decay
    SPRING = 0.18         # pull-back fraction once past an edge (gentle bounce)
    IMPULSE = 14          # px of velocity added per wheel tick
    VEL_MAX = 70          # cap so rapid trackpad events don't build runaway speed

    def __init__(self, canvas, inner, win, scrollbar):
        self.c, self.inner, self.win, self.sb = canvas, inner, win, scrollbar
        self.offset = 0.0
        self.vel = 0.0
        self._job = None
        scrollbar.config(command=self._on_scrollbar)
        canvas.bind("<Configure>", self._on_resize)
        inner.bind("<Configure>", self._on_resize)
        canvas.bind("<Enter>", lambda e: canvas.bind_all("<MouseWheel>", self._wheel))
        canvas.bind("<Leave>", lambda e: canvas.unbind_all("<MouseWheel>"))

    def _max(self):
        return max(0.0, self.inner.winfo_reqheight() - self.c.winfo_height())

    def _on_resize(self, _e=None):
        self.c.itemconfig(self.win, width=self.c.winfo_width())
        if self._job is None:
            self.offset = min(max(self.offset, 0.0), self._max())
        self._render()

    def _wheel(self, e):
        self.vel += (-1 if e.delta > 0 else 1) * self.IMPULSE
        self.vel = max(-self.VEL_MAX, min(self.VEL_MAX, self.vel))
        if self._job is None:
            self._job = self.c.after(16, self._tick)

    def _on_scrollbar(self, *args):
        m = self._max()
        if args[0] == "moveto":
            self.offset = float(args[1]) * m
        elif args[0] == "scroll":
            step = self.c.winfo_height() if args[2] == "pages" else self.IMPULSE
            self.offset += int(args[1]) * step
        self.offset = min(max(self.offset, 0.0), m)
        self.vel = 0.0
        self._render()

    def _tick(self):
        m = self._max()
        self.offset += self.vel
        self.vel *= self.FRICTION
        edge = 0 if self.offset < 0 else (m if self.offset > m else None)
        if edge is not None:
            self.offset += (edge - self.offset) * self.SPRING
            self.vel *= 0.6
            if abs(self.offset - edge) < 0.6 and abs(self.vel) < 0.6:
                self.offset, self.vel = float(edge), 0.0
        self._render()
        if abs(self.vel) > 0.5 or self.offset < 0 or self.offset > m:
            self._job = self.c.after(16, self._tick)
        else:
            self._job = None

    def _render(self):
        self.c.coords(self.win, 0, -self.offset)
        h = max(self.inner.winfo_reqheight(), 1)
        vh = self.c.winfo_height()
        top = max(self.offset, 0.0)
        self.sb.set(top / h, min(top + vh, h) / h)


def postprocess(P, *, kill_zero_weight=False, r_cut=None, z_slice=None):
    """Return a (possibly filtered) copy of P.

    kill_zero_weight : drop zero/negative-weight macroparticles.
    r_cut            : keep only r ≤ r_cut [mm] (transverse collimation preview).
    z_slice          : (center_mm, halfwidth_mm) — keep |z − center| ≤ halfwidth.
    """
    mask = np.ones(len(P.x), dtype=bool)
    if kill_zero_weight:
        mask &= P.weight > 0
    if r_cut is not None:
        mask &= P.r <= r_cut * 1e-3
    if z_slice is not None:
        c, hw = z_slice
        mask &= np.abs(P.z - c * 1e-3) <= hw * 1e-3
    if mask.all():
        return P
    if not mask.any():
        return None
    data = {k: getattr(P, k)[mask] for k in ("x", "y", "z", "px", "py", "pz", "t", "weight")}
    data["status"] = P.status[mask]
    data["species"] = P.species
    return ParticleGroup(data=data)


# ═════════════════════════════════════════════════════════════════════════════
# Subprocess runner: stream a stage script's output to the console (one at a time —
# pywarpx binds one geometry per interpreter, so stage sims must not overlap).
# ═════════════════════════════════════════════════════════════════════════════
class Runner:
    def __init__(self, console_write, on_state):
        self._write = console_write
        self._on_state = on_state          # called(busy: bool) on the main thread
        self.q = queue.Queue()
        self._proc = None
        self._pending = []                 # remaining (argv, title) jobs in the sequence
        self._stopped = False
        self.busy = False

    def run(self, argv, title):
        """Run a single job (skipped if a sequence/job is already live)."""
        self.run_many([(argv, title)])

    def run_many(self, jobs):
        """Run `jobs` = [(argv, title), …] sequentially, each starting when the prior ends."""
        if self.busy:
            self._write("\n[busy — a job is already running]\n")
            return
        if not jobs:
            return
        self._pending = list(jobs)
        self._stopped = False
        self.busy = True
        self._on_state(True)
        self._start_next()

    def _start_next(self):
        argv, title = self._pending.pop(0)
        self._write(f"\n$ python {' '.join(argv)}\n")
        env = dict(os.environ)
        env.setdefault("OMP_NUM_THREADS", "1")
        env["HDF5_USE_FILE_LOCKING"] = "FALSE"
        env["PYTHONPATH"] = _ROOT + os.pathsep + env.get("PYTHONPATH", "")
        env["PYTHONUNBUFFERED"] = "1"
        env["COLUMNS"] = "100"   # tqdm sizes its bar to this when there's no real terminal width

        def worker():
            # Drive the child through a pty so tqdm / lume-warpx keep their progress bars:
            # both gate the bar on stdout being a tty, which a plain pipe is not. The bar's
            # in-place \r updates are reflected by _console_write.
            try:
                master, slave = pty.openpty()
                self._proc = subprocess.Popen(
                    [sys.executable, *argv], cwd=_ROOT, env=env,
                    stdout=slave, stderr=slave, stdin=slave, close_fds=True)
                os.close(slave)
                while True:
                    try:
                        data = os.read(master, 4096)
                    except OSError:           # master closes when the child exits
                        break
                    if not data:
                        break
                    self.q.put(("line", data.decode("utf-8", "replace")))
                os.close(master)
                rc = self._proc.wait()
                self.q.put(("done", f"\n[{title} exited {rc}]\n"))
            except Exception as e:
                self.q.put(("done", f"\n[{title} failed: {e}]\n"))

        threading.Thread(target=worker, daemon=True).start()

    def stop(self):
        """Abort the current job and cancel any remaining queued ones."""
        self._stopped = True
        self._pending = []
        if self._proc is not None and self._proc.poll() is None:
            self._proc.terminate()

    def drain(self):
        """Pump queued output to the console; advance the sequence when a job ends. Main thread."""
        try:
            while True:
                kind, payload = self.q.get_nowait()
                self._write(payload)
                if kind == "done":
                    self._proc = None
                    if self._pending and not self._stopped:
                        self._start_next()        # next job in the sequence
                    else:
                        self.busy = False
                        self._on_state(False)
        except queue.Empty:
            pass


class BeamGUI:
    def __init__(self, root):
        self.root = root
        root.title("CESR Injector Linac Simulation (2026)")
        self.stage_data = {}          # name -> StageData (lazy)
        self.q = queue.Queue()        # worker-thread → main-thread results (explorer)
        self._busy = False
        self._gen = 0                 # monotonic token; a newer _run_async supersedes older
        self._progress_text = ""      # worker-written; reflected to Tk only on the main thread

        self.runner = Runner(self._console_write, self._on_run_state)
        self._run_buttons = []        # toggled while a subprocess is live

        self._set_app_icon()
        self._build_widgets()
        self._build_menubar()
        self._bind_keys()
        self._poll_runner()
        # Beam explorer starts on the last stage that actually has dumps on disk.
        self._on_stage_change()

    # ── chrome: app icon + menu bar ──────────────────────────────────────────
    def _set_app_icon(self):
        """Draw a small beam-line icon and set it as the window/dock icon.

        Rendered with matplotlib to an in-memory PNG so there's no committed asset; the
        PhotoImage is kept on `self` so Tk doesn't garbage-collect it.
        """
        try:
            fig = plt.figure(figsize=(1, 1), dpi=64)
            ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
            ax.add_patch(plt.Rectangle((0, 0), 1, 1, color="#101830"))
            t = np.linspace(0, 1, 200)
            ax.plot(t, 0.5 + 0.32 * np.sin(2 * np.pi * 2.2 * t) * np.exp(-1.6 * t),
                    color="#36d6ff", lw=4, solid_capstyle="round")
            ax.scatter([0.92], [0.5], s=180, color="#ffe14d", zorder=3)
            ax.set_xlim(0, 1); ax.set_ylim(0, 1)
            buf = BytesIO()
            fig.savefig(buf, format="png", transparent=False)
            plt.close(fig)
            self._icon = tk.PhotoImage(data=base64.b64encode(buf.getvalue()))
            self.root.iconphoto(True, self._icon)
        except Exception:
            pass

    def _build_menubar(self):
        menubar = tk.Menu(self.root)
        pipeline = tk.Menu(menubar, tearoff=0)
        pipeline.add_command(label="Run Selected", command=self._run_selected)
        pipeline.add_command(label="Stop", command=self.runner.stop)
        pipeline.add_separator()
        pipeline.add_command(label="Select All", command=lambda: self._select_all(True))
        pipeline.add_command(label="Select None", command=lambda: self._select_all(False))
        pipeline.add_separator()
        pipeline.add_command(label="Quit", command=self.root.destroy)
        menubar.add_cascade(label="Pipeline", menu=pipeline)

        helpm = tk.Menu(menubar, tearoff=0)
        helpm.add_command(label="About", command=self._about)
        menubar.add_cascade(label="Help", menu=helpm)
        self.root.config(menu=menubar)

    def _about(self):
        messagebox.showinfo(
            "About",
            "CESR Injector Linac Simulation (2026)\n\n"
            "Cornell CHESS electron-source beam-dynamics chain:\n"
            "cathode → gun → injector → linac 1–4 → converter → linac 5–8.\n\n"
            "Run stages from the Pipeline panel; inspect beams in Beam Explorer.")

    # ── layout ───────────────────────────────────────────────────────────────
    def _build_widgets(self):
        outer = ttk.PanedWindow(self.root, orient=tk.VERTICAL)
        outer.pack(fill=tk.BOTH, expand=True)
        work = ttk.PanedWindow(outer, orient=tk.HORIZONTAL)
        outer.add(work, weight=4)

        self._build_pipeline_panel(work)
        self._build_right_notebook(work)
        self._build_console(outer)

    # ── left: per-stage pipeline cards ───────────────────────────────────────
    def _build_pipeline_panel(self, parent):
        wrap = ttk.Frame(parent, width=320)
        wrap.pack_propagate(False)
        parent.add(wrap, weight=0)

        self.stage_selected = {}          # stage name -> BooleanVar (run-this-stage)
        self.stage_autophase = {}         # stage name -> BooleanVar (autophase before run)
        self._stage_dot = {}              # stage name -> readiness Checkbutton (●/○)

        ttk.Label(wrap, text="Pipeline", font=("", 12, "bold")).pack(anchor=tk.W, padx=8, pady=(8, 2))
        bar = ttk.Frame(wrap)
        bar.pack(fill=tk.X, padx=8)
        b = ttk.Button(bar, text="▶  Run Selected", command=self._run_selected)
        b.pack(side=tk.LEFT)
        self._run_buttons.append(b)
        ttk.Button(bar, text="■  Stop", command=self.runner.stop).pack(side=tk.LEFT, padx=4)
        bar2 = ttk.Frame(wrap)
        bar2.pack(fill=tk.X, padx=8, pady=(2, 0))
        ttk.Label(bar2, text="Select:").pack(side=tk.LEFT)
        ttk.Button(bar2, text="☑ All", width=6,
                   command=lambda: self._select_all(True)).pack(side=tk.LEFT, padx=2)
        ttk.Button(bar2, text="☐ None", width=7,
                   command=lambda: self._select_all(False)).pack(side=tk.LEFT, padx=2)

        # Scrollable, rubber-banding stack of full-width cards. BouncyScroll owns the
        # canvas/scrollbar/wheel and keeps the inner frame pinned to the canvas width.
        canvas = tk.Canvas(wrap, highlightthickness=0)
        sb = ttk.Scrollbar(wrap, orient=tk.VERTICAL)
        cards = ttk.Frame(canvas)
        win = canvas.create_window((0, 0), window=cards, anchor=tk.NW)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(8, 0), pady=8)
        sb.pack(side=tk.RIGHT, fill=tk.Y, pady=8)
        self._pipe_scroll = BouncyScroll(canvas, cards, win, sb)

        for st in STAGES:
            self._make_stage_card(cards, st)

    def _make_stage_card(self, parent, st):
        ready = has_dumps(st["diag"])
        card = ttk.Frame(parent, relief=tk.GROOVE, borderwidth=2)
        card.pack(fill=tk.X, pady=4, padx=2)

        # Header: selection checkbox + stage name (● = dumps on disk, ○ = none yet).
        sel = tk.BooleanVar(value=True)
        self.stage_selected[st["name"]] = sel
        cb = ttk.Checkbutton(card, variable=sel,
                             text=("● " if ready else "○ ") + st["name"])
        cb.pack(anchor=tk.W, padx=4, pady=(4, 0))
        self._stage_dot[st["name"]] = cb   # reconfigured by _refresh_after_run
        grid = ttk.Frame(card)
        grid.pack(fill=tk.X, padx=4, pady=4)

        # (text, command, is_run_button). Wrap into a 2-column grid so the buttons
        # always fit the narrow card regardless of how many a stage has.
        buttons = [("✎  Edit Config", lambda s=st: self._edit_config(s), False),
                   ("▶  Run", lambda s=st: self.runner.run(s["sim"], f"{s['name']} run"), True),
                   ("📊  Plot", lambda s=st: self.runner.run(s["plot"], f"{s['name']} plot"), True)]

        ncol = 2
        grid.columnconfigure(tuple(range(ncol)), weight=1, uniform="b")
        for i, (text, cmd, is_run) in enumerate(buttons):
            b = ttk.Button(grid, text=text, command=cmd)
            b.grid(row=i // ncol, column=i % ncol, sticky="ew", padx=1, pady=1)
            if is_run:
                self._run_buttons.append(b)

        # Autophase (linac stages) is a toggle, not a one-shot: when checked it runs as
        # a pre-step before this stage's sim in a "Run Selected" pass (as sim/main.py does).
        if st["autophase"]:
            ap = tk.BooleanVar(value=True)
            self.stage_autophase[st["name"]] = ap
            ttk.Checkbutton(card, variable=ap,
                            text="⚡  Autophase before run").pack(anchor=tk.W, padx=6, pady=(0, 4))

    # ── right: notebook with Beam Explorer + Plots tabs ──────────────────────
    def _build_right_notebook(self, parent):
        nb = ttk.Notebook(parent)
        parent.add(nb, weight=4)
        explorer = ttk.Frame(nb)
        plots = ttk.Frame(nb)
        nb.add(explorer, text="Beam Explorer")
        nb.add(plots, text="Plots")
        self._build_explorer(explorer)
        self._build_plots_tab(plots)

    def _build_explorer(self, parent):
        controls = ttk.Frame(parent, padding=8, width=340)
        controls.pack_propagate(False)
        controls.pack(side=tk.RIGHT, fill=tk.Y)
        figframe = ttk.Frame(parent)
        figframe.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        def row(parent, label):
            f = ttk.Frame(parent)
            f.pack(fill=tk.X, pady=2)
            ttk.Label(f, text=label, width=14).pack(side=tk.LEFT)
            return f

        # Only stages with openPMD dumps on disk are browsable.
        self.available = [s for s in STAGES if has_dumps(s["diag"])]
        names = [s["name"] for s in self.available] or ["(no data)"]
        default = "Gun" if "Gun" in names else names[-1]   # KE-vs-z phase space on the gun

        f = row(controls, "Stage")
        self.stage_var = tk.StringVar(value=default)
        self._stage_menu = ttk.OptionMenu(f, self.stage_var, default, *names,
                                          command=lambda _: self._on_stage_change())
        self._stage_menu.pack(side=tk.LEFT)

        f = row(controls, "Plot type")
        self.mode_var = tk.StringVar(value="2D Distribution")
        ttk.OptionMenu(f, self.mode_var, "2D Distribution",
                       "Trends", "1D Distribution", "2D Distribution",
                       command=lambda _: self._on_mode_change()).pack(side=tk.LEFT)

        ttk.Separator(controls, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=6)

        # Screen selector (1D / 2D modes)
        self.screen_frame = ttk.Frame(controls)
        self.screen_frame.pack(fill=tk.X)
        ttk.Label(self.screen_frame, text="Screen (by ⟨z⟩)").pack(anchor=tk.W)
        self.screen_scale = ttk.Scale(self.screen_frame, from_=0, to=0,
                                      orient=tk.HORIZONTAL, command=self._on_screen_slide)
        self.screen_scale.pack(fill=tk.X)
        self.screen_label = ttk.Label(self.screen_frame, text="—")
        self.screen_label.pack(anchor=tk.W)

        pf = ttk.Frame(self.screen_frame)
        pf.pack(fill=tk.X, pady=(4, 0))
        self._playing = False
        self.play_btn = ttk.Button(pf, text="▶ Play", width=8, command=self._toggle_play)
        self.play_btn.pack(side=tk.LEFT)
        self.loop_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(pf, text="loop", variable=self.loop_var).pack(side=tk.LEFT)
        ttk.Label(pf, text="ms").pack(side=tk.RIGHT)
        self.play_delay = tk.IntVar(value=200)
        ttk.Entry(pf, textvariable=self.play_delay, width=6).pack(side=tk.RIGHT)

        ttk.Separator(controls, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=6)

        # Variable / option controls (rebuilt per mode in _refresh_controls)
        self.ctl = ttk.Frame(controls)
        self.ctl.pack(fill=tk.X)

        ttk.Separator(controls, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=6)

        ttk.Label(controls, text="Postprocessing", font=("", 10, "bold")).pack(anchor=tk.W)
        self.kill_zero = tk.BooleanVar(value=False)
        ttk.Checkbutton(controls, text="Drop zero-weight", variable=self.kill_zero,
                        command=self.replot).pack(anchor=tk.W)
        f = row(controls, "r cut [mm]")
        self.rcut_on = tk.BooleanVar(value=False)
        ttk.Checkbutton(f, variable=self.rcut_on, command=self.replot).pack(side=tk.LEFT)
        self.rcut_val = tk.DoubleVar(value=9.547)
        ttk.Entry(f, textvariable=self.rcut_val, width=8).pack(side=tk.LEFT)
        f = row(controls, "z slice ±[mm]")
        self.zslice_on = tk.BooleanVar(value=False)
        ttk.Checkbutton(f, variable=self.zslice_on, command=self.replot).pack(side=tk.LEFT)
        self.zslice_hw = tk.DoubleVar(value=1.0)
        ttk.Entry(f, textvariable=self.zslice_hw, width=8).pack(side=tk.LEFT)

        ttk.Button(controls, text="Redraw", command=self.replot).pack(fill=tk.X, pady=(8, 2))
        self.status = ttk.Label(controls, text="", foreground="#555", wraplength=240)
        self.status.pack(fill=tk.X)

        ttk.Separator(controls, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=6)
        ttk.Label(controls, text="Beam statistics", font=("", 10, "bold")).pack(anchor=tk.W)
        self.stats = tk.Text(controls, width=30, height=12, font=("Menlo", 9),
                             relief=tk.FLAT, background="#f4f4f4", foreground="#111111")
        self.stats.pack(fill=tk.X)

        # Figure + matplotlib toolbar
        self.fig, self.ax = plt.subplots(figsize=(7.0, 5.5))
        self.cbar = None
        self.canvas = FigureCanvasTkAgg(self.fig, master=figframe)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        NavigationToolbar2Tk(self.canvas, figframe).update()

        if not self.available:
            self.status.config(text="No diagnostics on disk. Run a stage from the Pipeline panel.")
        self._refresh_controls()

    def _build_plots_tab(self, parent):
        left = ttk.Frame(parent, width=240)
        left.pack_propagate(False)
        left.pack(side=tk.LEFT, fill=tk.Y)
        ttk.Label(left, text="Figures", font=("", 10, "bold")).pack(anchor=tk.W, padx=6, pady=4)
        ttk.Button(left, text="Refresh", command=self._refresh_plot_list).pack(fill=tk.X, padx=6)
        self.plot_list = tk.Listbox(left, activestyle="none")
        self.plot_list.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)
        self.plot_list.bind("<<ListboxSelect>>", self._show_plot)

        self.pfig, self.pax = plt.subplots(figsize=(7.0, 5.5))
        self.pax.axis("off")
        self.pcanvas = FigureCanvasTkAgg(self.pfig, master=parent)
        self.pcanvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self._plot_files = []
        self._refresh_plot_list()

    def _refresh_plot_list(self):
        self.plot_list.delete(0, tk.END)
        self._plot_files = []
        for st in STAGES:
            for path in sorted(glob.glob(st["plots"])):
                self._plot_files.append(path)
                self.plot_list.insert(tk.END, f"{st['name']}: {os.path.basename(path)}")
        if not self._plot_files:
            self.plot_list.insert(tk.END, "(no figures — run a stage's plot)")

    def _show_plot(self, _evt=None):
        sel = self.plot_list.curselection()
        if not sel or sel[0] >= len(self._plot_files):
            return
        self.pax.clear()
        self.pax.axis("off")
        try:
            self.pax.imshow(plt.imread(self._plot_files[sel[0]]))
        except Exception as e:
            self.pax.text(0.5, 0.5, f"Cannot load:\n{e}", ha="center", va="center")
        self.pfig.tight_layout()
        self.pcanvas.draw()

    # ── console ──────────────────────────────────────────────────────────────
    def _build_console(self, parent):
        frame = ttk.Frame(parent)
        parent.add(frame, weight=1)
        bar = ttk.Frame(frame)
        bar.pack(fill=tk.X)
        ttk.Label(bar, text="Console", font=("", 10, "bold")).pack(side=tk.LEFT, padx=6, pady=2)
        ttk.Button(bar, text="Clear", command=lambda: self.console.delete("1.0", tk.END)
                   ).pack(side=tk.RIGHT, padx=6)
        self.console = tk.Text(frame, height=10, font=("Menlo", 9),
                               background="#101010", foreground="#d0d0d0",
                               insertbackground="#d0d0d0", wrap=tk.NONE)
        self.console.pack(fill=tk.BOTH, expand=True)

    def _console_write(self, text):
        # Honor \r as "overwrite the current line" so tqdm/WarpX progress bars animate in
        # place instead of stacking. `insert` stays at END (this is the only writer).
        c = self.console
        i, n = 0, len(text)
        while i < n:
            ch = text[i]
            if ch == "\r":
                c.delete("insert linestart", "insert")
                i += 1
            elif ch == "\n":
                c.insert("insert", "\n")
                i += 1
            else:
                j = i
                while j < n and text[j] not in "\r\n":
                    j += 1
                c.insert("insert", text[i:j])
                i = j
        c.see(tk.END)

    def _set_run_state(self, busy):
        state = tk.DISABLED if busy else tk.NORMAL
        for b in self._run_buttons:
            try:
                b.config(state=state)
            except Exception:
                pass

    def _on_run_state(self, busy):
        """Runner lifecycle hook: toggle the run buttons, and on the busy→idle edge
        resync the UI to whatever dumps/figures the just-finished jobs wrote."""
        self._set_run_state(busy)
        if not busy:
            self._refresh_after_run()

    def _refresh_after_run(self):
        """A run just finished: drop cached beam data and resync readiness dots, the stage
        menu, the Plots list, and (if a stage is loaded) the explorer to what's now on disk."""
        self.stage_data.clear()
        for name, cb in self._stage_dot.items():
            st = next(s for s in STAGES if s["name"] == name)
            cb.config(text=("● " if has_dumps(st["diag"]) else "○ ") + name)

        prev = self.stage_var.get()
        self.available = [s for s in STAGES if has_dumps(s["diag"])]
        names = [s["name"] for s in self.available] or ["(no data)"]
        menu = self._stage_menu["menu"]
        menu.delete(0, "end")
        for nm in names:
            menu.add_command(label=nm, command=lambda v=nm: self._select_stage(v))
        if prev not in names:
            self.stage_var.set(names[-1])

        self._refresh_plot_list()
        if self.available:
            self._on_stage_change()   # re-index dumps + replot the current stage

    def _select_stage(self, name):
        self.stage_var.set(name)
        self._on_stage_change()

    def _poll_runner(self):
        self.runner.drain()
        self.root.after(100, self._poll_runner)

    # ── keyboard: ←/→ step the screen slider, space toggles play ──────────────
    _KBD_TYPING = ("Entry", "TEntry", "Text", "Spinbox", "TSpinbox", "TCombobox")
    _KBD_CLICKY = ("TButton", "Button", "TCheckbutton", "Checkbutton",
                   "TMenubutton", "Menubutton", "TRadiobutton")

    def _bind_keys(self):
        self.root.bind("<Left>", lambda e: self._step_screen(-1))
        self.root.bind("<Right>", lambda e: self._step_screen(1))
        self.root.bind("<space>", self._space_toggle)

    def _focus_cls(self):
        w = self.root.focus_get()
        return w.winfo_class() if w is not None else ""

    def _step_screen(self, d):
        # Let the Scale handle arrows itself when it has focus (avoids double-stepping).
        if (self._focus_cls() in self._KBD_TYPING + self._KBD_CLICKY + ("TScale", "Scale")
                or self.mode_var.get() == "Trends" or not self.available):
            return
        self._stop_play()
        n = len(self._data().screens or [])
        if n == 0:
            return
        i = int(float(self.screen_scale.get())) + d
        self.screen_scale.set(max(0, min(i, n - 1)))   # fires _on_screen_slide → redraw

    def _space_toggle(self, _e=None):
        if (self._focus_cls() in self._KBD_TYPING + self._KBD_CLICKY
                or self.mode_var.get() == "Trends" or not self.available):
            return
        self._toggle_play()
        return "break"

    def _select_all(self, value):
        for var in self.stage_selected.values():
            var.set(value)

    def _run_selected(self):
        """Run the checked stages in chain order — autophase (if any) → sim → plot each,
        mirroring sim/main.py but limited to the selection."""
        jobs = []
        for st in STAGES:
            if not self.stage_selected[st["name"]].get():
                continue
            if st["autophase"] and self.stage_autophase[st["name"]].get():
                jobs.append((st["autophase"], f"{st['name']} autophase"))
            jobs.append((st["sim"], f"{st['name']} run"))
            jobs.append((st["plot"], f"{st['name']} plot"))
        if not jobs:
            self._console_write("\n[no stages selected]\n")
            return
        self.runner.run_many(jobs)

    # ── config editor ────────────────────────────────────────────────────────
    # YAML token → (color); a dark VS-Code-ish palette so highlighting reads well.
    _YAML_COLORS = {"comment": "#6a9955", "key": "#9cdcfe", "string": "#ce9178",
                    "number": "#b5cea8", "const": "#569cd6"}

    def _edit_config(self, st):
        path = os.path.join(_ROOT, st["config"])
        try:
            text = open(path, encoding="utf-8").read()
        except Exception as e:
            messagebox.showerror("Cannot open", f"{path}\n\n{e}")
            return
        top = tk.Toplevel(self.root)
        top.title(st["config"])
        ed = tk.Text(top, width=100, height=40, font=("Menlo", 12), wrap=tk.NONE, undo=True,
                     background="#1e1e1e", foreground="#d4d4d4", insertbackground="#d4d4d4",
                     selectbackground="#264f78", tabs="2c")
        for tag, color in self._YAML_COLORS.items():
            ed.tag_config(tag, foreground=color)
        ed.pack(fill=tk.BOTH, expand=True)
        ed.insert("1.0", text)

        def highlight(_=None):
            self._highlight_yaml(ed)

        def save(_=None):
            try:
                with open(path, "w", encoding="utf-8") as fh:
                    fh.write(ed.get("1.0", "end-1c"))
                self._console_write(f"[saved {st['config']}]\n")
                top.title(f"{st['config']} — saved ✓")
                top.after(1500, lambda: top.winfo_exists() and top.title(st["config"]))
            except Exception as e:
                messagebox.showerror("Cannot save", str(e))
            return "break"   # swallow the keystroke so no literal char is inserted

        def close(_=None):
            top.destroy()
            return "break"

        highlight()
        ed.bind("<KeyRelease>", highlight)
        for seq in ("<Command-s>", "<Control-s>"):
            ed.bind(seq, save)
            top.bind(seq, save)
        for seq in ("<Command-w>", "<Control-w>"):
            ed.bind(seq, close)
            top.bind(seq, close)

        bar = ttk.Frame(top)
        bar.pack(fill=tk.X)
        ttk.Label(bar, text="⌘/Ctrl+S save · ⌘/Ctrl+W close", foreground="#888").pack(side=tk.LEFT, padx=8)
        ttk.Button(bar, text="Save", command=save).pack(side=tk.RIGHT, padx=6, pady=4)
        ttk.Button(bar, text="Close", command=top.destroy).pack(side=tk.RIGHT, pady=4)
        ed.focus_set()

    def _highlight_yaml(self, ed):
        """Re-tag the whole (small) YAML buffer: comments, keys, strings, numbers, consts.

        Comments and strings are raised last so a '#' or digits inside them keep their color.
        """
        for tag in self._YAML_COLORS:
            ed.tag_remove(tag, "1.0", "end")
        for n, line in enumerate(ed.get("1.0", "end-1c").split("\n"), start=1):
            code = line
            m = re.search(r"(?:^|\s)#", line)         # comment: '#' at line start or after space
            if m:
                col = m.start() if line[m.start()] == "#" else m.start() + 1
                ed.tag_add("comment", f"{n}.{col}", f"{n}.end")
                code = line[:col]
            mk = re.match(r"(\s*)([A-Za-z0-9_.\-]+)(?=\s*:)", code)   # mapping key before ':'
            if mk:
                ed.tag_add("key", f"{n}.{mk.start(2)}", f"{n}.{mk.end(2)}")
            for sm in re.finditer(r"\"[^\"]*\"|'[^']*'", code):       # quoted strings
                ed.tag_add("string", f"{n}.{sm.start()}", f"{n}.{sm.end()}")
            for cm in re.finditer(r"\b(true|false|null|yes|no|on|off)\b", code, re.I):
                ed.tag_add("const", f"{n}.{cm.start()}", f"{n}.{cm.end()}")
            for nm in re.finditer(r"(?<![\w.])-?\d+\.?\d*(?:[eE][-+]?\d+)?", code):
                ed.tag_add("number", f"{n}.{nm.start()}", f"{n}.{nm.end()}")
        ed.tag_raise("string")
        ed.tag_raise("comment")

    # ── per-mode variable controls ───────────────────────────────────────────
    def _var_list(self):
        keys = [k for k in VARS if not (self._stage()["geom"] == "2d" and k in VARS_2D_ONLY)]
        return keys, [VARS[k][0] for k in keys]

    def _refresh_controls(self):
        for w in self.ctl.winfo_children():
            w.destroy()
        if not self.available:
            return
        mode = self.mode_var.get()
        keys, labels = self._var_list()
        self._key_by_label = {VARS[k][0]: k for k in keys}

        def var_row(label, default_label):
            f = ttk.Frame(self.ctl)
            f.pack(fill=tk.X, pady=2)
            ttk.Label(f, text=label, width=14).pack(side=tk.LEFT)
            v = tk.StringVar(value=default_label)
            ttk.OptionMenu(f, v, default_label, *labels,
                           command=lambda _: self._on_var_change()).pack(side=tk.LEFT)
            return v

        if mode == "Trends":
            f = ttk.Frame(self.ctl); f.pack(fill=tk.X, pady=2)
            ttk.Label(f, text="Y quantity", width=14).pack(side=tk.LEFT)
            self.trend_var = tk.StringVar(value="Norm. emittance x, y")
            ttk.OptionMenu(f, self.trend_var, "Norm. emittance x, y", *TRENDS.keys(),
                           command=lambda _: self.replot()).pack(side=tk.LEFT)
        elif mode == "1D Distribution":
            self.x_var = var_row("X variable", VARS["z"][0])
            f = ttk.Frame(self.ctl); f.pack(fill=tk.X, pady=2)
            ttk.Label(f, text="Bins", width=14).pack(side=tk.LEFT)
            self.nbins_var = tk.IntVar(value=80)
            ttk.Entry(f, textvariable=self.nbins_var, width=8).pack(side=tk.LEFT)
        else:  # 2D Distribution — default to the z–KE longitudinal phase space
            self.x_var = var_row("X variable", VARS["z"][0])
            self.y_var = var_row("Y variable", VARS["kinetic_energy"][0])
            f = ttk.Frame(self.ctl); f.pack(fill=tk.X, pady=2)
            ttk.Label(f, text="Method", width=14).pack(side=tk.LEFT)
            self.method_var = tk.StringVar(value="histogram")
            ttk.OptionMenu(f, self.method_var, "histogram", "histogram", "scatter",
                           command=lambda _: self.replot()).pack(side=tk.LEFT)
            f = ttk.Frame(self.ctl); f.pack(fill=tk.X, pady=2)
            ttk.Label(f, text="Bins X, Y", width=14).pack(side=tk.LEFT)
            self.nbins_x = tk.IntVar(value=120)
            self.nbins_y = tk.IntVar(value=120)
            ttk.Entry(f, textvariable=self.nbins_x, width=6).pack(side=tk.LEFT, padx=(0, 2))
            ttk.Entry(f, textvariable=self.nbins_y, width=6).pack(side=tk.LEFT)
            self.fixaxes_x = tk.BooleanVar(value=False)
            self.fixaxes_y = tk.BooleanVar(value=False)
            f = ttk.Frame(self.ctl); f.pack(fill=tk.X, pady=2)
            ttk.Label(f, text="Fix axis", width=14).pack(side=tk.LEFT)
            ttk.Checkbutton(f, text="x", variable=self.fixaxes_x,
                            command=self._on_fixaxes).pack(side=tk.LEFT)
            ttk.Checkbutton(f, text="y", variable=self.fixaxes_y,
                            command=self._on_fixaxes).pack(side=tk.LEFT, padx=(8, 0))

    # ── fixed-axis-range machinery (lock the 2D window across animation frames) ─
    def _fixaxes_on(self):
        return (self.mode_var.get() == "2D Distribution"
                and getattr(self, "fixaxes_x", None) is not None
                and (self.fixaxes_x.get() or self.fixaxes_y.get()))

    def _needed_range_keys(self):
        """The axis keys whose global range a locked-axis 2D plot needs (locked only)."""
        keys = []
        if self.fixaxes_x.get():
            keys.append(self._key_by_label[self.x_var.get()])
        if self.fixaxes_y.get():
            keys.append(self._key_by_label[self.y_var.get()])
        return keys

    def _on_var_change(self):
        if self._fixaxes_on():
            self._compute_ranges_async()
        else:
            self.replot()

    def _on_fixaxes(self):
        if self._fixaxes_on():
            self._compute_ranges_async()
        else:
            self.replot()

    def _compute_ranges_async(self):
        d = self._data()
        keys = list(dict.fromkeys(self._needed_range_keys()))

        def work():
            for k in keys:
                d.var_range(k, progress=lambda i, n, _k=k: setattr(
                    self, "_progress_text", f"axis range {VARS[_k][0]}: {i}/{n}…"))
            return d
        self._run_async(work, lambda _d: self.replot())

    def _apply_fixed_range(self, d, kx, sx, ky, sy):
        if self.fixaxes_x.get():
            rx = d.cached_range(kx)
            if rx is not None:
                self.ax.set_xlim(*self._pad(rx[0] * sx, rx[1] * sx))
        if self.fixaxes_y.get():
            ry = d.cached_range(ky)
            if ry is not None:
                self.ax.set_ylim(*self._pad(ry[0] * sy, ry[1] * sy))

    @staticmethod
    def _pad(lo, hi, frac=0.05):
        if not (np.isfinite(lo) and np.isfinite(hi)):
            return lo, hi
        if hi <= lo:
            d = abs(lo) * 0.05 or 1.0
            return lo - d, hi + d
        m = (hi - lo) * frac
        return lo - m, hi + m

    # ── helpers ──────────────────────────────────────────────────────────────
    def _stage(self):
        return next(s for s in self.available if s["name"] == self.stage_var.get())

    def _data(self):
        name = self.stage_var.get()
        if name not in self.stage_data:
            self.stage_data[name] = StageData(self._stage())
        return self.stage_data[name]

    def _set_status(self, msg):
        self.status.config(text=msg)
        self.root.update_idletasks()

    # ── events ───────────────────────────────────────────────────────────────
    def _on_stage_change(self):
        if not self.available:
            return
        self._stop_play()
        self._refresh_controls()
        d = self._data()
        self._run_async(lambda: self._load_screens(d), self._screens_ready)

    def _on_mode_change(self):
        self._stop_play()
        self._refresh_controls()
        is_screen_mode = self.mode_var.get() != "Trends"
        state = tk.NORMAL if is_screen_mode else tk.DISABLED
        self.screen_scale.config(state=state)
        self.play_btn.config(state=state)
        self.replot()

    def _on_screen_slide(self, _val):
        d = self._data()
        if not d.screens:
            return
        i = int(float(_val))
        i = max(0, min(i, len(d.screens) - 1))
        _it, zmean = d.screens[i]
        self.screen_label.config(text=f"#{i}   ⟨z⟩ = {zmean * 1e3:.2f} mm")
        if not self._busy:
            self.replot()

    def _toggle_play(self):
        if self._playing:
            self._stop_play()
        else:
            self._start_play()

    def _start_play(self):
        if self.mode_var.get() == "Trends":
            return
        d = self._data()
        if not d.screens:
            return
        self._playing = True
        self.play_btn.config(text="⏸ Pause")
        self._play_tick()

    def _stop_play(self):
        self._playing = False
        if getattr(self, "_play_job", None) is not None:
            try:
                self.root.after_cancel(self._play_job)
            except Exception:
                pass
            self._play_job = None
        if hasattr(self, "play_btn"):
            self.play_btn.config(text="▶ Play")

    def _play_tick(self):
        if not self._playing:
            return
        d = self._data()
        n = len(d.screens) if d.screens else 0
        if n == 0:
            self._stop_play()
            return
        i = int(float(self.screen_scale.get())) + 1
        if i >= n:
            if self.loop_var.get():
                i = 0
            else:
                self._stop_play()
                return
        self.screen_scale.set(i)
        delay = max(20, int(self.play_delay.get()))
        self._play_job = self.root.after(delay, self._play_tick)

    def _run_async(self, work, done):
        """Run `work()` off-thread; call `done(result)` on the main thread when finished.

        `work` reports progress only via `self._progress_text` (NEVER touch Tk from the
        worker). Re-entrant: each call bumps `self._gen` so a stale worker's result is
        never delivered to a newer request's `done` (see _drain).
        """
        self._gen += 1
        gen = self._gen
        self._busy = True
        self._progress_text = "Loading…"
        self._set_status(self._progress_text)

        def runner():
            try:
                self.q.put((gen, "ok", work()))
            except Exception as e:
                self.q.put((gen, "err", e))
        threading.Thread(target=runner, daemon=True).start()
        self._drain(done, gen)

    def _drain(self, done, gen):
        if gen != self._gen:
            return
        try:
            item_gen, kind, payload = self.q.get_nowait()
        except queue.Empty:
            self._set_status(self._progress_text)
            self.root.after(60, lambda: self._drain(done, gen))
            return
        if item_gen != gen:
            self.root.after(0, lambda: self._drain(done, gen))
            return
        self._busy = False
        if kind == "err":
            self._set_status(f"Error: {payload}")
        else:
            done(payload)

    def _load_screens(self, d):
        d.build_screen_list(progress=lambda i, n: setattr(self, "_progress_text",
                                                          f"Indexing dumps {i}/{n}…"))
        return d

    def _screens_ready(self, d):
        n = len(d.screens)
        self._set_status(f"{d.name}: {n} screens, species '{d.species}'")
        self.screen_scale.config(from_=0, to=max(n - 1, 0))
        self.screen_scale.set(n - 1)
        self._on_screen_slide(n - 1)
        if self._fixaxes_on():
            self._compute_ranges_async()
        else:
            self.replot()

    # ── the plot ─────────────────────────────────────────────────────────────
    def replot(self, *_):
        if self._busy or not self.available:
            return
        try:
            mode = self.mode_var.get()
            if mode == "Trends":
                self._plot_trends()
            elif mode == "1D Distribution":
                self._plot_1d()
            else:
                self._plot_2d()
        except Exception as e:
            self._set_status(f"Plot error: {e}")

    def _current_pg(self):
        d = self._data()
        if not d.screens:
            return None
        i = int(float(self.screen_scale.get()))
        i = max(0, min(i, len(d.screens) - 1))
        it, _z = d.screens[i]
        P = d.particle_group(it)
        return postprocess(
            P,
            kill_zero_weight=self.kill_zero.get(),
            r_cut=self.rcut_val.get() if self.rcut_on.get() else None,
            z_slice=(P["mean_z"] * 1e3, self.zslice_hw.get()) if self.zslice_on.get() else None,
        )

    def _reset_axes(self):
        if self.cbar is not None:
            try:
                self.cbar.remove()
            except Exception:
                pass
            self.cbar = None
        self.ax.cla()

    def _plot_trends(self):
        d = self._data()
        label = self.trend_var.get()
        keys, ylabel, scale = TRENDS[label]
        self._set_status(f"Computing '{label}' over {d.name}…")
        self._run_async(
            lambda: d.trend(label, progress=lambda i, n: setattr(self, "_progress_text",
                                                                 f"{label}: {i}/{n}")),
            lambda res: self._draw_trends(res, label, keys, ylabel, scale))

    def _draw_trends(self, res, label, keys, ylabel, scale):
        z, series = res
        self._reset_axes()
        for k in keys:
            self.ax.plot(z * 1e3, series[k] * scale, "-o", ms=3, label=k)
        self.ax.set_xlabel("dump ⟨z⟩ [mm] (stage-local)")
        self.ax.set_ylabel(ylabel)
        self.ax.set_title(f"{self._data().name}: {label}")
        if len(keys) > 1:
            self.ax.legend(fontsize=8)
        self.ax.grid(alpha=0.3)
        self.fig.tight_layout()
        self.canvas.draw()
        self._set_status(f"{self._data().name}: {len(z)} screens")

    def _plot_1d(self):
        P = self._current_pg()
        self._reset_axes()
        if P is None:
            self._set_status("No particles after cuts.")
            self.canvas.draw()
            return
        k = self._key_by_label[self.x_var.get()]
        lbl, sc = VARS[k]
        vals = P[k] * sc
        self.ax.hist(vals, bins=self.nbins_var.get(), weights=P.weight * 1e9,
                     color="C0", alpha=0.85)
        self.ax.set_xlabel(lbl)
        self.ax.set_ylabel("charge / bin [nC]")
        self.ax.set_title(f"{self._data().name}: {lbl} distribution")
        self.ax.grid(alpha=0.3)
        self.fig.tight_layout()
        self.canvas.draw()
        self._update_stats(P)

    def _plot_2d(self):
        P = self._current_pg()
        self._reset_axes()
        if P is None:
            self._set_status("No particles after cuts.")
            self.canvas.draw()
            return
        kx = self._key_by_label[self.x_var.get()]
        ky = self._key_by_label[self.y_var.get()]
        lx, sx = VARS[kx]
        ly, sy = VARS[ky]
        xv, yv = P[kx] * sx, P[ky] * sy
        if self.method_var.get() == "histogram":
            bins = [max(1, self.nbins_x.get()), max(1, self.nbins_y.get())]
            h = self.ax.hist2d(xv, yv, bins=bins, weights=P.weight * 1e9, cmap="viridis")
            self.cbar = self.fig.colorbar(h[3], ax=self.ax, label="charge [nC]")
        else:
            order = np.argsort(P.weight)   # heavy macroparticles drawn on top
            sc = self.ax.scatter(xv[order], yv[order], c=P.weight[order] * 1e9,
                                 s=4, cmap="viridis")
            self.cbar = self.fig.colorbar(sc, ax=self.ax, label="charge [nC]")
        if self._fixaxes_on():
            self._apply_fixed_range(self._data(), kx, sx, ky, sy)
        self.ax.set_xlabel(lx)
        self.ax.set_ylabel(ly)
        self.ax.set_title(f"{self._data().name}: {ly} vs {lx}")
        self.fig.tight_layout()
        self.canvas.draw()
        self._update_stats(P)

    # ── stats readout ────────────────────────────────────────────────────────
    def _update_stats(self, P):
        def g(k, default=np.nan):
            try:
                return P[k]
            except Exception:
                return default
        lines = [
            f"screen ⟨z⟩ : {g('mean_z')*1e3:8.2f} mm",
            f"macroparts: {len(P.x):8d}",
            f"charge    : {g('charge')*1e9:8.4f} nC",
            f"⟨KE⟩      : {g('mean_kinetic_energy')*1e-6:8.4f} MeV",
            f"σ_E       : {g('sigma_energy')*1e-3:8.3f} keV",
            f"σ_x       : {g('sigma_x')*1e3:8.4f} mm",
            f"σ_y       : {g('sigma_y')*1e3:8.4f} mm",
            f"σ_z       : {g('sigma_z')*1e3:8.4f} mm",
            f"ε_n,x     : {g('norm_emit_x')*1e6:8.4f} mm·mrad",
            f"ε_n,y     : {g('norm_emit_y')*1e6:8.4f} mm·mrad",
            f"⟨x⟩       : {g('mean_x')*1e3:8.4f} mm",
            f"⟨γ⟩       : {g('mean_gamma'):8.3f}",
        ]
        self.stats.delete("1.0", tk.END)
        self.stats.insert(tk.END, "\n".join(lines))


def main():
    root = tk.Tk()
    try:
        root.tk.call("tk", "scaling", 1.3)
    except Exception:
        pass
    root.geometry("1280x980")
    BeamGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
