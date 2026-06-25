"""Tk control panel + beam explorer for the Cornell2 linac chain.

One window over every stage's openPMD dumps, generated figures, and config, with per-stage
actions wired to the real drivers:
  • Edit Config     — raw YAML editor for config/<stage>.yaml (Load/Save)
  • Run Section      — python sim/<driver> [args]
  • Run From Here    — this stage + all downstream stages, in chain order
  • Autophase        — sim/autophase.py N (linac1-4) / sim/autophase_impact.py (linac5-8)
  • Generate Plots   — python sim/plot/<driver> [args]

Right pane is a notebook: Plots (PNG gallery from logs/plots/<stage>), Beam Explorer
(Trends / 1D / 2D over the dumps), Config, and a live Run Log. Subprocesses run in the
active env via sys.executable from the repo root; output streams to the Run Log.

Nothing is re-simulated by the explorer — it reads existing logs/diags dumps. See the
per-stage docs in docs/ for the physics.
"""

import os
import sys
import queue
import signal
import threading
import subprocess
import warnings

warnings.filterwarnings("ignore")

# Repo root so the stage-relative diagnostic / config / plot paths resolve, and prepare_env's
# OMP + fd-limit setup applies to this process (the explorer leaks an fd per get_particle).
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from sim.helpers.tools import prepare_env  # noqa: E402

prepare_env()
_ROOT = os.getcwd()

import numpy as np  # noqa: E402
import tkinter as tk  # noqa: E402
from tkinter import ttk, messagebox  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk  # noqa: E402

from openpmd_viewer import OpenPMDTimeSeries  # noqa: E402
from pmd_beamphysics import ParticleGroup  # noqa: E402

from sim.helpers.loadparticles import make_particle_group  # noqa: E402

# Optional Pillow for clean PNG scaling; without it we fall back to integer subsample.
try:
    from PIL import Image, ImageTk
    _HAVE_PIL = True
except Exception:
    _HAVE_PIL = False


# ── Stages, in chain order (mirrors sim/main.py's STAGES) ─────────────────────
# particles : openPMD series dir (a list = first existing wins, e.g. gun handoff→particles)
# plots_dir : where the stage plotter writes PNGs; plot_prefix filters shared dirs (linac1-4)
# autophase : argv (after sys.executable) for the stage's autophase tool, or None
def _stage(name, driver, args, config, particles, plots_dir,
           geom="rz", plot_prefix="", autophase=None):
    return dict(name=name, driver=driver, plot="sim/plot/" + driver.split("/")[-1],
                args=[str(a) for a in args], config=config,
                particles=particles if isinstance(particles, list) else [particles],
                plots_dir=plots_dir, geom=geom, plot_prefix=plot_prefix, autophase=autophase)


STAGES = [
    _stage("cathode", "sim/cathode.py", [], "config/cathode.yaml",
           "logs/diags/cathode/particles", "logs/plots/cathode", geom="2d"),
    _stage("gun", "sim/gun.py", [], "config/gun.yaml",
           ["logs/diags/gun/handoff", "logs/diags/gun/particles"], "logs/plots/gun"),
    _stage("injector", "sim/injector.py", [], "config/injector.yaml",
           "logs/diags/injector/main/particles", "logs/plots/injector"),
    _stage("linac1", "sim/linac1-4.py", [1], "config/linac1.yaml",
           "logs/diags/linac1-4/sec1/main/particles", "logs/plots/linac1-4",
           plot_prefix="sec1_", autophase=["sim/autophase.py", "1"]),
    _stage("linac2", "sim/linac1-4.py", [2], "config/linac2.yaml",
           "logs/diags/linac1-4/sec2/main/particles", "logs/plots/linac1-4",
           plot_prefix="sec2_", autophase=["sim/autophase.py", "2"]),
    _stage("linac3", "sim/linac1-4.py", [3], "config/linac3.yaml",
           "logs/diags/linac1-4/sec3/main/particles", "logs/plots/linac1-4",
           plot_prefix="sec3_", autophase=["sim/autophase.py", "3"]),
    _stage("linac4", "sim/linac1-4.py", [4], "config/linac4.yaml",
           "logs/diags/linac1-4/sec4/main/particles", "logs/plots/linac1-4",
           plot_prefix="sec4_", autophase=["sim/autophase.py", "4"]),
    _stage("converter", "sim/converter.py", [], "config/converter.yaml",
           "logs/diags/converter/main/particles", "logs/plots/converter"),
    _stage("linac5-8", "sim/linac5-8.py", [], "config/linac5-8.yaml",
           "logs/diags/linac5-8/main/particles", "logs/plots/linac5-8",
           autophase=["sim/autophase_impact.py"]),
]

# ── Per-particle variables: ParticleGroup key → (label, SI→display scale) ─────
VARS = {
    "x": ("x [mm]", 1e3), "y": ("y [mm]", 1e3), "z": ("z [mm]", 1e3), "r": ("r [mm]", 1e3),
    "px": ("px [keV/c]", 1e-3), "py": ("py [keV/c]", 1e-3),
    "pz": ("pz [keV/c]", 1e-3), "pr": ("pr [keV/c]", 1e-3),
    "xp": ("x' [mrad]", 1e3), "yp": ("y' [mrad]", 1e3),
    "energy": ("energy [MeV]", 1e-6), "kinetic_energy": ("KE [MeV]", 1e-6),
    "gamma": ("gamma", 1.0),
}
VARS_2D_ONLY = {"y", "py", "yp"}   # hidden when the active stage is the 2D cathode

# Bunch length is σ_z (NOT σ_t — WarpX dumps are time snapshots).
TRENDS = {
    "Beam size σ_x, σ_y": (["sigma_x", "sigma_y"], "σ [mm]", 1e3),
    "Bunch length σ_z": (["sigma_z"], "σ_z [mm]", 1e3),
    "Norm. emittance x, y": (["norm_emit_x", "norm_emit_y"], "ε_n [mm·mrad]", 1e6),
    "Mean kinetic energy": (["mean_kinetic_energy"], "⟨KE⟩ [MeV]", 1e-6),
    "Energy spread σ_E": (["sigma_energy"], "σ_E [keV]", 1e-3),
    "Charge": (["charge"], "q [nC]", 1e9),
    "Trajectory ⟨x⟩, ⟨y⟩": (["mean_x", "mean_y"], "⟨pos⟩ [mm]", 1e3),
}


# ═════════════════════════════════════════════════════════════════════════════
# Data layer: lazy per-stage loader with caching.
# ═════════════════════════════════════════════════════════════════════════════
class StageData:
    """One stage's openPMD series: cached ParticleGroups and a ⟨z⟩-ordered screen list."""

    def __init__(self, stage):
        self.name = stage["name"]
        self.geom = stage["geom"]
        self.ts = None
        self.path = None
        for p in stage["particles"]:                      # first dir with a valid openPMD series wins
            if not os.path.isdir(p):
                continue
            try:
                self.ts = OpenPMDTimeSeries(p)             # raises if the dir holds no valid files
                self.path = p
                break
            except Exception:
                continue
        self.species = (self.ts.avail_species[0] if self.ts and self.ts.avail_species
                        else "electrons")
        self.iterations = list(self.ts.iterations) if self.ts else []
        self.screens = None          # [(iteration, mean_z)] sorted by ⟨z⟩
        self._pg_cache = {}
        self._trend_cache = {}
        self._range_cache = {}

    def build_screen_list(self, progress=None):
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
        if len(self._pg_cache) > 16:                          # bounded LRU
            self._pg_cache.pop(next(iter(self._pg_cache)))
        self._pg_cache[iteration] = P
        return P

    def trend(self, label, progress=None):
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
        return self._range_cache.get(key)


def postprocess(P, *, kill_zero_weight=False, r_cut=None, z_slice=None):
    """Return a (possibly filtered) copy of P (drop zero-weight / r-cut / z-slice)."""
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
# Subprocess runner: a serial command queue streaming output to a callback.
# ═════════════════════════════════════════════════════════════════════════════
class ProcessRunner:
    """Run argv lists one at a time in the repo-root env, streaming combined stdout/stderr.

    on_line(text), on_done(label, returncode), on_all_done() are all invoked on the worker
    thread — the GUI marshals them back to Tk via its own after()-polled queue.
    """

    def __init__(self, on_line, on_done, on_all_done):
        self._on_line, self._on_done, self._on_all_done = on_line, on_done, on_all_done
        self._proc = None
        self._thread = None
        self._stop = False

    @property
    def busy(self):
        return self._thread is not None and self._thread.is_alive()

    def start(self, commands):
        """commands : list of (label, argv-after-python). Runs them sequentially; a non-zero
        return code aborts the rest."""
        if self.busy:
            return False
        self._stop = False
        self._thread = threading.Thread(target=self._run, args=(commands,), daemon=True)
        self._thread.start()
        return True

    def stop(self):
        self._stop = True
        p = self._proc
        if p and p.poll() is None:
            try:                                  # kill the whole group (WarpX/Impact spawn children)
                os.killpg(os.getpgid(p.pid), signal.SIGTERM)
            except Exception:
                try:
                    p.terminate()
                except Exception:
                    pass

    def _run(self, commands):
        for label, argv in commands:
            if self._stop:
                self._on_line(f"\n■ stopped before: {label}\n")
                break
            self._on_line(f"\n$ {' '.join([os.path.basename(sys.executable)] + argv)}\n")
            try:
                self._proc = subprocess.Popen(
                    [sys.executable, *argv], cwd=_ROOT,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, bufsize=1, start_new_session=True)
            except Exception as e:
                self._on_line(f"failed to launch: {e}\n")
                self._on_done(label, -1)
                break
            for line in self._proc.stdout:
                self._on_line(line)
            rc = self._proc.wait()
            self._on_done(label, rc)
            if rc != 0:
                if not self._stop:
                    self._on_line(f"\n■ {label} exited {rc} — aborting remaining steps.\n")
                break
        self._proc = None
        self._on_all_done()


# ═════════════════════════════════════════════════════════════════════════════
# Main application.
# ═════════════════════════════════════════════════════════════════════════════
class App:
    def __init__(self, root):
        self.root = root
        root.title("Cornell Linac — Control & Beam Explorer")
        self.current = STAGES[0]["name"]   # set before widgets (controls read it during build)
        self.stage_data = {}          # name -> StageData (lazy)
        self.q = queue.Queue()        # explorer worker-thread → main-thread results
        self._busy = False
        self._gen = 0
        self._progress_text = ""
        self._png_imgref = None       # keep a ref so Tk doesn't GC the shown image

        # Marshal ProcessRunner callbacks (worker thread) → main thread.
        self._proc_q = queue.Queue()
        self.runner = ProcessRunner(
            on_line=lambda s: self._proc_q.put(("line", s)),
            on_done=lambda lbl, rc: self._proc_q.put(("done", (lbl, rc))),
            on_all_done=lambda: self._proc_q.put(("all_done", None)))

        self._build_widgets()
        self.root.after(120, self._poll_proc_q)
        self._select_stage(STAGES[0]["name"])

    def _stage(self, name=None):
        name = name or self.current
        return next(s for s in STAGES if s["name"] == name)

    def _data(self):
        name = self.current
        if name not in self.stage_data:
            self.stage_data[name] = StageData(self._stage(name))
        return self.stage_data[name]

    # ── layout ────────────────────────────────────────────────────────────────
    def _build_widgets(self):
        outer = ttk.Frame(self.root, padding=6)
        outer.pack(fill=tk.BOTH, expand=True)

        left = ttk.Frame(outer)
        left.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 6))
        right = ttk.Frame(outer)
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # ── Stage navigator ──
        ttk.Label(left, text="Stage", font=("", 11, "bold")).pack(anchor=tk.W)
        self.stage_list = tk.Listbox(left, height=len(STAGES), exportselection=False,
                                     activestyle="none", width=18)
        for s in STAGES:
            self.stage_list.insert(tk.END, s["name"])
        self.stage_list.pack(fill=tk.X, pady=(2, 6))
        self.stage_list.bind("<<ListboxSelect>>", self._on_stage_pick)

        # ── Per-stage actions ──
        af = ttk.LabelFrame(left, text="Actions", padding=6)
        af.pack(fill=tk.X)
        self.btn_edit = ttk.Button(af, text="Edit Config", command=self._action_edit_config)
        self.btn_run = ttk.Button(af, text="▶ Run Section", command=self._action_run_section)
        self.btn_runfrom = ttk.Button(af, text="▶▶ Run From Here", command=self._action_run_from)
        self.btn_autophase = ttk.Button(af, text="◴ Autophase", command=self._action_autophase)
        self.btn_plots = ttk.Button(af, text="📈 Generate Plots", command=self._action_plots)
        for b in (self.btn_edit, self.btn_run, self.btn_runfrom, self.btn_autophase, self.btn_plots):
            b.pack(fill=tk.X, pady=2)
        self.btn_stop = ttk.Button(af, text="■ Stop", command=self.runner.stop, state=tk.DISABLED)
        self.btn_stop.pack(fill=tk.X, pady=(8, 2))

        self.run_status = ttk.Label(left, text="idle", foreground="#357", wraplength=170)
        self.run_status.pack(fill=tk.X, pady=(6, 0))

        # ── Beam statistics (explorer) ──
        ttk.Separator(left, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=6)
        ttk.Label(left, text="Beam statistics", font=("", 10, "bold")).pack(anchor=tk.W)
        self.stats = tk.Text(left, width=28, height=12, font=("Menlo", 9),
                             relief=tk.FLAT, background="#f4f4f4")
        self.stats.pack(fill=tk.X)

        # ── Notebook ──
        self.nb = ttk.Notebook(right)
        self.nb.pack(fill=tk.BOTH, expand=True)
        self._build_plots_tab()
        self._build_explorer_tab()
        self._build_config_tab()
        self._build_log_tab()

    # ── Plots tab ──
    def _build_plots_tab(self):
        tab = ttk.Frame(self.nb)
        self.nb.add(tab, text="Plots")
        top = ttk.Frame(tab)
        top.pack(fill=tk.X, pady=2)
        ttk.Button(top, text="⟳ Refresh", command=self._refresh_png_list).pack(side=tk.LEFT)
        self.png_caption = ttk.Label(top, text="")
        self.png_caption.pack(side=tk.LEFT, padx=8)

        body = ttk.Frame(tab)
        body.pack(fill=tk.BOTH, expand=True)
        self.png_list = tk.Listbox(body, width=28, exportselection=False)
        self.png_list.pack(side=tk.LEFT, fill=tk.Y)
        self.png_list.bind("<<ListboxSelect>>", lambda _e: self._show_png())
        # Scrollable image canvas (PNGs can exceed the pane).
        self.png_canvas = tk.Canvas(body, background="#222")
        self.png_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.png_canvas.bind("<Configure>", lambda _e: self._show_png())

    def _refresh_png_list(self):
        st = self._stage()
        self.png_list.delete(0, tk.END)
        d = st["plots_dir"]
        pref = st["plot_prefix"]
        files = []
        if os.path.isdir(d):
            files = sorted(f for f in os.listdir(d)
                           if f.endswith(".png") and f.startswith(pref))
        for f in files:
            self.png_list.insert(tk.END, f)
        self.png_caption.config(text=f"{d}  ({len(files)} png)")
        self._png_files = [os.path.join(d, f) for f in files]
        if files:
            self.png_list.selection_set(0)
            self._show_png()
        else:
            self.png_canvas.delete("all")
            self._png_imgref = None

    def _show_png(self):
        sel = self.png_list.curselection()
        if not sel:
            return
        path = self._png_files[sel[0]]
        cw = max(self.png_canvas.winfo_width(), 50)
        ch = max(self.png_canvas.winfo_height(), 50)
        self.png_canvas.delete("all")
        try:
            if _HAVE_PIL:
                im = Image.open(path)
                scale = min(cw / im.width, ch / im.height, 1.0)
                if scale < 1.0:
                    im = im.resize((max(1, int(im.width * scale)),
                                    max(1, int(im.height * scale))), Image.LANCZOS)
                img = ImageTk.PhotoImage(im)
            else:
                img = tk.PhotoImage(file=path)
                factor = max(1, int(max(img.width() / cw, img.height() / ch) + 0.999))
                if factor > 1:
                    img = img.subsample(factor, factor)
        except Exception as e:
            self.png_canvas.create_text(cw // 2, ch // 2, fill="#ccc",
                                        text=f"cannot display\n{os.path.basename(path)}\n{e}")
            self._png_imgref = None
            return
        self._png_imgref = img        # prevent GC
        self.png_canvas.create_image(cw // 2, ch // 2, image=img, anchor=tk.CENTER)

    # ── Explorer tab ──
    def _build_explorer_tab(self):
        tab = ttk.Frame(self.nb)
        self.nb.add(tab, text="Beam Explorer")
        ctlcol = ttk.Frame(tab, padding=4)
        ctlcol.pack(side=tk.LEFT, fill=tk.Y)
        figcol = ttk.Frame(tab)
        figcol.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        def row(parent, label):
            f = ttk.Frame(parent)
            f.pack(fill=tk.X, pady=2)
            ttk.Label(f, text=label, width=14).pack(side=tk.LEFT)
            return f

        f = row(ctlcol, "Plot type")
        self.mode_var = tk.StringVar(value="2D Distribution")
        ttk.OptionMenu(f, self.mode_var, "2D Distribution",
                       "Trends", "1D Distribution", "2D Distribution",
                       command=lambda _: self._on_mode_change()).pack(side=tk.LEFT)

        # Screen selector
        self.screen_frame = ttk.Frame(ctlcol)
        self.screen_frame.pack(fill=tk.X, pady=(4, 0))
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

        ttk.Separator(ctlcol, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=6)
        self.ctl = ttk.Frame(ctlcol)
        self.ctl.pack(fill=tk.X)

        ttk.Separator(ctlcol, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=6)
        ttk.Label(ctlcol, text="Postprocessing", font=("", 10, "bold")).pack(anchor=tk.W)
        self.kill_zero = tk.BooleanVar(value=False)
        ttk.Checkbutton(ctlcol, text="Drop zero-weight", variable=self.kill_zero,
                        command=self.replot).pack(anchor=tk.W)
        f = row(ctlcol, "r cut [mm]")
        self.rcut_on = tk.BooleanVar(value=False)
        ttk.Checkbutton(f, variable=self.rcut_on, command=self.replot).pack(side=tk.LEFT)
        self.rcut_val = tk.DoubleVar(value=9.547)
        ttk.Entry(f, textvariable=self.rcut_val, width=8).pack(side=tk.LEFT)
        f = row(ctlcol, "z slice ±[mm]")
        self.zslice_on = tk.BooleanVar(value=False)
        ttk.Checkbutton(f, variable=self.zslice_on, command=self.replot).pack(side=tk.LEFT)
        self.zslice_hw = tk.DoubleVar(value=1.0)
        ttk.Entry(f, textvariable=self.zslice_hw, width=8).pack(side=tk.LEFT)
        ttk.Button(ctlcol, text="Redraw", command=self.replot).pack(fill=tk.X, pady=(8, 2))
        self.status = ttk.Label(ctlcol, text="", foreground="#555", wraplength=200)
        self.status.pack(fill=tk.X)

        self.fig, self.ax = plt.subplots(figsize=(7.0, 5.6))
        self.cbar = None
        self.canvas = FigureCanvasTkAgg(self.fig, master=figcol)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        NavigationToolbar2Tk(self.canvas, figcol).update()
        self._refresh_controls()

    # ── Config tab ──
    def _build_config_tab(self):
        tab = ttk.Frame(self.nb)
        self.nb.add(tab, text="Config")
        bar = ttk.Frame(tab)
        bar.pack(fill=tk.X, pady=2)
        self.config_path_lbl = ttk.Label(bar, text="")
        self.config_path_lbl.pack(side=tk.LEFT)
        ttk.Button(bar, text="Reload", command=self._load_config_text).pack(side=tk.RIGHT)
        ttk.Button(bar, text="💾 Save", command=self._save_config_text).pack(side=tk.RIGHT, padx=4)
        wrap = ttk.Frame(tab)
        wrap.pack(fill=tk.BOTH, expand=True)
        self.config_text = tk.Text(wrap, wrap=tk.NONE, font=("Menlo", 11), undo=True)
        ysb = ttk.Scrollbar(wrap, orient=tk.VERTICAL, command=self.config_text.yview)
        self.config_text.configure(yscrollcommand=ysb.set)
        ysb.pack(side=tk.RIGHT, fill=tk.Y)
        self.config_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

    def _load_config_text(self):
        path = self._stage()["config"]
        self.config_path_lbl.config(text=path)
        self.config_text.delete("1.0", tk.END)
        try:
            with open(path) as fh:
                self.config_text.insert(tk.END, fh.read())
        except Exception as e:
            self.config_text.insert(tk.END, f"# cannot read {path}: {e}")

    def _save_config_text(self):
        path = self._stage()["config"]
        try:
            with open(path, "w") as fh:
                fh.write(self.config_text.get("1.0", "end-1c"))
            self.run_status.config(text=f"saved {path}")
        except Exception as e:
            messagebox.showerror("Save failed", str(e))

    # ── Log tab ──
    def _build_log_tab(self):
        tab = ttk.Frame(self.nb)
        self.nb.add(tab, text="Run Log")
        bar = ttk.Frame(tab)
        bar.pack(fill=tk.X, pady=2)
        ttk.Button(bar, text="Clear", command=lambda: self.log_text.delete("1.0", tk.END)).pack(
            side=tk.RIGHT)
        self.log_text = tk.Text(tab, wrap=tk.NONE, font=("Menlo", 10),
                                background="#111", foreground="#ddd", insertbackground="#ddd")
        ysb = ttk.Scrollbar(tab, orient=tk.VERTICAL, command=self.log_text.yview)
        self.log_text.configure(yscrollcommand=ysb.set)
        ysb.pack(side=tk.RIGHT, fill=tk.Y)
        self.log_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

    def _log(self, text):
        self.log_text.insert(tk.END, text)
        self.log_text.see(tk.END)

    # ── stage selection ─────────────────────────────────────────────────────
    def _on_stage_pick(self, _e):
        sel = self.stage_list.curselection()
        if sel:
            self._select_stage(STAGES[sel[0]]["name"])

    def _select_stage(self, name):
        self._stop_play()
        self.current = name
        idx = next(i for i, s in enumerate(STAGES) if s["name"] == name)
        self.stage_list.selection_clear(0, tk.END)
        self.stage_list.selection_set(idx)

        st = self._stage()
        self.btn_autophase.config(state=(tk.NORMAL if st["autophase"] else tk.DISABLED))
        self._refresh_png_list()
        self._load_config_text()
        self._refresh_controls()

        # Explorer: index dumps if this stage has any.
        d = self._data()
        if d.ts is None:
            self.status.config(text=f"{name}: no dumps yet (run the stage)")
            self.screen_scale.config(from_=0, to=0)
            self.screen_label.config(text="—")
            self._reset_axes()
            self.canvas.draw()
            self.stats.delete("1.0", tk.END)
            return
        self._run_async(lambda: self._load_screens(d), self._screens_ready)

    # ── actions ───────────────────────────────────────────────────────────────
    def _set_running(self, on):
        state = tk.DISABLED if on else tk.NORMAL
        for b in (self.btn_run, self.btn_runfrom, self.btn_plots, self.btn_edit):
            b.config(state=state)
        self.btn_autophase.config(
            state=(tk.DISABLED if on or not self._stage()["autophase"] else tk.NORMAL))
        self.btn_stop.config(state=(tk.NORMAL if on else tk.DISABLED))

    def _launch(self, commands, what):
        if self.runner.busy:
            messagebox.showinfo("Busy", "A run is already in progress.")
            return
        self.nb.select(3)             # Run Log tab
        self.run_status.config(text=f"running: {what}")
        self._set_running(True)
        self.runner.start(commands)

    def _action_run_section(self):
        st = self._stage()
        self._launch([(st["name"], [st["driver"], *st["args"]])], f"{st['name']} sim")

    def _action_run_from(self):
        start = next(i for i, s in enumerate(STAGES) if s["name"] == self.current)
        cmds = [(s["name"], [s["driver"], *s["args"]]) for s in STAGES[start:]]
        self._launch(cmds, f"{self.current} → {STAGES[-1]['name']}")

    def _action_plots(self):
        st = self._stage()
        self._launch([(f"{st['name']} plots", [st["plot"], *st["args"]])], f"{st['name']} plots")

    def _action_autophase(self):
        st = self._stage()
        if not st["autophase"]:
            return
        self._launch([(f"{st['name']} autophase", list(st["autophase"]))],
                     f"{st['name']} autophase")

    def _action_edit_config(self):
        self.nb.select(2)             # Config tab
        self._load_config_text()

    def _poll_proc_q(self):
        try:
            while True:
                kind, payload = self._proc_q.get_nowait()
                if kind == "line":
                    self._log(payload)
                elif kind == "done":
                    lbl, rc = payload
                    self.run_status.config(text=f"{lbl}: {'ok' if rc == 0 else f'exit {rc}'}")
                elif kind == "all_done":
                    self._set_running(False)
                    self.run_status.config(text="idle")
                    self._on_run_finished()
        except queue.Empty:
            pass
        self.root.after(120, self._poll_proc_q)

    def _on_run_finished(self):
        """A run finished — refresh derived views for the current stage."""
        self._refresh_png_list()
        self._load_config_text()      # autophase rewrites the YAML
        self.stage_data.pop(self.current, None)   # dumps may have changed → drop cache
        d = self._data()
        if d.ts is not None:
            self._run_async(lambda: self._load_screens(d), self._screens_ready)

    # ── per-mode variable controls ─────────────────────────────────────────────
    def _var_list(self):
        keys = [k for k in VARS if not (self._stage()["geom"] == "2d" and k in VARS_2D_ONLY)]
        return keys, [VARS[k][0] for k in keys]

    def _refresh_controls(self):
        if not hasattr(self, "ctl"):
            return
        for w in self.ctl.winfo_children():
            w.destroy()
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
        else:  # 2D Distribution
            self.x_var = var_row("X variable", VARS["x"][0])
            self.y_var = var_row("Y variable", VARS["xp"][0])
            f = ttk.Frame(self.ctl); f.pack(fill=tk.X, pady=2)
            ttk.Label(f, text="Method", width=14).pack(side=tk.LEFT)
            self.method_var = tk.StringVar(value="histogram")
            ttk.OptionMenu(f, self.method_var, "histogram", "histogram", "scatter",
                           command=lambda _: self.replot()).pack(side=tk.LEFT)
            f = ttk.Frame(self.ctl); f.pack(fill=tk.X, pady=2)
            ttk.Label(f, text="Bins", width=14).pack(side=tk.LEFT)
            self.nbins_var = tk.IntVar(value=120)
            ttk.Entry(f, textvariable=self.nbins_var, width=8).pack(side=tk.LEFT)
            self.fixaxes_var = tk.BooleanVar(value=False)
            ttk.Checkbutton(self.ctl, text="Fixed axis range",
                            variable=self.fixaxes_var,
                            command=self._on_fixaxes).pack(anchor=tk.W)

    # ── fixed-axis-range machinery ─────────────────────────────────────────────
    def _fixaxes_on(self):
        return (self.mode_var.get() == "2D Distribution"
                and getattr(self, "fixaxes_var", None) is not None
                and self.fixaxes_var.get())

    def _needed_range_keys(self):
        return [self._key_by_label[self.x_var.get()], self._key_by_label[self.y_var.get()]]

    def _on_var_change(self):
        if self._fixaxes_on():
            self._compute_ranges_async()
        else:
            self.replot()

    def _on_fixaxes(self):
        if self.fixaxes_var.get():
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
        rx, ry = d.cached_range(kx), d.cached_range(ky)
        if rx is not None:
            self.ax.set_xlim(*self._pad(rx[0] * sx, rx[1] * sx))
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

    # ── explorer async loading ─────────────────────────────────────────────────
    def _set_status(self, msg):
        self.status.config(text=msg)
        self.root.update_idletasks()

    def _run_async(self, work, done):
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
        d.build_screen_list(progress=lambda i, n: setattr(
            self, "_progress_text", f"Indexing dumps {i}/{n}…"))
        return d

    def _screens_ready(self, d):
        n = len(d.screens)
        self._set_status(f"{d.name}: {n} screens, species '{d.species}'")
        self.screen_scale.config(from_=0, to=max(n - 1, 0))
        if n:
            self.screen_scale.set(n - 1)
            self._on_screen_slide(n - 1)
        if self._fixaxes_on():
            self._compute_ranges_async()
        else:
            self.replot()

    # ── explorer events ─────────────────────────────────────────────────────────
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
        self._stop_play() if self._playing else self._start_play()

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

    # ── the plot ─────────────────────────────────────────────────────────────
    def replot(self, *_):
        if self._busy:
            return
        d = self._data()
        if d.ts is None:
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
            lambda: d.trend(label, progress=lambda i, n: setattr(
                self, "_progress_text", f"{label}: {i}/{n}")),
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
        self.ax.hist(P[k] * sc, bins=self.nbins_var.get(), weights=P.weight * 1e9,
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
        nb = self.nbins_var.get()
        if self.method_var.get() == "histogram":
            h = self.ax.hist2d(xv, yv, bins=nb, weights=P.weight * 1e9, cmap="viridis")
            self.cbar = self.fig.colorbar(h[3], ax=self.ax, label="charge [nC]")
        else:
            order = np.argsort(P.weight)
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
    root.geometry("1280x820")
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
