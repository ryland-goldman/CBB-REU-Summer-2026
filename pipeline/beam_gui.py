"""Interactive beam-properties explorer for the Cornell Linac pipeline.

A standalone GUI — modelled on GPT_tools' ``gpt_plot_gui.py`` / ``gpt_plot.py``
(https://github.com/AdamCBartnik/GPT_tools) — for manually investigating the beam
produced by ``run_pipeline.py``. It is a **separate** program: run the pipeline first
to generate the openPMD dumps, then launch this to browse them.

    conda activate CBB
    python pipeline/beam_gui.py

Unlike GPT_tools (ipywidgets, Jupyter-only) this is a Tkinter desktop window with an
embedded matplotlib canvas, so it runs from a plain terminal after the pipeline. It
reads each stage's existing openPMD particle series — nothing is re-simulated — and
loads every dump into a ``pmd_beamphysics.ParticleGroup`` so the rich derived
quantities (emittance, σ, energy, …) come for free, exactly as GPT_tools leans on its
ParticleGroup.

Three plot modes mirror the GPT GUI:
  * Trends         — a beam statistic vs dump ⟨z⟩ across every dump of a stage.

NOTE ON THE z AXIS: each dump is placed at its raw openPMD charge-weighted ⟨z⟩. For cathode/gun/
injector that ⟨z⟩ is already lab-frame; for linac_sec1 and linac_rest the openPMD z is STAGE-LOCAL
(linac_sec1 writes z−z.min()+Z_INJECT; linac_rest zeroes z at injection), so their ⟨z⟩ is offset
from lab z by ~1.9 m / ~5.1 m respectively. `plot_chain` applies those per-stage z0 shifts to stitch
the chain in true lab z; this explorer does NOT — it shows each stage one at a time, so the absolute
axis is labelled "dump ⟨z⟩ (stage-local for the two linac sections)". Every derived quantity (σ,
emittance, ⟨KE⟩, σ_E, charge) is z-offset-invariant, so only the absolute z position is affected.
  * 1D Distribution — a weighted histogram of one coordinate on one dump (screen).
  * 2D Distribution — a phase-space view (hist2d or scatter) of two coordinates.

ADAPTATION TO WARPX DATA (vs GPT screens): the pipeline's WarpX dumps are time
snapshots — every particle in a dump shares one simulation time t, so σ_t ≈ 0 and the
longitudinal spread lives in z. GPT "screens" are z-plane crossings (σ_t meaningful).
So here a "screen" is one dump, ordered by its charge-weighted ⟨z⟩, and bunch length
is reported as σ_z (not σ_t). Momenta are openPMD u = γβ → converted to eV/c for the
ParticleGroup; positions are lab-frame metres.
"""

import os
import sys
import threading
import queue
import warnings

warnings.filterwarnings("ignore")

# Run from the repo root so the stage-relative diagnostic paths below resolve, exactly
# like run_pipeline.py does. (This file lives in pipeline/, so the root is its parent's
# parent.)
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(_ROOT)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

# fd-limit raise: the injector series is ~290 dumps and openpmd-viewer leaks an fd per
# get_particle, so browsing a full stage would exhaust macOS's 256-fd soft limit. Reuse
# the pipeline's helper (same reason as plot_chain.py). Best-effort — skip if unavailable.
try:
    from pipeline._runner import _raise_fd_limit
    _raise_fd_limit()
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

# ── Physical constants ───────────────────────────────────────────────────────
MC2_EV = 0.51099895e6          # electron rest energy [eV]
Q_E = 1.602176634e-19          # elementary charge [C]

# ── Stages, in chain order (mirrors pipeline/plot_chain.py STAGES) ───────────
# Each WarpX dump in these series stores positions [m] and momenta u = γβ for the
# `electrons` species. The cathode is 2D (x–z, no y); the rest are RZ.
STAGES = [
    {"name": "cathode",     "path": "cathode/diags/particles",          "geom": "2d"},
    {"name": "gun",         "path": "gun/diags/particles",              "geom": "rz"},
    {"name": "injector",    "path": "injector/diags/main/particles",    "geom": "rz"},
    {"name": "linac_sec1",  "path": "linac_sec1/diags/main/particles",  "geom": "rz"},
    {"name": "linac_rest",  "path": "linac_rest/diags/main/particles",  "geom": "rz"},
]

# ── Per-particle variables: ParticleGroup key → (label, display scale) ───────
# The 2D cathode has no y/py, so y-derived keys are hidden for that stage (see
# BeamGUI._var_list). Scales convert SI ParticleGroup units to readable display units.
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

# ── Trend Y options: label → (ParticleGroup stat key(s), axis label, scale) ──
# Multiple keys (e.g. σ_x and σ_y) plot as several lines on one axis. Bunch length is
# σ_z (NOT σ_t — WarpX dumps are time snapshots; see module docstring).
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
# Data layer: lazy per-stage loader with caching.
# ═════════════════════════════════════════════════════════════════════════════
class StageData:
    """One stage's openPMD series, with cached ParticleGroups and a z-ordered screen list.

    `screens` is the list of (iteration, mean_z) sorted by mean_z — the "screens" the GUI
    browses. Building it reads only z/w per dump (cheap-ish); full ParticleGroups are built
    on demand and cached by iteration.
    """

    def __init__(self, stage):
        self.name = stage["name"]
        self.path = stage["path"]
        self.geom = stage["geom"]
        self.ts = OpenPMDTimeSeries(self.path)
        self.species = self.ts.avail_species[0] if self.ts.avail_species else "electrons"
        self.iterations = list(self.ts.iterations)
        self.screens = None          # filled by build_screen_list(); list of (it, mean_z)
        self._pg_cache = {}          # iteration -> ParticleGroup
        self._trend_cache = {}       # trend-label -> (z[N], {key: vals[N]})

    # ── screen (dump) index, ordered by ⟨z⟩ (stage-local for linac_sec1/linac_rest) ──
    def build_screen_list(self, progress=None):
        """Populate `self.screens` = [(iteration, mean_z), …] sorted by ⟨z⟩.

        Reads only z and w per dump. `progress(done, total)` is called for a status bar.
        Dumps with <2 particles are skipped (boundary/empty dumps).
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

    # ── full ParticleGroup for one dump (cached) ─────────────────────────────
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
        P = ParticleGroup(data=dict(
            x=x, y=y, z=z,
            px=ux * MC2_EV, py=uy * MC2_EV, pz=uz * MC2_EV,   # γβ → eV/c
            t=np.zeros_like(x),
            status=np.ones_like(x, dtype=int),
            weight=w * Q_E,                                   # macro-weight → charge [C]
            species="electron",
        ))
        # Bounded cache: keep the 16 most-recently-used dumps so trend sweeps over a
        # 290-dump stage don't pin all of them in RAM.
        if len(self._pg_cache) > 16:
            self._pg_cache.pop(next(iter(self._pg_cache)))
        self._pg_cache[iteration] = P
        return P

    # ── trend: one stat (or several) vs ⟨z⟩ across every screen ──────────────
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


# ═════════════════════════════════════════════════════════════════════════════
# Postprocessing on a ParticleGroup (a faithful subset of GPT_tools' options).
# ═════════════════════════════════════════════════════════════════════════════
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
# The GUI.
# ═════════════════════════════════════════════════════════════════════════════
class BeamGUI:
    def __init__(self, root):
        self.root = root
        root.title("Cornell Linac — Beam Properties Explorer")
        self.stage_data = {}          # name -> StageData (lazy)
        self.q = queue.Queue()        # worker-thread → main-thread results
        self._busy = False
        self._gen = 0                 # monotonic token: a newer _run_async supersedes older ones,
                                      # so an in-flight worker's payload can't be mis-delivered to a
                                      # later request's `done` callback (see _run_async / _drain).
        self._progress_text = ""      # written by worker threads, reflected to the
                                      # status label only on the main thread (Tk is not
                                      # thread-safe — see _run_async / _drain).

        # Discover which stages actually have data on disk.
        self.available = [st for st in STAGES if os.path.isdir(st["path"])]
        if not self.available:
            messagebox.showerror(
                "No data",
                "No stage diagnostic series found under the repo root.\n\n"
                "Run the pipeline first:\n    python pipeline/run_pipeline.py")
            root.destroy()
            return

        self._build_widgets()
        self._on_stage_change()

    # ── layout ───────────────────────────────────────────────────────────────
    def _build_widgets(self):
        left = ttk.Frame(self.root, padding=8)
        left.pack(side=tk.LEFT, fill=tk.Y)
        right = ttk.Frame(self.root)
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        def row(parent, label):
            f = ttk.Frame(parent)
            f.pack(fill=tk.X, pady=2)
            ttk.Label(f, text=label, width=16).pack(side=tk.LEFT)
            return f

        # Stage + plot type
        f = row(left, "Stage")
        self.stage_var = tk.StringVar(value=self.available[0]["name"])
        ttk.OptionMenu(f, self.stage_var, self.available[0]["name"],
                       *[s["name"] for s in self.available],
                       command=lambda _: self._on_stage_change()).pack(side=tk.LEFT)

        f = row(left, "Plot type")
        self.mode_var = tk.StringVar(value="2D Distribution")
        ttk.OptionMenu(f, self.mode_var, "2D Distribution",
                       "Trends", "1D Distribution", "2D Distribution",
                       command=lambda _: self._on_mode_change()).pack(side=tk.LEFT)

        ttk.Separator(left, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=6)

        # Screen selector (1D / 2D modes)
        self.screen_frame = ttk.Frame(left)
        self.screen_frame.pack(fill=tk.X)
        ttk.Label(self.screen_frame, text="Screen (by ⟨z⟩)").pack(anchor=tk.W)
        self.screen_scale = ttk.Scale(self.screen_frame, from_=0, to=0,
                                      orient=tk.HORIZONTAL, command=self._on_screen_slide)
        self.screen_scale.pack(fill=tk.X)
        self.screen_label = ttk.Label(self.screen_frame, text="—")
        self.screen_label.pack(anchor=tk.W)

        ttk.Separator(left, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=6)

        # Variable / option controls (rebuilt per mode in _refresh_controls)
        self.ctl = ttk.Frame(left)
        self.ctl.pack(fill=tk.X)

        ttk.Separator(left, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=6)

        # Postprocessing
        ttk.Label(left, text="Postprocessing", font=("", 10, "bold")).pack(anchor=tk.W)
        self.kill_zero = tk.BooleanVar(value=False)
        ttk.Checkbutton(left, text="Drop zero-weight", variable=self.kill_zero,
                        command=self.replot).pack(anchor=tk.W)
        f = row(left, "r cut [mm]")
        self.rcut_on = tk.BooleanVar(value=False)
        ttk.Checkbutton(f, variable=self.rcut_on, command=self.replot).pack(side=tk.LEFT)
        self.rcut_val = tk.DoubleVar(value=9.547)
        ttk.Entry(f, textvariable=self.rcut_val, width=8).pack(side=tk.LEFT)
        f = row(left, "z slice ±[mm]")
        self.zslice_on = tk.BooleanVar(value=False)
        ttk.Checkbutton(f, variable=self.zslice_on, command=self.replot).pack(side=tk.LEFT)
        self.zslice_hw = tk.DoubleVar(value=1.0)
        ttk.Entry(f, textvariable=self.zslice_hw, width=8).pack(side=tk.LEFT)

        ttk.Button(left, text="Redraw", command=self.replot).pack(fill=tk.X, pady=(8, 2))
        self.status = ttk.Label(left, text="", foreground="#555", wraplength=220)
        self.status.pack(fill=tk.X)

        # Stats readout
        ttk.Separator(left, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=6)
        ttk.Label(left, text="Beam statistics", font=("", 10, "bold")).pack(anchor=tk.W)
        self.stats = tk.Text(left, width=30, height=12, font=("Menlo", 9),
                             relief=tk.FLAT, background="#f4f4f4")
        self.stats.pack(fill=tk.X)

        # Figure + matplotlib toolbar
        self.fig, self.ax = plt.subplots(figsize=(7.5, 6))
        self.cbar = None
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        NavigationToolbar2Tk(self.canvas, right).update()

        self._refresh_controls()

    # ── per-mode variable controls ───────────────────────────────────────────
    def _var_list(self):
        keys = [k for k in VARS if not (self._stage()["geom"] == "2d" and k in VARS_2D_ONLY)]
        return keys, [VARS[k][0] for k in keys]

    def _refresh_controls(self):
        for w in self.ctl.winfo_children():
            w.destroy()
        mode = self.mode_var.get()
        keys, labels = self._var_list()
        self._key_by_label = {VARS[k][0]: k for k in keys}

        def var_row(label, default_label):
            f = ttk.Frame(self.ctl)
            f.pack(fill=tk.X, pady=2)
            ttk.Label(f, text=label, width=16).pack(side=tk.LEFT)
            v = tk.StringVar(value=default_label)
            ttk.OptionMenu(f, v, default_label, *labels,
                           command=lambda _: self.replot()).pack(side=tk.LEFT)
            return v

        if mode == "Trends":
            f = ttk.Frame(self.ctl); f.pack(fill=tk.X, pady=2)
            ttk.Label(f, text="Y quantity", width=16).pack(side=tk.LEFT)
            self.trend_var = tk.StringVar(value="Norm. emittance x, y")
            ttk.OptionMenu(f, self.trend_var, "Norm. emittance x, y", *TRENDS.keys(),
                           command=lambda _: self.replot()).pack(side=tk.LEFT)
        elif mode == "1D Distribution":
            self.x_var = var_row("X variable", VARS["z"][0])
            f = ttk.Frame(self.ctl); f.pack(fill=tk.X, pady=2)
            ttk.Label(f, text="Bins", width=16).pack(side=tk.LEFT)
            self.nbins_var = tk.IntVar(value=80)
            ttk.Entry(f, textvariable=self.nbins_var, width=8).pack(side=tk.LEFT)
        else:  # 2D Distribution
            self.x_var = var_row("X variable", VARS["x"][0])
            self.y_var = var_row("Y variable", VARS["xp"][0])
            f = ttk.Frame(self.ctl); f.pack(fill=tk.X, pady=2)
            ttk.Label(f, text="Method", width=16).pack(side=tk.LEFT)
            self.method_var = tk.StringVar(value="histogram")
            ttk.OptionMenu(f, self.method_var, "histogram", "histogram", "scatter",
                           command=lambda _: self.replot()).pack(side=tk.LEFT)
            f = ttk.Frame(self.ctl); f.pack(fill=tk.X, pady=2)
            ttk.Label(f, text="Bins", width=16).pack(side=tk.LEFT)
            self.nbins_var = tk.IntVar(value=120)
            ttk.Entry(f, textvariable=self.nbins_var, width=8).pack(side=tk.LEFT)
            self.equal_var = tk.BooleanVar(value=False)
            ttk.Checkbutton(self.ctl, text="Equal axis scale", variable=self.equal_var,
                            command=self.replot).pack(anchor=tk.W)

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
        self._refresh_controls()
        # Resolve the StageData on the MAIN thread (it reads the Tk stage_var and is cheap —
        # just opens the series). Only the heavy screen indexing (reading z/w over ~hundreds
        # of dumps) goes off-thread, so the worker never touches Tk.
        d = self._data()
        self._run_async(lambda: self._load_screens(d), self._screens_ready)

    def _on_mode_change(self):
        self._refresh_controls()
        is_screen_mode = self.mode_var.get() != "Trends"
        state = tk.NORMAL if is_screen_mode else tk.DISABLED
        self.screen_scale.config(state=state)
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

    # ── async machinery (worker thread + main-thread queue drain) ────────────
    def _run_async(self, work, done):
        """Run `work()` off-thread; call `done(result)` on the main thread when finished.

        `work` may report progress by assigning to `self._progress_text` (a plain string —
        NEVER touch Tk widgets from the worker). `_drain`, which runs on the main thread via
        `root.after`, reflects that string to the status label.

        Re-entrant by design: each call bumps `self._gen` and stamps its worker payload with that
        token. If the user triggers a second `_run_async` (e.g. switches stage while a Trends
        computation is still running), the older drain loop sees `gen != self._gen` and stops
        WITHOUT delivering — so a stale worker's result is never routed to the newer request's
        `done`, and the orphaned 60 ms poll loop is not left rescheduling forever.
        """
        self._gen += 1
        gen = self._gen
        self._busy = True
        self._progress_text = "Loading…"
        self._set_status(self._progress_text)

        def runner():
            try:
                self.q.put((gen, "ok", work()))
            except Exception as e:        # surface loader errors to the status line
                self.q.put((gen, "err", e))
        threading.Thread(target=runner, daemon=True).start()
        self._drain(done, gen)

    def _drain(self, done, gen):
        if gen != self._gen:             # superseded by a newer _run_async — stop, deliver nothing
            return
        try:
            item_gen, kind, payload = self.q.get_nowait()
        except queue.Empty:
            self._set_status(self._progress_text)   # main-thread progress reflection
            self.root.after(60, lambda: self._drain(done, gen))
            return
        if item_gen != gen:              # a stale leftover payload — discard, keep draining for `gen`
            self.root.after(0, lambda: self._drain(done, gen))
            return
        self._busy = False
        if kind == "err":
            self._set_status(f"Error: {payload}")
        else:
            done(payload)

    def _load_screens(self, d):
        # Runs on the worker thread — `d` was resolved on the main thread. Progress callback
        # sets the plain string only (no Tk).
        d.build_screen_list(progress=lambda i, n: setattr(self, "_progress_text",
                                                           f"Indexing dumps {i}/{n}…"))
        return d

    def _screens_ready(self, d):
        n = len(d.screens)
        self._set_status(f"{d.name}: {n} screens, species '{d.species}'")
        self.screen_scale.config(from_=0, to=max(n - 1, 0))
        # Default to the last screen (stage exit) — usually the most interesting.
        self.screen_scale.set(n - 1)
        self._on_screen_slide(n - 1)
        self.replot()

    # ── the plot ─────────────────────────────────────────────────────────────
    def replot(self, *_):
        if self._busy:
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
        """ParticleGroup for the selected screen, with postprocessing applied."""
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
        # trend() runs on the worker thread; its progress callback must not touch Tk.
        self._run_async(
            lambda: d.trend(label, progress=lambda i, n: setattr(self, "_progress_text",
                                                                 f"{label}: {i}/{n}")),
            lambda res: self._draw_trends(res, label, keys, ylabel, scale))

    def _draw_trends(self, res, label, keys, ylabel, scale):
        z, series = res
        self._reset_axes()
        for k in keys:
            self.ax.plot(z * 1e3, series[k] * scale, "-o", ms=3, label=k)
        self.ax.set_xlabel("dump ⟨z⟩ [mm] (stage-local for linac_sec1 / linac_rest)")
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
        nb = self.nbins_var.get()
        if self.method_var.get() == "histogram":
            h = self.ax.hist2d(xv, yv, bins=nb, weights=P.weight * 1e9, cmap="viridis")
            self.cbar = self.fig.colorbar(h[3], ax=self.ax, label="charge [nC]")
        else:
            order = np.argsort(P.weight)   # heavy macroparticles drawn on top
            sc = self.ax.scatter(xv[order], yv[order], c=P.weight[order] * 1e9,
                                 s=4, cmap="viridis")
            self.cbar = self.fig.colorbar(sc, ax=self.ax, label="charge [nC]")
        if self.equal_var.get():
            self.ax.set_aspect("equal", adjustable="datalim")
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
    BeamGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
