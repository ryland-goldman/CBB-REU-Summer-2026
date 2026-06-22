"""Standalone Tk + matplotlib beam-properties explorer over a run's openPMD dumps.

Trends / 1D / 2D distribution modes over each stage's existing particle series
(nothing re-simulated). A "screen" is one WarpX dump ordered by ⟨z⟩.
See pipeline/README.md -> "Interactive beam explorer" for usage and the
stage-local-z / σ_z-vs-σ_t conventions.
"""

import os
import sys
import threading
import queue
import warnings

warnings.filterwarnings("ignore")

# Run from the repo root so the stage-relative diagnostic paths resolve.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(_ROOT)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

# openpmd-viewer leaks an fd per get_particle; raise RLIMIT_NOFILE so a full-stage browse
# doesn't hit the 256-fd wall. Best-effort.
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

# ── Shared beam-handoff helpers (γβ→eV/c ParticleGroup build) ─────────────────
from pipeline.beam_io import make_particle_group

# ── Stages, in chain order ───────────────────────────────────────────────────
# Dumps store positions [m] and momenta u = γβ. The cathode is 2D (x–z, no y); rest RZ.
STAGES = [
    {"name": "cathode",     "path": "cathode/diags/particles",          "geom": "2d"},
    {"name": "gun",         "path": "gun/diags/particles",              "geom": "rz"},
    {"name": "injector",    "path": "injector/diags/main/particles",    "geom": "rz"},
    {"name": "linac_sec1",  "path": "linac_sec1/diags/main/particles",  "geom": "rz"},
    {"name": "linac_rest",  "path": "linac_rest/diags/main/particles",  "geom": "rz"},
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
# Bunch length is σ_z (NOT σ_t — WarpX dumps are time snapshots).
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
    """One stage's openPMD series: cached ParticleGroups and a ⟨z⟩-ordered screen list."""

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


class BeamGUI:
    def __init__(self, root):
        self.root = root
        root.title("Cornell Linac — Beam Properties Explorer")
        self.stage_data = {}          # name -> StageData (lazy)
        self.q = queue.Queue()        # worker-thread → main-thread results
        self._busy = False
        self._gen = 0                 # monotonic token; a newer _run_async supersedes older (see _drain)
        self._progress_text = ""      # worker-written; reflected to Tk only on the main thread (Tk not thread-safe)

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

        # Play/pause: a root.after timer advances the screen slider, whose own callback redraws.
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
                           command=lambda _: self._on_var_change()).pack(side=tk.LEFT)
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
            self.fixaxes_var = tk.BooleanVar(value=False)
            ttk.Checkbutton(self.ctl, text="Fixed axis range (steady animation)",
                            variable=self.fixaxes_var,
                            command=self._on_fixaxes).pack(anchor=tk.W)

    # ── fixed-axis-range machinery (lock the 2D window across animation frames) ─
    def _fixaxes_on(self):
        return (self.mode_var.get() == "2D Distribution"
                and getattr(self, "fixaxes_var", None) is not None
                and self.fixaxes_var.get())

    def _needed_range_keys(self):
        """The variable keys whose global range the fixed-axis 2D plot needs."""
        return [self._key_by_label[self.x_var.get()],
                self._key_by_label[self.y_var.get()]]

    def _on_var_change(self):
        # New variable ⇒ new locked-window ranges, so recompute off-thread before redraw.
        if self._fixaxes_on():
            self._compute_ranges_async()
        else:
            self.replot()

    def _on_fixaxes(self):
        if self.fixaxes_var.get():
            self._compute_ranges_async()   # precompute the global window, then draw
        else:
            self.replot()                  # back to per-frame autoscale

    def _compute_ranges_async(self):
        """Compute & cache the global 2D-axis range over every screen (off-thread, once)."""
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
        self._stop_play()
        self._refresh_controls()
        # Resolve StageData on the MAIN thread (it reads Tk stage_var); only the heavy
        # screen indexing goes off-thread, so the worker never touches Tk.
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
            return                       # no per-screen frames to animate
        d = self._data()
        if not d.screens:
            return
        self._playing = True
        self.play_btn.config(text="⏸ Pause")
        self._play_tick()

    def _stop_play(self):
        self._playing = False
        if hasattr(self, "_play_job") and self._play_job is not None:
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
        # Setting the slider fires _on_screen_slide, which updates the label and redraws.
        self.screen_scale.set(i)
        delay = max(20, int(self.play_delay.get()))
        self._play_job = self.root.after(delay, self._play_tick)

    def _run_async(self, work, done):
        """Run `work()` off-thread; call `done(result)` on the main thread when finished.

        `work` reports progress only via `self._progress_text` (NEVER touch Tk from the
        worker). Re-entrant: each call bumps `self._gen` and stamps its payload; a later
        call supersedes earlier ones so a stale worker's result is never delivered to the
        newer request's `done` (see _drain).
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
        # Worker thread — `d` resolved on the main thread; progress sets the string only.
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
        # A new stage has its own per-dump ranges, so recompute the locked window.
        if self._fixaxes_on():
            self._compute_ranges_async()
        else:
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
    BeamGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
