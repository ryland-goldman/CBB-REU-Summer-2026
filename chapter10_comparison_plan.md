# Chapter 10 (PARMELA) Figure Reproduction — Implementation Plan

## 0. Summary

Build a new cross-stage plotting module, `pipeline/plot_chapter10.py`, that reads four existing pipeline snapshots (gun exit, injector pre-Preb2, injector handoff, linac_sec1 exit), synthesizes a Tan-style RF-phase axis from each particle's **arrival time** at a common reference plane, and emits four 2-panel figures (KE-vs-phase scatter + peak-normalized charge-vs-phase histogram) plus a comparison table against Tan's Table 10.2.

**Critical scoping correction (reviewer C1 — physics):** The repo default operating point is **NOT** Tan condition (i). The repo default (`injector_sim.py:100,105,102,107`) is the LinacSim "8 kW / 10 kW two-cavity" point — **Preb1 8 kW @ −70°, Preb2 10 kW @ −45°**, both at intermediate crest-referenced offsets. Tan condition (i) is **Preb1 50 kV @ phase-null, Preb2 150 kV @ crest** — different powers *and* different phases, which drive different bunching dynamics and hence every shape Tan describes. Therefore there are **two distinct deliverables**, and the plan must not conflate them:

- **Deliverable A — "repo default" figures (primary, no re-run):** Generate the four 2-panel figures and the table from the existing default diags, labeled **"repo default operating point (8 kW/10 kW two-cavity), qualitative comparison to Tan Fig 10.2–10.5."** Make **no** claim of numeric agreement with Table 10.2; the Tan numbers appear only as a labeled reference column with deltas, annotated as a *different operating point*.
- **Deliverable B — true Tan condition (i) (secondary, requires an injector re-run):** Re-run the injector with Preb-1/Preb-2 powers and phases set to Tan's 50 kV @ phase-null / 150 kV @ crest (with the kV→kW conversion of §1.1), into `injector/diags/cond_i`, then re-point the plotter. Only *these* figures may be compared numerically to Table 10.2. This is the physics-correct cond (i) and is the one that gates any "Table 10.2 agreement" statement.

The plotter takes the four diag paths as arguments, so it serves both deliverables and the optional conditions (ii)–(iv) without code changes.

---

## 1. Scope & Conditions

### 1.1 kV → kW conversion (load-bearing — reviewer C2)

Tan specifies cavities in **kV** (peak gap voltage); the repo parametrizes in **kW** (dissipated power). The relation is set in `injector_sim.py:19`:

```
V_peak = sqrt(1e3 · Q_L · P_kW / (2π f))      # volts, per cavity
⇒  P_kW = (2π f) · V_peak² / (1e3 · Q_L)
```

Use each cavity's loaded Q (`Q_L_1`, `Q_L_2`) and `F_RF=214.18 MHz` to invert. For every condition that specifies a cavity in kV, compute the equivalent `PREB1_KW`/`PREB2_KW`, set it, and **state the achieved kV in the figure caption** (mirror of the Sec1 MV/m↔MW handling in §1.3). Do not write `40 kV → PREB1_KW=…` with no conversion.

### 1.2 Phase-null resolution (reviewer #5 — completeness)

The repo's prebunchers are crest-referenced: crest ≡ `PHI_OFF=0` nominally. Tan's "phase null / not accelerating" = the zero-net-energy-gain (zero-crossing) phase = **±90° from crest**. Concretely set `PREB2_PHI_OFF = -90.0` (or `+90.0`), and **verify** by the printed Preb-2 kick-sign / net-ΔE diagnostic being ≈0. Do not leave this as `<null>`.

### 1.3 Conditions table

Condition (i) is split into the two deliverables of §0. Conditions (ii)–(iv) are **optional follow-on re-runs**, each into a distinct `OUTDIR`, then re-point the plotter.

| Cond | Tan figs | Gun | Preb1 | Preb2 | Sec1 | Knobs / where |
|------|----------|-----|-------|-------|------|---------------|
| (i) [Deliv. B] | 10.2–10.5, Table 10.2 | 150 kV (default) | 50 kV @ phase-null | 150 kV @ crest | 11 MV/m | `injector.config(PREB1_KW=<kV→kW>, PREB1_PHI_OFF=-90, PREB2_KW=<kV→kW>, PREB2_PHI_OFF=0)`; **single injector re-run** |
| (ii) | 10.6/7/8 | 150 kV (default) | 40 kV @ zero | phase-null | 11 MV/m | **two-pass — see §1.4 caveat** (`PREB1_KW`/`PREB2_PHI_OFF=-90`) |
| (iii)| 10.9/10/11 | **175 kV** | 50 kV @ zero | phase-null | 11 MV/m | `gun.config(GUN_VOLTAGE=175e3)` + re-run gun → injector → sec1; `PREB2_PHI_OFF=-90`; single-pass per stage |
| (iv) | 10.12 | **175 kV** | 50 kV @ zero | phase-null | **15 MV/m** | same as (iii) **plus** `linac_sec1.config(POWER_MW=<solve for 15 MV/m>)` |

**Knob-mapping notes (verified against code):**

- **Gun 175 kV** lives in `gun/gun_sim.py:52` (`GUN_VOLTAGE=150e3`), **not** an injector knob. Conditions (iii)/(iv) require `gun.config(GUN_VOLTAGE=175e3); gun.run()` then re-run injector + linac_sec1. This is a full partial-chain re-run.
- **Gun OUTDIR caveat (reviewer — software):** the gun writes to the **fixed** `gun/diags/particles` with no per-condition OUTDIR. Re-running cond (iii)/(iv) **clobbers the cond-(i)/default gun baseline**. The plotter's "diag-path argument" design cannot help, because the path is identical. Therefore: run each gun-voltage condition's full chain as a self-contained unit and plot it **before** moving on, and re-run the default gun afterward to restore the baseline (or add a gun `OUTDIR` knob as a follow-on).
- **Sec1 by power, not gradient (reviewer #6).** The repo parametrizes Sec1 RF by **power (MW)** via `POWER_MW` (`linac_sec1_sim.py:80`), scaled `sc = √(POWER_MW / RF_NORM_MW)`. The build report prints peak gradient = `env.max()·scale` and on-crest gain = `scale·V1KW_KEV` (`build_linac_sec1_field.py:185-187`, with `V1KW_KEV=331.2 keV` at line 66). To hit Tan's literal **15 MV/m**: run a trial `POWER_MW`, read the printed peak gradient, and scale `POWER_MW` by `(15 / printed_gradient)²` (since gradient ∝ √POWER_MW); confirm by re-reading the print. **Document the achieved gradient in the caption.** Do not assume 15 MW.

**Per-condition OUTDIR convention** (so cases don't clobber `diags/main`):
- injector: `injector.config(OUTDIR=f"injector/diags/cond{n}")` (matches `injector.resolve_outdir()`)
- linac_sec1: `linac_sec1.config(OUTDIR=f"linac_sec1/diags/cond{n}")`
- gun: fixed `gun/diags/particles` — see the clobber caveat above.

### 1.4 Condition (ii) two-pass requirement (reviewer #4)

A **Preb-1 power** change (cond ii, 50→40 kV) desyncs the Preb-2 phase reference: Preb-2 is timed analytically from the injection β + Preb-1's kick, so a hardened-Preb-1 study needs a **two-pass run** (see `injector/README.md`). The §1.3 table marks cond (ii) "two-pass — see §1.4 caveat" so no one runs it single-pass and gets a desynced Preb-2. Nulling Preb-2 *itself* (no Preb-1 power change) is safe; cond (i) Deliverable B changes Preb-1 power **and** phase, so it must also follow the two-pass timing fix when matching Tan exactly — apply the README recipe and verify the printed Preb-2 net-ΔE diagnostic.

---

## 2. New Code

### 2.1 New module: `pipeline/plot_chapter10.py`

Cross-stage, runs in-process (no subprocess isolation needed for plotting), from repo root, writes to repo-root `results/`.

**Module-level constants:**
```python
import os
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from openpmd_viewer import OpenPMDTimeSeries
from pipeline._runner import _raise_fd_limit                     # sig: target=16384
from injector.build_injector_field import F_RF as F_214         # 214.1847 MHz — import, don't re-derive
from injector.build_injector_field import Z_GAP_CENTER_2 as Z_PREB2   # 1.318 m
from injector.build_injector_field import Z_HANDOFF                   # 2.03 m
from linac_sec1.linac_sec1_sim import F_RF as F_2856            # 2856.0e6 — import, stays in sync with config()

MC2 = 0.51099895          # electron rest energy [MeV]
C   = 299792458.0
Q_E = 1.602176634e-19
RESULTS = "results"
SPECIES = "electrons"
```
(Import the stage constants rather than hard-coding — extended to `F_2856` per reviewer L3. The `HDF5_USE_FILE_LOCKING` env line is **dropped**: per `project_hdf5_file_locking` memory, locking was a red herring — fd exhaustion was the real cause and is handled by `_raise_fd_limit`.)

**Snapshot reader** — single source for velocity, used by *both* φ and σ so a sign fix can't desync them (reviewers G1/#1):

```python
def _read_snapshot(series_path, iteration=None, target_z=None):
    """Return (z, ux, uy, uz, w, zbar) for one dump (RZ stages only).
       iteration='last' -> last dump with >= 50 live particles (walk reversed; match
                           plot_chain.py's `< 50` skip threshold — one canonical threshold).
       target_z set     -> dump whose charge-weighted <z> is nearest target_z.
       Records and returns the ACTUAL zbar of the chosen dump for the caption.
    """
```
Read idiom (from `plot_chain.py:178-188`):
```python
ts = OpenPMDTimeSeries(series_path)
x,y,z, ux,uy,uz, w = ts.get_particle(
    ["x","y","z","ux","uy","uz","w"], species=SPECIES, iteration=it)
```
For `target_z`: loop `ts.iterations`, skip `len(z) < 50` (wrap in try/except), compute `zbar = np.average(z, weights=w)`, pick `min(|zbar - target_z|)`. All four Tan locations are RZ — no 2D branch.

```python
def ke_mev(ux, uy, uz):
    gamma = np.sqrt(1.0 + ux**2 + uy**2 + uz**2)   # ux/uy/uz are gamma*beta
    return (gamma - 1.0) * MC2
```

**Arrival-time phase (reviewers C3/G1/#1 — the central physics fix).** Tan's x-axis is the **arrival phase at a fixed reference plane**, i.e. drift each particle to the common centroid plane and read the *time* it crosses — NOT a per-particle spatial-to-phase conversion. A particle ahead (larger z) arrives *earlier* (t<0). Use **one** velocity definition shared by φ and σ_t:

```python
def _vz(ux, uy, uz):
    gamma = np.sqrt(1.0 + ux**2 + uy**2 + uz**2)
    return (uz / gamma) * C

def arrival_time(z, ux, uy, uz, w, z_ref=None):
    """t_i = -(z_i - z_ref)/v_z,i  (ahead -> earlier -> t<0). Returns (t, z_ref, good).
       GUARD (reviewer G1): mask |v_z| below a floor (e.g. 1e-3*c) or non-positive
       v_z — at the gun, space-charge tails can have small/negative uz, giving
       |t|->inf. Excluded particles are dropped from the phase axis AND reported
       as a count in the caption."""
    if z_ref is None:
        z_ref = np.average(z, weights=w)
    vz = _vz(ux, uy, uz)
    good = vz > (1e-3 * C)
    t = np.full_like(z, np.nan)
    t[good] = -(z[good] - z_ref) / vz[good]
    return t, z_ref, good

def rf_phase_deg(t, f):
    return 360.0 * f * t      # head (t<0) -> phi<0
```
**Definitional note (reviewer C3):** for the per-*particle* arrival time we keep `v_z,i` (physically correct time-of-flight, well-defined once the small-/negative-`v_z` guard removes the pathological tail). This differs from a *common-β* z→φ map only for the energy-spread populations (the gun sine and before-Sec1 parabola) — and that difference *is* the nonlinear KE(φ) correlation Tan's figures display, so retain it. The local-frame z-reset in `linac_sec1_sim.py:178` (`z = z - z.min() + Z_INJECT`) is **offset-only and monotonic**, so `(z - z_ref)` is frame-invariant; add a one-line comment asserting this so a future editor does not "fix" it.

```python
def sigma_t(t, w, good):
    """Charge-weighted RMS arrival-time spread [s] on the SAME t/good as the phase axis."""
    tt, ww = t[good], w[good]
    m = np.average(tt, weights=ww)
    return np.sqrt(np.average((tt - m)**2, weights=ww))

def wstats(v, w):
    """Canonical charge-weighted (mean, std) — the plot_chain.py idiom."""
    m = np.average(v, weights=w)
    return m, np.sqrt(np.average((v - m)**2, weights=w))
```

```python
def location_panel(ax_top, ax_bot, phi, ke, w, title="", phi_label="@214 MHz",
                   ke_window=None, tan_overlay=None):
    """Top: KE[MeV] vs phi[deg] scatter. Bottom: PEAK-normalized charge histogram vs phi
       (sharex). ke_window=(lo,hi) optional KE gate (after-Sec1 captured core)."""
```
- Top: `ax_top.scatter(phi, ke, s=2, c="C3", alpha=...)`; ylabel `"Kinetic Energy (MeV)"`.
- Bottom (**peak**-normalized, not `density=True`): `counts,edges = np.histogram(phi, bins=120, weights=w); counts = counts/counts.max()`; `ax_bot.step(...)`. ylabel `"normalized charge distribution"`.
- xlabel `f"RF phase φ (deg) {phi_label}"`.

```python
def main(case,                    # e.g. "repo_default" or "cond_i"
         gun_diag="gun/diags/particles",
         inj_diag="injector/diags/main/particles",
         sec1_diag="linac_sec1/diags/main/particles",
         sec1_summary="linac_sec1/diags/main/injection_summary.json"):
    _raise_fd_limit()
    os.makedirs(RESULTS, exist_ok=True)
    # read 4 snapshots, build 4 figures, build table, write all.
```

### 2.2 The four location reads (exact)

| Tan loc | Source | Selection | φ freq | z_ref / σ population |
|---------|--------|-----------|--------|----------------------|
| **A at Gun** (10.2) | `gun/diags/particles` | `'last'` (last live dump, `< 50` skip) | `F_214` | full beam; report gun σ_t and ⟨KE⟩ in caption (repo ~146–148 keV) |
| **B before Preb2** (10.3) | `injector/diags/main/particles` | nearest `zbar` to `Z_PREB2` | `F_214` | full beam |
| **C before Sec1** (10.4) | `injector/diags/main/particles` | nearest `zbar` to `Z_HANDOFF` | `F_214` | **pre-scrape full injector population** (Tan's "before Sec1" is pre-capture — do NOT apply the iris scrape here) |
| **D after Sec1** (10.5) | `linac_sec1/diags/main/particles` | `'last'` (`its[-1]`) | `F_2856` | **captured core**: KE ≥ 12 MeV (`MIN_KE_MEV`) or within a few σ of ⟨KE⟩≈25 MeV; read **before** any downstream 12 MeV handoff cut |

**Plane-selection caveat (reviewer C2):** injector `N_DIAGS=60` over `ZMAX=2.10 m` ⇒ ~35 mm/dump spacing, so the nearest dump to `Z_HANDOFF` can be ±17 mm off, and the beam is converging hard through 1.922→2.03 m. Pick the nearest dump, **report its actual `zbar` in the caption**, treat the plane as approximate. Likewise confirm the gun dump being read sits at the same longitudinal location Tan samples, or σ_z@214 differs by the drift-bunching factor.

**Capture % for the table (reviewers G2/#3/#7):** report all three summary denominators:
- `injection_summary.json` exposes `q_injected_C` (all-r, pre-iris injected — honest denominator), `q_in_bore_C` (post-iris), `q_in_domain_C`.
- **all-buckets** = `w.sum()·Q_E` at Sec1 exit / `q_injected_C`.
- **in-bucket** = charge within ±180°@2856 of the captured-core centroid AND inside the KE window / `q_injected_C`.
- Also report **in-bucket / q_in_bore_C** so iris loss (~68%) is separable from capture loss.
- Report as `"in-bucket / all-buckets"` like Tan, but **explicitly footnote that the denominators differ from Tan's**: Tan's 100% (at C) and 89.4/96.8 (at D) reference the **gun-emitted** bunch with no loss yet, whereas the repo's `q_injected_C` already sits past the converging halo and already includes the real 9.547 mm iris scrape. So the repo pair (~7%/~32%) is **much lower than Tan's by construction** — γ² self-field pessimism **plus** real iris loss **plus** (for repo-default) a different operating point. Annotate; do **not** "fix."
- **Location C capture column = "n/a (pre-iris)"**: C is pre-scrape; C's charge population differs from D's (injector vs iris survivors) — different populations, C→D ratio is **not** a clean repo fraction.

---

## 3. Figure Spec

Four 2-panel figures, `figsize=(7.0, 7.0)`, `constrained_layout=True`, `dpi=140`, two rows sharing x (top scatter, bottom histogram). Each saved to repo-root `results/`. Filenames carry the `case` so Deliverable-A and Deliverable-B figures coexist.

| Output PNG (case-tagged) | Tan fig | Location | Top y | Bottom y | x-axis |
|--------------------------|---------|----------|-------|----------|--------|
| `results/tan_fig10p2_at_gun_<case>.png` | 10.2 | A Gun | KE (MeV) | normalized charge | RF phase (deg) @214 MHz |
| `results/tan_fig10p3_before_preb2_<case>.png` | 10.3 | B | KE (MeV) | normalized charge | @214 MHz |
| `results/tan_fig10p4_before_sec1_<case>.png` | 10.4 | C | KE (MeV) | normalized charge | @214 MHz |
| `results/tan_fig10p5_after_sec1_<case>.png` | 10.5 | D | KE (MeV) | normalized charge | @2856 MHz |

(Optional combined `results/tan_chapter10_panels_<case>.png` — 4-wide scorecard analogous to the existing `chain_*` figures.)

**Reference-overlay:** Tan's figures are dissertation scans, not digitized data, so a true overlay is not feasible without digitizing. **Deliver clearly-labeled standalone figures** whose axes/units/shape match Tan's, with a caption stating which Tan figure to compare and the expected shape (sine → S-curve → left-opening parabola → falling-KE spike). `location_panel` carries an optional `tan_overlay=(phi,ke)` hook for a future digitization — leave it unused, do not block on it.

**Caption per figure (mandatory):**
- **Operating point** — Deliverable A: "repo default 8 kW/10 kW two-cavity (Preb1 −70°, Preb2 −45°), **qualitative** comparison only — NOT Tan cond (i)." Deliverable B: "Tan cond (i): Preb1 50 kV @ phase-null, Preb2 150 kV @ crest; achieved kV = …."
- Sign convention (head = φ<0); φ=0 = charge-weighted centroid (captured-core centroid for D).
- Actual ⟨z⟩ (`zbar`) of the chosen snapshot; count of particles dropped by the `v_z` guard.
- Reference frequency (214.18 MHz, not 214.0 / 2856 MHz).
- Model caveats: charge **peak-normalized** (Tan 1e11 e⁻ ≈16 nC vs repo ~0.8 nC — shapes compare, not absolute pC); repo gun σ_t (state explicitly — Tan's gun base is 3.7 ns); non-relativistic self-field (~γ² pessimism → conservative capture); real 9.547 mm iris scrape (~32%); after-Sec1 from the Sec1 exit dump (before the 12 MeV `linac_rest` handoff cut).

---

## 4. Comparison Table Artifact

`main` writes, per `case`:
- `results/tan_comparison_<case>.md` (human-readable) and `results/tan_comparison_<case>.csv` (machine-readable).

Columns (one row per location): `location, zbar_m, Ebar_MeV, sigE_MeV, sigz_deg@214, sigz_deg@2856, capture_in_bucket_pct, capture_all_buckets_pct` — plus a parallel **Tan published** row and a **delta/ratio** column.

`sigz@214 = 360·F_214·σ_t`, `sigz@2856 = 360·F_2856·σ_t` (same `σ_t`). After Sec1, leave @214 blank (Tan does) and compute @2856 on the captured core. Location C capture = "n/a (pre-iris)".

**σ_z internal-consistency vs Tan (reviewer C3 — important):** because both σ_z columns derive from one `σ_t`, the repo @2856/@214 ratio is *by construction* `F_2856/F_214 = 13.335`. **Tan's own published ratio is 13.52** (428.8/31.7, 194.8/14.4, 75.3/5.57 all ≈13.52) — a consistent ~1.4% offset, meaning PARMELA's two σ_z columns are NOT one σ_t scaled by two frequencies (likely a spatial σ_z converted with a per-location β). Therefore the 13.335 check validates **our code's internal consistency only**; it is **not** agreement with Tan and **must not** be presented as such. Footnote the table.

Hard-code Tan's Table 10.2 (condition i) reference block:
```
Gun:        Ebar 0.150, sigE 0.000, 31.7/428.8,  100%
prePreb2:   Ebar 0.139, sigE 0.023, 14.4/194.8,  100%
preSec1:    Ebar 0.253, sigE 0.043, 5.57/75.3,   100%
postSec1:   Ebar 27.2,  sigE 3.5,   -/11.1,      89.4 (in-bucket) / 96.8 (all)
```
(Tables 10.3–10.5 for conditions ii–iv: leave a dict keyed by condition; fill from the dissertation when run.)

---

## 5. Documentation Updates (each exact edit)

1. **`FIGURES.md`** — under `## 0. Cross-stage — results/`, add one `###` subsection per new PNG (`### tan_fig10p2_at_gun_<case>.png` … `### tan_fig10p5_after_sec1_<case>.png`, plus the panel scorecard if added), each with the `![alt](results/<name>.png)` embed and 1–3 sentences: location, Tan figure compared, expected shape, the **operating-point** label (default vs cond i), the peak-normalized / sign / caveat notes. Add `### tan_comparison_<case>.md` describing the table artifact.
2. **Root `README.md`** — add a component/script row for `pipeline/plot_chapter10.py`. Note explicitly that default-diags figures are a **qualitative** comparison and that true Tan cond (i) needs an injector re-run (Deliverable B).
3. **`.claude/CLAUDE.md`** — under **Commands**, add `python pipeline/plot_chapter10.py`. Note the **default ≠ Tan cond (i)** distinction, the kV→kW conversion, the (ii)–(iv) re-run recipe, the OUTDIR convention, the gun-diag clobber caveat, and the cond-(ii)/Deliverable-B two-pass Preb phasing requirement. Recommend keeping it a **separate manual command**, not in the default chain.
4. **Stage READMEs** — one-line cross-reference in `injector/README.md` and `linac_sec1/README.md`.
5. **`reference/Papers/README.md`** — if the Tan dissertation (Cheng-Yang Tan, "The CESR Injector", 1997) is not already indexed, add it (the comparison target).

---

## 6. Validation

**Scope note:** quantitative Table-10.2 comparison applies **only to Deliverable B** (true cond i). Deliverable A (repo default) is validated on **shape and internal consistency only**.

**Quantitative (Deliverable B vs Table 10.2):**
| Quantity | Tan (i) | Acceptance |
|----------|---------|------------|
| Gun ⟨E⟩ | 0.150 MeV | repo ~0.146 MeV; within ~5% ✓ |
| Gun σ_E | 0.000 | repo small but nonzero; ≤ ~0.01 MeV |
| Gun σ_z@214 / @2856 | 31.7 / 428.8 | **first verify the repo cathode emission window**; internal @2856/@214 ratio must be 13.335 (our code), **not** Tan's 13.52 |
| prePreb2 ⟨E⟩ | 0.139 | within ~10% |
| preSec1 ⟨E⟩ | 0.253 | within ~15% |
| postSec1 ⟨E⟩ | 27.2 MeV | repo ⟨KE⟩≈25 MeV; within ~15% ✓ |
| postSec1 σ_E | 3.5 | repo σ_KE larger — **expected**; flag, don't fail |
| postSec1 σ_z@2856 | 11.1° | within ~factor 2 |
| capture in/all | 89.4 / 96.8 | repo ~7% / ~32% — **expected far lower**; annotate, NOT a failure |

**Qualitative (shape — the real test, both deliverables):**
- **A Gun:** sine-like KE-vs-φ, φ<0 side at higher KE; histogram ~gaussian. Also validates sign convention.
- **B before Preb2:** monotonic S-curve (KE vs φ); histogram peakier (bunched).
- **C before Sec1:** left-opening parabola / "C"; sharp charge spike near φ≈0 with a tail.
- **D after Sec1:** steeply falling KE-vs-φ near 25–30 MeV; captured core a narrow spike.

**Programmatic cross-checks (in code):**
- **Sign convention:** at the gun, assert `corr(φ, KE) < 0`. Fail loudly if positive — replaces a manual "if mirrored, negate" loop with an assertion.
- **σ_z internal ratio:** `σ_z@2856 / σ_z@214 == 13.335` at every upstream location (shared `σ_t`). Internal check — do not compare to Tan's 13.52.
- **φ=0 on the histogram peak** (centroid reference).
- **`v_z`-guard count** is small (≪1% of macroparticles).

**Tolerances rationale:** energies/shapes should match closely for Deliverable B. σ_E and capture% diverge in a **documented direction** (repo = conservative lower bound on capture, broader σ_E). Any figure where the **shape** is wrong (parabola opens the wrong way, no bunching at Preb2) is a real bug.

---

## 7. Task Breakdown (ordered, self-contained)

1. **Confirm baseline diags exist.** If `gun/diags/particles`, `injector/diags/main/particles`, `linac_sec1/diags/main/particles`, or `linac_sec1/diags/main/injection_summary.json` are absent/stale, run `python pipeline/run_pipeline.py`. Verify each opens with `OpenPMDTimeSeries`.
2. **Create `pipeline/plot_chapter10.py` skeleton + helpers.** Constants (imported from stage modules), `_read_snapshot`, `ke_mev`, `_vz`, `arrival_time` (with the `v_z` guard), `rf_phase_deg`, `sigma_t`, `wstats`. Call `_raise_fd_limit()` first in `main`. Unit-test helpers on the gun exit dump.
3. **Implement `location_panel` + the four reads (Deliverable A — repo default).** Wire the four snapshots (§2.2) with correct per-location frequency, peak-normalized histograms, captured-core KE window for D, actual-`zbar` captions. Build the four 2-panel figures (§3) with Tan axis labels and the **"repo default, qualitative"** caption.
4. **Verify the sign convention at the Gun (gates all figures).** Render the gun figure; assert `corr(φ, KE) < 0`. If positive, negate in `rf_phase_deg` and update the caption.
5. **Implement capture bookkeeping + comparison table.** Read `q_injected_C`, `q_in_bore_C`, `q_in_domain_C`; compute in-bucket, all-buckets, in-bucket/in-bore. Write `tan_comparison_repo_default.md`/`.csv` with repo, Tan-reference, delta rows; C capture = "n/a (pre-iris)"; footnote σ_z and denominator caveats.
6. **Deliverable B — true Tan cond (i).** Compute kV→kW (§1.1); re-run injector with Preb1 50 kV @ phase-null (`PREB1_PHI_OFF=-90`, two-pass per §1.4) and Preb2 150 kV @ crest into `injector/diags/cond_i`; re-run linac_sec1 from it into `linac_sec1/diags/cond_i`. Re-point `main(case="cond_i", ...)`. Generate the four `_cond_i` figures + `tan_comparison_cond_i.*`. **Only these may be compared to Table 10.2.**
7. **Validation pass (§6).** Deliverable A: shapes + internal cross-checks. Deliverable B: shapes + quantitative table within tolerances (after verifying the cathode emission window). Fix wrong-shape bugs.
8. **Documentation sync (§5).** Update `FIGURES.md`, root `README.md`, `.claude/CLAUDE.md`, the two stage READMEs, `reference/Papers/README.md` — same change as the code.
9. **Commit.** Commit `pipeline/plot_chapter10.py` + doc edits; `git add -f results/tan_*.png results/tan_comparison_*.{md,csv}`. Do not commit `diags/`, `.h5`, or logs.
10. **(Optional) Conditions (ii)–(iv).** `main(case="cond{n}", ...)`. Re-run recipe (§1.3); gun-voltage conditions respect the clobber caveat; cond (ii) two-pass. Follow-on, non-blocking.

---

**Key files referenced:** `pipeline/plot_chapter10.py` (new), `pipeline/_runner.py` (`_raise_fd_limit`, `target=16384`), `pipeline/beam_metrics.py` (weighted-moment idiom), `pipeline/plot_chain.py:178-188` (read idiom), `gun/plot_gun.py:52,112`, `gun/gun_sim.py:52` (`GUN_VOLTAGE=150e3`), `injector/injector_sim.py:19,100,102,105,107` (kV↔kW, default 8 kW/−70°/10 kW/−45°), `injector/build_injector_field.py` (`Z_GAP_CENTER_2=1.318`, `Z_HANDOFF=2.03`, `F_RF=214.18 MHz`), `linac_sec1/linac_sec1_sim.py:56,80,178` (`F_RF=2856 MHz`, `POWER_MW=11`, local-frame z-reset), `linac_sec1/build_linac_sec1_field.py:66,185-187`, `linac_sec1/diags/main/injection_summary.json` (`q_injected_C`, `q_in_bore_C`, `q_in_domain_C`).
