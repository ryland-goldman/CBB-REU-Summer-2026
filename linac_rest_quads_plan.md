# linac_rest — Implement the Quads (FODO Focusing) Plan

> **⚠️ SUPERSEDED IN IMPLEMENTATION.** This is the pre-revision design. The single-quad,
> thin-lens, μ=72° approach below was found unable to focus both transverse planes and was
> **replaced by a thick-lens H/V doublet at μ=50°** during implementation. The gradient table
> here (~0.1–0.35 T/m) is NOT the as-built value (shipped ~1–2.7 T/m). For the as-built recipe
> see `linac_rest/README.md` → *Space charge & quads* and `fodo_quad_gradients`/`_solve_doublet_k1`
> in `linac_rest/build_linac_rest_lattice.py`. Kept as a historical record of the planning phase.

## Goal & scope

Turn the currently-inert quadrupoles (present at real tabulated lengths but `b1_gradient = 0`) into a working FODO focusing lattice that actually contains the beam through Cornell linac sections 2–8. The headline run **stays quads-OFF** (energy result ≈308 MeV is the robust, quad-independent deliverable; transmission ≈78.5% remains the no-focusing LOWER BOUND). `QUADS_ON=True` becomes a real, contained-beam **exploratory** path: per-quad `b1_gradient` values are derived from accelerator optics (constant-phase-advance FODO design, scaled by local beam momentum) — NOT from measured quad current, since the A→T (current→field) calibration is undocumented (`details.md`, A→T-unknown / "correct order of magnitude" note ≈ line 188–194). Every QUADS_ON output stays labeled "placeholder optics — guessed-strength FODO, A→T undocumented." The deliverable is a focused, bounded σ_x(z) / σ_y(z) envelope and a higher transmission against the **real tapered bore** (the box `XYRAD_M` is never widened — un-gameable-transmission ruling).

**State of the code this plan builds on (verified, corrects the draft):**
- The **real tapered bore aperture is ALREADY ON** for the quads-OFF headline. `BORE_APERTURE_ON = True` (`build_linac_rest_lattice.py:247`), `bore_aperture_on = bool(BORE_APERTURE_ON or quads_on)` (`:342`), and the solrf body already takes `radius = section_bore_radii(index)[0]` (`:263,:284`). The published 78.5% lower bound is produced *against this bore today*. This plan therefore does **NOT** "enable the bore" — that work is already done. The **only genuinely-new aperture** is the **quad element `radius`** (`:371`, hard-coded `0.0`) and the inter-section drift `radius` (`:365,:374`).
- There is exactly **one** `build_impact(...)` call in the sim (`linac_rest_sim.py:284`), and that single `I` is **both** calibrated (`calibration.calibrate_sections(I, …)` `:299`) **and** run (`_run_with_progress(I, …)` `:304`). `calibrate_sections` does `Ic = I.copy()` per section (`calibration.py:224`) and writes fitted scale/phase back onto the live `I`. So "calibrate quads-OFF then enable" requires a real restructure (Task 4), not a post-hoc tweak.

## Physics: the chosen K1 / b1_gradient recipe

### Units & conversion (load-bearing)

Impact-T quadrupole V2 = `b1_gradient` = field gradient **Bg in T/m** (confirmed in `impact/parsers.py:603-619` `quadrupole_v` writes `v=[ele,zedge,b1_gradient,file_id,radius]`; `8.10_quadrupole.md:3-9`: `By = Bg·x`, `Bx = Bg·y`), NOT geometric K1. lume-impact passes it verbatim. Convert from optics:

```
K1 [1/m^2] = Bg / (Bρ)
Bρ [T·m]   = p/e = √(KE·(KE + 2·mc²)) / c                       (electron, mc²=0.511 MeV)
         ⇒ Bρ [T·m] ≈ 3.3356 · p[GeV/c],   p[GeV/c] = √(KE_GeV·(KE_GeV + 0.001022))
⇒  Bg [T/m] = b1_gradient = K1 · Bρ
```

### Chosen approach — constant-phase-advance FODO, energy-scaled gradient

Two recipes were proposed: (A) constant per-cell phase advance μ with energy-scaled gradient, vs. holding a fixed Bg. **Choose (A).** Justification: `b1_gradient` is a fixed T/m per element, but Bρ rises **~5×** across the placed quads (0.18→0.88 T·m at the quad locations). A fixed Bg therefore gives K1 that *falls* as the beam accelerates, driving the downstream cells toward the under-focused end. The standard linac-FODO rule holds the per-cell phase advance μ roughly constant by raising Bg with the local beam momentum. This is robust (needs no injected-Twiss match) and is the design-intent reading of "BMAD generic quad function" (`details.md`).

**Modeling honesty (do NOT over-claim).** Two approximations are baked in and must be stated on every QUADS_ON output:
- **RF-section-as-drift.** The thin-lens FODO phase-advance relation below is derived for a *periodic* cell with field-free drifts. Here the "drift" between quads is a **multi-meter accelerating RF section** (2.94–5.15 m), the half-cells are unequal, and the beam accelerates (adiabatic damping) through them. So the derived Bg gives a **nominal** μ that the real lattice will **not** exactly reproduce. The acceptance is **boundedness** of σ_x(z)/σ_y(z), NOT a measured 72° per cell.
- **Single collapsed quad ≠ H/V doublet.** Every real machine quad in `details.md` (≈ line 184–191) is an **H/V doublet** — the tabulated length (e.g. "4+3+4" in) is the *whole assembly* (focus-pole + spacer + defocus-pole), which nets focusing in *both* planes. This plan models each gap with a **single thick quad of one sign** spanning the full tabulated length, alternating sign gap-to-gap. That is a real FODO of thick lenses, but it is a **coarser approximation than a doublet** and it focuses one plane / defocuses the other *per element*. The collapse, on top of the missing A→T, is an explicit caveat on every figure and summary.

**Formulas** (thin-lens FODO, **field-free half-cell drift** `L_h`):

```
1/f      = K1 · L_q                              (thin-lens integrated focal)
1/f      = (2 / L_h) · sin(μ/2)                  (FODO per-cell, μ = phase advance/cell)
⇒ K1_i  = 2·sin(μ/2) / (L_h,i · L_q,i)          [1/m²]
⇒ Bg_i  = K1_i · (Bρ)_i                          [T/m]   (energy-scaling enters here)
```

`L_h` is the **field-free half-cell** between lenses, so the quad's own length is **excluded** (thin-lens drift is *between* lenses): `L_h,i = 0.5·(L_section(i+1) + gap)`, `gap = DRIFT_M = 0.4 m` (corrects the draft, which included `L_q` and ~halved `L_h`, doubling the realized μ). Bρ_i is computed from the **per-section exit energy** `KE_in + Σ_{j≤i} ΔE_target,j` (deterministic from `section_de_target`, `build_linac_rest_lattice.py:113-119`).

`KE_in` is the **measured** sec-1 handoff energy (the sim reads it at runtime as `ke_in` in `load_sec1_core`, `linac_rest_sim.py:150` — there is **no module constant**). The sim therefore passes the measured `ke_in` into `fodo_quad_gradients(ke_in_mev=…)`; the table below uses the ≈25 MeV nominal and is approximate (the helper's standalone default).

Choose **μ = 72°/cell** (nominal design value, inside the 0–180° stability window with margin).

**Sign pattern (FODO):** for an electron, `By = Bg·x` gives `Fx = +e·vz·Bg·x` ⇒ **negative Bg focuses in x, positive Bg focuses in y** (matches lume-impact's own convention `impact/lattice.py:208,216-218`, `q_sign=-1`). FODO alternates sign across the 6 placed quads: `[+, −, +, −, +, −]` (global start sign is a free 1-bit choice; the alternation is what matters). Encode the sign directly in `QUAD_K[i]` (signed T/m). **Verify the sign empirically in the smoke run** (σ_x should *decrease* through the first negative-Bg quad), do not trust the hand-derivation alone.

### Per-quad values (μ = 72°, **field-free L_h**, recomputed)

There are **6 placed quads** (`quad2…quad7`, one after every section except the last — `build_linac_rest_lattice.py:354,367-368`). The quad after section *i* uses `section_quad_length_m(i)`. The table is regenerated by the helper (Task 1); the **acceptance is the helper's own printed output**, not hand-tuned magic numbers. With `L_q` dropped from `L_h`, `L_h` ≈ doubles vs. the draft, so K1 (∝ 1/L_h) ≈ halves and Bg lands lower (~0.1–0.35 T/m). Indicative values:

| deck quad | after sec | L_q [m] | exit KE [MeV] | Bρ [T·m] | b1_gradient [T/m] (sign) |
|-----------|-----------|---------|---------------|----------|---------------------------|
| quad2 | sec2 (CEA 2) | 0.2794 | ~53.3 | ~0.179 | **+** |
| quad3 | sec3 (CEA 3) | 0.4572 | ~81.5 | ~0.274 | **−** |
| quad4 | sec4 (CU 5) | 0.6350 | ~125.2 | ~0.419 | **+** |
| quad5 | sec5 (CEA 4) | 0.4115 | ~172.3 | ~0.576 | **−** |
| quad6 | sec6 (CEA 5) | 0.5588 | ~219.4 | ~0.734 | **+** |
| quad7 | sec7 (CU 3) | 0.5588 | ~263.1 | ~0.879 | **−** |

Exit-KE column verified self-consistent with `section_de_target` @11 MW (cumulative ΣΔE after sec2 ≈ 28.3 ⇒ ~53.3; after sec7 ⇒ ~263.2). All gradients stay 0.1–0.35 T/m — modest, physically reasonable pole-tip fields for a ~cm-bore S-band linac quad ("correct order of magnitude" per `details.md`). If a quad over-pinches in a smoke run, cap `|K1|` (≈1.5) or lower μ for that quad only (tuning knob, document it). The 7th SECTIONS entry (Q8, after the final section 8) is **never placed** — no quad follows the last section.

> **Drop the contradictory "~12× Bρ" sentence.** The energy-scaling is applied at the per-section exit (0.18→0.88 T·m ≈ 4.9×). Because Bg is fixed per element while Bρ keeps rising ~20–40% *across* each downstream RF section, the realized μ sags below nominal between quads — this is part of the "nominal, not realized" caveat above, not a separate 12× claim.

## Implementation tasks (ordered)

> The deck keeps **one** collapsed quad per gap (the real machine H/V doublet is approximated by a single thick quad with one sign per element, alternating gap-to-gap — see "Single collapsed quad ≠ H/V doublet" above; the existing plumbing supports this with no structural change). `file_id=0` hard-edge is kept (effective length = geometric L); do NOT set `0<V3<100` (that reinterprets V3 as Enge effective length — `6_lattice_elements.md` V3 rule, `parsers.py:603-619` maps the dict keys `b1_gradient`/`file_id`/`radius` → V2/V3/V4 correctly).

**1. Add the FODO gradient helper** — `build_linac_rest_lattice.py`, new function after `section_quad_length_m` (**line 143**, i.e. immediately after its `return` and before `total_rf_length_m`).
- Add `def fodo_quad_gradients(*, phase_adv_deg=72.0, mc2_mev=0.510998950, ke_in_mev=25.0):` returning a **length-`N_SECTIONS` (7)** list of signed `b1_gradient` [T/m], with the **last (Q8) entry `0.0`** so the shape matches the existing `[0.0]*N_SECTIONS` default and no length asymmetry can cause an IndexError downstream. The first 6 entries are the placed quads. Internally: cumulative section exit energies from `section_de_target` (lines 113-119), `Bρ_i = √(KE·(KE+2·mc2))/c`, **`L_h,i = 0.5·(L_section(i+1)+DRIFT_M)`** (field-free, quad length excluded), `K1_i = 2·sin(μ/2)/(L_h,i·L_q,i)`, `Bg_i = K1_i·Bρ_i`, sign `(-1)**i` (start +). `ke_in_mev` defaults to the ≈25 MeV nominal for standalone calls; the sim passes the **measured** `ke_in`.
- **Acceptance:** `python -c "import build_linac_rest_lattice as L; print(L.fodo_quad_gradients())"` prints 7 values — 6 alternating-sign T/m in the ~0.1–0.35 range plus a trailing `0.0`.

**2. Wire the FODO fallback into `build_impact`** — `build_linac_rest_lattice.py:337`.
- Replace `quad_k = quad_k or [0.0] * N_SECTIONS` with explicit `is None` logic: if `quad_k is not None`, use it; elif `quads_on`, `quad_k = fodo_quad_gradients()`; else `quad_k = [0.0]*N_SECTIONS`. Use `is None` (not `or`) so an explicit all-zero list isn't silently discarded; add `assert len(quad_k) >= N_SECTIONS - 1`. Keep `quad_k[i]` indexing at `:369` (only `i<N_SECTIONS-1` consumed). Ensure the `quads_on=False` branch is byte-for-byte unchanged.
- **Acceptance:** with `quads_on=False`, `quad_k` stays all-zeros and `b1_gradient` is 0.0 on every quad; with `quads_on=True, quad_k=None`, the 6 derived values land in `I.ele["quad2..7"]["b1_gradient"]`.

**3. Add the NEW quad + drift element bore aperture — gated on `quads_on`, NOT `bore_aperture_on`** — `build_linac_rest_lattice.py:365,371,374`.
- Replace the hard-coded `"radius": 0.0` on the quad (`:371`) and the two inter-section drifts (`:365,:374`) with `section_bore_radii(i)[1] if quads_on else 0.0` — the section **exit** bore radius, because the quad sits *downstream* of the section exit at the narrower taper (the solrf body already uses the **entrance** radius `[0]` at `:263`; entrance-on-body / exit-on-quad is an intentional taper, not an inconsistency — state this in the code comment).
- **Gate on `quads_on`, NOT `bore_aperture_on`** (corrects the draft). `bore_aperture_on` is already `True` for the quads-OFF headline, so gating the new quad/drift radius on it would add a never-before-present scrape plane at every inter-section quad and **change the published 78.5%** — a silent headline regression. Gating on `quads_on` keeps the quads-OFF path byte-identical (`radius` stays `0.0`) and adds the new aperture only on the focused exploratory path. (Impact-T ignores drift `radius` per `parse_drift`, so the drift change is harmless/forward-looking; the **quad** radius is the load-bearing one — element `radius>0` is a real loss aperture, `2_recently_added_features.md:9`.)
- **Acceptance:** with `quads_on=True`, `I.ele["quad3"]["radius"]` equals `section_bore_radii(1)[1]` (>0); with `quads_on=False`, `I.ele["quad3"]["radius"]` is `0.0` (unchanged from today) and the quads-OFF transmission is **bit-for-bit the published 78.5%**.

**4. Restructure build: calibrate on a quads-OFF deck, then build a FRESH quads-ON final deck** — `linac_rest_sim.py:284-304`.
- The sim currently builds one `I` with `quads_on=QUADS_ON` and both calibrates and runs it. Restructure to:
  1. Build the **calibration deck** with `quads_on=False` unconditionally: `I_cal = L.build_impact(..., quads_on=False, quad_k=None)`. Calibrate on it (`cal.calibrate_sections(I_cal, …)`) — energy gain is transverse-independent **only on-axis**, so a quads-OFF, no-bore-scrape calib bunch keeps the `mean_energy` fit clean and prevents the decimated `Np_calib=400` bunch from being focused/scraped mid-fit (which would remove particles from the mean and shift gates 1/2). Collect the fitted `scales` (and phases) from the calibration result.
  2. If `QUADS_ON`: build a **fresh** final deck `I = L.build_impact(..., quads_on=True, quad_k=(QUAD_K or fodo_quad_gradients(ke_in_mev=ke_in)), scales=calibrated_scales)` (note `build_impact` already accepts `scales=`). Pass the **measured** `ke_in` into `fodo_quad_gradients`. Set `I.initial_particles = P_in`, `I.configure()`, and **assert `I.ele["quad3"]["b1_gradient"]` is non-zero** before `_run_with_progress`. If not `QUADS_ON`, the run deck is just the quads-OFF `I_cal` (headline unchanged).
  - This avoids the in-place-mutate-then-`configure()` ambiguity (whether `configure()` re-reads mutated `ele` dicts into the deck is not relied upon): the final deck is assembled fresh from `build_impact` with the calibrated `scales`, so what runs is exactly what `build_impact` produced.
  - **The phase reference must carry over.** Calibration finds each section's absolute crest phase; the fresh final deck must apply those same per-section phases (not just scales). Apply the calibrated phases to the final `I` (the calib result carries both scale and phase per section) so the fresh deck is on-crest identically to the calibrated one.
- **Acceptance:** gate 1 (per-section ΔE ±3%) and gate 2 (exit ≈308 MeV ±3%) PASS **identically** for QUADS_ON and QUADS_OFF (a *consequence* of calibrating on the quads-physically-absent deck — verify by running both); the QUADS_ON final track shows non-zero `b1_gradient` on the quads and non-zero quad `radius`.

**5. Record optics in the summary** — `linac_rest_sim.py:357-381`.
- Add `quad_k` (the applied per-quad T/m list actually placed), `quad_phase_adv_deg` (72.0), and keep existing `quads_on`/`bore_aperture_on`/`transmission_core`. Preserve the count-before-recharge ordering (`linac_rest_sim.py:309-330`) — do NOT move charge re-imposition before the `n_out/n_in` count.
- **Acceptance:** `injection_summary.json` (QUADS_ON) contains `quad_k` (the applied values) and `quad_phase_adv_deg=72.0`.

**6. Add `sigma_y` to the stat table, then the envelope-in-bore soft gate** — `linac_rest_sim.py:187-192` (stat) and `calibration.py` (~line 405 add, ~422 print; NOT in the `426-439` asserts).
- First add `"sigma_y_m": I.stat("sigma_y")[idx].tolist()` to `_stat_vs_z` (`linac_rest_sim.py:187-192`) — needed both for the gate and for the plane-asymmetry plot (Task 7). Confirm `I.stat("sigma_y")` is a valid key on first run (it is the expected counterpart to `sigma_x`); if it raises, fall back to `sigma_x` and note it.
- Compute `gates["envelope_in_bore"] = (BEAM_EDGE_SIGMA * max(σ_x, σ_y)) < min(section bore radius)` with `BEAM_EDGE_SIGMA = 3.0` (σ is RMS, not the beam edge; a 3σ multiple is the meaningful edge-vs-bore check). Print as diagnostic PASS/FAIL. Keep **soft** (K1 unvalidated — must not gate the energy headline).
- **Do NOT claim "FAIL when quads OFF" as a liveness test** (corrects the draft): the quads-OFF beam already transmits ~78.5% with adiabatic damping, so its *RMS* σ may stay below the bore even while distribution tails scrape — the gate could legitimately print PASS quads-OFF. Liveness is instead confirmed by the figure (Task 7) showing bounded oscillation with quads ON vs. monotonic RMS rise with quads OFF.
- **Acceptance:** with the derived FODO, `envelope_in_bore` prints PASS and the gate is exercised (printed, not asserted); the value is recorded in the summary.

**7. Add σ_y and retitle plots for the derived-FODO path** — `plot_linac_rest.py:156,185,188-194`.
- Add a **σ_y(z)** curve to the `fodo_optics.png` panel (currently plots only `sx` at `:185`) using the new `sigma_y_m` — the "x & y out of phase" / plane-asymmetry argument is unobservable without it.
- `fodo_optics.png` title (already branches on `summ["quads_on"]`, `:188-191`): change "(QUADS ON — exploratory)" → "(QUADS ON — derived energy-scaled FODO, nominal μ=72°; single-quad collapse of H/V doublet; A→T calib missing)".
- `emittance.png` title (`:156`) "quads OFF ⇒ RF + drift only" → make it conditional on `quads_on` (quad chromaticity grows εn — see gates below).
- **Acceptance:** `linac_rest.plot()` on a QUADS_ON run renders σ_x(z) and σ_y(z) (bounded, oscillating) with the updated title; a quads-OFF run renders the unchanged headline figure.

## Calibration & gates

- **Calibration runs on a quads-OFF, no-new-aperture deck** (Task 4): per-section crest-find + `rf_field_scale` `brentq` fit operate on `mean_energy` (longitudinal, transverse-independent **on-axis only**). With the quads physically absent (and the new quad/drift `radius` therefore `0.0`), the decimated `Np_calib=400` bunch is neither focused nor bore-scraped mid-fit, so no particle is removed from the `mean_energy` and the fit is clean and identical to today's headline. Quads (+ the new quad bore radius) appear only on the **fresh final-run deck** assembled with the calibrated `scales`+phases.
- **Hard gates unchanged and must still PASS** with QUADS_ON: gate 1 (per-section ΔE ±3%, `:387-389,429`), gate 2 (exit ⟨KE⟩ ±3% ≈308 MeV, `:392,431`), gate 5 (mean_z reached final zedge, `:401,434`), gate 6 (β>0.999, `:404,437`). Gate identity QUADS_ON↔OFF is a *consequence* of calibrating on the quads-absent deck (Task 4 acceptance), not an assumption.
- **Gate 4 (εn)** stays soft/diagnostic — now expect **real, possibly sizeable chromatic growth**, not ~0. Off-energy particles (the ±~20% injection spread) get different K1·(1/Bρ) focusing ⇒ chromatic emittance growth that can be **tens of %**, not "a few %." Treat a *runaway* (filamentation/over-focus) as the failure, not any growth.
- **Gate 7 (transmission)** stays soft/print-only — with QUADS_ON it becomes a real focused number but is NOT promoted to a hard gate (K1 is guessed; a hard target would invite reverse-fitting).
- **NEW soft gate `envelope_in_bore`** (Task 6) — `3σ < min bore`; the meaningful new check that the FODO actually contains the beam (vs. a few core particles passing). Print-only, RMS-level indicator (documented).

## Docs to update

- **`.claude/CLAUDE.md`** — long linac_rest paragraph, deviation **(3)**: scope "transmission ≈78.5%/capture ≈4.78%" explicitly to the **quads-OFF** headline; note the real bore is *already* the binding aperture for that headline. Add the QUADS_ON path (derived energy-scaled FODO, nominal μ=72°/cell, alternating sign, single-quad collapse of the H/V doublet, contained higher transmission). Change "(Meaningful predictive transmission needs `QUADS_ON` — exploratory, guessed K1.)" → "energy-scaled constant-phase-advance FODO derived from optics (nominal μ; RF-section-as-drift + single-quad-collapse approximations); magnitude un-validated (A→T missing)."
- **`linac_rest/README.md`** — *Space charge & quads* §: document `fodo_quad_gradients` (μ=72° nominal, field-free `L_h`, energy-scaling, alternating sign, single-quad collapse caveat), the calibrate-quads-OFF-then-build-fresh-quads-ON flow, the **new** quad/drift bore radii (gated on `quads_on`); clarify the solrf bore is already on for the headline; scope the "78.5%" sentence to quads-OFF and add the QUADS_ON contained number; *validation gates* list: add `envelope_in_bore` soft gate; *performance knobs*: `QUADS_ON`/`QUAD_K` now auto-derive K1 (QUAD_K override optional), document the Q8-inert off-by-one (only QUAD_K[0..5]/Q2–Q7 act, the 7th entry is unplaced).
- **`FIGURES.md`** — `fodo_optics.png` entry: update "exploratory/placeholder FODO envelope" → "derived energy-scaled FODO, contained oscillating σ_x/σ_y envelope (nominal μ; single-quad collapse; A→T still missing)."
- **`pipeline/run_pipeline.py`** — PERFORMANCE-KNOBS block (~89-91): add a **commented** `# linac_rest.config(QUADS_ON=True)  # exploratory FODO` preset (headline stays OFF).

## Verification

From repo root in the `CBB` env:

```bash
conda activate CBB

# Helper sanity (no run):
python -c "import build_linac_rest_lattice as L; print(L.fodo_quad_gradients())"
#   expect 7 values: 6 alternating-sign T/m in ~0.1–0.35 plus a trailing 0.0

# Baseline regression (headline, quads OFF) — MUST be unchanged:
python -c "import linac_rest; linac_rest.run(plots=False)"
#   exit ⟨KE⟩ ≈ 308 MeV (307.97 survivors / 309.2 full), transmission ≈ 78.5%, all hard gates PASS
#   (quad radius stays 0.0 ⇒ no new scrape plane ⇒ 78.5% bit-for-bit)

# Exploratory FODO (derived K1, default; or pass QUAD_K to override):
python -c "import linac_rest; linac_rest.config(QUADS_ON=True); linac_rest.run()"

# Figures:
python -c "import linac_rest; linac_rest.plot()"
```

Runtime ≈3 min (calibration ~140 s at `Np_calib=400` dominates; SC-off serial `ImpactTexe` track is fast).

**Numbers proving focusing works** (QUADS_ON vs. the 78.5% no-focusing baseline):

| Quantity | Source | Expected with working FODO |
|---|---|---|
| `transmission_core` (n_out/n_in, real bore) | `injection_summary.json` / gate 7 | **> 78.5%** (do NOT pin ~1.0 — a chromatic beam against the 9.9 mm exit bore may land ~85–92%; the deliverable is ">78.5% and bounded σ") |
| `σ_x(z)`, `σ_y(z)` | `stat_vs_z` → `fodo_optics.png` | **bounded / oscillating** (FODO beating); plane asymmetry expected (single-quad collapse) — NOT monotonic blow-up |
| exit ⟨KE⟩ (gate 2) | summary / `:392` | **still ≈308 MeV, PASS ±3%** |
| gates 1/5/6 | `validate_run` | **all PASS** |
| εn in→out (gate 4) | `:397-398,416` | **visible chromatic growth (tens of % possible)** — call a *runaway* the failure, not any growth |
| `envelope_in_bore` (new) | gate print | **PASS** (`3σ < min bore`; RMS-level indicator, soft) |

Most diagnostic single plot: `σ_x(z)` / `σ_y(z)` bounded along the full 36 m vs. the quads-OFF monotonic RMS rise.

## Risks & guards

| # | Risk | Symptom | Guard |
|---|---|---|---|
| R1 | Unstable FODO (realized μ/cell >180° → resonant blow-up) | transmission *drops below* 78.5%; σ grows faster than no-focus; εn explodes | μ=72° nominal is mid-band; verify σ_x/σ_y bounded; focusing that *lowers* transmission vs 78.5% is wrong. Cap `|K1|`≈1.5 if a quad over-pinches. |
| R2 | b1_gradient unit confusion (K1 [1/m²] vs T/m) | off by Bρ (~tens–hundreds): no effect or instant loss | V2 = T/m (`parsers.py:603-619`, `6_lattice_elements.md`). Always `Bg = K1·Bρ`; **per-quad** Bg scales with local Bρ (Task 1) — one Bg for all quads silently weakens focusing ~5× by the exit. |
| R3 | Sign error → net defocus | both planes diverge; transmission collapses | Alternate sign across the 6 quads (`[+,−,+,−,+,−]`); verify empirically σ_x *decreases* through the first negative-Bg quad. (Note: single-quad collapse means planes are asymmetric, not clean anti-phase doublet behavior.) |
| R4 | Aperture double-count / box-gaming | transmission too low (box+bore both bind) or faked ~1.0 (box widened) | Keep `XYRAD_M=0.02` (`build_linac_rest_lattice.py:236`, already 20 mm > 14.7 mm widest bore) so **bore is binding**; never widen the box (locked ruling, README:131-147). `section_bore_radii` is radius (halves diameter, `:138`). |
| R5 | Calibration breaks / shifts gates with quads on | gate 1 fails; `brentq` bracket error; mean shifts from scraped calib bunch | Calibrate on the quads-OFF, zero-quad-radius deck (Task 4); build the fresh quads-ON deck only after. Confirm gate 1/2 PASS identically. |
| R6 | Ntstep truncation | gate 5 FAIL (mean_z short) | `Ntstep=200000` sized for ~36 m; quads don't change path length. Gate 5 assert catches it. |
| R7 | QUAD_K length / Q8-inert off-by-one | user sets 7 values expecting all applied | Only QUAD_K[0..5] (Q2–Q7) act; the 7th (Q8, after last section) is never placed (`:354`). Helper returns length-7 with a trailing 0.0 to match the default shape (no IndexError). Document in README + helper docstring. |
| R8 | Default-path regression (the big one) | quads-OFF exit/transmission numbers change | The new quad/drift `radius` is gated on **`quads_on`**, not `bore_aperture_on` (Task 3) ⇒ headline quad radius stays `0.0`. All new FODO logic under `if quads_on:`; `quad_k[i] if quads_on else 0.0` at `:369` untouched. Re-run baseline, diff exit ⟨KE⟩ + transmission (must be bit-for-bit 78.5%). |
| R9 | Reverse-fitting K1 to a transmission target | K1 not traceable to optics | K1 must come from the phase-advance formula (Task 1), documented as placeholder; reverse-fitting is the same gaming class as widening the box — flag in review. |
| R10 | Chromatic εn / σ_KE growth misread as failure | gate 4 / σ_KE grows tens of % with quads ON | Expected (off-energy particles focus differently). Only a *runaway* is a failure; record growth, don't gate it. |
| R11 | Single-quad collapse vs H/V doublet | one plane over-defocuses; "anti-phase" acceptance fails | Modeled deliberately as a single thick quad (caveat on every output). Acceptance is **boundedness**, NOT clean doublet anti-phase. If one plane is lost, note it as a collapse artifact, not a bug. |

## Acceptance checklist

- [ ] `QUADS_ON` defaults `False` (`linac_rest_sim.py:86`); baseline reproduces exit ≈308 MeV (307.97/309.2) and transmission ≈78.5% **bit-for-bit** — **no regression** (R8); quad `radius` stays `0.0` quads-OFF.
- [ ] `fodo_quad_gradients()` returns **7** values (6 alternating-sign T/m ≈ 0.1–0.35 + trailing 0.0), derived from nominal μ=72° + **field-free `L_h`** + per-quad Bρ at measured `ke_in` (R2), documented as placeholder optics citing `details.md`.
- [ ] `QUAD_K` override → `build_impact(quad_k=…)` via `is None` logic (not `or`) → `b1_gradient` in T/m, applied only when `quads_on` (R8); `len(quad_k) >= N_SECTIONS-1` asserted.
- [ ] NEW quad (and inter-section drift) `radius` set to the real **exit** section bore when `quads_on` (Task 3) — gated on `quads_on`, NOT `bore_aperture_on`; solrf entrance-bore already on (unchanged).
- [ ] Calibration runs on the quads-OFF zero-quad-radius deck; the quads-ON run deck is built **fresh** with calibrated scales+phases (Task 4); `I.ele["quad3"]["b1_gradient"]` asserted non-zero pre-run; hard gates 1/2/5/6 **PASS** identically QUADS_ON vs OFF (R5).
- [ ] `XYRAD_M` NOT widened; transmission measured count-based vs the real bore **before** charge re-imposition (R4).
- [ ] QUADS_ON: transmission **> 78.5%** (not pinned ~1.0); σ_x(z) and σ_y(z) **bounded/oscillating** (R3, R11); εn growth recorded (tens of % OK, runaway = fail, R10); `envelope_in_bore` (`3σ<min bore`) PASS.
- [ ] New soft gate `envelope_in_bore` added (print-only, not in hard asserts); no "FAIL-when-OFF" liveness claim.
- [ ] `_stat_vs_z` records `sigma_y_m`; `fodo_optics.png` plots σ_x AND σ_y; summary records `quad_k` + `quad_phase_adv_deg`; figure + emittance titles updated for derived-FODO path with collapse/nominal-μ caveats.
- [ ] Q8-inert / QUAD_K length off-by-one documented; helper returns length-7 (R7).
- [ ] Docs synced: `.claude/CLAUDE.md` (deviation 3, bore-already-on), `linac_rest/README.md` (Space-charge-&-quads, gates, knobs, single-quad caveat), `FIGURES.md`, `pipeline/run_pipeline.py` commented preset.
- [ ] No "predictive transmission" claim; every QUADS_ON output labeled "placeholder optics — guessed K1, A→T undocumented, single-quad collapse, nominal μ"; K1 not reverse-fit (R9).

## Review resolutions

**Physics review**
- **B1 (bore already ON):** Adopted — Goal, Task 3, Task 4 rewritten; the bore is already binding for the headline, the *only* new aperture is the quad/drift element `radius`.
- **B2 (quad uses exit bore `[1]`, solrf uses entrance `[0]`):** Resolved decisively — exit radius `[1]` for the quad is correct (quad is downstream of the section exit taper); the entrance/exit split is an intentional taper, stated in Task 3.
- **B3 (`quad_k or` truthiness + length assert):** Adopted — Task 2 uses `is None`, adds `assert len(quad_k) >= N_SECTIONS-1`.
- **B4/B5 (calibration must be genuinely quads-OFF; gate identity is a consequence):** Adopted — Task 4 builds a separate quads-OFF calib deck and a fresh quads-ON run deck; gate identity framed as a consequence, not an assumption.
- **S1 (`L_h` excludes `L_q`):** Adopted — `L_h = 0.5·(L_section(i+1)+gap)`; table flagged as recomputed (~halved Bg), helper output is the acceptance.
- **S2 (H/V doublet):** Resolved — keep the single collapsed quad (no doublet split), documented loudly as an approximation; the "anti-phase" hard acceptance is dropped in favor of boundedness.
- **S3 (cite `parsers.py`):** Adopted — V2/V3/V4 mapping cited via `parsers.py:603-619`.
- **S4 (chromatic εn/σ_KE growth):** Adopted — gate 4 expects tens-of-% growth; runaway is the failure (R10).
- **S5 (RMS vs edge):** Adopted — `envelope_in_bore` uses `3σ < min bore`, documented as RMS-level.

**Integration review**
- **B1 (build/calibrate restructure):** Adopted the cleaner option — fresh `build_impact(quads_on=True, scales=…)` final deck (no in-place-mutate-then-configure reliance); calibrated phases carried over explicitly.
- **B2 (bore already on):** Same as physics B1 — adopted.
- **B3 (single-quad vs doublet + length misuse):** Resolved as physics S2 — single collapsed quad kept, loudly caveated, no anti-phase acceptance.
- **B4 (measured `ke_in`, drop "12×"):** Adopted — sim passes measured `ke_in` into the helper; contradictory "~12×" sentence removed; exit-KE consistency with `section_de_target` noted.
- **Regression (gate quad radius on `quads_on` not `bore_aperture_on`):** Adopted as the load-bearing fix (Task 3, R8) — prevents the silent 78.5% headline regression.
- **S1 (non-periodic RF-as-drift):** Adopted — "nominal μ, not realized" caveat; boundedness is the acceptance.
- **S2 (envelope liveness):** Adopted — `3σ` multiple; dropped the "FAIL-when-OFF" claim.
- **S3/S4 (add `sigma_y` to stat + plot):** Adopted — Task 6 adds `sigma_y_m`, Task 7 plots σ_y.
- **S5 (`is None` not `or`):** Adopted in Task 2.
- **S6 (bore vs box double-bind):** Confirmed — bore < box so bore binds; entrance/exit taper documented.
- **NITs:** N1 line citations corrected (helper after line 143; `details.md` A→T ≈188-194; `parsers.py:603-619`). N4 helper returns length-7 with trailing 0.0 (avoids list-length asymmetry). N2 sign verified empirically in smoke run (not hand-derivation alone). N3/N5 (file_id hard-edge, title strings) confirmed accurate, kept.
