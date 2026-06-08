# Implementation Plan: Cornell Linac Sections 2–8 via Impact-T (`linac_rest/` stage)

## 1. SCOPE & DECISIONS

| Decision | Choice | Rationale |
|---|---|---|
| **Sections** | 2, 3, 4, 5, 6, 7, 8 (seven electron TW sections) | The straight electron line to CHESS. |
| **e+ compressor (CU 2)** | **DEFER, out of scope** | Scope is the electron line, sections 2–8. The e+ comp section is omitted; its role/placement in the lattice is **not established in `details.md`**. Document as a one-line future-work note — do NOT assert a specific converter topology (no QWT/branch claim) the source does not support. |
| **Packaging** | **ONE combined Impact-T stage**: `linac_rest/` (single `ImpactT.in` deck, all 7 sections chained by increasing `zedge`) | Impact-T integrates one beam through one time-ordered lattice; per-section packages would force 7 exe round-trips and 7 openPMD conversions. One deck mirrors how Impact-T is meant to be driven and how LinacSim treats 2–8 as one BMAD generic linac. |
| **Field model** | **Generic constant-gradient TW, NO field maps** | Confirmed: no GPT/CST maps exist for sec 2–8; LinacSim/BMAD use the generic constant-gradient linac function (`details.md` lines 151–153). **Reuse the shipped `rfdata4–7` template field shape verbatim**; per section only rescale `rf_field_scale` / length / `zedge` (see Field-shape decision below). The rfdata is geometry-neutral (shape only); all per-section physics lives in the calibrated scale — do NOT encode R/τ/shunt impedance into the field profile (that would double-count, since the table gain already embeds R/τ). |
| **Field-shape recipe** | **Reuse shipped `rfdata4–7` + the template's 4-line `solrf` TW superposition (entrance + body_A + body_B·90° + exit), rescaled per section** | Resolves the draft's internal contradiction (it asserted BOTH a hand-synthesized flat-top AND the 4-line `1/sin(β₀d)` decomposition — these are inconsistent; the `1/sin` body factor and +30°/+90° inter-line phases are specific to the SLAC standing-wave decomposition encoded in `rfdata4–7`). Cloning the template verbatim and only rescaling is the lowest-risk path and de-risks the rfdata task. |
| **Phase** | **On-crest, `theta0_deg = 0`** per section (config knob, defaults 0); β>0.999 ⇒ fix θ=0, do NOT per-section re-autophase | Matches `linac_sec1` `PHASE_DEG=0` and BMAD generic. At β>0.999 the synchronous-phase reference is trivial; per-section β-dependent re-phasing is over-engineering and risks `autophase` non-convergence on a non-cathode beam. Use `autophase_and_scale` only for the **scale** calibration (Task 5). |
| **Space charge** | **OFF** (`Bcurr = 0`, no `-5` 3D line) | Negligible at >25 MeV (transverse SC ∝ 1/γ²; γ>49 at entry). Measured: SC-off cuts the template run 87.6 s → 5.4 s (16×). Coarse SC grid kept only as an optional sanity knob. |
| **Operating power** | **`POWER_MW = 11`** default (section-1 faithful klystron point), `√P` scaled per section from the @15 MW table; per-section config knob | Chain-consistency with `linac_sec1`. **One power convention for the whole linac.** Print the assumed power loudly. The @15 MW table column corresponds to a *different* (15 MW) section-1 run with a *different, higher* captured input — NOT a co-equal target for the same beam. |
| **Quads** | Include as `quadrupole` elements at **real tabulated lengths**, but **OFF (K1=0) for the headline beam**; the FODO line is a separate, clearly-labeled exploratory figure | A→T calibrations are genuinely undocumented (`details.md` line 178). A guessed-K1 emittance/transmission must NOT silently become "the linac exit emittance/transmission" in `_beam_summary`/`plot_chain`. Headline number = drift+RF only; quad-on is exploratory and every quad-on output is labeled "placeholder optics — not predictive." |
| **Inter-section drifts** | **Placeholder `DRIFT_M ≈ 0.4 m`** config knob, flagged as placeholder | Girder gaps not in `details.md`. |
| **Subprocess isolation** | **Run in-process** (do NOT reuse `Stage`'s subprocess launcher / `_launch_sim` / `run_step` / `afterstep`) | Impact-T is an external exe (`ImpactTexe`, serial) via lume-impact; no pywarpx global-geometry binding. BUT the in-process `ImpactStage` MUST still replicate `_prepare_environment()` (chdir to repo root + raise fd limit) and `setup_logging()`, and must capture/redirect `ImpactTexe` stdout into the pipeline log (lume-impact runs its own subprocess — verify it doesn't spew to the parent terminal mid-pipeline, wrap if it does). |

**On-crest scaling + calibration (per section, ONE power convention):**
```
ΔE_target(P_op) = ΔE_table × sqrt(P_op / 15)          # per-section, P_op = POWER_MW (default 11)
ΔKE_in          = measured ⟨KE⟩ from the read-in sec-1 exit dump (NOT hardcoded 25)
exit ⟨KE⟩       = ΔKE_in + Σ_i ΔE_target,i(P_op)       # validation gate, computed from actuals
```
`rf_field_scale` per section is calibrated (Task 5) via `autophase_and_scale` to hit `ΔE_target` — NOT computed analytically (solrf has no scalar-gradient input).

---

## 2. ARCHITECTURE

### New package `linac_rest/`

| File | Role |
|---|---|
| `linac_rest/__init__.py` | Facade: `config(**kwargs)`, `run(plots=True)`, `plot()`, `resolve_outdir()`, `DEFAULT_OUTDIR = "linac_rest/diags/main"`. Uses the new in-process `ImpactStage` (NOT `Stage`) which replicates cumulative-dict `config()` semantics and `_warn_unknown_params`. |
| `linac_rest/build_linac_rest_lattice.py` | Reuses shipped `rfdata4–7` field shapes (copies them into the run workdir, rescaled per section via `rf_field_scale`/length/`zedge`). Builds the `Impact` lattice (`add_ele` per the 4-element TW pattern + `drift` gaps + `quadrupole`). Module-level section-table constants (lengths, gains, bores, quad lengths) overridable by `config()`. Analogue of `build_*_field.py`. |
| `linac_rest/linac_rest_sim.py` | `main()`: read upstream beam → `ParticleGroup`; build deck; calibrate per-section scale (Task 5); size `Ntstep` from lattice length; `I.run()`; **assert final `mean_z` reached the last `zedge`**; convert `ParticleGroup` outputs → WarpX-style openPMD into `diags/main/particles/`; write `injection_summary.json`. Module-level **performance knobs**: `Np`, `Nx/Ny/Nz`, `Ntstep`, `Dt`, `POWER_MW`, `PHASE_DEG`, `DRIFT_M`, `QUAD_K`, `QUADS_ON`, `OUTDIR`, `N_DIAGS`. |
| `linac_rest/plot_linac_rest.py` | `main()`: read `diags/main`, write PNGs to `linac_rest/results/`. Reuse `pipeline.beam_metrics.rms_emit`. Use Impact-T `I.stat(...)` arrays (energy/emit/σ vs z) for per-section evolution panels. |
| `linac_rest/README.md` | Physics, generic-CG-TW rationale, template-reuse decision, on-crest/β>0.999 assumption (with the min-captured-KE check), SC-off, per-section table, quad-OFF-headline + FODO-exploratory caveat, drift placeholder, e+ defer wording, Impact-T deck description, perf knobs, gotchas, figure list (must match FIGURES.md + plot script). |
| `pipeline/impact_io.py` (new shared) | Adapter: `ParticleGroup → WarpX openPMD layout` — species literally **`"electrons"`** (plural), records `x,y,z,ux,uy,uz,w` with `ux=γβ` (dimensionless momenta) and `w` = macro-weight **count** (NOT charge). Flag the asymmetry explicitly: the `ParticleGroup` uses `species="electron"` (singular) but the openPMD output MUST be `"electrons"` (plural) for `plot_chain`/`_beam_summary`. Also the reverse reader helper for the handoff IN. Do NOT use `ParticleGroup.write()` (verified: emits openPMD 2.0 with a STRING `openPMDextension` attr that `OpenPMDTimeSeries` rejects). |

### New in-process runner: `ImpactStage` (in `pipeline/_runner.py` or new `pipeline/_impact_runner.py`)
Mirrors `Stage`'s public surface (`config`, `run`, `plot`, `_params`, `_warn_unknown_params`) but `run()` calls `build.main()`, `sim.main()`, `plot.main()` **in-process** (no `subprocess.run`, no `_launch_sim`, no `run_step`/tqdm `afterstep`). It MUST:
- call `_prepare_environment()` (chdir repo root + raise RLIMIT_NOFILE — needed: `impact_io.py` loops openPMD dumps and `plot_chain` already hit the 256-fd wall) and `setup_logging()` (reuse, don't reimplement);
- capture/redirect the `ImpactTexe` subprocess stdout into the pipeline log;
- `_warn_unknown_params` must **AST-introspect the sim module** AND **live-check build+plot modules** (mirror `Stage` exactly), so a `config()` key targeting a build-module section-table constant doesn't spuriously warn.

### Handoff IN (read linac_sec1 captured beam → ParticleGroup)
- Select `linac_sec1`'s **last (exit) dump** — the captured ~25 MeV coasting beam (mirror `_exit_row`'s `rows[-1]` for the linac, NOT a handoff-plane nearest-⟨z⟩ pick).
- Read via `ParticleGroup`; `P_in.drift_to_t()`; zero z (`P_in.z -= P_in["mean_z"]`); Impact wants `z==0`, t-coordinates, `species=="electron"`.
- Carry actual surviving charge/weights via `I.initial_particles = P_in` (canonical lume-impact pattern; carries charge — confirm it propagates to the header, don't rely on a bare `total_charge` assignment). Do NOT renormalize.
- Measure ⟨KE⟩_in and **min captured KE** from this dump; assert min-KE β>0.999 (justifies the rigid-crest no-slip assumption across 30 m). Record the true-injected denominator from `linac_sec1` so end-to-end capture stays ~7%.

### Handoff OUT (Impact-T → openPMD for plot_chain + summary)
- `P_out = I.particles["final_particles"]` (a `ParticleGroup`).
- `pipeline/impact_io.py` writes the explicit WarpX layout (species `"electrons"`, records `x,y,z,ux,uy,uz,w`, `ux=γβ`) — NOT `ParticleGroup.write()`.
- Optionally dump several iterations along z (from Impact-T slice stats) for the `chain_evolution` vs-z panels; a single exit dump is tolerated (plot_chain skips <50-macroparticle dumps, sorts by ⟨z⟩).
- Write `injection_summary.json` recording: `q_injected_C` = charge read in from linac_sec1's exit (so `_beam_summary` reports ~100% within-stage and the chain capture narrative stays coherent); `z_inject_lab_m` = the lab-z the beam was injected at; note Impact-T output z is local-frame.

### Orchestration wiring
- **`pipeline/run_pipeline.py`**: `import linac_rest`; add `linac_rest.config(POWER_MW=11.0)` to PHYSICS block; add a perf-knob line (`linac_rest.config(Np=..., Ntstep=...)`) to PERFORMANCE KNOBS; add `linac_rest.run()` after `linac_sec1.run()`; add `_beam_summary(linac_rest.resolve_outdir(), "linac exit (8 sections)", "MeV")`; update banner `... → linac_sec1 → linac_rest` and the "Figures:" line.
- **`pipeline/plot_chain.py`**: add `STAGES` entry `{"name": "linac_rest", "path": "linac_rest/diags/main/particles", "z0": <computed>, "geom": "rz", "color": "C5"}` (C5/C6 — NOT C4; C4 is the "passes iris" waterfall bar). Add `_apply_linac_rest_z0()` analogous to `_apply_linac_z0()` that computes the lab offset from **recorded values** (`linac_sec1_exit_lab_z − linac_rest_local_inject_z`), NOT a literal 5.1 m. Confirm `_exit_row` needs no special case (falls through to `rows[-1]`, correct for the linac exit). `render_*` pick it up once the openPMD layout matches.
- **`pipeline/__init__.py`**: confirm no edit needed (adding a STAGES segment doesn't require it) — check only.

---

## 3. IMPLEMENTATION STEPS (ordered, each a reviewable unit)

**Task 1 — Section-table + scaling module.**
Files: `linac_rest/build_linac_rest_lattice.py` (constants only).
Encode the per-section table (type, L, ΔE@15MW, bore taper, quad length) as module-level constants + the `√P` scaling helper `section_gradient(power_mw)`.
Acceptance: importing the module and calling the helper at 11 MW reproduces the table below within rounding.
Validation: assert `G_section(15) ≈ ΔE_table/L` for all 7 sections; assert `√P` scaling gives the @11 MW column.

**Task 2 — rfdata field-shape reuse.**
Files: `linac_rest/build_linac_rest_lattice.py`.
**Reuse the shipped `rfdata4–7`** from `/Users/rylandgoldman/Downloads/lume-impact-master/docs/examples/templates/traveling_wave_cavity/` as the field shape (copy into the run workdir; do NOT hand-synthesize a flat-top). Provide a helper that places the shape per section with rescaled length/`zedge`.
Acceptance: `rfdata` files load in a minimal single-section `Impact` run without error; integrated on-axis voltage ≈ G·L (±2%) after scale calibration.
Validation: a single-section smoke run reaches ~the template's energy; numerically integrate reconstructed Ez via `process_fieldmap_solrf_fourier`.

**Task 3 — Build the chained lattice (`ImpactT.in`).**
Files: `linac_rest/build_linac_rest_lattice.py`.
Per section: 4 `solrf` sub-elements (entrance, body_A, body_B with the template's +90° `theta0`, exit; body scale = entrance·1/sin(β₀d), sin(β₀d)≈0.8657 — reuse the template's recipe verbatim), placed at increasing `zedge`; **`drift`** (NOT `DriftTube`) gaps (`DRIFT_M`); `quadrupole` elements at real lengths with `QUAD_K` placeholders (default K1=0, `QUADS_ON=False`). Header: `Npcol=Nprow=1`, `Bcurr=0`, `Nemission<0`, `Flagimg=0`, `Perdlen` > total length (~30 m), `Ntstep` sized from total lattice length / (c·Dt) with margin (~50k for ~30 m at Dt≈2e-12), power-of-2 SC mesh (unused with SC off).
Acceptance: `I = Impact(...)` loads the generated deck; `len([e for e in I.lattice if e['type']=='solrf'])==28` (7×4); elements sorted by `zedge`; all element `type` strings parse (`drift`, `quadrupole`, `solrf`).
Validation: total lattice length ≈ Σ(L+drift) within the section table.

**Task 4 — Handoff IN reader + `pipeline/impact_io.py`.**
Files: `pipeline/impact_io.py`, `linac_rest/linac_rest_sim.py`.
Read `linac_sec1/diags/main/particles` last dump → `ParticleGroup`; `drift_to_t()`; zero z. Write the reverse adapter (`ParticleGroup → WarpX openPMD`, species `"electrons"`, `ux=γβ`, `w`=count). Set `I.initial_particles = P_in`.
Acceptance: `ParticleGroup` loads `len>0`, `species=="electron"`, `mean_energy≈25 MeV`; min-KE β>0.999 asserted/printed; the openPMD writer output is readable by **both** `OpenPMDTimeSeries(...).get_particle(["x","y","z","ux","uy","uz","w"], species="electrons")` AND a dry-run of `_beam_summary`/`build_moment_table` over the directory (no error); `final_particles.charge == P_in.charge`.
Validation: round-trip a known PG and assert `⟨KE⟩`, charge, count preserved within 0.1%.

**Task 5 — Per-section scale calibration via `autophase_and_scale`.**
Files: `linac_rest/linac_rest_sim.py`.
For each section in lab order: use `impact.autophase.autophase_and_scale(...)` scaling a `ControlGroup` (via `add_group`) over the 4 sub-elements to `metric=mean_energy` target `ΔE_target` with `isolate=True` (NOT `I.autophase()` no-arg — verified it throws on a non-cathode beam; NOT a hand-rolled linear ratio loop — gain is not exactly linear in scale once phase couples). At β>0.999 fix θ=0 (crest); calibrate scale only. Preserve the 1/sin body ratio across the group.
Acceptance: after calibration, each section's measured ΔKE matches its target within ±3%.
Validation: print a per-section `target vs achieved` table.

**Task 6 — Full sim run + output.** *(depends on Task 7 facade existing — see ordering note)*
Files: `linac_rest/linac_rest_sim.py`.
Wire: build deck → `initial_particles`/`Np` → calibrate (Task 5) → `I.run()` → assert `I.finished and not I.error` **AND `I.stat("mean_z")[-1] ≈ Σ L` (beam reached final zedge — Ntstep not truncated)** → convert `final_particles` + along-z slices to openPMD in `diags/main/particles/` → write `injection_summary.json` (`q_injected_C`, `z_inject_lab_m`).
Acceptance: `linac_rest_sim.main()` completes; `diags/main/particles/` contains readable dumps; exit ⟨KE⟩ ≈ `⟨KE⟩_in + Σ ΔE_target` (≈307 MeV @11 MW from a ~25 MeV input); `_beam_summary` reads the dir reporting ~100% within-stage capture.
Validation: `mean_z` reached; per-section target-vs-achieved table within ±3%.

**Task 7 — Facade + in-process `ImpactStage`.** *(land before/with Task 6)*
Files: `linac_rest/__init__.py`, `pipeline/_runner.py` (or `pipeline/_impact_runner.py`).
Implement `ImpactStage` (in-process build/sim/plot; cumulative `config`; `_warn_unknown_params` AST-on-sim + live-on-build/plot; `_prepare_environment()` + `setup_logging()` + `ImpactTexe` stdout redirect); wire the facade + `resolve_outdir()`.
Acceptance: `config()` overrides reach module constants in BOTH build and sim modules; unknown keys warn, valid build-module keys do NOT warn; `resolve_outdir()` returns the OUTDIR override or default; `import linac_rest; linac_rest.run(plots=False)` runs end-to-end.
Validation: `linac_rest.config(POWER_MW=15); resolve_outdir()` correct; a typo key logs a warning; a valid section-table key does not.

**Task 8 — Plots.**
Files: `linac_rest/plot_linac_rest.py`.
Energy-vs-z (per-section gains annotated), σ_KE/⟨KE⟩ vs z, ε_n,x/ε_n,y vs z, transmission, and a **separate clearly-labeled FODO exploratory** β/σ_x vs z (only meaningful with `QUADS_ON`). Write PNGs to `linac_rest/results/`.
Acceptance: `linac_rest.plot()` writes the exact PNG filenames listed in README/FIGURES.md.
Validation: PNGs non-empty; energy plot endpoint matches Task 6.

**Task 9 — Orchestrator + cross-stage wiring.**
Files: `pipeline/run_pipeline.py`, `pipeline/plot_chain.py` (+ check `pipeline/__init__.py`).
Add import, config lines, `linac_rest.run()`, `_beam_summary`, banner/Figures strings; add `STAGES` entry (color C5/C6) + `_apply_linac_rest_z0()` computing offset from recorded `z_inject_lab_m` and linac_sec1 exit lab-z.
Acceptance: `python pipeline/run_pipeline.py` runs the full 5-stage chain end-to-end; `plot_chain` renders the new segment in lab z with NO overlap with `linac_sec1`.
Validation: chain figures show monotonic energy rise to ~307 MeV @11 MW; end-to-end capture stays ~7%.

**Task 10 — Doc-sync (see §6).**

**Task 11 — Memory entry** for the Impact-T conventions (no subprocess isolation but still `_prepare_environment`/fd-raise; generic CG TW via reused `rfdata4–7`; `autophase_and_scale` not `autophase()`; Ntstep-truncation-reports-success gotcha → assert `mean_z`; `ParticleGroup.write()` STRING-attr incompatibility → custom writer; species `"electron"`/`"electrons"` asymmetry).

---

## 4. PHYSICS PARAMETERS (per section)

Frequency = **2856 MHz**, phase = **0° (on-crest)**, SC **off**, quads **OFF for headline** for all. `sin(β₀d)≈0.8657` (S-band, d=3.5 cm).

| Sec | Type | L (m) | ΔE@15MW | G@15MW (MV/m) | ΔE@11MW | G@11MW (MV/m) | Bore radius (m) | Drift after |
|----|------|------|--------|------|--------|------|------|------|
| 2 | CEA 2 | 2.94 | 33 | 11.2 | 28.3 | 9.6 | 0.0126→0.0099 | DRIFT_M |
| 3 | CEA 3 | 2.94 | 33 | 11.2 | 28.3 | 9.6 | 0.0126→0.0099 | DRIFT_M |
| 4 | CU 5 | 4.97 | 51 | 10.3 | 43.7 | 8.8 | 0.0147→0.0117 | DRIFT_M |
| 5 | CEA 4 | 5.15 | 55 | 10.7 | 47.1 | 9.1 | 0.0147→0.0117 | DRIFT_M |
| 6 | CEA 5 | 5.15 | 55 | 10.7 | 47.1 | 9.1 | 0.0147→0.0117 | DRIFT_M |
| 7 | CU 3 | 4.97 | 51 | 10.3 | 43.7 | 8.8 | 0.0147→0.0117 | DRIFT_M |
| 8 | CU 4 | 4.97 | 51 | 10.3 | 43.7 | 8.8 | 0.0147→0.0117 | — |

- **`rf_field_scale`** per section calibrated to G (Task 5 via `autophase_and_scale`), not analytic.
- **rfdata** is shape-only (reused `rfdata4–7`); R/τ/shunt impedance are NOT encoded (already embedded in the table gain).
- **Drifts:** `DRIFT_M ≈ 0.4 m` placeholder between sections (config knob, flagged).
- **Quads:** `quadrupole` at real lengths (Q1=10.5″, Q2=11″, Q3=18″, Q4=25″, Q5=16.2″, Q6/Q7=22″, Q8=20.9″); **K1 defaults to 0 (`QUADS_ON=False`)** for the headline. `QUAD_K` placeholder list is exploratory only; A→T unknown, flagged loudly. With K1=0 the headline emittance/transmission are NOT distorted by guessed optics.
- **Energy budget (ONE power convention):** exit ⟨KE⟩ = measured `⟨KE⟩_in` (from sec-1 exit dump, ~25 MeV @11 MW) + Σ ΔE_target,i(P_op). At 11 MW: ≈ 25 + 282 ≈ **307 MeV**. The @15 MW column is a *different* beam (15 MW sec-1, higher captured input) and is NOT a co-equal target — do not present "354/300+ MeV bracket" as a faithfulness check. Faithfulness gate = measured-in + Σ-table-scaled, computed from actuals.
- **σ_KE:** input σ_KE≈8 MeV; min captured KE ~10–15 MeV ⇒ β>0.999 ⇒ rigid-crest no-slip holds (asserted in Task 4).

---

## 5. VALIDATION

Each run prints (per section + cumulative), from `I.stat(...)`/`fort.18`/`ParticleGroup`:

1. **Per-section energy gain** (⟨KE⟩_out − ⟨KE⟩_in) vs `ΔE_table×√(P/15)` — within **±3%** (primary; Task 5 enforces).
2. **Cumulative ⟨KE⟩ and γ** after each section — exit ≈ measured `⟨KE⟩_in` + Σ ΔE_target (≈307 MeV @11 MW), computed from actuals not hardcoded.
3. **σ_KE/⟨KE⟩** — absolute σ_KE conserved to first order (second-order crest curvature adds a small correlated term, negligible vs the 8 MeV input); relative spread adiabatically shrinks as ⟨KE⟩ grows (8 MeV/25 MeV ≈ 32% → ≈2.6% at 307 MeV).
4. **Normalized emittance ε_n,x/ε_n,y** in vs out — with quads OFF, near-conserved (RF + drift only); flags any numerical growth. Quad-ON ε is exploratory only, never the headline.
5. **Beam reached final zedge** — `I.stat("mean_z")[-1] ≈ Σ L` (catches Ntstep truncation, which reports `finished=True` falsely).
6. **min captured KE β>0.999** — justifies no-slip (Task 4).
7. **Transmission** — with no aperture scraping, ~100% is a tautology not a result; either add tapered bore as `radius` scraping (then meaningful) or state transmission is assumed and drop it from the gates. End-to-end stays ~7% of true-injected (all loss is sec-1 capture).
8. **(exploratory)** β-function / σ_x through the FODO line — sanity only, labeled placeholder optics.

**Sanity vs `details.md`:** Σ sections 2–8 = 329 MeV @15 MW (282 MeV @11 MW); printed per-section achieved-vs-table.

---

## 6. DOC-SYNC (part of "done")

- **`.claude/CLAUDE.md`**: add inter-stage-contract row (`linac_rest/linac_rest_sim.py` | reads `linac_sec1/diags/main/particles` (last/exit dump) | writes `linac_rest/diags/main`); update chain diagram (`... → linac_sec1 → linac_rest`); add architecture prose for the Impact-T-not-pywarpx → in-process (no `_launch_sim`) but still `_prepare_environment`/fd-raise exception, generic-CG-TW via reused `rfdata4–7` (no field maps), `autophase_and_scale` calibration, Ntstep-assert-`mean_z` gotcha, and the `ParticleGroup`-not-`OpenPMDTimeSeries`/`ParticleGroup.write()` round-trip; add `linac_rest` to Commands/Stage-API + PERFORMANCE KNOBS (serial-only `ImpactTexe`; SC-off makes it cheap, ~tens of s).
- **`README.md`**: component-table row; update both ASCII chain diagrams.
- **`linac_sec1/README.md`** and **`linac_sec1/__init__.py:3` docstring** ("later sections → `linac_sec2`, …"): point at the IMPACT-T `linac_rest` stage.
- **`FIGURES.md`**: add `## 5. Linac sections 2–8 (IMPACT-T)` listing every PNG `plot_linac_rest.py` writes (names must match); update chain diagram + any waterfall/scorecard column.
- **`requirements.txt`**: pin `lume-impact==0.11.0`, `openpmd-beamphysics` (+ `distgen` if used); note `impact-t`/`ImpactTexe` installs via conda-forge (comment, not a pip line) like pywarpx/openpmd-api.
- **`.gitignore`**: `fort.*`, `*.h5`, `diags/` already covered. **Add `rfdata*` and the generated `ImpactT.in`** (only `ImpactT.in.bak` is currently ignored). Decide+state the workdir: if `Impact(workdir=...)` is in-tree (so diags persist for `plot()`), these matter; if `use_temp_dir=True` (default) they're belt-and-suspenders. Resolve before implementation.
- **Memory**: new entry (Task 11).
- **Commit**: `linac_rest/*.py` + `README.md`, edited `pipeline/*.py`, docs, `requirements.txt`; **plain `git add linac_rest/results/*.png`** — `results/` is NOT git-ignored in this repo (the docs' "git-ignored, use `git add -f`" claim is a pre-existing doc bug; do not propagate it, do not fix repo-wide here). Do NOT commit `diags/`, `ImpactT.in`, `rfdata*`, `.h5`, `fort.*`, logs.

---

## 7. RISKS & OPEN QUESTIONS

1. **rfdata fidelity** — reusing `rfdata4–7` + the template's 4-line `1/sin`/+90° recipe and only rescaling avoids the synthesize-vs-decompose inconsistency; verify integrated voltage post-calibration (Task 2). *Mitigation: clone verbatim, rescale only.*
2. **`ParticleGroup.write()` ↔ `OpenPMDTimeSeries`** (verified `openPMDextension STRING`) — custom `impact_io.py` writer, species `"electrons"`, `ux=γβ`, `w`=count. Highest-risk integration; Task 4 acceptance covers via dual-reader assert.
3. **Ntstep truncation reports `finished=True`** (verified: stopped at z=2.33 m, returned success) — size from lattice length and **assert `mean_z` reached** (Task 6 gate).
4. **`autophase()` no-arg throws** on a non-cathode beam — use `autophase_and_scale` + `ControlGroup` (Task 5).
5. **Operating power** — 11 MW, one convention; print loudly. @15 MW column is a different beam, not a bracket.
6. **Quad A→T calibrations UNKNOWN** — quads OFF for headline; FODO is exploratory, every quad-on output labeled placeholder.
7. **Inter-section drift lengths undocumented** — `DRIFT_M` placeholder, flagged.
8. **lab-z offset** — `_apply_linac_rest_z0()` from recorded `z_inject_lab_m`, not literal 5.1 m, else overlap.
9. **Sec 5/6 (155*) and bore (†)** are source-table guesses — propagate uncertainty in docs.
10. **e+ comp** — deferred; its lattice role is not established in `details.md` (no topology claim).

**Open questions to resolve before/during implementation:** (a) confirmed 11 MW (one convention). (b) aperture scraping (tapered bore as `radius`) wanted, or transmission assumed? — decide in Task 8/validation. (c) `autophase_and_scale` is the entry point (no-arg form does not converge). (d) lume-impact workdir location — decide in-tree vs temp before implementation (drives `.gitignore`).

Key reference paths: TW template `/Users/rylandgoldman/Downloads/lume-impact-master/docs/examples/templates/traveling_wave_cavity/` (`rfdata4–7`); `ImpactTexe` `~/miniforge3/envs/CBB/bin/ImpactTexe` (serial); handoff pattern `linac_sec1/linac_sec1_sim.py` (`load_injector_bunch`, `injection_summary.json`, `_exit_row` rows[-1]); orchestrator `pipeline/run_pipeline.py` (`_beam_summary`); cross-stage `pipeline/plot_chain.py` (`STAGES`, `_apply_linac_z0`, `build_moment_table`, `get_particle species="electrons"`); section table `reference/Linac Simulation Documentation/details.md:159-167`; Impact-T TW model `reference/Impact-T Documentation/physical_models/8.7_traveling_wave_structures.md`; lume-impact API `impact/autophase.py` (`autophase_and_scale`), `impact/fieldmaps.py`, `impact/parsers.py`.

---

## Changelog vs draft (review fixes applied)

**review:physics**
- HIGH#1 (energy-budget power-point conflation): adopted ONE power convention (11 MW); exit = measured `⟨KE⟩_in` + Σ ΔE_target computed from actuals, not hardcoded 25; dropped the "brackets ~300+ MeV / 354 co-equal" framing; @15 MW column reframed as a different beam.
- HIGH#2 (low-energy tail no-slip): added an explicit min-captured-KE β>0.999 check (Task 4) as the actual no-slip justification.
- MED: rfdata is shape-only, R/τ NOT encoded (avoid double-count); e+ comp deferral reworded to "role not established in details.md" (no QWT/branch claim); quads OFF for headline so guessed-K1 ε/transmission can't leak into the summary; `DriftTube`→`drift`; Ntstep sized + `mean_z` assert.
- LOW: at β>0.999 fix θ=0, skip per-section autophase; softened σ_KE conservation to first-order (+second-order curvature note); transmission tautology flagged (add bore `radius` or drop from gates); confirmed SC-off/β≈1/generic-CG/calibrate-scale/one-deck/in-process as sound.

**review:impact-tooling**
- HIGH (Ntstep truncation reports success): assert `I.stat("mean_z")[-1] ≈ Σ L` as a hard gate (Task 6).
- MED (`autophase()` throws): switched Task 5 to `autophase_and_scale` + `ControlGroup`/`add_group`.
- MED (synthesized-flat-top vs 4-line `1/sin` inconsistency): resolved by reusing shipped `rfdata4–7` verbatim and only rescaling (Task 2 rewritten); recommendation (a) adopted.
- LOW: runtime de-budgeted (SC-off ≈ 5 s, not 90 s — dropped "large speedup" as load-bearing); `drift` not `DriftTube`; `I.initial_particles = P_in` canonical pattern, confirm charge propagation; species `"electron"`/`"electrons"` asymmetry flagged in `impact_io.py`; workdir/`use_temp_dir` resolution noted.

**review:integration**
- HIGH (openPMD writer contract): writer MUST emit species literally `"electrons"`, records `x,y,z,ux,uy,uz,w`, `ux=γβ`, `w`=count; Task 4 acceptance now asserts both `_beam_summary` and `build_moment_table` read it.
- HIGH (`_apply_linac_rest_z0` under-specified): computes offset from recorded `z_inject_lab_m` + linac_sec1 exit lab-z, not literal 5.1 m; sim must record those z fields.
- HIGH (in-process drops fd-raise/log/stdout contract): `ImpactStage` must reuse `_prepare_environment()` + `setup_logging()` and redirect `ImpactTexe` stdout (Task 7).
- HIGH (workdir/`.gitignore`): decide+state workdir; add `rfdata*` and generated `ImpactT.in` to `.gitignore`.
- MED: `injection_summary.json` records `q_injected_C` = sec-1 captured charge (Task 6 acceptance); STAGES color C5/C6 not C4 (C4 collides with waterfall bar); `_exit_row` confirmed needs no special-case; `_warn_unknown_params` must AST-introspect sim + live-check build/plot.
- MED/LOW: dropped `git add -f` (results/ is NOT git-ignored — pre-existing doc bug, flagged not propagated); `impact-t` binary as a comment not a pip line; **task ordering fixed — Task 7 facade lands before/with Task 6** (Task 6 acceptance otherwise couldn't pass); `pipeline/__init__.py` confirm-only.
