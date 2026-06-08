# `linac_rest` — Cornell Linac Sections 2–8 (Impact-T)

The rest of the straight electron line to CHESS, after `linac_sec1`: seven S-band
(2856 MHz) traveling-wave (TW) accelerating sections (CEA 2/3/4/5 + CU 3/4/5),
chained into **one** Impact-T deck and integrated as one time-ordered beam. Reads
`linac_sec1`'s captured ~25 MeV exit beam and accelerates the captured core on-crest
to **≈308 MeV** at the default 11 MW klystron point (307.97 MeV survivor mean through the real
bore; 309.2 MeV full-beam mean before the bore scrapes the lower-energy off-axis tail; ≈307 MeV
expected from the table).

```
cathode → gun → injector → linac_sec1 → linac_rest (this, sections 2–8)
   (WarpX 2D)   (WarpX RZ ×3)            (Impact-T, 7 TW sections, ~36 m)
```

## Why Impact-T, not WarpX

The four upstream stages are WarpX/pywarpx runs. This stage is an **external serial
Impact-T run** (`ImpactTexe`) driven through **lume-impact**. Impact-T integrates one
beam through one time-ordered lattice — the natural fit for 7 sections that BMAD/LinacSim
also treat as one generic linac. Because there is no pywarpx global-geometry binding, the
stage runs **in-process** (`pipeline._impact_runner.ImpactStage`), not in a per-stage
subprocess like the WarpX stages — but it still reuses the pipeline's repo-root chdir +
fd-limit raise + shared log, and redirects `ImpactTexe` stdout into the pipeline log.

**Progress bars.** The WarpX stages drive a tqdm bar from a pywarpx `afterstep` callback;
Impact-T is an opaque external exe with no such hook, so this stage shows two
`pipeline._impact_runner.terminal_progress` bars instead: a **calibrate** bar (one tick per
TW section, `calibration.calibrate_sections(bar=…)`) and a **track** bar for the final
`I.run()` driven by a background thread watching the run's `fort.18` (column 2 = reference
`mean_z`) advance 0 → lattice length. Both write to a saved duplicate of the terminal fd, so
they survive the sim-phase stdout redirect (the same trick `_runner.run_step` uses).

## Field model — generic constant-gradient TW, no field maps

Sections 2–8 have **no GPT/CST field maps** (none exist; LinacSim/BMAD model them with the
generic constant-gradient linac function). We reuse the shipped lume-impact
`traveling_wave_cavity` template field **shape** — `rfdata4–7`, vendored into
`linac_rest/rfdata/` (committed) — verbatim, and put **all** per-section physics in the
calibrated field scale. Each section is the template's 4-line `solrf` superposition of two
standing-wave maps (G. A. Loew et al., SLAC-PUB-2295):

| line     | rfdata | length            | θ₀ offset | scale         |
|----------|--------|-------------------|-----------|---------------|
| entrance | 4      | short coupler cell| base + 0° | S             |
| body_1   | 5      | bulk (L−couplers) | base + 30°| S / sin(β₀d)  |
| body_2   | 6      | bulk (L−couplers) | base + 90°| S / sin(β₀d)  |
| exit     | 7      | short coupler cell| base + 0° | S             |

`sin(β₀d) ≈ 0.8657` (S-band, d = 3.5 cm). The rfdata carries **no** R/τ/shunt impedance —
those are already embedded in the per-section ΔE table, so encoding them in the field
profile would double-count.

**Key field-reuse fact:** the rfdata Fourier reconstruction uses the fundamental period stored
*inside* the file (~0.105 m, the 3-cell block) as its wavelength — NOT the lattice element `L`.
(The body maps `rfdata5`/`rfdata6` store ~2.94 m of fitted Fourier data, but the field is
*periodic* at that ~0.105 m wavelength, so it tiles over any element length.)
The element `L` only sets the active z-range `[zedge, zedge+L]` that Impact-T integrates the
periodic field over, so a longer section is simply *more cells of the same per-cell field*.
That is why "rescale length" just means setting the body element `L` per section; the field
shape is reused unchanged.

## Operating point & energy budget (one power convention)

Default `POWER_MW = 11` — the section-1 faithful klystron point — for the **whole** linac,
`√P`-scaled per section from the @15 MW `details.md` table:

```
ΔE_target,i(P_op) = ΔE_table,i × sqrt(P_op / 15)
exit ⟨KE⟩         = measured ⟨KE⟩_in + Σ_i ΔE_target,i(P_op)   # computed from actuals
```

| Sec | Type | L (m) | ΔE@15MW | ΔE@11MW | bore r (mm) | quad after |
|-----|------|-------|---------|---------|-------------|------------|
| 2 | CEA 2 | 2.94 | 33 | 28.3 | 12.6→9.9  | Q2 |
| 3 | CEA 3 | 2.94 | 33 | 28.3 | 12.6→9.9  | Q3 |
| 4 | CU 5  | 4.97 | 51 | 43.7 | 14.7→11.7 | Q4 |
| 5 | CEA 4 | 5.15 | 55 | 47.1 | 14.7→11.7 | Q5 |
| 6 | CEA 5 | 5.15 | 55 | 47.1 | 14.7→11.7 | Q6 |
| 7 | CU 3  | 4.97 | 51 | 43.7 | 14.7→11.7 | Q7 |
| 8 | CU 4  | 4.97 | 51 | 43.7 | 14.7→11.7 | — |

Σ ΔE @11 MW ≈ 282 MeV ⇒ exit ≈ 307 MeV expected from a ~25 MeV input; the **achieved**
calibrated exit is **309.2 MeV full-beam / 307.97 MeV for the survivors through the real bore**
(every per-section ΔE within <0.05 % of its target; the ~1.2 MeV survivor-vs-full-beam gap is the
real bore scraping the lower-energy, more-divergent off-axis particles — see *Space charge &
quads*). The @15 MW column is a *different* (15 MW) section-1 beam — NOT a co-equal target. Sec 5/6
cell counts and the bore taper are `details.md` source-table guesses (propagated as uncertainty).

The per-section `rf_field_scale` is **calibrated**, not analytic (the `solrf` element has no
scalar-gradient input). Each section is set to its **local on-crest base phase** first —
Impact-T `theta0_deg` is *absolute*, so the chained-deck phase walk (drifts + finite β shift
the bunch arrival phase hundreds of degrees per section) means θ₀ = 0 is on-crest **only for
section 2**; downstream sections find their crest via a coarse phase scan + parabolic refine
(`calibration._find_crest_phase`), then the scale is fit to ΔE_target (`brentq`). Calibration
runs on copies via lume-impact `track`, not full runs.

## Captured-core cut (handoff IN)

The `linac_sec1` exit dump carries a low-energy **un-captured tail** (min KE ~0.2 MeV; ~16%
of charge below 10 MeV) whose β < 0.999 would break the rigid-crest no-slip assumption across
the 36 m line. The handoff keeps only the **captured core** (KE ≥ `MIN_KE_MEV`, default
12 MeV ⇒ β > 0.99917, ~88% of the exit charge) — a **model-validity cut**, since the rigid-
crest TW model is only valid for the relativistic core.

**Honest capture denominator:** `injection_summary.json` records `q_injected_C` = the **full**
`linac_sec1` captured charge that arrived at the handoff (NOT the post-cut core). So
`_beam_summary`'s within-stage capture (`q_out / q_injected`) counts both the dropped tail
(the model-validity cut) and any in-run loss as loss — normalizing to the core would overstate
capture. The core charge actually tracked is recorded separately (`q_core_injected_C`,
`core_charge_frac_of_sec1_exit`). The dominant real loss is the upstream sec-1 capture; the
within-stage `q_out/q_injected` additionally includes the **quads-OFF transverse-divergence loss**
(see *Space charge & quads*), which is a no-focusing model artifact — NOT real-machine loss — so
the **end-to-end capture ≈ 4.78 %** of the true-injected (cathode) charge is a **lower bound**, not
a prediction of the real linac's sections-2–8 transmission (the real FODO contains the beam, so the
true value is higher). The robust deliverable is the energy gain, not the transmission.

## Space charge & quads

- **SC OFF by default** (`SPACE_CHARGE=False` ⇒ `Bcurr = 0`): the headline. SC is negligible at
  >25 MeV (transverse SC ∝ 1/γ², γ > 49 at entry). `SPACE_CHARGE=True` (a `config()`-overridable
  module constant) sets `Bcurr = |q_injected|·Bfreq` so Impact-T's single-bunch SC solve carries
  `Q = Bcurr/Bfreq = q_injected` on the `Nxyz³` mesh. This path is **exploratory/unvalidated** (like
  QUADS_ON): the per-section ΔE gates were validated SC-off, and `Nxyz=16` is an order-of-magnitude
  mesh (~1 cell per σ), not converged — re-confirm gates 1/2/5/6 and raise `Nxyz` before relying on
  SC-on numbers. The per-section field-scale **calibration always runs SC-free** (the on-crest ΔE it
  fits is SC-independent at γ>49); `Bcurr` is applied only to the final run deck.
- **Quads present at real tabulated lengths but OFF (`K1 = 0`)** for the headline beam — the
  A→T (current→field) calibrations are undocumented. Each inter-section spacing is a
  `DRIFT_M` (0.4 m placeholder) field-free margin split around the real-length quadrupole
  (`gap/2` drift, quad, `gap/2` drift — the quad length is NOT subtracted from the gap, since
  several real quads exceed 0.4 m). With `K1 = 0` a quad is optically a drift of its length.
- **`QUADS_ON=True` → derived energy-scaled FODO (exploratory).** Turns the inert quads into a
  real focusing lattice. `build_linac_rest_lattice.fodo_quad_gradients(ke_in_mev=…)` derives the
  per-quad `b1_gradient` [T/m] from **accelerator optics**, NOT measured quad current (A→T
  undocumented):
  - **constant per-cell phase advance** (nominal **μ = 50°**), with the gradient **energy-scaled by
    the local beam momentum** Bρ_i = √(KE·(KE+2·mc²))/c at the section *exit* energy (Bρ rises ~5×
    across the line, so a fixed Bg would under-focus downstream);
  - **each gap is a real H/V doublet** — the tabulated machine quad is a focus-pole + defocus-pole
    assembly, so the deck splits it into **two opposite-sign `L_q/2` halves back-to-back** (`quad{N}a`
    lead pole at the alternating sign `(-1)**i`, `quad{N}b` its negation; the halves sum to `qL`).
    A doublet **net-focuses in BOTH planes**, so the envelope stays bounded in x *and* y (a single
    thick quad of one sign would defocus one plane over the multi-metre half-cell and that plane would
    scrape — that was the over-pinched first attempt, 49–58 % transmission). The lead-pole sign
    alternates gap-to-gap.
  - K1 is solved from the **exact thick-lens cell matrix** (`_solve_doublet_k1`), NOT a thin-lens
    formula — the half-quad phase `√K1·(L_q/2) ≈ 0.46 rad (~27°)` is not small. The cell is
    `drift(gap/2)·(+K1 half-quad, L_q/2)·(−K1 half-quad, L_q/2)·drift(gap/2)·drift(L_section(i+1))`
    with the following **RF section treated as a field-free drift**; per gap, K1_i solves
    `cos μ = ½·Tr(cell_i)` by bisection (the symmetric ± doublet gives the same `½·Tr` in both
    planes ⇒ both get μ), then `g_i = (-1)**i · K1_i · Bρ_i` [T/m]. `k1_max` (=14) caps the bracket;
    a cell that can't reach μ within it falls back to `k1_max`. μ=50° is reachable by every gap
    (the weakest, gap 2 / short CEA-2, tops out near 67°) and sits mid-band (0–180° stable).

  Returns **length-`N_SECTIONS` (7)**: the first 6 are the placed quads (`quad2…quad7`, one after
  every section except the last, each placed as an `a`/`b` doublet); the **7th (Q8, after the final
  section 8) is never placed** and is fixed at `0.0` (`QUAD_K[6]` / the helper's last entry are inert
  — the Q8-inert off-by-one). A `QUAD_K=[…]` override (T/m, signed, the per-gap doublet strength)
  replaces the derived list. Two approximations are stamped on **every** QUADS_ON output and must stay
  there: (1) the inter-quad multi-metre RF section is treated as a **field-free drift** (the lattice is
  non-periodic and the beam accelerates through it), so **μ is nominal, not realized** — the acceptance
  is σ_x/σ_y *boundedness*, NOT a measured 50°/cell; (2) the K1
  **magnitude** is guessed (A→T undocumented), so the focusing strength is order-of-magnitude only.
  Label: "placeholder optics — guessed K1, A→T undocumented, H/V doublet (±g qL/2 halves), nominal μ."
- **Calibrate quads-OFF, then build a fresh quads-ON run deck.** The per-section crest-phase + field-
  scale calibration fits `mean_energy`, which is transverse-independent **only on-axis**, so it runs on
  an unconditionally **quads-OFF, zero-quad-radius** deck (the decimated `Np_calib` bunch is neither
  focused nor bore-scraped mid-fit ⇒ gates 1/2 PASS *identically* OFF/ON). The quads-ON path then
  builds a **fresh** deck with the calibrated `scales`, re-applies each section's absolute crest phase
  AND its `rf_field_scale` ControlGroup value (`cal._set_group_scale` — the group is `absolute=True`
  defaulting 0, so a naive scale carry-over would silently run **zero-field**), and asserts
  `quad3.b1_gradient ≠ 0` before tracking. The quads-OFF headline run deck is byte-identical to before.
- **New quad / inter-section-drift bore aperture — gated on `quads_on`, NOT `bore_aperture_on`.** When
  `quads_on`, the quad (and the two inter-section drifts) get `radius = section_bore_radii(i)[1]` — the
  section **exit** bore (the quad sits *downstream* of the section exit at the narrower taper; the solrf
  body already uses the **entrance** radius `[0]`, an intentional taper). Gating on `quads_on` (not the
  already-True `bore_aperture_on`) keeps the quads-OFF path byte-identical (`radius` stays `0.0`) so the
  published 78.5% can't silently regress; the box `XYRAD_M` is never widened.
- **Transverse confinement: the quads-OFF headline does NOT contain the beam — the energy gain
  is the result, transmission is only a no-focusing lower bound.** With **no quad focusing** over
  the 36 m line the beam genuinely diverges. Transmission is measured against the **real tapered
  bore** (`section_bore_radii`, 12.6→9.9 / 14.7→11.7 mm): `BORE_APERTURE_ON` **defaults `True`**,
  so the aperture is the actual machine beampipe — NOT a tunable numerical box. (Deliberate: an
  oversized containment box, e.g. `Xrad/Yrad` widened to 0.30 m, would fake transmission ≈ 1.0 by
  counting a 30-cm-radius beam as "transmitted" — physically meaningless; the real bore makes the
  number physically anchored and un-gameable.) The acceleration partly counteracts the divergence
  via adiabatic damping (σ_r ∝ 1/√(γβ), γβ ≈ 50 → 605 over the line ⇒ ~3.5× shrink), so a fraction
  passes the real bore, but the rest scrapes — measured **transmission ≈ 78.5 %** (1783/2271
  macroparticles through the real bore). **This is a no-focusing MODEL ARTIFACT, NOT real-machine
  loss** — the real FODO lattice (quad A→T calibrations undocumented, `details.md`) contains the
  beam, so the true transmission is higher. The reported transmission (count-based `n_out/n_in`
  against the real bore, measured *before* the openPMD charge re-imposition so it can never be
  masked) is therefore a **quads-OFF lower bound, not a prediction**; capture/transmission through
  sections 2–8 is not predicted by the quads-OFF
  headline. The robust, quad-independent deliverable is the **longitudinal physics** (exit ⟨KE⟩,
  per-section ΔE), which does not depend on transverse confinement.
- **`QUADS_ON` contains the beam ⇒ both transverse planes are BOUNDED — that is the deliverable, NOT
  transmission.** With the μ=50° doublet the measured envelope stays bounded and **out-of-phase
  oscillating** the full 36 m: σ_x ≈ 0.6–4.4 mm, σ_y ≈ 0.7–3.9 mm (FODO beating, no blow-up), vs the
  quads-OFF monotonic rise. The longitudinal headline is preserved (exit ⟨KE⟩ ≈ **309.0 MeV**, gates
  1/2/5/6 PASS; chromatic εn growth as expected, not a runaway). The soft `envelope_in_bore` 3σ is
  ~13 mm (vs ~19 mm no-focus / ~25 mm an over-pinched single-sign quad) — much improved, the win is
  the **bounded** envelope visible in `fodo_optics.png`, not the 3σ < narrowest-bore threshold.
- **Transmission lands ≈ the quads-OFF baseline (~78.2 %, 1775/2271), NOT above it.** This is expected
  and is **not** a regression of the focusing, for two structural reasons: (1) the doublet halves now
  carry the real **exit-bore `radius` aperture** that the no-quad baseline lacks — extra scrape planes
  the quads-OFF run simply does not have; and (2) the injected beam **expands σ from ~1.2 mm to ~4.4 mm
  over the first ~3.5 m** (through section 2) **before the first quad** (placed *after* section 2) can
  act — the real lattice has no quad ahead of section 2, so that initial divergence is already past the
  narrowing bore by the time focusing begins. So there is **no "> 78.5 %" or predictive-transmission
  claim** — transmission stays a soft/print-only number ≈ baseline; the contained, bounded σ_x/σ_y is
  the exploratory result. Still placeholder optics — guessed K1, A→T undocumented, nominal μ.

## e+ compressor (CU 2)

**Out of scope** — its lattice role is not established in `details.md`; no converter topology
is asserted. Future work.

## Files

| File | Role |
|------|------|
| `build_linac_rest_lattice.py` | Per-section table + `√P` scaling helpers; `fodo_quad_gradients(...)` (derived energy-scaled FODO K1, exploratory); `build_impact(...)` assembles the chained 7-section deck in-memory (reuses vendored `rfdata4–7`, drifts, quads). |
| `linac_rest_sim.py` | `main()`: handoff IN (captured core from `linac_sec1`) → build deck → calibrate → `I.run()` → §5 validation gates → openPMD handoff OUT + `injection_summary.json`. |
| `calibration.py` | `calibrate_sections` (per-section local-crest + scale fit) and `validate_run` (§5 gates). |
| `plot_linac_rest.py` | Energy/σ_KE/emittance vs z (from `I.stat`), per-section target-vs-achieved bars, FODO transverse envelope σ_x/σ_y. |
| `rfdata/rfdata4–7` | Vendored S-band TW field shapes (committed). |

## Performance knobs (via `linac_rest.config(...)`)

`POWER_MW`, `PHASE_DEG`, `MIN_KE_MEV` (capture cut), `Np` (tracked macroparticles), `Ntstep`
(step cap; the run asserts `mean_z` reached the final zedge so a truncated run fails loudly,
NOT silently), `Dt`, `Nxyz`, `DRIFT_M`, `QUADS_ON` (default `False`; `True` ⇒ the exploratory
derived-FODO path — `fodo_quad_gradients` auto-derives the per-quad K1 from optics, no need to
supply `QUAD_K`), `QUAD_K` (optional signed-`b1_gradient` [T/m] override; `None` ⇒ auto-derive
when `QUADS_ON`, else zeros — only `QUAD_K[0..5]` / Q2–Q7 are placed, the 7th entry/Q8 after the
last section is **never installed**, the Q8-inert off-by-one), `BORE_APERTURE_ON` (default
`True` — the real tapered bore is the aperture; set `False` only to disable bore scraping, e.g.
a divergence study), `Np_calib` (decimated bunch for the per-section calibration;
the final run uses full `Np`), `OUTDIR`. SC-off + the seeded/decimated calibration make the run
cheap (~tens of seconds; the 7-section calibration alone is ~140 s at `Np_calib=400`, down ~7×
from a full-bunch full-scan calibration); `ImpactTexe` is serial.

## Validation gates (§5, printed each run)

1. per-section ΔE within ±3% of target (hard);
2. cumulative exit ⟨KE⟩ ≈ measured ⟨KE⟩_in + Σ ΔE_target (hard, ±3%);
3. σ_KE absolute + relative — diagnostic. NOTE: absolute σ_KE is **not** conserved — it
   **grows ~3.9×** (5.42 → 21.21 MeV, measured) because the crest is an energy *maximum*, so a
   bunch of finite phase length accumulates a second-order correlated spread over the 7 sections
   (each adding ~40 MeV on a cosine; the curvature term the plan called "negligible" is sizeable
   at 282 MeV total gain). The **relative** spread still shrinks (20.0% → 6.9%) only because ⟨KE⟩
   grows faster (27 → ~308 MeV, ~11.4×) than σ_KE. This affects spread only — calibration uses
   ⟨KE⟩, so the energies/scales are unaffected;
4. normalized emittance εn,x/εn,y in vs out — diagnostic. NOTE: quads-OFF εn is **not**
   conserved; the recorded vs-z εn sawtooths ~2.4× — a fort.10N `norm_emit` artifact at
   bore/section crossings (σ_x stays smooth across the jumps ⇒ not physical growth);
5. beam reached the final zedge (`I.stat("mean_z")[-1] ≈ Σ L`) — hard, catches Ntstep
   truncation (which falsely reports `finished=True`);
6. min captured KE ⇒ β > 0.999 — hard, justifies no-slip;
7. transmission (count-based `n_out/n_in`) — diagnostic + **no-focusing LOWER BOUND, not a
   prediction**: quads-OFF the beam diverges and a fraction leaves the transverse domain (model
   artifact, not real-machine loss — the real FODO contains it). With `QUADS_ON` it lands **≈ the
   quads-OFF baseline (~78.2 %), NOT above it** (the doublet's win is the bounded envelope, not
   transmission — see *Space charge & quads* for the two structural reasons); stays **soft/print-only**,
   never a hard gate (K1 is guessed).
8. `envelope_in_bore` (**soft**, print-only, `QUADS_ON`-relevant) — `3σ·max(σ_x, σ_y) < min bore`,
   a **conservative 3σ** check on whether the derived FODO contains the envelope. NOT a hard gate (K1
   unvalidated — it must never gate the energy headline). It currently prints **FAIL both OFF and ON**
   because the **3σ multiple** exceeds the ~9.9 mm narrowest bore (~19 mm OFF, ~13 mm ON) — but the
   **RMS σ itself stays well inside** the bore quads-ON (≤ 4.4 mm vs 9.9 mm), so the FAIL is the
   conservative 3σ edge, not the beam hitting the pipe. It is therefore **not** a "FAIL-when-OFF /
   PASS-when-ON" liveness test — the meaningful improvement is the much smaller 3σ and the **bounded,
   oscillating** σ_x/σ_y with quads ON (vs the OFF monotonic rise), seen in `fodo_optics.png`.

## Gotchas (Impact-T / lume-impact)

- **In-process, no subprocess** — but still `_prepare_environment()` (repo-root chdir +
  RLIMIT_NOFILE raise) + `setup_logging()` + `ImpactTexe` stdout redirect.
- **Ntstep truncation reports `finished=True`** falsely (stops mid-line) — size `Ntstep` from
  the lattice length and **assert `mean_z` reached** (gate 5).
- **`autophase()` no-arg throws** on a non-cathode beam — calibration scales/phases per
  section explicitly (coarse crest scan + `brentq`), not `I.autophase()`.
- **`ParticleGroup.write()` is incompatible** with `openpmd-viewer` (emits openPMD 2.0 with a
  STRING `openPMDextension`; the viewer wants the integer ED-PIC extension). The handoff OUT
  uses `pipeline.impact_io.write_openpmd_particles` (replicates WarpX's byte layout).
- **Species name asymmetry**: `ParticleGroup.species` is `"electron"` (singular); the openPMD
  output and every cross-stage reader key on `"electrons"` (plural). `impact_io` translates.
- **`write_beam` elements break calibration** (`track_to_s`'s `load_many_fort` raises
  "Conflicting data for key:t" on the fort.10N stat columns), so per-section vs-z evolution
  comes from `I.stat(...)` (written to the summary), not particle slices. `write_beams=False`.
- **The generated `ImpactT.in`/`rfdata*` live in lume-impact's temp workdir** (not in-tree);
  the vendored source `linac_rest/rfdata/` is committed and is the single source of truth.

Run:
```bash
conda activate CBB
python -c "import linac_rest; linac_rest.run()"          # build deck + Impact-T + plots
python -c "import linac_rest; linac_rest.run(plots=False)"
python -c "import linac_rest; linac_rest.plot()"
```
