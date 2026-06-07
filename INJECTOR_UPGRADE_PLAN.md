# Cornell Linac Injector Upgrade — Implementation Plan

## Summary

This upgrade makes the WarpX rebuild faithful to Adam Bartnik's LinacSim across the entire low-energy front end by collapsing the unfocused, single-cavity `prebuncher/` stage into one self-consistent RZ **`injector`** stage that contains the full LinacSim injector: Lens 0A → Prebuncher 1 → Prebuncher 2 → Sol 0 / Lens 0E, handing off a focused beam to `linac_sec1` at the true linac entrance z ≈ 2.03 m.

What changes and why:
- **Task 1 (power):** Nothing changes. **8 kW is correct and faithful** — it is the GUI default `prebuncher1_input_power` and the WarpX stage already reproduces it (V_gap ≈ 58.6 kV, scale ≈ 0.133). It is *intentionally* weak: the gun bunch (σ_z ≈ 1 mm) sits within ~0.1% of one 214 MHz RF period, so a single cavity can only apply a small linear chirp; the single-cavity bunching threshold is ~95 kW (scale ≈ 0.46) because the cavity must first cancel the gun's intrinsic +1.40 keV/mm debunching chirp, and space charge at ~0.8–1 nC / β≈0.63 dominates the drift. The real design is a **two-cavity + solenoid distributed buncher**, not one strong kick. The earlier 160–800 kW scan was exploration, not a correction; do not raise the default. We only sharpen the README to state this verdict up front.
- **Task 2:** Add Prebuncher 2 (z=1.318 m, 10 kW, Q=4300), installed **reversed** via a mirrored field map.
- **Task 3:** Add the three energized solenoids (Lens 0A 6 A, Sol 0 40 A, Lens 0E 10 A) as static B-only maps in the injector, and **remove the `linac_sol.h5` in-linac hack** — this is the physical fix for the known 68% radial-scrape loss in `linac_sec1`.
- **Task 4:** Add a cross-stage `pipeline/plot_chain.py` writing `results/chain_evolution.png` (⟨KE⟩, ε_n, σ_x, σ_z, transmission, I_peak vs lab-z across the whole chain) plus an emittance-growth budget, transmission waterfall, and scorecard.

## Architecture decision

**Chosen approach: extend the prebuncher into a single RZ `injector` stage; do NOT add elements piecewise.** (Sections 2, 3, and 5 all converge on this; Sections 2 and 3 flagged it as a dependency and Section 5 resolves it.)

Rationale:
1. **LinacSim models all of these in ONE GPT drift with space charge** (z = 0.06–2.1 m carries both prebunchers and all lenses/Sol 0). The bunching, two-cavity phasing, and transverse focusing are coupled *through the self-field*. Piecewise WarpX stages would re-inject a snapshot and re-solve Poisson from a re-seeded beam at each boundary, destroying the continuous space-charge history that sets both bunch length and radial envelope.
2. **It is the physical fix for the 68% radial scrape.** The loss exists *because* there is no transverse focusing between gun and linac. Applying Lens 0A + Sol 0 / Lens 0E focusing *in the same drift as the bunching* keeps the beam inside the bore. The `linac_sol` hack (SOL_Z=−0.60, mis-located inside the linac) stands in for Sol 0 / Lens 0E that physically live at z ≈ 1.9 m; moving them upstream to their true z makes them physical and lets the linac stop scraping.
3. **The contract stays four stages** (`cathode → gun → injector → linac_sec1`); only the handoff plane moves from 1.30 m to 2.03 m, and the linac sheds a field it should not own.

Piecewise (`lens0a/`, `preb2/`, `sol0/` stages) is rejected: it multiplies subprocess handoffs, re-seeds the space-charge solve at each, and still cannot co-locate focusing with bunching.

**Implementation note on naming:** `injector/` is created by copying `prebuncher/`. The old `prebuncher/` package is removed in the same change (its logic migrates wholesale). All `prebuncher.*` facade calls and the `P8_zc` output dir disappear.

### Final inter-stage contract

| Stage | Reads | Writes |
|-------|-------|--------|
| `cathode/cathode_diode.py` | — | `cathode/diags/particles` |
| `gun/gun_sim.py` | `cathode/diags/particles` + `gun/gun_field/gun_E.h5` | `gun/diags` |
| **`injector/injector_sim.py`** | `gun/diags/particles` + `injector/injector_field/{preb1_EB,preb2_EB,lens0a,sol0,lens0e}.h5` | `injector/diags/main` |

> **As-built note:** Preb 2 reuses the **forward** `preb2_EB.h5` (same field as `preb1_EB.h5`, gap at z=1.318 m) — the reversal is `PREB2_REV_PHASE=0` at run time (the crest-referenced drive auto-absorbs the −1,0,0 install; see the Task 2 FINAL VERDICT). There is no `preb2_EB_rev.h5` mirrored map; the `_rev` suffix used earlier in this plan is superseded.
| `linac_sec1/linac_sec1_sim.py` | **`injector/diags/main/particles`** (snapshot at the z ≈ 2.03 m handoff plane, see Task 3) + `linac_sec1/linac_sec1_field/{linac_rf1,linac_rf2}.h5` | `linac_sec1/diags/main` |

- Default chain writes `injector/diags/main`; the power/phase scan facility survives as an optional `OUTDIR` override (`injector/diags/<case>`). The facade `resolve_outdir()` must return `injector/diags/main` by default (see "File-by-file changes / facade" below).
- `linac_sec1` no longer builds `linac_sol.h5` or carries `I_SOL`/`SOL_Z`/`SOL_MAP`; its only applied fields are the two SLAC quadrature maps.

### Domain & grid for the injector box

- **z-domain:** `ZMIN=0`, `ZMAX = 2.10 m` (LinacSim prebuncher-subsection ZSTOP). Handoff snapshot taken at the **z = 2.03 m** plane (Z_acc_1), with a field-free exit drift past 2.03 m so the handoff beam coasts (the "stop in the drift, not at the wall" rule). The snapshot must be selected *by z-proximity to 2.03 m*, NOT by min-σ_z or max-in-bore-charge — see Task 3.
- **r-domain:** keep `RMAX = 0.036 m`, `NR = 80` (dr = 0.45 mm). RF maps r-padded to RMAX as today; lens/sol maps already reach 40 mm. **Keep `NR=80` — do not copy the linac's `NR=16`** (the RF map reaches 36 mm and needs the radial resolution).
- **Cell aspect (binding constraint):** hold dz ≈ 1.26 mm to keep dz/dr ≈ 2.80:1 (the ≈3:1 rule) or MLMG diverges. `NZ = 2.10 / 0.00126 ≈ 1654 → 1664`. **Justification: at `NR=80` (dr=0.45 mm), `NZ=1664` gives dz=1.262 mm ⇒ dz/dr = 2.80:1, and 1664 is divisible by the blocking factor 8.** (1656 is also ÷8; 1664 is chosen as the next ÷8 value that keeps dz slightly under 1.27 mm with margin.) **This is unrelated to the linac's NZ=1664** — the linac reaches 2.8:1 via `NR=16`/dr=0.75 mm, a different dr; the coincidence of the same NZ value is not the rationale. **Do not coarsen NZ.**
- **Cost:** cells and steps each rise ≈ 2.10/1.30 = 1.6× vs the current prebuncher. This long-thin box is **convergence-bound, not cell-bound** (the MLMG solve at NZ=1024 is already near its iteration cap; per CLAUDE.md), so iterations-to-tolerance over the 1.6× longer box rise super-linearly and the wall-time penalty likely **exceeds** a naive 1.6–2.0×. Since the prebuncher was ~75% of the chain, the full chain more than doubles unless knobs loosen. Mitigate via the documented levers (`CFL` up, `MAX_ITERS`/`REQUIRED_PRECISION` looser) — never by coarsening NZ. **Take an actual wall-time + `MLMG failed` divergence checkpoint at build step 2 (Preb-1-only on the extended box) before committing to two cavities + 3 solenoids.** The linac gets cheaper (loses the solenoid add + heavy radial scrape), partly offsetting.

## Task 1: Prebuncher power

**Verdict: 8 kW is faithful to LinacSim — do NOT change it.** No code/constant changes. The only deliverable is documentation that pre-empts the "are you sure?" question.

Source-of-truth constants are all correct and stay as-is (now in the `injector` build, see Task 2): `V1J_KEV = 438.6`, `F_RF = 214.18 MHz`, preb-1 `Q_L_1 = 3000`, `Z_GAP_CENTER_1 = 0.534`, scale formula `sqrt(1e3·Q·P_kW/(2π f_RF))` ⇒ scale ≈ 0.1335, V_gap ≈ 58.6 kV.

Documentation changes (in the new `injector/README.md`):
1. **State the verdict up front:** "8 kW is the faithful LinacSim default and is intentionally weak — single-cavity bunching is not the design; the injector is a two-prebuncher + solenoid distributed buncher (each cavity ~12× below the ~95 kW single-cavity threshold in power, ~4× in voltage)."
2. **Keep the stale-chirp caveat:** the −3.05 keV/mm cavity coefficient and +1.40 keV/mm gun chirp were measured at 0.1 nC; at the reconciled ~1 nC charge, space charge is stronger and the threshold is *higher*, so 8 kW is even further below threshold — reinforcing the verdict.
3. **Cross-reference downstream:** the old 8 kW + single-cavity + no-injector-solenoid config is *why* `linac_sec1` received a radially-blown-up beam and dropped ~68% of charge — a modeling-incompleteness artifact (missing Preb 2 + Lens 0A/Sol 0/Lens 0E), not a power error. Tasks 2–3 fix it.
4. Keep `POWER_KW=8` (Preb 1) and `10` (Preb 2, Q=4300) as defaults; the 160–800 kW scan remains available as a non-default `OUTDIR` override, documented as exploration.

## Task 2: Second prebuncher

> **⚠ CORRECTION 2 (phase base convention — found at the step-1/2 boundary).**
> The GUI `phi_off` values (−70° Preb-1, −45° Preb-2) are **referenced to ON-CREST**, but the
> original `make_cavity` stacks them on `base = π/2` (the zero-crossing "zc" reference). Stacking a
> crest-referenced offset on a zc base puts Preb-1 at ~160° from crest (near anti-crest) ⇒ ~−62 keV
> deceleration + debunching (impl measured exactly this; physics predicted −55 keV). **Fix: when
> applying a GUI `phi_off`, use the CREST base:** `phi = −ω·t_gap + π + radians(phi_off_deg)`. Then
> −70° means 70° from crest = the faithful LinacSim point (mild +~20 keV kick, strong bunching slope).
> `PHI_OFF_1_DEG = −70` is the correct faithful value — only its application was wrong; 8 kW unchanged.
> The "`base = π/2` for zc / mirror-not-base" guidance below is amended accordingly: the faithful
> operating point uses `base = π` (crest) + GUI `phi_off`; the bare `zc`/`crest` PHASE knob remains only
> for the exploratory scan. Document `phi_off` as crest-referenced in `injector/README.md`.
>
> **⚠ CORRECTION (physics validation, supersedes the original Task-2 mechanism below).**
> The team measured the actual `prebuncher_25D.gdf` z-flip parity about the gap centre, off-axis at
> beam radii (corr ±0.9999): **Ez EVEN, Er ODD, Bφ EVEN** — a definite-parity TM0 standing wave
> (Bφ ∝ ∂Er/∂z − ∂Ez/∂r, both even, also required by Maxwell). *(An initial "Bφ odd" reading was
> on-axis numerical noise where Bφ ≈ 1e-8; corrected after re-measuring across all radii.)*
> This inverts the original plan's reversed-install mechanism:
> - The originally-prescribed build ("z-reverse all components + negate Er & Bφ, keep Ez") is WRONG
>   for this map and silently mis-models the cavity. The stated self-check ("preb2 Ez(z) mirrors
>   preb1") would falsely PASS because an even Ez mirrors to itself. **Do not build the asymmetric map.**
> - The geometric `-1,0,0` install is a **180° ROTATION**, not just a spatial z-mirror: it flips the
>   field VECTOR (ẑ, φ̂ → −) in addition to remapping z → −z. Folding the vector flip with the measured
>   parities (Ez even, Er odd, Bφ even), all three components flip sign ⇒ a global E,B sign flip ⇒ a
>   `+π` phase shift **in the ABSOLUTE drive frame**.
>
> **FINAL VERDICT (physics + empirical, LOCKED): `PREB2_REV_PHASE = 0` — NO extra π.** The sim does
> not phase in the absolute frame; it references phase to the **crest of the LOADED field** (`base = π`
> = max −cos of whatever map is loaded) plus the GUI `phi_off`. Since the loaded map *is* the reversed
> cavity's field, its crest already sits at the flipped phase — **crest-referencing auto-absorbs the
> reversal**, so adding +π would double-count. Equivalently: reversed(φ) ≡ forward(φ+π), and the GUI
> on-crest reference for Preb-2 (178.9°, vs Preb-1's 304.7°, Δ=125.8°) is defined for the
> already-reversed install, so the geometric +π and the GUI's built-in reversal cancel.
> **Verified empirically** (impl ran both ways; the discriminator is the chirp-slope change, not net ΔKE):
> `rev_phase=0` → dchirp −0.33 keV/mm (compressive = bunching ✓); `rev_phase=+π` → +0.04 keV/mm + strong
> decel (debunching ✗). Preb-2 is confirmed a faithful 2nd velocity-buncher at −45° from crest
> (V_gap 78.4 kV, Q=4300, 10 kW: bunching-slope ~55 kV-equiv, net +55 kV-equiv — more net accel than
> Preb-1's −70°, sensible for the cavity nearer the linac).
>
> **Corrected implementation:** reuse the SAME forward `preb1_EB.h5` map for BOTH cavities; keep
> `PREB2_REV_PHASE` as a named constant **set to 0**. No second `.h5` is built. The code comment +
> README MUST pre-empt a future reader re-adding the π ("the GUI on-crest phase absorbs the −1,0,0
> reversal; crest-referencing the loaded field cancels the geometric +π; verified empirically").
> (Quantified phasing: inter-cavity drift 0.534→1.318 m at β≈0.628 ⇒ ~321° advance at 214 MHz;
> constant-v error <1° at the sub-threshold 8 kW operating point — injection β is fine.)

### What "reversed install" (`-1,0,0`) means in WarpX — ORIGINAL (superseded by the correction above)

GPT's `Map25D_TM(..., -1,0,0, ...)` for Preb 2 is a 180° rotation of the cavity about a transverse axis — a z-mirror (parity) of the axisymmetric m=0 standing-wave TM map about the gap centre. Under z → −z: the **Ez longitudinal profile is preserved but Er and Bφ flip sign**. WarpX's `read_from_file` / `LoadAppliedField` reads meshes verbatim and only multiplies them by a scalar `warpx_{E,B}_time_function`; a single scalar scales Ez and Er together and *cannot* flip Er relative to Ez. ~~Therefore the reversed install must be baked into a separately-built mirrored field map.~~ *(This reasoning assumed a generic z-reversal; for this map's measured definite parity it is wrong — see the correction block above. A scalar π phase shift is correct.)*

**Polarity bookkeeping (reviewer-flagged silent-error risk).** Because the map's spatial parity is baked in, the time-function for Preb 2 must use the **same cos/sin convention as Preb 1** (`base = π/2` for `zc`, `+π` for crest). The mirror handles the spatial Er/Bφ parity; the time function handles only arrival phase. The reversed cavity's GUI on-crest reference (178.9° vs Preb-1's 304.7°) is absorbed entirely by `PHI_OFF_2_DEG = −45°` plus the arrival term `t_gap2`. **Do not add an extra sign to the Preb-2 time function and do not double-count the reversal in `base`.** Whether the geometric Ez polarity flip needs any compensation in `base`/`phi_off` is exactly what build step 3's validation must resolve empirically: validate that the −45° offset on the mirrored map produces the **intended kick sign and bunching**, not merely that Ez(z) is a geometric mirror.

**Highest-risk item:** validate by confirming `preb2`'s on-axis Ez(z) is the mirror of `preb1`'s AND that the per-cavity transit-time energy-kick sign at the −45° design phase is as intended (against a single-cavity GPT kick), before trusting two-cavity bunching.

### File-by-file changes

**`injector/build_injector_field.py`** (migrated from `build_prebuncher_field.py`):
- Constants: rename `Z_GAP_CENTER → Z_GAP_CENTER_1 = 0.534` and `Q_L → Q_L_1 = 3000` (keep back-compat aliases until the sim is migrated). Add `Z_GAP_CENTER_2 = 1.318`, `Q_L_2 = 4300`, `PHI_OFF_1_DEG = -70.0`, `PHI_OFF_2_DEG = -45.0`, `MAP_HALF_Z` (= 0.1524 m).
- Add `OUT_FILE_1 = .../preb1_EB.h5` and `OUT_FILE_REV = .../preb2_EB_rev.h5`.
- Refactor the inline `write_mesh` + writing block into `write_field(out_file, z_offset, mirror=False)`. When `mirror=True`: reverse every component along z (`arr[:, ::-1]` on the z axis) and **negate Er and Bφ**, set `z_offset = Z_GAP_CENTER_2 - MAP_HALF_Z`. **State explicitly that `z_offset = Z_GAP_CENTER_2 − MAP_HALF_Z` is correct only because the prebuncher map is symmetric in z-extent (±152.4 mm about the gap); add an assertion that the map's native z-extent is symmetric before applying it.**
- `main()` writes both: forward at `Z_GAP_CENTER_1`, mirrored at `Z_GAP_CENTER_2`. Keep the `V1J_KEV` assertion (mirroring preserves ∫|Ez|).

**`injector/injector_sim.py`** (migrated from `prebuncher_sim.py`):
- Import `Z_GAP_CENTER_1, Z_GAP_CENTER_2, V1J_KEV, F_RF, Q_L_1, Q_L_2, PHI_OFF_1_DEG, PHI_OFF_2_DEG`.
- Module constants (tunable via `injector.config(...)`): `PREB1_KW=8.0`, `PREB1_Q=3000`, `PREB1_PHI_OFF=-70`, `PREB2_KW=10.0`, `PREB2_Q=4300`, `PREB2_PHI_OFF=-45`, `PREB2_REVERSED=True`, `PREB1_FIELD=".../preb1_EB.h5"`, `PREB2_FIELD=".../preb2_EB_rev.h5"`. Set `ZMAX=2.10`, `NZ=1664`.
- Factor RF setup into a helper (**keep `.10e` precision on every term — `ω·t` truncation accumulates over the ~5 ns transit at 214 MHz**):
  ```python
  def make_cavity(field_path, power, q_l, z_gap, v_at_gap, phi_off_deg, phase, omega):
      scale = sqrt(1e3 * q_l * power / (2*pi*F_RF))
      t_gap = (z_gap - Z_INJECT) / v_at_gap     # caller passes accumulated arrival for preb2
      base  = pi/2 if phase == "zc" else pi
      phi   = -omega*t_gap + base + radians(phi_off_deg)
      e_time = f"{scale:.10e}*cos({omega:.10e}*t + ({phi:.10e}))"
      b_time = f"{scale:.10e}*sin({omega:.10e}*t + ({phi:.10e}))"
      return picmi.LoadAppliedField(read_fields_from_path=field_path, load_E=True, load_B=True,
                                    warpx_E_time_function=e_time, warpx_B_time_function=b_time)
  ```
  The new `phi_off` term reproduces LinacSim's 304.7° / 178.9° on-crest definitions (the current code hardcodes `−ωt_gap + π/2` with no offset). `base` is `π/2`/`π` for **both** cavities (the mirror, not `base`, encodes Preb-2's reversal).
- **Two-cavity phasing (reviewer-flagged):** cavity 1 arrival uses `v_beam`. Cavity 2 arrival = `t_gap1 + (Z2 − Z1)/v_after_preb1`. `v_after_preb1` is the **post-cavity mean β of the actual injected bunch at Preb-1 exit**, *not* re-derived from an analytic on-crest energy. **Mechanism constraint:** the Preb-2 time-function string is baked at field-construction time, before WarpX integrates Preb 1, so the true post-Preb-1 β is not yet known. At the faithful 8 kW `zc` operating point the mean kick is ~0 ⇒ `v_after_preb1 ≈ v_beam`, which is what we use. **Document explicitly in the README and in a code comment: Preb-2 phasing uses the injection β and is valid only while both cavities are sub-threshold (the design case); a hard Preb-1 power scan desyncs the Preb-2 phase reference.** If a future hardened-Preb-1 study is needed, it requires a two-pass run (read post-Preb-1 β from a diagnostic, rebuild Preb-2 time function). Quantify the constant-v phase error in degrees at 214 MHz in the README.
- `main()` builds cavity 1 (forward map); if `PREB2_KW > 0`, builds cavity 2 (mirrored map). Both added with `sim.add_applied_field(...)` — **after** the solenoid B-only maps (Task 3 ordering).
- Transit / `n_steps` uses the longer ZMAX; for the default `zc` case `transit = (zmax − Z_INJECT)/v_beam` still holds (both kicks sub-threshold).
- Default OUTDIR → `injector/diags/main`.

**`injector/plot_injector.py`** (migrated from `plot_prebuncher.py`):
- Extend the cavity figure to overlay both on-axis Ez(z) lobes (at Z1, Z2) and both RF waveforms at their arrival times; re-derive Preb-2's scale/phase the same way the sim does (single source of truth). Add vertical gap markers at Z1 and Z2 to the σ_z(z) and bunch-profile panels.

**`injector/__init__.py` (facade) + `resolve_outdir`:** reimplement `resolve_outdir()` so it returns `injector/diags/main` by default (the old `prebuncher` facade derived `P{kw}_{phase}`; that derivation is gone). `_beam_summary` in `run_pipeline.py` must read `injector/diags/main`. When an `OUTDIR` override is set for a scan, `resolve_outdir()` returns that override.

## Task 3: Solenoid lenses

### Active solenoids and placement (Section 3, corrected per reviews)

Three energized maps (GUI defaults); 0B/0C/0D are 0 A and omitted; Sol 1A/1B/1C are downstream of Section 1 and omitted.

| Map | Current | GUI lab-z | Native peak z | Native grid z-range | `grid_global_offset` |
|-----|---------|-----------|---------------|---------------------|----------------------|
| LENS_0A | 6 A | 0.225 m | **0.2333 m** | 0.0–0.5 m | **0.0** (accept ~8 mm peak offset, see below) |
| SOL_0 | 40 A | 1.897 m | **0.8209 m** | 0.0–2.5 m | **+1.0761 m** (= 1.897 − 0.8209) |
| LENS_0E | 10 A | 1.914 m | **1.9147 m** | 0.8–2.4 m | **0.0** |

(All three peak locations and grid ranges measured directly from the `.gdf` files: SOL_0 peak 0.8209 m on a 0–2.5 m grid; LENS_0A peak 0.2333 m on a 0–0.5 m grid, nr=189/nz=16; LENS_0E peak 1.9147 m on a 0.8–2.4 m grid.)

**Conventions, corrected:**
- **LENS_0A and LENS_0E store absolute lab z in their native grid** (LENS_0E index-0 at native 0.8 m = lab 0.8 m; peak native 1.9147 = lab 1.9147 ≈ GUI 1.914). So `grid_global_offset = 0` is correct *because the map stores absolute lab z*, and the offset is applied to **grid index 0** (lab-z of index 0), not the peak.
- **LENS_0A has an ~8 mm discrepancy** (native peak 0.2333 m vs GUI 0.225 m). The map is also **coarse/transposed in z (nr=189, nz=16 ⇒ ~31 mm/cell over 500 mm)**, far coarser than the 1.26 mm sim grid, so WarpX will interpolate a very coarse axial profile. **Accept offset = 0 and document the 8 mm peak offset and the coarse-axial-interpolation caveat** (the 8 mm is below the lens's own ~30 mm cell size, so a correcting shift is meaningless at this map's resolution).
- **SOL_0 is the exception:** its native grid stores a *local* 2.5 m window (index-0 at native 0, peak at native ~0.813 m — the long solenoid + its 2.5 m steel frame) but the solenoid belongs at lab 1.897 m. So `grid_global_offset = 1.897 − 0.813 ≈ +1.084 m`. (This is why the old linac hack used SOL_Z=−0.60 to drag the peak into the linac-local frame.) **Derive the offset programmatically** (`offset = GUI_z − Z[argmax|Bz|]`) rather than hard-coding, so it cannot drift.
- **PLACEMENT RULE LOCKED (physics, after a native-Z detour):** place each solenoid's field PEAK at its **GUI element z** (0.225 / 1.897 / 1.914) via the programmatic offset — the unifying faithful rule. A "preserve native Z" rule was considered and **rejected**: the three maps use different z-conventions (LENS_0E's grid is pre-shifted into absolute lab z; LENS_0A is a local window near lab 0, a near-coincidence; SOL_0 is a *local* 2.5 m window). Decisive evidence SOL_0's native is local: the old `linac_sec1` hack used `SOL_Z=−0.60` to slide SOL_0's 0.813 m peak, its code comment stating "SOL_0 peaks 813 mm into ITS OWN GRID." Native-Z would place the 40 A solenoid at 0.813 m (mid-injector) — a ~1 m error vs its iris-matching role at 1.897 m. **LENS_0A: ship GUI 0.225** (the native-vs-GUI 8 mm is below the map's ~31 mm/cell axial resolution, so neither is "more accurate"; GUI is the consistent rule). Note the **~7× capture sensitivity to LENS_0A's 8 mm** — capture is a tune-sensitive ~order-1% conservative lower bound, not a precise number; the step-7 current-scan characterizes it.

**z-overrun (reviewer-flagged, BLOCKING-VERIFY).** With these offsets the *native* grids exceed the 2.10 m box: SOL_0 spans lab 1.076–3.576 m, LENS_0E spans lab 0.8–2.4 m. WarpX `read_from_file` tolerates a map larger than the domain (it reads the in-domain window), so this runs — but **the build must add, per solenoid, an assertion that the field's physical lab-z peak falls inside [0, 2.10] m** (`assert 0 <= lab_peak <= 2.10` and `assert abs(lab_peak − GUI_z) < 0.002` for SOL_0/LENS_0E; for LENS_0A use a `< 0.009` tolerance to accommodate the documented 8 mm). The per-solenoid sanity report must print the **in-domain physical lab-z peak**, not the native peak.

**Peaks are upstream of the handoff plane.** SOL_0 peak (1.897 m) and LENS_0E peak (1.914 m) are upstream of the z≈2.03 m handoff and so focus the injected beam; their map *tails* extend past 2.03 m but the beam is extracted at 2.03 m so the tail is harmless in the injector. **Add an assertion that every focusing peak is upstream of the snapshot plane** (`assert lab_peak < Z_HANDOFF`) so the linac never inherits a beam still inside a solenoid the linac does not model.

**Grid caveat (three separate files):** LENS_0A is (nr=189, nz=16); the others are (nr=16, nz=601). Different grids/offsets ⇒ write **three separate single-mesh openPMD files**, each with its own `LoadAppliedField`. Do not combine.

**⚠ ADDED SCOPE — 9.547 mm injector collimator (physics-found faithfulness gap, greenlit).** `gpt_master.in`
collimates the injector at a **9.547 mm iris at z=1.922 m** (`scatteriris` + a 9.547 mm pipe to 2.1 m), just
downstream of the Sol 0 / Lens 0E peaks — the matching solenoids exist to squeeze the beam THROUGH this iris.
The rebuild had only the 36 mm domain wall, so the uncollimated handoff charge (0.766 nC) overstated
transmission against the wrong aperture. **Decision: enforce the 9.547 mm collimation.** Implementation
(impl picks the cleanest reliable option): (1) preferred — an explicit collimator absorbing r>9.547 mm for
z≥1.922 m in the injector (embedded boundary / z-dependent radial scrape) so the injector's own handoff charge
is the real collimated number; (2) acceptable minimal — apply the 9.547 mm radial cut at the 2.03 m handoff
snapshot AND keep the linac `BORE_R=9.55 mm` injection cut, documenting the 1.922 m iris as the basis (same
physical beam since the envelope only grows past 1.922 m). **RMAX/aperture = 9.547 mm — do NOT widen the linac
RMAX to the 16 mm re-expanded envelope** (would accept charge the real machine scrapes). Faithful 6/40/10 A
currents are KEPT (no retune — the 1.45 m-waist→2.03 m re-expansion is partly the γ²≈1.66× ES space-charge
overestimate, so retuning would bake the artifact into unfaithful currents). Capture is reported THROUGH
9.547 mm as a **conservative lower bound** (γ² pessimism). The current-scan (#7) is optional exploration, now
an honest "focus through 9.547 mm" optimization rather than against a 36 mm wall.

**Current → scale:** maps are 1-A-normalised, scale linearly. In WarpX: static constant `warpx_B_time_function=f"{current:.8e}"`, `load_E=False, load_B=True`. No cos/sin/ω.

**Ordering gotcha:** WarpX forces the global `E_ext_particle_init_style` to "none" if the last-added `LoadAppliedField` has `load_E=False`. So **append all B-only solenoid maps BEFORE the RF cavity maps** (which have `load_E=True`). Reuse the linac's **unconditional** `assert getattr(applied[-1], "load_E", False)` guard (there is always an RF map, so an `if any(...)` wrapper just adds noise — keep the unconditional form).

Because the injector domain now spans to 2.10 m, all three solenoids are in-domain (the old 1.30 m box would have left Sol 0 / Lens 0E inert — the gating dependency that motivated the domain extension).

### File-by-file changes

**`injector/build_injector_field.py`** — add solenoid conversion (reuse the linac builder's `to_grid`/`pad_r`/`write_series` pattern); compute each offset programmatically from the loaded peak:
```python
SOL_GUI_Z = {"LENS_0A": 0.225, "SOL_0": 1.897, "LENS_0E": 1.914}
# offset = GUI_z - Z[argmax|Bz|]  (derived per map, not hard-coded)
SOL_FILES = {"LENS_0A": ".../lens0a.h5", "SOL_0": ".../sol0.h5", "LENS_0E": ".../lens0e.h5"}
```
For each, load `R,Z,Br,Bz`, `to_grid`, optional `pad_r` to RMAX=36 mm (no-op but robust), compute `offset = GUI_z − Z[argmax|Bz|]`, write one single-mesh file with `grid_global_offset = [0.0, offset]`. Print a per-solenoid report (peak |Bz|·current in mT, **physical lab-z of peak**) for the sanity log, and assert `0 ≤ lab_peak ≤ 2.10`, `lab_peak < Z_HANDOFF`, and `abs(lab_peak − GUI_z) < tol`. (Currents live on the sim side for `config()` tunability.)

**`injector/injector_sim.py`** — apply static solenoid maps:
- Import `SOL_FILES` paths from the build module. Module constants (tunable via `config`): `I_LENS0A=6.0, I_SOL0=40.0, I_LENS0E=10.0`.
- Build the applied-field list **solenoids first, RF last**:
  ```python
  applied = []
  for path, cur in [(SOL_FILES["LENS_0A"], I_LENS0A),
                    (SOL_FILES["SOL_0"],   I_SOL0),
                    (SOL_FILES["LENS_0E"], I_LENS0E)]:
      if cur != 0.0:
          applied.append(picmi.LoadAppliedField(read_fields_from_path=path,
                          load_E=False, load_B=True, warpx_B_time_function=f"{cur:.8e}"))
  applied.append(make_cavity(PREB1_FIELD, ...))         # load_E=True
  if PREB2_KW > 0:
      applied.append(make_cavity(PREB2_FIELD, ...))     # load_E=True
  assert getattr(applied[-1], "load_E", False), "B-only field must not be last (E init-style)"
  for fld in applied:
      sim.add_applied_field(fld)
  ```

**`linac_sec1/build_linac_sec1_field.py`** — remove the in-linac solenoid build: drop `SOL_FILE`/`SOL_MAP`/`SOL_Z` and the `linac_sol.h5` write.

**`linac_sec1/linac_sec1_sim.py`** — drop `SOL_FIELD`, `I_SOL`, `SOL_Z`, and the `if I_SOL != 0.0: applied.append(...)` block. **Concrete import edit (hard breakage if missed):** the `from .build_linac_sec1_field import ...` line currently pulls `SOL_Z` (and `SOL_MAP`/`SOL_FILE`) — **remove `SOL_Z`/`SOL_MAP`/`SOL_FILE` from that import** or the module fails to compile the moment the build drops them. Repoint the bunch reader to `injector/diags/main/particles`.
- **Snapshot selector (reviewer-flagged, highest silent-failure risk):** rename `load_prebuncher_bunch → load_injector_bunch` **and change its selection logic.** The current selector maximizes in-bore charge (it warns it "optimizes the TRANSVERSE bore fit, not the LONGITUDINAL bunch" and "lands on the earliest/least-expanded post-gate snapshot"). With two cavities + focusing the bunch now forms a real longitudinal waist near 2.03 m, so max-q_bore would silently pick an early debunched snapshot and discard the bunching. **The new selector must target the handoff plane: pick the diagnostic dump whose bunch ⟨z⟩ is closest to `Z_HANDOFF = 2.03 m`, not min-σ_z and not max-q_bore.** This is a logic change, not a rename.
- **Handoff diagnostic resolution (resolves former Open Q#8, now a step-5 blocker):** time-spaced dumps land at irregular ⟨z⟩, so a nearest-⟨z⟩-to-2.03 selection may miss the plane. **The injector must place a diagnostic at/near 2.03 m** — preferred: add a fixed z-station diagnostic (a `BoundaryScrapingDiagnostic`/station plane at z=2.03 m if available in this pywarpx build) or a targeted dense `DIAG_PERIOD` window around 2.03 m so at least one dump has ⟨z⟩ within a few mm of the plane. Decide and implement this in step 4/5; the README must state which mechanism is used.
- Re-evaluate `RMAX` (currently 12 mm, sized for the old blown-up beam) once the focused beam exists (defer the exact value — likely toward the SLAC ~9.5 mm bore — to step 5). The RF maps stay; the `applied[-1].load_E` guard still holds.

## Task 4: Figures

### Headline: `pipeline/plot_chain.py` → `results/chain_evolution.png`

New in-process plotter (no pywarpx), exposed as `pipeline.plot_chain()` and called at the end of `run_pipeline.py:main()` after the four `*.run()` calls (next to `_beam_summary`). Output to a **new repo-root `results/`** — already matched by the existing `.gitignore` line `results/` (verified), so PNGs need `git add -f results/*.png` to commit.

**Input contract** (reads each stage's existing openPMD particle series; the injector path is a top-of-module constant so it tracks the contract):

| Stage | particles path |
|-------|----------------|
| cathode | `cathode/diags/particles` (2D x–z; no `y`/`uy`) |
| gun | `gun/diags/particles` |
| injector | `injector/diags/main/particles` |
| linac_sec1 | `linac_sec1/diags/main/particles` |

For each stage, iterate every dump, compute a row of beam moments, and place each segment at its **lab z**. **CORRECTION (software-designer): the openPMD `z` is already lab-frame, so segments abut naturally on each dump's own ⟨z⟩ — do NOT add a per-stage injection-z offset (that double-counts). Use `z0=0` everywhere** (keep it only as an inert escape-hatch param). Segments land at cathode ~0, gun 0–0.05, injector 0.06–2.1, linac to ~3.5 m. Draw stage dividers + name bands. **The plotter must read `injector/diags/main` (not stale `prebuncher/diags/P8_zc`, which may persist on disk).**

**Per-stage axis handling (reviewer-flagged):**
- The **cathode is 2D x–z** — it has no `y`/`uy`. The plotter must request only `x`/`ux` from the cathode series; never request `y`/`uy` (it raises). Transverse emittance there is the slab x-emittance.
- The **gun/injector/linac are RZ** — `get_particle` returns Cartesian x,y reconstructed from r, so ε_n,x is the projected emittance. This is only *approximately* comparable to the cathode's 2D-Cartesian x-emittance. **Annotate the cathode→gun seam on the cross-stage ε_n panel as a 2D→RZ definitional discontinuity** so the emittance-budget does not read a spurious jump there.

**Shared transverse emittance helper** `rms_emit(q, uq, w)` (openPMD `u = γβ`, already includes γ — no extra γ multiply; matches `gun/plot_gun.py`):
```
⟨q²⟩  = Σw q²/Σw − (Σw q/Σw)²
⟨uq²⟩ = Σw uq²/Σw − (Σw uq/Σw)²
⟨q·uq⟩= Σw q·uq/Σw − ⟨q⟩⟨uq⟩
ε_n   = sqrt(max(⟨q²⟩⟨uq²⟩ − ⟨q·uq⟩², 0)) * 1e6   # mm·mrad
```
Use this for ε_n,x; ideally the stage plotters import it too (removes duplicated inline math — a doc-sync win).

**Longitudinal emittance (reviewer-flagged units pitfall):** ε_n,z = z·(γβ_z) has units mm·(dimensionless), **not** mm·mrad. **Do not apply the ×1e6 mm·mrad scaling to the longitudinal plane.** Define ε_n,z explicitly (e.g. the z–(γβ_z) longitudinal emittance in mm, or a z–δ form) and label its axis accordingly; state the definition in the figure caption / README.

**Panels (3×2, shared x = ⟨z⟩):** (1) ⟨KE⟩(z) ±σ_KE band, **log-y** (0.06 keV → 146 keV → ~137 keV → ~15 MeV); (2) ε_n,x(z) with the cathode→gun discontinuity annotation; (3) σ_x(z) / σ_r; (4) σ_z(z) **log-y**; (5) charge fraction q(z)/q_emitted(z) using each stage's renorm + the linac `injection_summary.json` (annotate where renorm resets absolute charge); (6) I_peak(z). Per-stage colored segments + entry→exit annotations; mirror the numbers to stdout/log.

### Ship in the same change (all are views of the one moment table — marginal effort small)
- **#1 `emittance_budget.png`** — ε_n,x entry vs exit per stage (waterfall), showing *which* stage degrades quality; carry the cathode→gun 2D→RZ caveat as a footnote.
- **#2 `transmission_waterfall.png`** — cathode → gun exit → injector exit → enters bore → captured, each loss annotated; makes the scrape + capture losses one picture (central to the solenoid motivation).
- **#10 `chain_scorecard.png`/stdout** — per-stage entry/exit ⟨KE⟩, σ_KE, ε_n,x, ε_n,z, σ_x, σ_z, q, transmission; numeric backbone validated against each stage's own prints (cathode ε_n,x≈2.29 mm·mrad, gun ⟨KE⟩≈146 keV, linac capture %). **Capture denominator = true injected charge** (the linac already reports vs `q_inj`, not post-scrape); keep that.

### Second tranche (after the headline lands)
- **#3 `lps_montage.png`** (z−⟨z⟩ vs KE hexbin per stage exit) and **#4 `espec_waterfall.png`** (KE-histogram pcolormesh vs z) — both effort M, high value.

### Sweeps (after Tasks 2 & 3 exist)
- **#6 `preb_phase_scan.png`** — capture / σ_z,min / ε_n vs Preb-2 phase offset around the −45° design point (relative phasing is the key two-cavity knob).
- **#7 `sol_capture_scan.png`** — capture vs Lens0A / Sol0 / Lens0E current (quantifies the radial-blow-up fix).

### File plan
- New: `pipeline/plot_chain.py`, repo-root `results/`.
- Edit `pipeline/run_pipeline.py` (call `plot_chain.main()` in `main()`; update the closing "Figures:" banner). Expose `plot_chain` so `python -c "import pipeline; pipeline.plot_chain()"` works.
- No new deps (numpy/matplotlib/openpmd-viewer already pinned).

## Build / run / validate order

1. **Scaffold `injector/`** by copying `prebuncher/` (build, sim, plot, README, `__init__.py` facade); wire the `Stage` facade + `_launch_sim` path; reimplement `resolve_outdir()` to default to `injector/diags/main`. Keep Preb 1 working at z=0.534 m as baseline. Remove the old `prebuncher/` package once parity is confirmed.
2. **Extend domain** to ZMAX=2.10 m, NZ=1664 (NR=80); run **Preb 1 only** (8 kW). *Sanity (hard gate):* MLMG converges (no `MLMG failed`); **record wall-time here and decide whether knobs must loosen before adding more fields**; the z=0.534 m single-cavity result (V_gap≈58.6 kV, scale≈0.133, energy kick, σ_z) matches the old prebuncher.
3. **Add Preb 2** (mirrored `preb2_EB_rev.h5` at z=1.318 m, Q=4300, 10 kW). *Sanity (highest-risk):* on-axis Ez(z) of preb2 is the mirror of preb1 AND the **per-cavity transit-time energy-kick sign at the −45° design phase is correct** (validate against a single-cavity GPT kick — geometric mirror alone is insufficient); two-cavity bunching ratio improves over the drift-only baseline; include a per-cavity *phase* check (degrees at 214 MHz) for the constant-v Preb-2 timing.
4. **Add solenoids** (programmatic offsets: lens0a≈0, sol0≈+1.0761, lens0e≈0; B-only, current scale; added **before** RF). *Sanity:* per-solenoid **physical lab-z peak** report (assert in [0, 2.10] and < Z_HANDOFF, and within tol of GUI z); σ_r(z) shows focusing; beam stays inside RMAX through z=2.03 m. Implement the 2.03 m handoff diagnostic (z-station or dense DIAG window) here.
5. **Repoint `linac_sec1`** to `injector/diags/main`; remove its solenoid build/config and the `SOL_Z` import; switch the selector to **nearest-⟨z⟩-to-2.03 m**; re-evaluate RMAX. *Sanity:* injection scrape collapses; **measure and report capture fraction (vs true injected charge) — expect substantial improvement; do not assert a hard target number until measured**; ⟨KE⟩≈15 MeV.
6. **Add cross-stage figures** (`pipeline/plot_chain.py` + #1/#2/#10). *Sanity:* scorecard numbers match each stage's own prints; cathode→gun ε_n seam annotated.
7. **Sweeps** (#6, #7) once the above is stable.

## Doc-sync checklist (same change, per CLAUDE.md)

- **`CLAUDE.md`** — replace prebuncher rows in the inter-stage contract table with the `injector` row; update architecture prose (4-stage chain `cathode → gun → injector → linac_sec1`, drop "we model one prebuncher", drop the linac-solenoid-hack note); update field-map list (two prebuncher maps `preb1_EB`/`preb2_EB_rev` + `lens0a`/`sol0`/`lens0e`; linac no longer builds `linac_sol.h5`); update Commands / PERFORMANCE-KNOBS (rename `prebuncher.config`→`injector.config`, NZ=1664, **>2× cost / convergence-bound penalty**, linac cheaper); note repo-root `results/` in the git-ignore/commit Conventions.
- **Root `README.md`** — component table row prebuncher → injector; chain diagram; note the new repo-root `results/`.
- **`pipeline/README.md`** (reviewer-flagged GAP) — document the new cross-stage `pipeline/plot_chain.py`, the `pipeline.plot_chain()` facade export, the `run_pipeline.py` call, and the repo-root `results/`; update the inter-stage contract description (handoff plane moves to 2.03 m; `injector/diags/main`).
- **`injector/README.md`** (new, replacing `prebuncher/README.md`) — Task 1 verdict up front; two cavities + phi_off offsets (304.7°/178.9° → `phi_off` + arrival); reversed-install mirrored-map rationale (z-reverse + Er/Bφ sign-flip, not a time-function sign; symmetric-extent assumption); **Preb-2 phasing uses injection β, valid only sub-threshold, with the constant-v phase error quoted in degrees**; three solenoids with the **measured peak z's and corrected offsets** (LENS_0A 0.2333 m/offset 0 with the accepted 8 mm + coarse-axial caveat; SOL_0 0.8209 m local frame → +1.0761 m; LENS_0E 1.9147 m absolute/offset 0), programmatic offset derivation, z-overrun-tolerated-by-read_from_file note, 1-A normalisation, B-before-RF ordering gotcha; domain/grid (ZMAX=2.10, NZ=1664 justified by dz=1.26 mm at NR=80 → 2.80:1, ÷8 — *not* "same as linac"); handoff plane z≈2.03 m and the diagnostic mechanism used to land a snapshot there; γ²-overestimate self-field caveat over the longer ~2 m drift.
- **`linac_sec1/README.md`** — remove the in-stage solenoid sections; "reads `injector/diags/main`"; document the new nearest-⟨z⟩-to-2.03 m selector; revisit capture/RMAX now that focusing is upstream (no hard capture target until measured).
- **`FIGURES.md`** — replace prebuncher figure entries with injector figures (two-cavity bunching, envelope with lens/sol markers); add a new top section **"0. Cross-stage — `results/`"** documenting `chain_evolution.png`, `emittance_budget.png`, `transmission_waterfall.png`, scorecard, with the cathode→gun ε_n 2D→RZ caveat and the longitudinal-emittance-units note, regenerate line `python -c "import pipeline; pipeline.plot_chain()"`; update the "Regenerate everything" block (note repo-root `results/` is `git add -f`'d, matched by the existing `.gitignore` `results/` line — no new ignore pattern needed); update linac envelope/capture entries.
- **`requirements.txt`** — no new pinned dependency (same easygdf/openPMD/pywarpx/numpy/matplotlib stack); add only if a new analysis package is introduced.
- **Commit convention** — commit `injector/*.py` + `injector/README.md`, modified `linac_sec1/*`, `pipeline/*`, and doc files; `git add -f injector/results/*.png linac_sec1/results/*.png results/*.png`; never commit `diags/`, `.h5`, or logs. Delete the old `prebuncher/` package files in the same commit.

## Risks & open questions

1. **Two-cavity reversed install (highest risk).** A wrong Er/Bφ parity or wrong kick sign is a *silent* wrong-physics result. Mitigation: a scalar time-function sign cannot do it — must use the mirrored map (z-reverse + Er/Bφ negate, symmetric-extent assumption asserted). Validate **both** the Ez(z) mirror **and** the per-cavity energy-kick sign at the −45° phase against a single-cavity GPT kick before trusting two-cavity bunching. (Section 5's "flip Ez sign / π phase" suggestion is insufficient and overridden by Section 2.)
2. **Handoff-plane snapshot selection (high silent-failure risk).** The current max-q_bore selector would discard the bunching by landing on an early debunched dump. Mitigation: switch to nearest-⟨z⟩-to-2.03 m selection AND place a diagnostic at/near 2.03 m (z-station or dense DIAG window); decided in step 4/5.
3. **MLMG divergence / cost on the longer thin box (convergence-bound).** Mitigation: hold dz≈1.26 mm (NZ=1664) for 2.80:1 aspect; re-confirm convergence at looser Balanced precision over more cells; take the step-2 wall-time + `MLMG failed` checkpoint before adding more fields. Cost likely **exceeds** the naive 1.6–2.0×. Never coarsen NZ.
4. **Focusing may still be insufficient.** Sol 0 at 40 A is weak; LENS_0E (10 A, peak ~28.6 mT at z=1.914, just before the 2.03 m handoff) is the element most likely to control the envelope into the bore, with LENS_0A (6 A) setting the early envelope. The knobs are now physical (lens currents at true z); `sol_capture_scan.png` (#7) quantifies it. Watch σ_r(z) at 2.03 m against the SLAC ~9.5 mm bore.
5. **Cost roughly doubles (or more).** Accept and loosen the documented solve levers (`CFL`, `MAX_ITERS`, `REQUIRED_PRECISION`); flag in `run_pipeline.py` comments. The cheaper linac partly offsets.
6. **Lab-frame electrostatic self-field overestimates transverse SC by ~γ²** at β≈0.62 — now over a longer ~2 m drift, so the focusing solution is pessimistic and the capture number conservative. Note in README; reconsider whether the longitudinal-only ES approximation is acceptable over 2 m before trusting the absolute capture fraction.
7. **Cross-stage emittance discontinuity.** The cathode is 2D-Cartesian, downstream is RZ-projected; the ε_n panel/budget must annotate the cathode→gun seam as definitional, not physical.
8. **Open: handoff RMAX for the linac.** Re-evaluate (likely → SLAC ~9.5 mm bore) after step 5 produces the focused envelope; defer the exact value.
9. **Open: charge reconciliation (~1 nC vs 0.1 nC).** Per the reconciliation backlog (memory), absolute charge is still being reconciled; the chirp/threshold coefficients in Task 1 were measured at 0.1 nC. Transmission/scorecard denominators in Task 4 must stay consistent with whatever charge the chain actually injects; flag if reconciliation lands mid-upgrade.
10. **Open: hardened-Preb-1 power scans desync Preb-2 phasing.** The baked-at-setup Preb-2 time function uses injection β; a hard Preb-1 power scan would need a two-pass run (read post-Preb-1 β, rebuild Preb-2). Documented as a known limitation of the scan facility.
