# WarpX CESR Injector (RZ)

Third stage of the Cornell Linac chain modelled in WarpX:

```
cathode (cathode/)  ->  gun (gun/)  ->  injector (this)  ->  linac_sec1 (linac_sec1/)
```

The injector is the **full LinacSim injector subsection in one self-consistent RZ
space-charge run** (it replaced the earlier single-cavity `prebuncher/` stage):

```
Lens 0A  ->  Prebuncher 1  ->  Prebuncher 2 (reversed)  ->  Sol 0 / Lens 0E  ->  9.547 mm collimator
 6 A          8 kW                10 kW                       40 A / 10 A          iris @ 1.922 m
 @0.225 m     @0.534 m            @1.318 m                    @1.897 / 1.914 m
```

It reads the gun's exit beam (~146 keV, β ≈ 0.63, ~0.83 nC, already RZ), velocity-bunches
it with two 214 MHz prebuncher cavities while focusing it with three static solenoid
lenses, and hands a focused, collimated beam to `linac_sec1` at the true linac entrance
**z ≈ 2.03 m** (Z_acc_1). Modelling all elements in ONE drift is essential: the bunching,
two-cavity phasing, and transverse focusing are coupled through the self-field.

## Running

```python
# from repo root, in the CBB env:
import injector
injector.run()        # build fields + sim + plots  (writes injector/diags/main)
# injector.plot()     # re-generate figures from existing diags/
```

`injector.run()` runs the faithful default operating point and writes `injector/diags/main`.
The build reads `fieldmaps/{prebuncher_25D,SOL_0,LENS_0A,LENS_0E}.gdf`; the sim reads the
gun output from `gun/diags/particles`. Run the whole chain with `pipeline/run_pipeline.py`.

## Task 1 — Prebuncher power: 8 kW is faithful, and intentionally weak

**8 kW is the faithful LinacSim default** (`prebuncher1_input_power`) and is intentionally
weak — single-cavity bunching is NOT the design. The injector is a **two-prebuncher +
solenoid distributed buncher**; each cavity sits ~12× below the ~95 kW single-cavity
bunching threshold in power (~4× in voltage). At 8 kW the gap voltage is V_gap ≈ 58.6 kV
(scale ≈ 0.133 from `scale = sqrt(1e3·Q·P/(2π f_RF))`, V1J = 438.6 kV, f_RF = 214.18 MHz,
Q = 3000). Do not raise it to "get bunching" — that misreads the architecture. (The old
160–800 kW single-cavity scan was exploration; the prior −3.05 keV/mm cavity / +1.40 keV/mm
gun chirp coefficients were measured at 0.1 nC, and at the reconciled ~1 nC charge space
charge is stronger and the threshold higher, so 8 kW is even further below threshold.)

## Field maps

`build_injector_field.py` writes five openPMD maps from two `.gdf` sources:

- `preb1_EB.h5` — Prebuncher 1, forward 1-J cavity field, gap at `Z_GAP_CENTER_1 = 0.534 m`.
- `preb2_EB.h5` — the **same forward field**, gap at `Z_GAP_CENTER_2 = 1.318 m`. The two
  cavities differ only in lab-z placement (`grid_global_offset`) and the run-time reversal
  phase; there is **no mirrored/negated map**.
- `lens0a.h5` / `sol0.h5` / `lens0e.h5` — the three static, per-Ampere B-only solenoid maps.

### RF drive and the two-cavity phasing

Each cavity drives its 1-J map as a standing-wave TM mode (GPT's `Map25D_TM`):
`Er,Ez(t) = map·scale·cos(ωt+φ)`, `Bφ(t) = map·scale·sin(ωt+φ)`. The drive phase is

```
φ = -ω·t_gap + base + radians(phi_off_deg) + rev_phase
```

- `base = π` ("crest") is the faithful reference: the GUI `phi_off` values are
  **crest-referenced**, so the operating point is crest + `phi_off` (Preb-1 −70°, Preb-2 −45°,
  reproducing the GUI's 304.7° / 178.9° on-crest definitions). `base = π/2` ("zc") is the
  bare zero-crossing kept only for the exploratory scan (use with `phi_off=0`).
- **Preb-1** (8 kW, Q=3000, −70° from crest): a mild **+20 keV** net kick with a strong
  bunching slope (on-axis kick fraction −cos(π − 70°) ≈ +0.34; slope sin ≈ +0.94 ⇒ the tail
  gains energy and the bunch compresses downstream).
- **Preb-2** (10 kW, Q=4300, reversed, −45° from crest): the second velocity buncher. Its
  mean kick is ~+43.5 keV; it drives the head→tail chirp compressive (dchirp ≈ −0.33 keV/mm
  in the Preb-2-only test).

**Two-cavity bunching (vs the P=0 drift baseline AND vs Preb-1 alone):** σ_drift/σ_2cav ≈
**4.4×** near the focus and ≈ **2.4× vs Preb-1 alone** — Preb-2 genuinely *adds* bunching.
σ_z monotonically tightens drift ~80 mm → Preb-1 ~37 mm → two-cavity ~12–18 mm.

### Reversed install (`PREB2_REV_PHASE = 0`) — the subtle part

GPT installs Preb 2 with `-1,0,0` (a 180° rotation). For this map's **definite parity** —
Ez EVEN, Er ODD, Bφ EVEN about the gap (measured corr ±0.9999; also forced by Maxwell for a
TM0 mode, Bφ ~ dEr/dz − dEz/dr with both even ⇒ Bφ even, asserted in the build) — the 180°
rotation flips ALL THREE lab components, i.e. a global E,B sign flip **≡ +π in ABSOLUTE
drive phase**.

BUT we do **not** phase in absolute phase: we phase to bunch arrival + a CREST base + the GUI
`phi_off`, and the GUI's 178.9° on-crest reference for Preb-2 was measured for the
**as-installed (already-reversed)** cavity. So the crest reference **already contains the
reversal**; adding a separate +π double-counts and lands Preb-2 on the debunching slope. In
the (forward-map + crest-base + GUI φ_off) parametrization, **`PREB2_REV_PHASE = 0` IS the
geometric reversal** (the +π and the crest-reference's built-in reversal cancel) — not its
absence. Verified empirically by a Preb-2-only kick-sign run: rev=0 bunches (dchirp
−0.33 keV/mm, tail gains); rev=+π decelerates (−67 keV, no bunching). **Do not "fix" this to
+π.** (The knob is retained for a future map whose loaded drive phase is not the as-installed
crest.)

### Preb-2 timing caveat (constant-v phase error)

Preb-2's arrival is timed in two segments — `v_beam` to Z1, then the post-Preb-1 speed over
Z1→Z2 — using an analytic estimate of Preb-1's +20 keV kick. The faithful crest-base Preb-1
raises β over the inter-cavity drift, so timing Preb-2 with the bare injection β would
mis-time it by ~−10° at 214 MHz; the two-segment estimate cuts this to a **~−4° residual**
(the analytic thin-gap kick slightly overshoots the space-charge-loaded ensemble mean, so it
does not reach <1°). This is below both the constant-v-per-segment approximation and the
γ²-scale self-field pessimism, and the σ-ratio bunching gate passes — accepted at the design
point, documented honestly as ~4° (NOT <1°). **The exact fix is a two-pass run** (read the
post-Preb-1 β from a diagnostic, rebuild the Preb-2 timing); needed only if a future study
changes Preb-1's mean kick (a Preb-1 power scan, or moving its phase toward crest) — the scan
facility carries this note.

## Solenoid lenses (the transverse focusing / radial-scrape fix)

Three energized LinacSim lenses (GUI defaults; 0B/0C/0D are 0 A and omitted). Each is a
separate per-Ampere B-only openPMD map (the grids differ — LENS_0A is nr=189/nz=16, the
others nr=16/nz~601 — so they cannot be combined), placed in the lab frame via
`grid_global_offset`. The 1-A maps scale linearly with current.

| Map | Current | GUI / lab-z peak | native peak z | programmatic offset |
|-----|---------|------------------|---------------|---------------------|
| LENS_0A | 6 A | 0.225 m | 0.2333 m | −0.0083 m |
| SOL_0   | 40 A | 1.897 m | 0.8129 m | +1.0841 m |
| LENS_0E | 10 A | 1.914 m | 1.9147 m | −0.0007 m |

- **Offsets are derived programmatically** per map: `offset = GUI_z − Z[argmax|Bz|]`, landing
  each peak dead-on at its GUI lab-z. This self-corrects against stale plan literals (the plan
  table's SOL_0 +1.0761 / native 0.8209 is stale — the actual file gives +1.0841 / 0.8129).
  Do NOT hard-code the offsets.
- **LENS_0A 8 mm placement (note):** native peak 0.2333 m vs GUI 0.225 m differ by 8 mm —
  below the map's own ~31 mm axial cell, so neither is "more accurate." We ship the
  GUI-position (programmatic) placement for consistency across all three. Capture is
  **tune-sensitive** to the upstream lens placement/currents — treat the default as a
  conservative (γ²) lower bound, not a precise value, and use the optional current scan to
  characterize it. *(An earlier ~7× LENS_0A sensitivity figure — 0.21% vs 1.6% — was measured
  before the LENS_0E grid_global_offset bug was fixed and is superseded; with all three lenses
  correctly placed the default captures ~18%.)*
- **Ordering gotcha:** picmi forces the global `E_ext_particle_init_style` to "none" if the
  last-added `LoadAppliedField` has `load_E=False`, so the B-only solenoids are added **before**
  the RF cavities; an unconditional `assert applied[-1].load_E` guards it (a pure-drift baseline
  with no RF field legitimately skips the guard).
- The build asserts each solenoid's in-domain lab-z peak is in [0, 2.10 m], **upstream of the
  2.03 m handoff** (so the linac never inherits a beam still inside a lens), and within tol of
  the GUI z.

## The 9.547 mm collimator

LinacSim's prebuncher subsection carries a `scatteriris` of radius **9.547 mm at z = 1.922 m**
followed by a 9.547 mm pipe to 2.1 m — the injector→linac aperture. Past 1.922 m the
restriction is the SLAC ~9.55 mm bore, and Sol 0 / Lens 0E peak just upstream of it precisely
to squeeze the beam through. So the faithful success metric is **transmission through the
9.547 mm iris**, NOT "contained within the 36 mm domain."

It is applied as a **radial cut at the handoff snapshot** (and the linac's `RMAX=9.547 mm`
injection cut), not an in-run scrape: this pywarpx RZ build's particle-position SoA accessors
raise *"Component x does not exist"* (the radial position is the AMReX particle position, not a
named real component), so an afterstep weight-zeroing callback is not available here. The cut
is physically equivalent to the continuous pipe: the pipe holds 9.547 mm from 1.922 m to 2.1 m
and the envelope grows monotonically over that 0.1 m tail, so a particle inside 9.547 mm at
2.03 m was inside the whole pipe and one outside hit the wall before 2.03 m. (The only
approximation is the scraped-halo self-field over that 0.1 m tail — small, late, β≈0.7.) **A
future current scan that focuses a waist INTO 1.922–2.03 m would break the monotonic-divergence
assumption and need a true multi-plane scrape.**

## Domain / grid

- **z:** `ZMIN=0`, `ZMAX=2.10 m` (LinacSim prebuncher-subsection ZSTOP), with a field-free
  exit drift so the handoff beam coasts. Handoff snapshot at **z ≈ 2.03 m** (Z_acc_1).
- **r:** `RMAX=0.036 m`, `NR=80` (dr=0.45 mm) — keep NR=80 (the RF map reaches 36 mm and needs
  the radial resolution; do NOT copy the linac's NR=16).
- **Cell aspect (binding):** `NZ=1664` gives dz=1.262 mm ⇒ dz/dr = 2.80:1 (the ≈3:1 rule) and
  is ÷8 (blocking factor). **Do not coarsen NZ** — this long-thin box is convergence-bound, so
  coarsening NZ slows the per-step MLMG solve faster than it removes cells AND under-resolves
  the ~1 mm bunch. Speed it via `CFL`, `MAX_ITERS`, `REQUIRED_PRECISION`. The coincidence with
  the linac's NZ=1664 is not the rationale (the linac reaches 2.8:1 via NR=16/dr=0.75 mm). The
  injector run is convergence-bound, so its cost over the 2.10 m box is >2× the old 1.30 m
  prebuncher (~60 s vs ~24 s).
- **Handoff diagnostic:** the dump cadence (`period`) is sized so the spacing near 2.03 m is
  ≤8 mm, landing a snapshot within ~1 mm of the plane (picmi exposes only a uniform `period`;
  a true z-station / multi-interval diagnostic isn't available in this build — two same-name
  diagnostics trip "Diagnostic attributes not consistent" and `warpx_intervals` is rejected).

## Capture / handoff result (the headline, with caveats)

At the faithful currents (6/40/10 A) the three lenses focus the beam through the injector:
Lens 0A (z ≈ 0.225 m) sets the early envelope, and the Sol 0 / Lens 0E matching telescope at
z ≈ 1.9 m — just upstream of the 1.922 m iris — squeezes it through the 9.547 mm aperture.
**~91% of the handoff charge passes the iris** (≈0.69 / 0.76 nC), and the downstream linac
captures **~18% of the true injected charge** into the RF bucket (⟨KE⟩ ≈ 26 MeV, σ_KE ≈ 8 MeV,
max ~32 MeV). The two prebunchers net-accelerate the beam **146 → ~220 keV** at the handoff
while velocity-bunching it. This is **faithful machine behavior** — real lenses at real z, real
iris, real handoff — and a large improvement over the old mislocated-solenoid hack.

> **Fixed (physics-review):** an earlier version placed LENS_0E ~800 mm out of position (a
> `grid_global_offset` bug that omitted the native grid origin `z[0]`, putting its peak at
> 1.114 m instead of 1.914 m — so no lens focused at the iris). That gave only ~8% iris
> transmission / ~1% capture. The fix (corrected offset + a read-back assertion that validates
> the *stored* peak, not an input recompute) restores the matching lens to the iris and yields
> the 91% / 18% above. The capture is **~18×** what the buggy build reported.

Three caveats frame this number:
1. **Conservative lower bound:** the lab-frame electrostatic self-field omits the 1/γ² magnetic
   pinch cancellation, overestimating transverse space charge by ~γ² (≈1.66× at β≈0.7) over the
   ~2 m drift, so the real machine focuses tighter and captures more.
2. **Tune-sensitive:** capture responds strongly to the upstream lens currents/placement; the
   faithful 6/40/10 A currents are not tuned to a capture optimum, so the optional current scan
   is the right tool to map the achievable capture (an earlier ~7×-from-8 mm-LENS_0A figure was
   a pre-fix artifact of the LENS_0E mislocation — superseded).
3. **Charge recovery is the real win:** the solenoids recover in-domain charge from ~0.04 nC
   (no focusing) to ~0.77 nC — the radial-scrape fix works; the iris then sets the true
   transmission. The optional current/phase scans characterize the achievable capture.
4. **Handoff-seam residual field (bounded approximation):** Sol 0 and Lens 0E peak just upstream
   of the 2.03 m handoff but their field *tails* are still substantial AT the plane (Sol 0 ≈ 97 %
   of peak, Lens 0E ≈ 15 % — together ~9 mT on-axis at 40/10 A). The linac stage models no
   solenoid, so this continuing transverse focus is dropped at the seam — an unphysical
   discontinuity. The build only asserts each lens *peak* is upstream of the handoff (not that the
   field has decayed there), so the beam is handed off while still being focused. The dropped focus
   makes capture **more** conservative (the real continuous field would hold the envelope tighter
   into the structure). A fully faithful treatment would carry the Sol 0 / Lens 0E tails into the
   linac stage as applied fields, or move the handoff downstream of the tails — a documented
   follow-up, not done here.

## Outputs

`injector.plot()` reads `injector/diags/main` (and any `injector/diags/P*` scan dirs) and
writes to `injector/results/`:

- `injector_line.png` — σ_z(z) (vs drift baseline) and peak current / mean energy, with both
  prebuncher-gap markers (Z1, Z2).
- `injector_phasespace.png` — z–KE at injection / cavity exit / best focus.
- `injector_cavity.png` — the RF drive: both on-axis Ez(z) lobes (Preb 1 @ 534 mm, Preb 2 @
  1318 mm) and both RF waveforms at their gap arrivals (scale/phase re-derived as the sim does).
- `injector_bunch_profile.png` — the longitudinal line-charge density λ(z).
- `compare_power_phase.png` — cross-case scan summary (when scan dirs are present).

## Notes / caveats

- The lab-frame electrostatic self-field is non-relativistic (omits 1/γ², overestimates
  transverse SC by ~γ²); the focusing solution is pessimistic and the capture a lower bound.
- Preb-2 phasing uses the injection β (+ analytic Preb-1 kick): valid only while both cavities
  are sub-threshold (the design point); a hardened-Preb-1 scan needs a two-pass run.
- The 9.547 mm collimation is a post-hoc handoff cut (the in-run scrape isn't available in this
  pywarpx RZ build) — physically equivalent for the monotonically-diverging design beam.
