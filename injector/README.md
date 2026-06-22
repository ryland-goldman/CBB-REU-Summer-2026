# WarpX CESR Injector (RZ)

Third stage of the Cornell Linac chain modelled in WarpX:

```
cathode (cathode/)  ->  gun (gun/)  ->  injector (this)  ->  linac_sec1 (linac_sec1/)
```

Driven through **lume-warpx**: every constant lives in `injector.yaml` (operating point + the
solenoid/prebuncher knobs in its `params:` block) and `injector_sim.py` reads them back,
imports the gun handoff via `WarpX(initial_particles=...)`, and overrides only the
runtime-computed values (the per-field RF/solenoid time functions, step count, `dt`, diagnostic
period). Edit `injector.yaml` to retune (the `config()` knob API is bypassed for this stage).

The injector is the **full LinacSim injector subsection in one self-consistent RZ
space-charge run** (it replaced the earlier single-cavity `prebuncher/` stage):

```
Lens 0A  ->  Prebuncher 1  ->  Prebuncher 2 (reversed)  ->  Sol 0 / Lens 0E  ->  9.547 mm collimator
 6 A          8 kW                10 kW                       40 A / 10 A          iris @ 1.922 m
 @0.225 m     @0.534 m            @1.318 m                    @1.897 / 1.914 m
```

It reads the gun's exit beam (~149 keV, β ≈ 0.63, ~1.0 nC time-release / ~0.83 nC legacy
snapshot, already RZ), velocity-bunches
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
gun output from **`gun/diags/handoff`** when the gun ran in its default time-release mode (the
full reconstructed ~2 ns exit beam — see `gun/README.md` → *Beam source*), else from
`gun/diags/particles` (the legacy snapshot). Run the whole chain with
`pipeline/run_pipeline.py`.

> **Caveat (time-release handoff).** The prebuncher phases tuned below were established against
> the compact *snapshot* gun handoff. The time-release beam is longer and lower-density (the
> physical 2 ns grid pulse), so its operating point should be re-validated — the handoff wiring
> is in place, the phase re-tuning is a follow-up (part of the LinacSim input-reconciliation
> backlog). Force the legacy input with `gun.config(BEAM_RELEASE="snapshot")` before `gun.run()`.

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

> **Gotcha (descending-z GDF → Er sign).** `prebuncher_25D.gdf` stores its z column
> **descending** (+152.4 → −152.4 mm), unlike the solenoid maps (ascending). The flat GDF
> rows must be reversed to match the ascending `np.unique(Z)` axis *before* the `(nz,nr).T`
> transpose (`load_prebuncher_map`). Skipping the reversal z-flips the map relative to its
> own axis, which negates the **odd** `Er` (Ez/Bφ are even and unaffected) — a silent bug:
> on-axis longitudinal bunching and the σ_z/energy headline are unchanged, but every off-axis
> transverse RF force is wrong-signed in both prebunchers. The gap-parity asserts do **not**
> catch it (parity is invariant under the flip), so `main()` adds a raw-GDF orientation check
> comparing the stored `Er` against the raw GDF column at a fixed off-axis +z point.

### RF drive and the two-cavity phasing (zero-crossing, centroid-referenced)

Each cavity drives its 1-J map as a standing-wave TM mode (GPT's `Map25D_TM`):
`Er,Ez(t) = map·scale·cos(ωt+φ)`, `Bφ(t) = map·scale·sin(ωt+φ)`. The drive phase is

```
φ = -ω·t_gap + base + radians(phi_off_deg) + rev_phase
```

where `t_gap` is the arrival time of the **bunch centroid** (`z_ref = z_centroid`, not the
tail) at the gap. The faithful default is **`PHASE="zc"`** (`base = π/2`) with
`phi_off = 0`: this lands the **centroid on the RF zero-crossing**, so the net mean-energy
kick is zero and each cavity acts as a **pure velocity buncher** — bunching does not change the
mean energy (gun-exit 149 keV → 152 keV at the handoff, ≈ flat; the transient ±18 keV swings
on the `mean KE` line near each gap are the long bunch *straddling* the cavity and cancel out
as it clears). `base = π` (`"crest"`) is the legacy *net-accelerating* reference (149 → 157 keV
ramp) kept for comparison.

> **Phase-reference caveat (which is "the" LinacSim operating point).** The GPT master deck is
> internally inconsistent about the prebuncher phase, so "faithful" here means matching the
> deck **body**, not the GUI annotation: the executed deck body (`gpt_master.in:213-219`) runs
> `phi_off = −90°` — the zero-crossing, i.e. exactly this `zc` default — while the `@GUI`
> annotation block (`:22-29`) *saves* −70°/−45° from crest (net-accelerating, `cos = 0.34/0.71`).
> `input_files.md` warns the committed GUI values are defaults, "not Adam's tuned working point."
> So neither phase is unambiguous machine truth; this repo follows the deck body (zc), Cornell2
> follows the GUI (−70/−45). Confirm the intended working point with the LinacSim author before
> treating either as canonical.

To reproduce that GUI working point in one call, the facade exposes
**`injector.faithful_gpt_deck()`** — it `config()`s `PHASE="crest"`, `PREB1_PHI_OFF=-70`,
`PREB2_PHI_OFF=-45`, and `PREB2_REV_PHASE=0` (the reversal is absorbed into the crest
reference). It does **not** change the shipped `zc` default, and its non-zero Preb-1 net kick
desyncs the analytic Preb-2 inter-cavity timing (see *Reversed install* below) — a hardened
study still needs the two-pass timing fix. Call it before `injector.run()`.

> **Why centroid-referenced.** The physical beam is the ~2 ns / ~380 mm time-release stream,
> which spans **~154° of RF** (the centroid sits **77° behind the tail**). Phasing the *tail*
> to the zero-crossing would leave the bulk on the decelerating slope (measured: −53 keV net),
> so the reference plane is the centroid. The cavity is energy-neutral by the symmetry of
> −cos(φ) about the zero-crossing over ±77°; the head (φ ≈ 13°) is decelerated and the tail
> (φ ≈ 167°) accelerated ⇒ compressive chirp.

- **Preb-1** (8 kW, Q=3000): centroid on the zero-crossing, full bunching slope (sin ≈ 1).
- **Preb-2** (10 kW, Q=4300, reversed): the second velocity buncher, also centroid-on-zero-crossing
  (`PREB2_REV_PHASE = π`, see below).

**Bunching.** σ_z monotonically tightens from ~108 mm (injection) to its ~33 mm **waist, which
now lands at the 2.03 m handoff plane** (min σ_z at ≈1962 mm) — the exit is the focus, not past
it. The waist is fold-limited (~33 mm): a single zero-crossing imparts a *sinusoidal* (not
linear) chirp across the 154°-wide bunch, so the longitudinal phase space folds into the
checkmark seen in `injector_phasespace.png` rather than collapsing to a point.

### Reversed install (`PREB2_REV_PHASE = π`) — the subtle part

GPT installs Preb 2 with `-1,0,0` (a 180° rotation). For this map's **definite parity** —
Ez EVEN, Er ODD, Bφ EVEN about the gap (measured corr ±0.9999; also forced by Maxwell for a
TM0 mode, Bφ ~ dEr/dz − dEz/dr with both even ⇒ Bφ even, asserted in the build) — the 180°
rotation flips ALL THREE lab components, i.e. a global E,B sign flip **≡ +π in ABSOLUTE
drive phase**.

In the **zc + `phi_off = 0`** parametrization, `phi_off` carries **no** reversal information
(unlike the old crest+GUI convention, where the GUI's 178.9° Preb-2 reference was measured for
the already-reversed cavity and so absorbed the +π). So Preb-2 must carry the genuine
geometric +π itself: **`PREB2_REV_PHASE = π`**.

Both `rev = 0` and `rev = π` are **energy-flat** (the centroid sits on a zero-crossing either
way); they differ only in the chirp **slope**, so the distinction is *which* of the two
bunching outcomes you get, not bunching-vs-debunching. A decisive SC-off A/B (same input)
measures it (σ_z [mm], injection → post-Preb-1 → post-Preb-2 → waist → 2.03 m handoff):

| `PREB2_REV_PHASE` | post-Preb-2 (1.45 m) | waist | σ_z at handoff |
|-------------------|----------------------|-------|----------------|
| **π** (faithful)  | 36.7 | **31.2 mm @ 2013 mm** | **31.3 mm** |
| 0 (forward)       | 29.9 | 20.3 mm @ 1644 mm | 48.6 mm |

So `rev = 0` actually bunches **harder** (the more-aggressive of the two energy-flat
zero-crossings) — it **over-compresses** to an earlier ~20 mm waist at ~1.64 m that then
**re-expands** to ~49 mm by the handoff. The faithful `rev = π` is the *less-aggressive*
slope: σ_z tightens monotonically through Preb-2 to its ~31 mm waist landing **at** the 2.03 m
handoff. `rev = π` is therefore both the geometrically-correct reversed install **and** the
operationally-correct default — do **not** "fix" it to `rev = 0`. **Note:** the value is
convention-dependent — under the old crest+GUI convention the reversal was absorbed into the
crest reference, so `PREB2_REV_PHASE = 0` there; if you switch `PHASE` back to `"crest"` with
the GUI `phi_off`, restore `rev = 0`.

### Preb-2 timing caveat (constant-v phase error)

Preb-2's arrival is timed in two segments — `v_beam` to Z1, then the post-Preb-1 speed over
Z1→Z2 — using an analytic estimate of Preb-1's net kick (`-cos(base+phi_off)·V_gap`). At the
**zc/centroid** default Preb-1's net centroid kick is ≈0 (energy-flat), so `v_after_preb1 ≈
v_beam` and the inter-cavity timing residual is small (the logged Δφ vs bare-injection-β timing
is set by the two-segment estimate; the σ-ratio bunching gate passes). **The exact fix is a
two-pass run** (read the post-Preb-1 β from a diagnostic, rebuild the Preb-2 timing); needed
only if a future study gives Preb-1 a non-zero net kick (a Preb-1 power scan, or `PHASE="crest"`
which net-accelerates ~+20 keV) — the scan facility carries this note.

## Solenoid lenses (the transverse focusing / radial-scrape fix)

All six LinacSim injector solenoids are now built and wired (in z-order). LENS_0A / SOL_0 /
LENS_0E carry the GUI default currents (6 / 40 / 10 A); **LENS_0B / 0C / 0D default to 0 A**
— the faithful LinacSim GUI value, so they are inert at the default operating point (a 0-A
lens is skipped) — but are real magnets, built and `config()`-overridable (`I_LENS0B=…`) for
matching/transport studies. Sol 1A-C sit downstream of the 2.03 m handoff (Section 1) and
remain omitted. Each lens is a separate per-Ampere B-only openPMD map (the grids differ —
LENS_0A is nr=189/nz=16, the others nr=16/nz~601 — so they cannot be combined), placed in the
lab frame via `grid_global_offset`. The 1-A maps scale linearly with current.

| Map | Current | GUI annotation | native peak z | lab-z offset (= z[0]) | lab-z span | peak \|Bz\| |
|-----|---------|----------------|---------------|-----------------------|------------|-------------|
| LENS_0A | 6 A | 0.225 m | 0.2333 m | 0.0 m | [0.00, 0.50] | 4.03 mT/A |
| LENS_0B | 0 A (default) | 1.603 m | 1.6108 m | +0.8 m | [0.80, 2.40] | 0.40 mT/A |
| LENS_0C | 0 A (default) | 1.692 m | 1.7014 m | +0.8 m | [0.80, 2.40] | 0.32 mT/A |
| LENS_0D | 0 A (default) | 1.838 m | 1.8241 m | +0.8 m | [0.80, 2.40] | 0.34 mT/A |
| SOL_0   | 40 A | 1.897 m | 0.8129 m | 0.0 m | [0.00, 2.50] | 0.15 mT/A |
| LENS_0E | 10 A | 1.914 m | 1.9148 m | +0.8 m | [0.80, 2.40] | 2.21 mT/A |

> **Note (SOL_0 is a long channel solenoid).** SOL_0's strong-field region (FWHM ≈ [0.35,
> 1.87] m) spans the **entire** prebuncher section, focusing PB1 (0.534 m) and PB2 (1.318 m).
> LENS_0B/0C/0D cluster near Lens 0E around 1.9 m. So the early line is **not** unfocused — the
> SOL_0 channel covers it; the Lens 0E telescope adds the final squeeze into the iris.

- **Placement is NATIVE absolute machine-z**, matching `gpt_master.in`, which installs every
  solenoid with `Map2D_B("wcs", "z", 0.0, …)` — the GDF's stored Z column *is* absolute machine
  z, so `grid_global_offset` is simply the native origin `z[0]` (0 for LENS_0A/SOL_0, 0.8 m for
  the pre-shifted LENS_0B…0E grids). The thin lenses' native peaks already match their GUI
  annotation (LENS_0A 0.233≈0.225, LENS_0E 1.915≈1.914), confirming native z is absolute.
  **Do NOT align argmax→GUI z:** that is correct only for the narrow lenses (argmax ≈ GUI) but
  **wrong for the flat-top SOL_0** — its argmax (0.813 m) is an arbitrary point on the [0.35,
  1.87] m plateau, and GUI "1.897" is a center/edge label, not the peak. Forcing argmax→1.897
  shifts the whole channel **+1.08 m** (to [1.08, 3.58], mostly past the 2.03 m handoff) and
  leaves PB1/PB2 unfocused — the placement bug fixed here (see *Capture / handoff result*).
- **LENS_0A 8 mm placement (note):** native peak 0.2333 m vs GUI 0.225 m differ by 8 mm —
  below the map's own ~31 mm axial cell, so neither is "more accurate." We ship the
  GUI-position (programmatic) placement for consistency across all three. Capture is
  **tune-sensitive** to the upstream lens placement/currents — treat the default as one operating
  point, not a precise capture optimum, and use the optional current scan to
  characterize it. *(An earlier ~7× LENS_0A sensitivity figure — 0.21% vs 1.6% — was measured
  before the LENS_0E grid_global_offset bug was fixed and is superseded; with all three lenses
  correctly placed, and with the multi-plane iris scrape, the default captures ~10%.)*
- **Ordering gotcha:** picmi forces the global `E_ext_particle_init_style` to "none" if the
  last-added `LoadAppliedField` has `load_E=False`, so the B-only solenoids are added **before**
  the RF cavities; an unconditional `assert applied[-1].load_E` guards it (a pure-drift baseline
  with no RF field legitimately skips the guard).
- The build asserts each solenoid's in-domain lab-z peak is in [0, 2.10 m] and **upstream of the
  2.03 m handoff** (so the linac never inherits a beam still inside a lens). A loose
  peak-within-half-FWHM-of-GUI-z cross-check runs **only for the narrow lenses** (the flat-top
  SOL_0 is exempt — its argmax is not its GUI annotation) to catch an 800 mm-class placement bug.

## The 9.547 mm collimator

LinacSim's prebuncher subsection carries a `scatteriris` of radius **9.547 mm at z = 1.922 m**
followed by a 9.547 mm pipe to 2.1 m — the injector→linac aperture. Past 1.922 m the
restriction is the SLAC ~9.55 mm bore, and Sol 0 / Lens 0E peak just upstream of it precisely
to squeeze the beam through. So the faithful success metric is **transmission through the
9.547 mm iris**, NOT "contained within the 36 mm domain."

It is applied as a **multi-plane id scrape** (`pipeline/collimator.py`), not an in-run scrape:
this pywarpx RZ build's particle-position SoA accessors raise *"Component x does not exist"*
(the radial position is the AMReX particle position, not a named real component), so an
afterstep weight-zeroing callback is not available here. **A single radial cut at the 2.03 m
handoff would be wrong**: the Sol 0 / Lens 0E telescope focuses the beam *hard* across the
1.922→2.03 m tail — it **converges**, not diverges (measured at the faithful 6/40/10 A tune:
in-iris ~38 % @1.92 m → ~93 % @2.03 m; σ_r 12.4 → 4.9 mm), so a particle outside the iris at
1.922 m — scraped by the real machine — can converge back inside it by 2.03 m and be wrongly
kept. Instead we emulate the continuous 9.547 mm pipe by tracking particle `id` across every
dump from z = 1.922 m on: a particle outside the aperture at *any* plane in the pipe hit the
wall and is removed; only the survivors are injected into the linac. (Approximations left: the
scraped-halo self-field over the tail — small, late, β≈0.7 — and the finite dump spacing
between planes.) This is exact in the dense-dump limit and reduces to the entrance-plane cut
when the envelope happens to be monotone.

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
- **Handoff diagnostic:** the dump cadence (`period`) is sized from the post-Preb-2 speed —
  the speed *at* the handoff plane, downstream of both cavities — so the spacing near 2.03 m is
  ≤8 mm, landing a snapshot within ~1 mm of the plane (picmi exposes only a uniform `period`;
  a true z-station / multi-interval diagnostic isn't available in this build — two same-name
  diagnostics trip "Diagnostic attributes not consistent" and `warpx_intervals` is rejected).
- **`SPACE_CHARGE`** (`config()`-overridable, default `True`): `False` passes
  `warpx_do_not_deposit` (beam self-field off, only the applied prebuncher/solenoid maps act). This
  is a diagnostic-only mode — the velocity bunching, two-cavity phasing, and transverse focusing are
  coupled *through* the self-field, so the ~1 mm bunch is set with SC on; the faithful default is
  `True`.

## Self-field solver (relativistic electromagnetostatic)

The beam self-field uses WarpX's **electromagnetostatic** solver (`warpx_magnetostatic=True`), not
the plain lab-frame electrostatic Poisson solve. In addition to `∇²φ = -ρ/ε₀` it solves the
Coulomb-gauge vector potential from the beam current (`∇²A = -μ₀ j`, `B = ∇×A`), so the self
magnetic field is included and the relativistic magnetic-pinch term `qβ×B` partially cancels the
radial space-charge repulsion: the net transverse self-force is `qE_r/γ²` rather than the
pure-electrostatic `qE_r`. This removes the ≈γ² (≈1.6–1.7× at β≈0.6–0.65 — the beam mean γ stays
≈1.29–1.31 across the line) transverse-SC over-repulsion the lab-frame solver incurs — the
WarpX–GPT 150 kV-gun benchmark's *cause 4*, which that writeup flags **grows for a longer line**,
and the injector is the longest line in the chain (~2 m at γ≈1.3 *throughout*, unlike the
low-γ-weighted gun where it was only +2%). Matches the
gun's solver. Measured A/B on the same time-release input: iris transmission **34% (lab-frame) →
42% (EMS)**; cost ≈3.7× (the extra A_z vector-Poisson over the long-thin box, ~78 → ~287 s).

> **Gotcha (same as the gun):** the outer radial wall BC is `dirichlet`, not `neumann`. The
> magnetostatic vector-Poisson's dominant A_z component (driven by `j_z`) has an all-Neumann,
> *singular* operator — the MLMG bottom solve then diverges — unless the outer wall is grounded.
> At RMAX=36 mm (well outside the beam) the self-field has decayed, so the beam dynamics are
> unaffected; this only makes A_z well-posed.

The A_z solve reuses φ's knobs (`REQUIRED_PRECISION=1e-4`, `MAX_ITERS=500`) and is muted
(`warpx_magnetostatic_verbosity=0`). A_z is the harder-conditioned solve, so it could in principle
hit the iteration cap and proceed *under-converged* (under-counting the pinch) silently — but a
verbosity≥2 spot-check confirms it converges in **3–4 V-cycles** (resid/bnorm ~1e-5, far below the
500 cap) at the default tune. If you tighten `CFL` or relax the knobs, re-check with
`warpx_magnetostatic_verbosity=2` that A_z still converges rather than capping.

## Capture / handoff result (the headline, with caveats)

At the faithful currents (6/40/10 A) the lenses focus the beam through the injector: the SOL_0
channel (FWHM ≈ [0.35, 1.87] m) confines the whole prebuncher section, Lens 0A (z ≈ 0.225 m)
sets the early envelope, and the Lens 0E telescope at z ≈ 1.9 m — just upstream of the 1.922 m
iris — squeezes it through the 9.547 mm aperture. On the time-release gun beam with the
relativistic EMS self-field (see *Self-field solver*) and the **zero-crossing (energy-flat)
cavity phasing**, **~64% of the in-domain charge passes the iris** (0.599 / 0.933 nC, via the
multi-plane scrape at the real 1.922 m iris plane — see *The 9.547 mm collimator*).

> **Fixed (physics-review, 2026-06-19): SOL_0 native-placement.** The dominant focusing element,
> the long channel solenoid SOL_0, was previously mis-placed by **+1.08 m**: the build forced its
> on-axis `argmax(Bz)` (0.813 m, an arbitrary point on its flat [0.35, 1.87] m plateau) to land
> at the GUI annotation 1.897 m, shifting the channel to [1.08, 3.58] m — mostly **past** the
> 2.03 m handoff, leaving PB1/PB2 and the early line **unfocused**. `gpt_master.in` (the
> authority) installs every solenoid with `Map2D_B("wcs", "z", 0.0, …)` — native absolute z, no
> peak-alignment shift. Placing SOL_0 at its native z restores the channel over the prebuncher
> section and raises iris transmission from **~19% → ~64%** (the thin lenses, whose argmax ≈ GUI,
> were unaffected — ≤ 14 mm shift). See *Solenoid lenses → Placement*.

> **⚠ Transverse match — largely resolved by the SOL_0 fix.** With SOL_0 correctly placed,
> iris transmission is **~64%** on the energy-flat ~150 keV zc beam — *above* the old crest value
> (~42%). The earlier ~19% figure was an artifact of the misplaced SOL_0, not (as previously
> framed) the energy-flat beam mismatching the 6/40/10 A telescope. The longitudinal operating
> point is also correct: zc/centroid phasing keeps the mean energy flat and lands the σ_z waist
> at the 2.03 m handoff (see *RF drive*). A clean re-optimization of the lens currents to the
> ~150 keV beam remains an optional **transverse** follow-up (the LinacSim reconciliation
> backlog), but is no longer load-bearing for a sane transmission. Note `linac_sec1`'s
> `n ≥ 0.8·nmax` dump selector may still fall back off-plane on a low-population near-handoff
> dump; the downstream RF-bucket capture should be re-measured at this operating point.

> **Earlier fixes (physics-review):** (1) LENS_0E was once placed ~800 mm out of position (a
> `grid_global_offset` bug that omitted the native grid origin `z[0]`); the fix (carry `z[0]` +
> a read-back assertion on the *stored* peak) restored the matching lens — the SOL_0 fix above
> generalizes that to all solenoids (native z, no argmax→GUI shift). (2) The iris scrape was once
> a single radial cut at the 2.03 m handoff, but the beam **converges** across the 1.922→2.03 m
> tail, so that cut kept converged halo the real 1.922 m iris scrapes and overstated transmission
> ~3×. The multi-plane scrape at the true iris plane is the fix (see *The 9.547 mm collimator*).

Three caveats frame this number:
1. **Self-field solver (resolved):** the injector now uses the relativistic electromagnetostatic
   solver (`warpx_magnetostatic`), which includes the 1/γ² magnetic-pinch cancellation — so the
   earlier lab-frame "~γ² conservative lower bound" no longer applies. See *Self-field solver*.
2. **Tune-sensitive:** capture responds strongly to the upstream lens currents/placement; the
   faithful 6/40/10 A currents are not tuned to a capture optimum, so the optional current scan
   is the right tool to map the achievable capture (an earlier ~7×-from-8 mm-LENS_0A figure was
   a pre-fix artifact of the LENS_0E mislocation — superseded).
3. **Charge recovery is the real win:** the solenoids recover in-domain charge from ~0.04 nC
   (no focusing) to ~0.933 nC — the radial-scrape fix works; the iris then sets the true
   transmission. The optional current/phase scans characterize the achievable capture.
4. **Handoff-seam residual field (bounded approximation):** Lens 0E peaks just upstream of the
   2.03 m handoff and its field *tail* is still substantial AT the plane (≈ 15 % of peak,
   ~3.4 mT on-axis at 10 A); the native-placed Sol 0 peaks far upstream (0.813 m) so its tail at
   the plane is small (≈ 8 % of peak, ~0.5 mT at 40 A). The linac stage models no
   solenoid, so this continuing transverse focus is dropped at the seam — an unphysical
   discontinuity. The build only asserts each lens *peak* is upstream of the handoff (not that the
   field has decayed there), so the beam is handed off while still being focused. The dropped focus
   makes capture **more** conservative (the real continuous field would hold the envelope tighter
   into the structure). A fully faithful treatment would carry the Sol 0 / Lens 0E tails into the
   linac stage as applied fields, or move the handoff downstream of the tails — a documented
   follow-up, not done here.

## Outputs

`injector.plot()` reads `injector/diags/main/particles` and writes to `injector/results/`,
generating three layers — lume-warpx's plotting helpers (particle diagnostics only — no field
diagnostic is dumped, so no `plot_fields`), the shared custom figures in `pipeline/plot_extras.py`,
and the stage-specific rich figures (raw openPMD over the whole run; the prebuncher field lobes
read from the `preb1/2_EB.h5` maps):

- `injector_phase_space_z_KE.png` — `plot2D("z","kinetic_energy")`: the energy-flat (~150 keV)
  velocity-bunched exit beam at the 2.03 m handoff.
- `injector_transverse_x_px.png` — `plot2D("x","px")`: the solenoid-focused exit transverse phase space.
- `injector_centroid_vs_t.png` — `plot1D("t","mean_z")`: the bunch traversing the ~2 m line.
- `injector_bunch_length_vs_t.png` — `plot1D("t","sigma_z")`: σ_z compressing to its waist at the handoff.
- `injector_emittance_vs_t.png` — `plot1D("t","norm_emit_x")`: transverse emittance over the run.
- `injector_beamsize_vs_t.png` — `plot1D("t","sigma_x")`: the transverse envelope σ_x — the Sol 0 / Lens 0E focusing.
- `injector_energy_spectrum.png` — `plot_extras.energy_spectrum`: charge-weighted handoff KE histogram — the energy-flat (~150 keV, zero net kick) spectrum.
- `injector_current_profile.png` — `plot_extras.current_profile`: longitudinal current I(z) — the velocity-bunched current peak at the waist.
- `injector_energy_chirp.png` — `plot_extras.energy_chirp`: slice ⟨KE⟩ vs z (with density) — the longitudinal energy chirp that drives the velocity bunching.
- `injector_beam_spot_xy.png` — `plot_extras.beam_spot`: the transverse x–y spot (RZ-reconstructed y) with marginals.
- `injector_cavity.png` — the on-axis E_z field lobes of both prebunchers (scaled by their drive amplitude V_g), in lab z, with the 2.03 m handoff plane marked.
- `injector_line.png` — σ_z, peak current, and mean KE vs ⟨z⟩ along the line: the velocity-bunching story (σ_z 108→27 mm, peak current → ~16 A) with the energy modulation at each cavity gap.
- `injector_bunch_profile.png` — the longitudinal line-charge density λ(z) at four stations (injection / after Preb 2 / best focus / handoff).
- `injector_phasespace.png` — the charge-weighted longitudinal phase space (z−⟨z⟩, KE−⟨KE⟩) at the same four stations: the cavity chirp rotating to the bunched S-curve at the waist.

## Notes / caveats

- The self-field is the relativistic electromagnetostatic solver (includes the 1/γ² magnetic
  pinch — see *Self-field solver*), so the earlier lab-frame "~γ² pessimistic lower bound" framing
  no longer applies.
- Preb-2 phasing uses the injection β (+ analytic Preb-1 kick): valid only while both cavities
  are sub-threshold (the design point); a hardened-Preb-1 scan needs a two-pass run.
- The 9.547 mm collimation is a post-hoc **multi-plane id scrape** (the in-run scrape isn't
  available in this pywarpx RZ build) — applied at the real 1.922 m iris plane, because the beam
  *converges* through the 1.922→2.03 m tail, so a single 2.03 m cut would overstate transmission.
- **openPMD fd-leak gotcha:** openpmd-viewer leaks one file descriptor per `get_particle()`, so
  looping over this stage's ~280 diagnostic dumps exhausts macOS's default 256-fd soft limit and
  the read fails with `IO Task OPEN_FILE failed … Inaccessible` (HDF5's report of EMFILE) at a
  fixed dump *count*, not a specific file — the diag files are intact. The pipeline raises
  `RLIMIT_NOFILE` (`pipeline/_runner._raise_fd_limit`, called from `_prepare_environment`, plus
  `_launch_sim.py`, `plot_chain.main`, and `plot_injector.main`) to fix it; the `_retry_io`
  backoff in the sim/plot scripts is only a backstop for a transient open and does **not** rescue
  fd exhaustion. Any new code that loops reads over a long diag series needs the raised limit.
