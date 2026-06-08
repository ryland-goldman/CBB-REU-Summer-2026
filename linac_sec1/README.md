# WarpX SLAC Linac — Section 1 (RZ)

Fourth stage of the Cornell Linac chain modelled in WarpX, and the **first downstream
accelerating section** (later sections → `linac_sec2`, … as their field maps are added):

```
cathode (cathode/) -> gun (gun/) -> injector (injector/) -> linac_sec1 (this)
```

The **injector's** focused, velocity-bunched beam — read at the **z ≈ 2.03 m handoff plane**,
already collimated to the 9.547 mm iris — enters a 3 m, 86-cell, 2π/3 **traveling-wave** SLAC
accelerating structure with self-consistent space charge. **Transverse focusing is upstream now**
(the injector's three real lenses at their true lab z); this stage carries **no solenoid**. The
linac selects the injector dump whose ⟨z⟩ is nearest 2.03 m and applies the **multi-plane 9.547 mm
iris scrape at injection** (`pipeline/collimator.py`) — that scrape IS the physical
injector→linac iris collimation. At the faithful 11 MW point the captured charge is **~7 % of the
true injected charge** to **⟨KE⟩ ≈ 25 MeV** (max ~32 MeV, σ_KE ≈ 8 MeV). That capture is faithful
machine behavior — the Sol 0 / Lens 0E matching telescope focuses ~32 % of the handoff charge
through the 9.547 mm iris, and the linac then captures the in-bucket fraction — and is a
**conservative lower bound** (the lab-frame ES
self-field overestimates transverse space charge by ~γ²≈1.66×, so the real machine captures more),
tune-sensitive to the upstream lens currents (see `injector/README.md`). Capture is reported
against the **true injected charge** — see *Capture bookkeeping*.

## Running

```python
# from repo root, in the CBB env:
import linac_sec1
linac_sec1.run()         # build field + sim + plots -> diags/main, results/
# linac_sec1.run(plots=False)
# linac_sec1.plot()      # re-generate figures from existing diags/
```

`run()` runs **one case** at the default operating point — the original LinacSim values
(`PHASE_DEG=0`, `POWER_MW=11`). The operating point and grid are module-level constants at the top
of `linac_sec1/linac_sec1_sim.py` (`POWER_MW`, `PHASE_DEG`, `NZ`, …), overridable via `config()` —
e.g. a `PHASE_DEG` sweep in a Python loop for the acceptance curve. (There is no `I_SOL` — focusing
moved upstream to the injector.)

`build_linac_sec1_field` reads the two SLAC maps from `fieldmaps/`; the sim reads the injector
output from `injector/diags/main/particles/` (repo-root-relative), selecting the dump nearest the
z ≈ 2.03 m handoff. To run the whole chain (cathode → gun → injector → linac_sec1), use
**`pipeline/run_pipeline.py`**.

## Field maps

**This stage uses two GPT maps** (the in-linac solenoid was removed — focusing is upstream in
the injector). The two SLAC files are **not two sections** — they are the **real and imaginary
(quadrature) components of one** 3 m traveling-wave structure
(`reference/Linac Simulation Documentation/details.md`):

| file | columns (used) | grid | normalisation |
|------|----------------|------|---------------|
| `SLAC-3mLinac-field1.gdf` | `ErRe, EzRe, HphiIm` | 21 (r) × 6379 (z), r ≤ 9.55 mm, z ≈ 3.0 m | 1 kW input power |
| `SLAC-3mLinac-field2.gdf` | `ErIm, EzIm, HphiRe` | same grid | 1 kW input power |

`build_linac_sec1_field.py` writes two openPMD files (both `thetaMode` m = 0, components r/t/z,
shape `(1, nr, nz)`, the layout WarpX's `read_from_file` reader expects):

- `linac_rf1.h5` — `E = (ErRe, 0, EzRe)`, `B = (0, HphiIm, 0)` (the in-phase quadrature),
- `linac_rf2.h5` — `E = (ErIm, 0, EzIm)`, `B = (0, HphiRe, 0)` (the 90° quadrature).

The SLAC maps reach the **9.55 mm bore** in r; they are zero-padded in r out to the sim domain
`RMAX = 9.547 mm` (the SLAC bore / injector→linac collimator iris — *not* the old 12 mm). So the
radial domain IS the aperture: a particle outside it is scraped at injection, exactly as the real
iris does. Each map is placed in the lab frame via `grid_global_offset` (`Z_STRUCT = 0.10 m`
structure entrance).

## RF drive — synthesising the traveling wave

GPT builds the traveling wave as the sum of two standing waves 90° apart (details.md):

```
Map25D_TM(... field1, "ErRe","EzRe","HphiIm", scale, 0, φ,        2π·f);
Map25D_TM(... field2, "ErIm","EzIm","HphiRe", scale, 0, φ+0.5π,   2π·f);
```

i.e. each map is driven `E(t) = map·scale·cos(ωt+φ)`, `Bφ(t) = map·scale·sin(ωt+φ)`, with field2
offset by +π/2; the sum `Re[(Ẽ_re + iẼ_im)e^{i(ωt+φ)}]` is a forward traveling wave.
`linac_sec1_sim.py` reproduces this with **two** `picmi.LoadAppliedField` objects (one per map,
constants baked into the AMReX parser strings); WarpX sums the named external fields on the
particles. Parameters fixed by details.md:

- **Frequency** `f_RF = 2856 MHz` (S-band).
- **Scale** `scale = sqrt(P_MW / 0.001)` from the 1-kW field normalisation. The map's on-axis
  traveling-wave voltage `∫|Ez|dz = 331 kV` at 1 kW, so on-crest gain ≈ `scale × 331 keV`
  (≈ **35 MeV at the default P = 11 MW**, the LinacSim `sec1_input_power`; cf. details.md's
  "37 MeV @ 15 MW").
- **Phase** `PHASE_DEG` — the absolute synchronous phase is undocumented (capture from β ≈ 0.63
  into a β_phase = 1 wave), so it is an offset relative to the bunch arrival; the default
  `PHASE_DEG = 0` is on the crest (max energy / capture). Sweep it in a `config()` loop to map the
  acceptance: capture/energy peak in a broad plateau near 0° and collapse ~180° away.

## RF capture (focusing is upstream)

A ~146 keV beam injected into a phase-velocity-c wave **slips in phase** and must be *captured*.
Focusing is no longer this stage's job — the injector's three real lenses (Lens 0A / Sol 0 /
Lens 0E) at their true lab z focus the beam upstream and hand it across the 9.547 mm iris.

- **Injection / collimation:** the linac reads the injector dump nearest the **z ≈ 2.03 m
  handoff** and applies the **multi-plane 9.547 mm iris scrape** at injection — the physical iris
  cut at the real 1.922 m iris plane (the beam converges through the 1.922→2.03 m tail, so a single
  2.03 m radial cut would overstate transmission; see `injector/README.md` and
  `pipeline/collimator.py`). At the faithful currents the Sol 0 / Lens 0E matching telescope
  focuses the beam through the iris, so **~32 % of the handoff charge passes** the 9.547 mm aperture.
- **Capture + adiabatic damping:** the captured fraction locks to the wave within the first
  ~0.4 m (β → 1), after which it accelerates and the transverse size **damps** (σ_r ∝ 1/√(γβ)).
  At the faithful 11 MW point capture is **~7 % of the true injected charge** to **⟨KE⟩ ≈
  25 MeV** (max ~32 MeV, σ_KE ≈ 8 MeV). The Sol 0 / Lens 0E matching telescope at z ≈ 1.9 m
  focuses ~32 % of the handoff charge through the 9.547 mm iris; the linac then captures the
  fraction that lands in the RF bucket. It is a **conservative lower bound** (the lab-frame ES
  self-field overestimates transverse SC by ~γ²≈1.66×, so the real machine captures more) and is
  **tune-sensitive to the upstream lens currents**. The optional injector current/phase scans
  (`injector/README.md`) characterize the achievable capture; the faithful currents are not
  retuned to inflate the number. *(A pre-fix ~1 %/0.21 %-vs-1.6 % LENS_0A-sensitivity figure was
  an artifact of the now-fixed LENS_0E 800 mm mislocation — superseded.)*

## Simulation parameters

| parameter | value |
|-----------|-------|
| geometry | RZ, `n_azimuthal_modes = 1` |
| grid | `NR`=16 (r) × `NZ`=1664 (z), r ∈ [0, 9.547 mm], z ∈ [0, 3.5 m] |
| solver | electrostatic, lab frame, Multigrid (self-field only), `REQUIRED_PRECISION`=1e-4, `MAX_ITERS`≤200 |
| applied fields | `linac_rf1/rf2.h5` × `scale` × cos/sin(ωt+φ) (E+B) — two quadrature RF maps only (no solenoid) |
| beam | injector snapshot nearest the **z ≈ 2.03 m handoff**, downsampled to `MAX_PART`=50k (reweighted), z head at `Z_INJECT` = 5 mm |
| time step | `dt = CFL · Δz / v_inject` (`CFL`=0.5; ≈ 5 ps; RF period 0.35 ns) |
| duration | segmented transit estimate to a stop plane short of the absorbing exit (`TRANSIT_MARGIN`=1.0), or fixed via `MAX_STEPS`; `N_DIAGS`=60 dumps |

**Injection point.** `load_injector_bunch` selects the injector dump whose bunch **⟨z⟩ is nearest
`Z_HANDOFF` = 2.03 m** (Z_acc_1) — NOT min-σ_z and NOT max-in-bore charge. With the two-cavity +
solenoid injector the bunch now forms a real longitudinal waist near the handoff, so the old
max-q_bore selector would land on an early, debunched snapshot and silently discard the bunching;
min-σ_z would pick the upstream waist (~1.45 m), not the handoff. The injector places a dump within
~1 mm of 2.03 m (its fine-cadence handoff diagnostic), so nearest-⟨z⟩ lands on the plane. The
**multi-plane 9.547 mm iris scrape** in `load_injector_bunch` (`pipeline/collimator.py`) IS the
physical iris collimation — applied at the real 1.922 m iris plane (the beam *converges* through
the 1.922→2.03 m tail, so a single `r ≤ RMAX` cut at 2.03 m would keep converged halo the real
iris scrapes). Only the survivors are injected; the halo removed is what the real aperture
scrapes, not a numerical-domain artifact.

## Capture bookkeeping

The handoff beam is collimated at the 9.547 mm iris, so **two charge baselines matter** and the
figures/log report both:

- **Injected charge** — every macroparticle handed to the sim at the handoff, recorded by
  `load_injector_bunch` to `diags/main/injection_summary.json` (with the in-iris / in-bore
  breakdown: `q_injected_C`, `q_in_domain_C`, `q_in_bore_C`).
- **In-iris charge** (`q_in_domain_C`) — the **multi-plane iris survivors**: what passes the
  9.547 mm collimator. `load_injector_bunch` scrapes these *before* injection, so WarpX is handed
  only survivors and the first dump already shows this post-collimation charge.

The **capture fraction is reported against the true injected charge** (the honest denominator):
`plot_linac_sec1.py` and the `run_pipeline.py` final-beam summary read the sidecar and fall back
to the first dump only if it is absent (old runs). So the two-stage loss is legible: iris
transmission (in-iris / true-injected ≈ 32 %) × in-iris capture (≈ 22 %) = end-to-end capture vs
true injected (~7 %, a conservative γ² lower bound).

## Gotchas

- **Cell aspect ratio must stay near ≈ 3:1 or MLMG diverges.** The box is long and thin (3.5 m ×
  9.547 mm). `NR=16` over the 9.547 mm bore gives dr=0.597 mm ⇒ ~3.5:1 cells (the captured-beam
  self-field is small, so it still converges); coarsening NR raises the aspect and the self-field
  Multigrid solve aborts (`MLMG failed`). Same lesson as the injector: this solve is
  convergence-bound, not cell-bound — raise NR (÷ blocking factor) rather than coarsening NZ if it
  ever diverges.
- **The applied-field E-init-style guard is unconditional now.** picmi's `LoadAppliedField` forces
  the *global* `E_ext_particle_init_style = "none"` if the last-added field has `load_E=False`. The
  linac now carries only the two RF maps (both `load_E=True`), so the `assert applied[-1].load_E`
  guard is always satisfied — but it is kept so a future reorder/added field fails loudly. (The
  in-linac solenoid that used to require B-before-RF ordering was removed; focusing is upstream.)
- **Stop the run in the exit drift, not at the wall.** Once the bunch clears the absorbing
  z-boundary the domain empties and MLMG aborts. `ZMAX = 3.5 m` leaves a field-free drift past the
  3.12 m structure exit; the transit estimate targets a plane short of `ZMAX` so the beam coasts
  (not absorbed) at the last dump. The estimate uses the on-crest (max) gain — the fastest case —
  so a slower off-crest `PHASE_DEG` (run manually) stays in-domain too.
- Fresh diags per run: `linac_sec1_sim.py` clears the case `OUTDIR` before each run, because WarpX
  appends one openPMD file per dump and a rerun of the same case would otherwise mix iterations.
- **The first diagnostic dump is already post-collimation.** Collimation happens *before*
  injection now (the multi-plane scrape in `load_injector_bunch` hands WarpX only iris survivors),
  and WarpX additionally drops any particle still at r > `RMAX` at injection before the first dump.
  Either way the first dump's charge is *not* the true injected charge. Any capture/survival metric
  must use the true injected charge from
  `injection_summary.json` (see *Capture bookkeeping*) — normalising to the first dump hides the
  iris collimation loss and overstates capture.

## Outputs

`linac_sec1.run()` writes `diags/main/particles/` and `linac_sec1.plot()` reads it,
writing five figures to `results/`:

- `linac_field.png` — the on-axis traveling-wave `|Ez|` amplitude (× scale) and a fixed-t field
  snapshot showing the 2π/3 cell structure.
- `energy_gain.png` — ⟨KE⟩ and max KE vs ⟨z⟩ (~220 keV → ~25 MeV mean / ~32 MeV max for the
  captured slice) with β → 1; the structure shaded.
- `long_phase_space.png` — (z − ⟨z⟩) vs KE at injection / mid / exit: capture into the RF bucket.
- `beam_envelope.png` — σ_x and surviving charge vs ⟨z⟩ with the bore line. The survival panel is
  normalised to the **injected** charge and shows the iris-collimation drop at the first dump
  (r > RMAX = 9.547 mm scrape) followed by the RF-capture loss.
- `exit_spectrum_capture.png` — exit energy spectrum (pC/bin) and the captured fraction **of the
  true injected charge** (~7 % at the faithful 11 MW point), annotated with how much charge
  passed the 9.547 mm iris (~32 %).

## Notes / caveats

- We model **one** section — the single 3 m structure these two quadrature maps describe (~37 MeV
  at 15 MW). The reference linac has 8 sections (2–8 are different CEA/CU designs, modelled in
  BMAD); they would become `linac_sec2`, … with their own maps.
- **Focusing moved upstream to the injector.** The old in-linac solenoid (`linac_sol.h5`,
  `I_SOL`, `SOL_Z`) was a stand-in for the real Lens 0A / Sol 0 / Lens 0E that physically live at
  z ≈ 0.23 / 1.90 / 1.91 m. The injector now applies those three real lenses at their true lab z
  and hands a focused, 9.547 mm-collimated beam across the 2.03 m plane; the linac carries only the
  two RF maps. RF power (11 MW) is the original LinacSim `sec1_input_power`; the absolute RF phase
  is undocumented (`PHASE_DEG` scanned for the crest).
- **Capture is ~7 % of true injected, a conservative lower bound, tune-sensitive.** The Sol 0 /
  Lens 0E matching telescope focuses ~32 % of the handoff charge through the 9.547 mm iris, and the
  fraction landing in the RF bucket is captured. The lab-frame ES self-field overestimates
  transverse SC by ~γ²≈1.66× ⇒ the real machine captures more; capture also responds strongly to
  the upstream lens currents. Not a precision-tuned number — the injector current/phase scans
  characterize the achievable value. *(A pre-fix ~1 %/~7×-LENS_0A figure was an artifact of the
  now-corrected LENS_0E 800 mm mislocation — superseded.)*
- The r-domain `RMAX = 9.547 mm` IS the SLAC bore / collimator iris; the RF maps reach it exactly.
  A particle that reaches the wall is scraped on the iris (counted against the true-injected
  capture); one starting beyond `RMAX` is dropped at injection (the collimation loss) and likewise
  accounted for via the injected-charge denominator — see *Capture bookkeeping*.
- The lab-frame electrostatic self-field omits the `1/γ²` magnetic-pinch cancellation (it applies
  the rest-frame Coulomb force `qE_r`, not `qE_r/γ²`), so it overestimates the transverse
  space-charge force by ~γ² — largest at the low-energy injection (~220 keV handoff; cf. the gun's
  β≈0.63 → ~66 % overestimate in `gun/README.md`) and
  shrinking toward negligible once captured (γ ≫ 1). Space charge is a small perturbation here, so
  this is acceptable for the demonstration.
