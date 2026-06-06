# WarpX SLAC Linac — Section 1 (RZ)

Fourth stage of the Cornell Linac chain modelled in WarpX, and the **first downstream
accelerating section** (later sections → `linac_sec2`, … as their field maps are added):

```
cathode (cathode/) -> gun (gun/) -> prebuncher (prebuncher/) -> linac_sec1 (this)
```

The prebuncher's beam (~137 keV, β ≈ 0.63, ~0.17 nC, **largely unbunched** at the
faithful-to-LinacSim 8 kW prebuncher default) enters a 3 m, 86-cell, 2π/3 **traveling-wave**
SLAC accelerating structure with a solenoid focusing channel and self-consistent space charge.
At the original LinacSim operating point (Sol 0 = 40 A, 11 MW), the weak focusing captures only
**~2 %** of the injected beam (~5.7 pC); that captured slice reaches **⟨KE⟩ ≈ 15.5 MeV**
(max ≈ 29.9 MeV, σ_KE ≈ 7.9 MeV). See *Solenoid focusing and RF capture* for why this is far
below the strongly-focused (`I_SOL ≈ 1000 A`) capture the stage can reach.

## Running

```python
# from repo root, in the CBB env:
import linac_sec1
linac_sec1.run()         # build field + sim + plots -> diags/main, results/
# linac_sec1.run(plots=False)
# linac_sec1.plot()      # re-generate figures from existing diags/
```

`run()` runs **one case** at the default operating point — now the original LinacSim values
(`PHASE_DEG=0`, `I_SOL=40`, `POWER_MW=11`). The operating point and grid are module-level
constants at the top of `linac_sec1/linac_sec1_sim.py` (`POWER_MW`, `PHASE_DEG`, `I_SOL`, `NZ`, …),
overridable via `config()` — e.g. `linac_sec1.config(I_SOL=0); linac_sec1.run()` for the unfocused case, or a
`PHASE_DEG` sweep in a Python loop if you want the acceptance curve.

`build_linac_sec1_field` reads the maps from `fieldmaps/`; the sim reads the prebuncher output
from `prebuncher/diags/P8_zc/particles/` (repo-root-relative). To run the whole chain
(cathode → gun → prebuncher → linac_sec1), use **`pipeline/run_pipeline.py`**.

## Field maps

**This stage uses three GPT maps.** The two SLAC files are **not two sections** — they are the
**real and imaginary (quadrature) components of one** 3 m traveling-wave structure
(`reference/Linac Simulation Documentation/details.md`):

| file | columns (used) | grid | normalisation |
|------|----------------|------|---------------|
| `SLAC-3mLinac-field1.gdf` | `ErRe, EzRe, HphiIm` | 21 (r) × 6379 (z), r ≤ 9.55 mm, z ≈ 3.0 m | 1 kW input power |
| `SLAC-3mLinac-field2.gdf` | `ErIm, EzIm, HphiRe` | same grid | 1 kW input power |
| `SOL_0.gdf` (`SOL_MAP`) | `Br, Bz` | 16 (r) × 939 (z), r ≤ 40 mm | 1 A coil current |

`build_linac_sec1_field.py` writes three openPMD files (all `thetaMode` m = 0, components r/t/z,
shape `(1, nr, nz)`, the layout WarpX's `read_from_file` reader expects):

- `linac_rf1.h5` — `E = (ErRe, 0, EzRe)`, `B = (0, HphiIm, 0)` (the in-phase quadrature),
- `linac_rf2.h5` — `E = (ErIm, 0, EzIm)`, `B = (0, HphiRe, 0)` (the 90° quadrature),
- `linac_sol.h5` — `B = (Br, 0, Bz)` (static solenoid, no E).

The SLAC maps only reach the **9.55 mm bore** in r; they are **zero-padded in r** out to the sim
domain `RMAX = 12 mm` so every applied field explicitly covers the domain (a particle in the bore
shadow then feels an exact zero RF field). The solenoid already reaches 40 mm. Each map is placed
in the lab frame via `grid_global_offset` (`Z_STRUCT = 0.10 m` structure entrance; `SOL_Z` slides
the solenoid peak into the low-energy capture region).

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

## Solenoid focusing and RF capture

The prebuncher beam arrives **diverging** (it had no focusing over its drift, so it expands
toward the bore), and a ~137 keV beam injected into a phase-velocity-c wave **slips in phase** and
must be *captured*. Both effects make focusing essential:

- **Solenoid** (`I_SOL`, A): the per-Ampere `SOL_0` map × `I_SOL`. Capture rises sharply with
  current — **~2 % at the LinacSim Sol 0 default (`I_SOL = 40 A`) → ~95 % at `I_SOL ≈ 1000 A`
  (≈ 0.15 T peak)** — because only a strong channel holds the diverging beam inside the 9.55 mm
  bore long enough to be captured. The default is now the original LinacSim **`I_SOL = 40 A`**,
  which is far too weak for this standalone WarpX capture, so most of the beam is scraped on the
  iris. Run `config(I_SOL=1000); run()` to recover the strongly-focused, high-capture case.
- **Capture + adiabatic damping:** the captured fraction locks to the wave within the first
  ~0.4 m (β → 1), after which it accelerates and the transverse size **damps** (σ_r ∝ 1/√(γβ)).
  At the faithful-to-LinacSim default (40 A, 11 MW), only **~2 %** of the injected beam is captured
  (~5.7 pC of ~169 pC), reaching **⟨KE⟩ ≈ 15.5 MeV** (max ≈ 29.9 MeV, σ_KE ≈ 7.9 MeV) — the rest
  is lost on the bore. (For comparison, the earlier strongly-focused 15 MW / 1000 A point captured
  ~97 % to a core near ~37–38 MeV; rerun with those overrides to reproduce it.) The low capture
  here is the expected consequence of the un-tuned original solenoid/power — see *Notes / caveats*.

## Simulation parameters

| parameter | value |
|-----------|-------|
| geometry | RZ, `n_azimuthal_modes = 1` |
| grid | `NR`=16 (r) × `NZ`=1664 (z), r ∈ [0, 12 mm], z ∈ [0, 3.5 m] |
| solver | electrostatic, lab frame, Multigrid (self-field only), `REQUIRED_PRECISION`=1e-4, `MAX_ITERS`≤200 |
| applied fields | `linac_rf1/rf2.h5` × `scale` × cos/sin(ωt+φ) (E+B) **+** `linac_sol.h5` × `I_SOL` (static B) |
| beam | prebuncher snapshot at its min-σ_z **bunching focus**, downsampled to `MAX_PART`=50k (reweighted), z head at `Z_INJECT` = 5 mm |
| time step | `dt = CFL · Δz / v_inject` (`CFL`=0.5; ≈ 5.6 ps; RF period 0.35 ns) |
| duration | segmented transit estimate to a stop plane short of the absorbing exit (`TRANSIT_MARGIN`=1.0), or fixed via `MAX_STEPS`; `N_DIAGS`=60 dumps |

**Injection point.** The prebuncher beam is bunched only transiently in its drift; the sim
auto-selects the **minimum-σ_z snapshot past the cavity** (`Z_FOCUS_MIN`), the only point where
the beam is both bunched and still inside the 9.5 mm bore. (The global σ_z minimum is the
*pre-modulation* injection snapshot near z = 0, deliberately excluded.)

## Gotchas

- **Cell aspect ratio must stay ≈ 3:1 or MLMG diverges.** The box is long and thin (3.5 m ×
  12 mm). At `NR=32` the cells are ≈5.6:1 and the self-field Multigrid solve aborts
  (`MLMG failed`); `NR=16` gives ≈2.8:1 (matching the prebuncher) and converges. Same lesson as
  the prebuncher: this solve is convergence-bound, not cell-bound.
- **Add the solenoid applied field BEFORE the RF maps.** The (env-local) picmi
  `LoadAppliedField` sets the *global* `E_ext_particle_init_style = "none"` for any field with
  `load_E=False`, and **the field initialised last wins**. The solenoid is B-only (`load_E=False`),
  so if it is added last it disables the accelerating E field for the whole run (silent: the beam
  simply coasts at ~137 keV). Adding it first — so an RF map (`load_E=True`) is last — keeps E on.
  Verify with `grep E_ext_particle warpx_used_inputs` (must be `read_from_file`, not `none`).
- **Stop the run in the exit drift, not at the wall.** As with the prebuncher, once the bunch
  clears the absorbing z-boundary the domain empties and MLMG aborts. `ZMAX = 3.5 m` leaves a
  field-free drift past the 3.12 m structure exit; the transit estimate targets a plane short of
  `ZMAX` so the beam coasts (not absorbed) at the last dump. The estimate uses the on-crest (max)
  gain — the fastest case — as the binding constraint, so a slower off-crest `PHASE_DEG` (run
  manually) stays in-domain too.
- Fresh diags per run: `linac_sec1_sim.py` clears the case `OUTDIR` before each run, because WarpX
  appends one openPMD file per dump and a rerun of the same case would otherwise mix iterations.

## Outputs

`linac_sec1.run()` writes `diags/main/particles/` and `linac_sec1.plot()` reads it,
writing five figures to `results/`:

- `linac_field.png` — the on-axis traveling-wave `|Ez|` amplitude (× scale) and a fixed-t field
  snapshot showing the 2π/3 cell structure.
- `energy_gain.png` — ⟨KE⟩ and max KE vs ⟨z⟩ (~137 keV → ~15.5 MeV mean / ~30 MeV max for the
  captured slice) with β → 1; the structure shaded.
- `long_phase_space.png` — (z − ⟨z⟩) vs KE at injection / mid / exit: capture into the RF bucket.
- `beam_envelope.png` — σ_r and surviving charge vs ⟨z⟩ with the bore line: focusing + adiabatic
  damping (only ~2 % survives at the weak 40 A focusing).
- `exit_spectrum_capture.png` — exit energy spectrum (pC/bin) and the captured-charge fraction.

## Notes / caveats

- We model **one** section — the single 3 m structure these two quadrature maps describe (~37 MeV
  at 15 MW). The reference linac has 8 sections (2–8 are different CEA/CU designs, modelled in
  BMAD); they would become `linac_sec2`, … with their own maps.
- **The solenoid current (40 A) and RF power (11 MW) are now the original LinacSim values**
  (`current_sol0`, `sec1_input_power`); element positions and the absolute RF phase remain
  undocumented. We inject at the prebuncher snapshot, place the structure after a short drift, and
  scan the RF phase for the crest. Note the original Sol 0 = 40 A is the *prebuncher-region*
  focusing solenoid and is far too weak to capture this diverging standalone beam into the bore
  (hence ~2 % capture); the real injector relies on the full lens chain (0A–0E) plus two
  prebunchers, which this single-stage WarpX model does not include. Set `I_SOL ≈ 1000 A` for a
  capture-optimised demonstration.
- The r-zero-padding of the RF maps at the 9.55 mm bore is a sharp truncation (the metal iris); a
  beam focused well inside the bore rarely samples it, and particles that reach the domain wall are
  treated as lost on the iris (counted against the capture fraction).
- The lab-frame electrostatic self-field omits the `1/γ²` magnetic-pinch cancellation (it applies
  the rest-frame Coulomb force `qE_r`, not `qE_r/γ²`), so it overestimates the transverse
  space-charge force by ~γ² — largest at injection (~137 keV → ~66 %, as in `gun/README.md`) and
  shrinking toward negligible once captured (γ ≫ 1). Space charge is a small perturbation here, so
  this is acceptable for the demonstration.
