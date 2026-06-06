# WarpX SLAC Linac — Section 1 (RZ)

Fourth stage of the Cornell Linac chain modelled in WarpX, and the **first downstream
accelerating section** (later sections → `linac_sec2`, … as their field maps are added):

```
cathode (cathode/) -> gun (gun/) -> prebuncher (prebuncher/) -> linac_sec1 (this)
```

The prebuncher's bunched beam (~148 keV, β ≈ 0.63, 0.066 nC) is captured by a 3 m,
86-cell, 2π/3 **traveling-wave** SLAC accelerating structure and accelerated to **~37 MeV**
(⟨KE⟩ ≈ 34 MeV after capture, peak 38 MeV), with a solenoid focusing channel and
self-consistent space charge.

## Running

```python
# from repo root, in the CBB env:
import linac_sec1
linac_sec1.demo()        # RF-phase acceptance scan + headline + focus-off + all 6 figures
# linac_sec1.run()       # a single case (the default operating point) -> diags/main
# linac_sec1.plot()      # re-generate figures from existing diags/
```

- **`linac_sec1.demo()`** is the full demonstration: a coarse RF-phase acceptance scan
  (`diags/scan_phi<deg>`), then a full-resolution headline run at the crest (`diags/main`) and a
  focus-off comparison (`diags/focusoff`), then the aggregate plots. Tune via
  `demo(phases=..., scan_nz=..., full_nz=..., i_sol=...)`.
- **`linac_sec1.run()`** runs one case; operating point + grid live at the top of
  `linac_sec1/linac_sec1_sim.py` (`POWER_MW`, `PHASE_DEG`, `I_SOL`, `NZ`, …), overridable via
  `config()`.

`build_linac_sec1_field` reads the maps from `fieldmaps/`; the sim reads the prebuncher output
from `prebuncher/diags/P800_zc/particles/` (repo-root-relative). To run the whole chain
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
  (≈ **40 MeV at the default P = 15 MW**; matches details.md's "37 MeV @ 15 MW").
- **Phase** `PHASE_DEG` — the absolute synchronous phase is undocumented (capture from β ≈ 0.63
  into a β_phase = 1 wave), so it is an offset **scanned** to find the crest; `demo()` picks the
  maximum-energy phase for the headline.

## Solenoid focusing and RF capture

The prebuncher beam arrives **diverging** (it had no focusing over its 1.3 m drift, so it expands
toward the bore), and a 148 keV beam injected into a phase-velocity-c wave **slips in phase** and
must be *captured*. Both effects make focusing essential:

- **Solenoid** (`I_SOL`, A): the per-Ampere `SOL_0` map × `I_SOL`. Capture rises sharply with
  current — **~4% (I = 0) → ~95% (I = 1000 A ≈ 0.15 T peak)** — because the strong channel holds
  the diverging beam inside the bore long enough to be captured. `I_SOL = 0` is the focus-off
  comparison.
- **Capture + adiabatic damping:** the captured fraction locks to the wave within the first
  ~0.4 m (β → 1), after which it accelerates linearly to ~37 MeV and the transverse size **damps**
  (σ_r ∝ 1/√(γβ)): the RMS σ_x falls to ~2 mm (r95 ~6 mm), well inside the 9.55 mm bore, by the
  structure exit. At the crest the headline run captures **~95 %** of the beam at ⟨KE⟩ ≈ 34 MeV
  (σ_KE ≈ 18 %); focus-off captures only the ~3 % on-axis core.

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
  simply coasts at 148 keV). Adding it first — so an RF map (`load_E=True`) is last — keeps E on.
  Verify with `grep E_ext_particle warpx_used_inputs` (must be `read_from_file`, not `none`).
- **Stop the run in the exit drift, not at the wall.** As with the prebuncher, once the bunch
  clears the absorbing z-boundary the domain empties and MLMG aborts. `ZMAX = 3.5 m` leaves a
  field-free drift past the 3.12 m structure exit; the transit estimate targets a plane short of
  `ZMAX` so the beam coasts (not absorbed) at the last dump. The estimate uses the on-crest (max)
  gain so the fastest case is the binding one — off-crest scan cases run the same length and may
  be measured slightly before full traversal (the acceptance curve's peak/shape is unaffected).
- Fresh diags per run: `linac_sec1_sim.py` clears the case `OUTDIR` before each run, because WarpX
  appends one openPMD file per dump and a rerun of the same case would otherwise mix iterations.

## Outputs

`linac_sec1.demo()` writes `diags/{scan_phi<deg>, main, focusoff}/{fields,particles}/` and
`linac_sec1.plot()` reads them, writing to `results/`:

- `linac_field.png` — the on-axis traveling-wave `|Ez|` amplitude (× scale) and a fixed-t field
  snapshot showing the 2π/3 cell structure.
- `energy_gain.png` — ⟨KE⟩ and max KE vs ⟨z⟩ (148 keV → ~37 MeV) with β → 1; the structure shaded.
- `long_phase_space.png` — (z − ⟨z⟩) vs KE at injection / mid / exit: capture into the RF bucket.
- `beam_envelope.png` — σ_r and surviving charge vs ⟨z⟩ with **focus ON vs OFF** and the bore line.
- `exit_spectrum_capture.png` — exit energy spectrum (pC/bin) and the captured-charge fraction.
- `phase_acceptance.png` — final ⟨KE⟩ and capture fraction vs injection RF phase (the scan).

## Notes / caveats

- We model **one** section — the single 3 m structure these two quadrature maps describe (~37 MeV
  at 15 MW). The reference linac has 8 sections (2–8 are different CEA/CU designs, modelled in
  BMAD); they would become `linac_sec2`, … with their own maps.
- **Element positions, the solenoid current, and the absolute RF phase are undocumented.** We
  inject at the prebuncher's bunching focus, place the structure after a short drift, scan the RF
  phase for the crest, and choose `I_SOL` for good capture — a physically sensible demonstration,
  not a calibrated reproduction of the real spacing.
- The r-zero-padding of the RF maps at the 9.55 mm bore is a sharp truncation (the metal iris); a
  beam focused well inside the bore rarely samples it, and particles that reach the domain wall are
  treated as lost on the iris (counted against the capture fraction).
- The lab-frame electrostatic self-field uses a single relativistic factor while the beam spans
  148 keV → 37 MeV; space charge is a small perturbation that becomes negligible once captured
  (γ ≫ 1), so this is acceptable for the demonstration.
