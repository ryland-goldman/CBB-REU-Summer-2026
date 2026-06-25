# Linac Sections 1–4 — SLAC Traveling-Wave Acceleration (WarpX RZ)

The four downstream accelerating sections of the Cornell linac chain, merged into ONE
parametrized WarpX driver. Each section is a 3 m, 86-cell, **2π/3 traveling-wave** S-band
(2856 MHz) accelerating structure modelled in WarpX RZ (single azimuthal mode) with
self-consistent space charge:

```
cathode → gun → injector → linac sec 1 (capture) → sec 2 → sec 3 → sec 4 → [4→5 boundary] → converter → linac sec 5-8
                            (WarpX RZ)               (WarpX RZ)  (WarpX RZ)  (WarpX RZ)        (G4beamline)  (Impact-T)
```

Section 1 **captures** the injector's ~150 keV velocity-bunched beam; sections 2, 3 and 4 each
**accelerate** the captured relativistic core by one more section's worth of energy. All four
share the same WarpX setup, the same SLAC field maps, and one driver (`sim/linac1-4.py`); they
differ only in the input beam and the frozen RF setpoints. Section 4's exit is the **4→5 boundary**:
the input to the e+/e- converter target (`sim/converter.py`), whose positron output then feeds the
Impact-T sections 5-8.

Built on `pywarpx`, driven through lume-warpx: every constant lives in `config/linac{1,2,3,4}.yaml`,
and `sim/linac1-4.py` reads them back, imports the upstream beam via `WarpX(initial_particles=...)`,
and overrides only the runtime-computed values (the two quadrature RF time functions, step count,
`dt`, diagnostic period).

Run (the section number `N ∈ {1,2,3,4}` is the sole argument):
```bash
conda activate CBB
python sim/linac1-4.py 1        # capture; reads the injector handoff, writes logs/diags/linac1-4/sec1/main
python sim/linac1-4.py 2        # accelerate; reads sec 1's exit
python sim/linac1-4.py 3        # accelerate; reads sec 2's exit
python sim/linac1-4.py 4        # accelerate; reads sec 3's exit (the 4→5 boundary handoff to the converter)
python sim/plot/linac1-4.py 1   # figures → logs/plots/linac1-4/sec1_*.png   (likewise 2, 3, 4)
```

The sections are a chain: section 2 reads section 1's output, …, section 4 reads section 3's, so
they must run **in order** (1 → 2 → 3 → 4). `sim/linac1-4.py main()` runs only the simulation;
`sim/plot/linac1-4.py main()` runs only the plotting (the section sim must have been run first).

---

## The field model — one shared SLAC traveling-wave structure, reused

Cornell linac sections 2–8 have **no dedicated GPT/CST field maps** (none exist; the reference
LinacSim/BMAD decks model them with a generic constant-gradient linac function). All four WarpX
sections therefore **reuse the SLAC 3 m traveling-wave quadrature maps** — the same spatial shape,
the same local entrance — and realise each section's energy gain by the field **amplitude scale**,
not by a different map.

The two SLAC files are **not two structures**: they are the **real and imaginary (quadrature)
components of one** 3 m structure. `sim.helpers.buildfields.build_linac_slac_fields()` parses the
two GPT maps and writes two openPMD `thetaMode` (m = 0) files, shared by all four sections:

| file | columns used | content |
|------|--------------|---------|
| `fieldmaps/h5/linac_rf1.h5` | `ErRe, EzRe, HphiIm` | `E = (ErRe, 0, EzRe)`, `B = (0, HphiIm, 0)` — the in-phase quadrature |
| `fieldmaps/h5/linac_rf2.h5` | `ErIm, EzIm, HphiRe` | `E = (ErIm, 0, EzIm)`, `B = (0, HphiRe, 0)` — the 90° quadrature |

`build_linac_slac_fields()` is **idempotent** and called at the start of `main()`, so the maps are
built once and re-used by every section (sections 2 and 3 do not rebuild them). The SLAC maps reach
the **9.55 mm bore** in r; they are zero-padded in r out to the sim domain `RMAX = 9.547 mm` (the
SLAC bore / injector→linac collimator iris). So the **radial domain IS the aperture**: a particle
that reaches the wall is scraped exactly as the real iris does. Each map is placed in the lab frame
via `grid_global_offset` (`Z_STRUCT = 0.10 m` structure entrance), which anchors the RF phase.

### Synthesising the traveling wave

GPT builds the forward traveling wave as the sum of two standing waves 90° apart: each map is
driven `E(t) = map·scale·cos(ωt+φ)`, `Bφ(t) = map·scale·sin(ωt+φ)`, with map 2 offset by **+π/2**.
The driver reproduces this with **two** applied-from-file fields (the two quadrature halves); WarpX
sums the named external fields on the particles, and the sum
`Re[(Ẽ_re + iẼ_im)·e^{i(ωt+φ)}]` is a forward traveling wave. The cos/sin amplitude+phase strings
are built by `sim.helpers.tools.rf_time_functions(scale, ω, φ)` and injected at runtime.

The structure's on-axis traveling-wave voltage is `∫|Ez|dz = V1KW_KEV = 331.2 keV` at the 1-kW
field normalisation, so the on-crest mean gain is ≈ `scale × 331.2 keV`.

---

## The RF setpoints — frozen, not found at runtime

Each section needs an RF **amplitude** (how hard to drive) and a **phase** (where on the RF wave the
bunch sits). How those are set differs between the capture section and the accelerating sections:

- **Section 1 (capture).** A ~150 keV beam injected into a phase-velocity-c wave **slips in phase**
  and must be *captured*. The amplitude is set from the klystron input power,
  `scale = sqrt(POWER_MW / RF_NORM_MW)` (`POWER_MW = 11 MW`, the original LinacSim `sec1_input_power`;
  `RF_NORM_MW = 1 kW` is the map normalisation). For section 1 `PHASE_DEG` is the **absolute**
  arrival-referenced base phase — the driver applies it directly as `base_deg`, NOT as a detune from
  a separate crest. It is the capture crest `sim/autophase.py` writes into `config/linac1.yaml`, so
  it is re-derived with the upstream beam and **`0` is not on-crest**. (Sections 2/3/4 are the other
  convention: there `PHASE_DEG` is a detune from the frozen `CREST_PHASE_DEG`, with `0` = on-crest.)

- **Sections 2, 3, 4 (accelerate).** The captured core is now **β ≈ 1 and locked** to the wave, and it
  is **micro-bunched** at λ_RF (≈105 mm): the charge sits in a narrow RF-phase band but spreads >1 λ
  in z, so the geometric z-centroid is offset from the charge's phase-centroid. The on-crest base
  phase therefore differs from section 1's literal 0° (at which a locked β ≈ 1 beam would
  *decelerate*) and from a bare single-particle crest. The original sims found this crest and the
  gain-per-scale **at runtime** (`tw_crest_phase`, folding the field integral against the bunch's
  phase distribution, then scaling to a per-section ΔE target).

  **This merged driver drops that runtime derivation.** The crest base phase and the field scale
  were derived **once** from the old `linac_sec1 → sec2 → sec3` chain and are **hardcoded** as
  `CREST_PHASE_DEG` and `FIELD_SCALE` in each section's `config/linac{2,3,4}.yaml`. Those YAML values
  are **authoritative** — they are re-derived (rewritten in place by `sim/autophase.py`) whenever the
  upstream beam changes, so the actual numbers are not reproduced here. `DE_TARGET_MEV` (the
  details.md CEA per-section ΔE @11 MW, √P-scaled) is the energy budget the scale was derived to hit;
  it is kept only as a comment/reference — the runtime no longer reads it for any field calculation.
  `PHASE_DEG` remains as a **detune** offset from the frozen crest (default 0 = on crest).

  > **Section 4 is new and not yet calibrated.** `config/linac4.yaml` ships `CREST_PHASE_DEG` /
  > `FIELD_SCALE` as **placeholders copied from section 3**. Once the chain through section 3 is
  > run, re-derive the crest with `python sim/autophase.py 4` and re-fit `FIELD_SCALE` to
  > `DE_TARGET_MEV` (the CU 5 √P-scaled budget). Until then section-4 output is not physical.

The RF block is otherwise **uniform** across all four sections:

```
omega = 2π·F_RF
t_in  = (Z_STRUCT − z_center) / v_beam            # z_center = weighted mean z of the injected bunch
phi   = −omega·t_in + deg2rad(base_deg)           # base_deg = PHASE_DEG (sec 1) | CREST_PHASE_DEG+PHASE_DEG (sec 2/3/4)
phi2  = phi + π/2                                  # the 90° quadrature half
e1,b1 = rf_time_functions(scale, omega, phi)
e2,b2 = rf_time_functions(scale, omega, phi2)
```

(`scale = sqrt(POWER_MW/RF_NORM_MW)` for section 1, `scale = FIELD_SCALE` for sections 2/3/4.)

---

## The input beam per section

### Section 1 — the injector handoff + the iris scrape

Section 1 reads the **injector** diagnostics (`logs/diags/injector/main/particles`), selects the
dump whose bunch ⟨z⟩ is nearest the `Z_HANDOFF = 2.03 m` handoff plane (NOT min-σ_z, NOT
max-in-bore charge — the injector forms a real longitudinal waist at the handoff), and applies the
**multi-plane 9.547 mm iris scrape** (`load_injector_bunch`, using `pipe_violator_ids` /
`survivor_mask`). The scrape is a union-of-ids cut over the dumps from `COLLIM_Z = 1.922 m` to just
past the handoff — NOT a single radial cut — because the beam **converges** through the 1.922 → 2.03 m
tail, so a single `r ≤ RMAX` cut at 2.03 m would keep converged halo the real iris scrapes. That
scrape IS the physical injector→linac iris collimation; only the survivors are injected.

Capture is reported against the **true injected charge** (every macroparticle at the handoff, all r),
so the iris transmission loss is legible. `load_injector_bunch` records `q_injected_C`,
`q_in_domain_C` (in-iris survivors), and `q_in_bore_C` (within the RF bore) to the sidecar.

### Sections 2, 3 — the previous section's captured core

Sections 2 and 3 read the previous section's exit dump
(`logs/diags/linac1-4/sec{N−1}/main/particles`) via `load_warpx_exit_bunch`, which picks the last
well-populated dump (the captured beam coasting in the field-free exit drift), keeps only the
**captured core** (`KE ≥ 0.5 · median KE`), downsamples it (reweighted), and shifts its tail to
`Z_INJECT`. The captured-core cut is essential: the section-exit dump trails a sparse slipping
low-energy tail that lags the relativistic core by ~metres in z and is **not in the RF bucket** —
genuinely lost between sections (the same physics as the linac sec 4–8 `MIN_KE_MEV` cut). There is
**no iris scrape** between sections — the 9.547 mm collimation is the one-time injector→linac event
at the section-1 entrance.

### Lab-z chaining

Each section runs in its own **local** frame (the bunch tail enters at `Z_INJECT = 5 mm`), but the
sections are physically contiguous, so each records where its injection sits in the **lab** frame.
Section 1 records `z_handoff_m`; sections 2 and 3 record `z_inject_lab_m`, computed by
`upstream_exit_lab_z(prev injection_summary.json, info["exit_zmean_local_m"])` — it reads the
upstream section's local→lab offset and adds the upstream exit dump's local ⟨z⟩, chaining the local
frames into one lab-frame z so a chain plotter can place each segment correctly.

---

## What the simulation does (`sim/linac1-4.py`)

- **Geometry / grid.** RZ, `n_azimuthal_modes = 1`, `NR = 16 × NZ = 1664`, r ∈ [0, 9.547 mm],
  z ∈ [0, 3.5 m] (≈ 3.5:1 cells). Identical across all four sections (the same reused SLAC map).
- **Solver.** Electrostatic, lab frame, Multigrid self-field only (`ES_MLMG_LF`,
  `required_precision = 1e-4`, `maximum_iterations ≤ 200`); space charge ON (a small perturbation at
  these energies).
- **Applied fields.** The two quadrature RF maps × `scale` × cos/sin(ωt+φ) (E+B); no solenoid
  (transverse focusing is upstream in the injector — these sections carry none).
- **Timestep.** `dt = CFL · (ZMAX / NZ) / v_inject` (`CFL = 0.5`).
- **Duration — the segmented-transit stop plane.** The run length is a segmented transit estimate to
  a stop plane **short of ZMAX** so the bunch finishes in the **field-free exit drift, NOT at the
  absorbing wall** (an emptied domain aborts the MLMG solve). The estimate sums the drift to the
  structure entrance, a capture/ramp length `L_cap = 0.40 m` at the mean β, and the rest at the final
  β; it uses the on-crest (max) gain — the fastest case — so a slower off-crest detune stays
  in-domain too. Sections 2 and 3 additionally subtract the bunch's z-extent (`z_span`) so the bunch
  **head** (not just the centroid) stops clear of the wall — the captured core is ~1 RF wavelength
  long. (`TRANSIT_MARGIN = 1.0`, `N_DIAGS = 60` dumps; `MAX_STEPS > 0` overrides with a fixed count.)
- **The load_E guard.** picmi's applied-from-file field forces the *global*
  `E_ext_particle_init_style = "none"` if the last-added field has `load_E = False`. Both RF maps
  `load_E = true`, so `assert fields[-1].load_E` always holds — kept so a future reorder/added field
  fails loudly.
- **Output.** openPMD particles to `logs/diags/linac1-4/sec{N}/main/particles`, plus
  `injection_summary.json` (charge in + the local→lab z chain) to `logs/diags/linac1-4/sec{N}/main/`.
  Stale diagnostics are removed before each run (the h5 backend appends one file per dump).

---

## The figures (`sim/plot/linac1-4.py` → `logs/plots/linac1-4/sec{N}_*.png`)

Four figures per section, all via lume-warpx's plotting helpers over the section's last populated
dump (no field diagnostic is dumped, so no `plot_fields`):

- **`sec{N}_phase_space_z_KE.png`** — `plot2D("z","kinetic_energy")`: the captured/accelerated
  longitudinal phase space (capture in sec 1; one more section of gain in sec 2, sec 3, sec 4).
- **`sec{N}_transverse_x_px.png`** — `plot2D("x","px")`: the exit transverse phase space within the bore.
- **`sec{N}_centroid_vs_t.png`** — `plot1D("t","mean_z")`: the bunch crossing the 3 m structure + drift.
- **`sec{N}_emittance_vs_t.png`** — `plot1D("t","norm_emit_x")`: transverse emittance over the run.

---

## Notes & gotchas

- **Cell aspect ratio must stay near ≈ 3:1 or MLMG diverges.** The box is long and thin (3.5 m ×
  9.547 mm). This solve is convergence-bound, not cell-bound — raise `NR` (÷ blocking factor) rather
  than coarsening `NZ` if the self-field solve ever aborts.
- **Stop the run in the exit drift, not at the wall.** `ZMAX = 3.5 m` leaves a field-free drift past
  the 3.12 m structure exit; the transit estimate targets a plane short of `ZMAX` so the beam coasts
  (not absorbed) at the last dump.
- **The first section-1 diagnostic dump is already post-collimation.** The iris scrape happens
  *before* injection (WarpX is handed only survivors), so the first dump's charge is *not* the true
  injected charge — any capture metric must use `q_injected_C` from `injection_summary.json`.
- **No in-section transverse focusing.** Like the real machine model, these sections carry no
  solenoid/lens; the long, unfocused bunch loses its phase/radial edges to the bore over 3 m, and
  off-axis particles sample a lower Ez than the on-axis ∫|Ez|dz integral — so the achieved mean gain
  runs below the 1-D on-axis target. This is a model artifact (the focusing reconciliation backlog),
  not a phasing error.
- **The lab-frame electrostatic self-field omits the 1/γ² magnetic-pinch cancellation**, so it
  overestimates the transverse space-charge force by ~γ² — largest at the low-energy section-1
  capture and negligible once the beam is relativistic (γ ≫ 1). Space charge is a small perturbation
  here, so this is acceptable for the demonstration.
- **The frozen setpoints were derived once and hardcoded.** They reproduce the old chain's operating
  point; re-deriving them (e.g. after a change to the upstream beam) is a deliberate, separate step,
  not a per-run computation — edit `CREST_PHASE_DEG` / `FIELD_SCALE` in the section yaml to retune.
