# Linac sections 4-8 (Impact-T)

The rest of the straight electron line to CHESS, after the WarpX linac sections 1-3: five S-band
(2856 MHz) traveling-wave (TW) accelerating sections (CU 5 + CEA 4/5 + CU 3/4), chained into **one**
Impact-T deck and integrated as one time-ordered beam. Reads the linac1-3 sec3 captured relativistic
exit beam (the **3->4 boundary**) and accelerates the captured core on-crest through the five
sections.

```
... linac1-3 (sec1/2/3, WarpX RZ) --[3->4 boundary]--> linac4-8 (this, Impact-T, 5 TW sections, ~29 m)
```

The 3->4 boundary (the sec3 exit) is the future slot for an e+/e- converter target (Geant4-class
physics, out of scope here).

## Why Impact-T, not WarpX

The upstream stages are WarpX/pywarpx runs. This stage is an external serial Impact-T run
(`ImpactTexe`) driven through **lume-impact**, in-process (no pywarpx global-geometry binding, so no
per-stage subprocess isolation is needed). Impact-T integrates one beam through one time-ordered
lattice -- the natural fit for the remaining sections that BMAD/LinacSim also treat as one generic
linac.

The run executes in `logs/diags/linac4-8/` (`use_temp_dir=False`, `workdir=` set), so its `fort.18`
lands at a known `logs/diags/linac4-8/fort.18`. The progress bar (`sim.helpers.tqdmwrapper.
impact_progress`) polls that file's column 1 (reference `z` [m]) and advances 0 -> total lattice
length while `I.run()` executes.

## Field model -- generic constant-gradient TW, no field maps

Sections 4-8 have **no GPT/CST field maps** (none exist; LinacSim/BMAD model them with the generic
constant-gradient linac function). We reuse the shipped lume-impact `traveling_wave_cavity` template
field **shape** -- `rfdata4-7`, vendored into `fieldmaps/rfdata/` -- verbatim, and put **all**
per-section physics in the per-section field scale. Each section is the template's 4-line `solrf`
superposition of two standing-wave maps (G. A. Loew et al., SLAC-PUB-2295):

| line     | rfdata | length             | theta0 offset | scale        |
|----------|--------|--------------------|---------------|--------------|
| entrance | 4      | short coupler cell | base + 0      | S            |
| body_1   | 5      | bulk (L - couplers)| base + 30     | S / sin(b0 d)|
| body_2   | 6      | bulk (L - couplers)| base + 90     | S / sin(b0 d)|
| exit     | 7      | short coupler cell | base + 0      | S            |

`sin(beta0 d) ~ 0.8657` (S-band, d = 3.5 cm; `beta0 d = 2*pi*f*d/c`). The rfdata carries **no**
R/tau/shunt impedance -- those are embedded in the per-section field scale already, so encoding them
in the field profile would double-count.

**Key field-reuse fact:** the rfdata Fourier reconstruction uses the fundamental period stored
*inside* the file (~0.105 m, the 3-cell block) as its wavelength -- NOT the lattice element `L`. The
element `L` only sets the active z-range `[zedge, zedge+L]` Impact-T integrates the periodic field
over, so a longer section is simply *more cells of the same per-cell field*. Setting the body
element `L` per section is all that "rescales length"; the field shape is reused unchanged.

## Frozen calibration (per-section field scale + crest phase)

The per-section `rf_field_scale` is **not analytic** (the `solrf` element has no scalar-gradient
input) and the chained-deck phase walk (drifts + finite beta shift the bunch arrival phase hundreds
of degrees per section) means `theta0 = 0` is on-crest only for the first section. The old stage
therefore searched, every run, for each section's local on-crest base phase (coarse phase scan +
parabolic refine) and fit its field scale to a ΔE target (`brentq`) -- expensive, plus a set of
§5 validation gates.

**Here that search is dropped.** The per-section `field_scale` [V/m] and `crest_phase_deg` were
**derived once** and are **hardcoded** in `config/linac4-8.yaml`, then applied directly to the deck:

- An `rf_field_scale` **ControlGroup** is created per section over its 4 solrf cells with factors
  `[1, 1/sin(b0 d), 1/sin(b0 d), 1]`, `absolute=True`. Its value `S` sets entrance/exit = `S`,
  body = `S/sin(b0 d)`, preserving the template body ratio. The group is set to the frozen
  `field_scale`. (The build-time element scales are also seeded, but the group apply is what makes
  the body ratio exact -- and since the group is `absolute=True` defaulting 0, it MUST be set or the
  deck silently runs zero-field.)
- Each section's absolute crest phase is pinned per solrf sub-element via `theta0_deg`
  (entrance +0, body_1 +30, body_2 +90, exit +0 on top of `crest_phase_deg`).

These setpoints were calibrated on-crest at the Balanced operating point (11 MW, sections 4-8); they
are an operating-point artifact, not first-principles values. There are **no per-run validation
gates**.

## Operating point & energy budget

`POWER_MW = 11` (the Balanced klystron point) for the whole linac. The per-section ΔE target
(recorded in the frozen-calibration table for the section-gains figure) is `sqrt(P_op/15)`-scaled
from the @15 MW `details.md` column:

```
dE_target,i(P_op) = dE_table,i * sqrt(P_op / 15)
```

| Sec | Type  | L (m) | dE@15MW | bore r (mm) |
|-----|-------|-------|---------|-------------|
| 4   | CU 5  | 4.97  | 51      | 14.7 -> 11.7 |
| 5   | CEA 4 | 5.15  | 55      | 14.7 -> 11.7 |
| 6   | CEA 5 | 5.15  | 55      | 14.7 -> 11.7 |
| 7   | CU 3  | 4.97  | 51      | 14.7 -> 11.7 |
| 8   | CU 4  | 4.97  | 51      | 14.7 -> 11.7 |

## Captured-core cut (handoff IN)

The sec3 exit dump (the 3->4 boundary) is already the relativistic captured core (sec1/2/3
dropped their slipping tails upstream), so the `MIN_KE_MEV` cut (default 12 MeV => beta > 0.99917)
now removes little -- it is kept as a **model-validity guard** (the rigid-crest no-slip TW model over
the ~29 m line is valid only for the relativistic core). The core is downsampled to `Np` (reweighted
to preserve core charge), drifted to mean `t`, and `z`-zeroed for Impact-T injection.

**Honest capture denominator:** `injection_summary.json` records `q_injected_C` = the **full** sec3
exit charge at the handoff (NOT the post-cut core), so within-stage capture (`q_out / q_injected`)
counts both the dropped tail and in-run loss. The tracked core charge is recorded separately
(`q_core_injected_C`).

## Space charge & quads

- **SC OFF** (`space_charge: false` => `Bcurr = 0`): the headline. Transverse SC is negligible at
  >25 MeV (`~ 1/gamma^2`, gamma > 49 at entry).
- **Quads OFF** (`K1 = 0`): the A->T (current->field) calibrations are undocumented. Each
  inter-section spacing is `gap/2` drift, a **real-length** zero-K1 quadrupole (its real tabulated
  drift-quad length, `quad_in`), `gap/2` drift, with `gap = DRIFT_M` (0.4 m). A K1=0 quad is
  optically a drift, but it keeps its **real length** so the cumulative path length -- and thus the
  bunch arrival time at sections 5-8 -- matches the deck the FROZEN ABSOLUTE crest phases were
  calibrated on (a shorter deck would shift the absolute `theta0` arrival phase by several RF periods
  per gap and throw sections 5-8 off-crest). Total lattice length ~28.97 m. The exploratory
  derived-FODO path of the old stage is dropped.
- **Transmission is a no-focusing LOWER BOUND, not a prediction.** With no quad focusing over the
  ~29 m line the beam diverges and a fraction scrapes the real tapered bore (`bore_aperture_on`,
  the solrf `radius` -- the binding aperture, deliberately not a widened numerical box). The robust,
  quad-independent deliverable is the **longitudinal physics** (exit `<KE>`, per-section ΔE), which
  does not depend on transverse confinement. Transmission is measured from the **macroparticle
  count** (`n_out / n_in`) BEFORE the openPMD charge re-imposition, so it can never be masked to 1.0.

## Output

`logs/diags/linac4-8/main/particles` -- the per-section and final beams as WarpX-layout openPMD
slices (sorted by `<z>`). `logs/diags/linac4-8/main/injection_summary.json` -- `q_injected_C`,
`q_core_injected_C`, `z_inject_lab_m` (the lab z of the 3->4 boundary, chained from the sec3
summary), `ke_in_mev`, `ke_out_mev`, `transmission_core`, `total_lattice_length_m`, `power_mw`, the
frozen calibration table, and `stat_vs_z` (`z_m`, `ke_mev`, `sigma_ke_mev`, `sigma_x_m`, `sigma_y_m`,
`norm_emit_x`, `norm_emit_y` from `I.stat`).

## Figures (`sim/plot/linac4-8.py` -> `logs/plots/linac4-8/`)

- `energy_gain` -- cumulative `<KE>` +/- sigma_KE vs z.
- `energy_spread` -- absolute sigma_KE and relative sigma_KE/`<KE>` vs z.
- `emittance` -- normalized emittance eps_n,x/eps_n,y vs z (the quads-OFF ~2.4x rise is a fort.10N
  diagnostic artifact, not physical growth).
- `section_gains` -- per-section achieved ΔE (from the vs-z KE curve) vs the frozen target ΔE.
- `fodo_optics` -- quads-OFF sigma_x / sigma_y vs z (placeholder optics, NOT predictive).

## Gotchas (Impact-T / lume-impact)

- **In-process, no subprocess** -- but `prepare_env()` (repo-root chdir + RLIMIT_NOFILE raise) still
  runs first.
- **Ntstep truncation reports `finished=True`** falsely (stops mid-line) -- `Ntstep` is sized from
  the lattice length; `mean_z_reached_m` is recorded so a truncated run is visible.
- **`autophase()` no-arg throws** on a non-cathode beam -- the crest phases are frozen, not searched.
- **`ParticleGroup.write()` is incompatible** with `openpmd-viewer` (openPMD 2.0 / STRING extension);
  the handoff OUT uses `sim.helpers.loadparticles.write_openpmd_particles` (WarpX byte layout).
- **Species name asymmetry**: `ParticleGroup.species` is `"electron"` (singular); the openPMD output
  and cross-stage readers key on `"electrons"` (plural). `loadparticles` translates.
- **No `write_beam` slice dumps** -- per-section vs-z evolution comes from `I.stat(...)` (written to
  the summary's `stat_vs_z`), not particle slices.
- **The frozen scale group is `absolute=True` defaulting 0** -- it MUST be set (via
  `_set_group_scale`) after `add_group`, or the deck silently runs zero-field.

---

## Positron mode -- accelerating the converter output

Sections 4-8 can be rewired to accelerate **positrons** from the converter that fills the 3->4
boundary (`docs/converter.md`) instead of the linac3 exit electrons. The same five-section deck and
the same frozen-calibration machinery are reused; only the input, the charge sign, and the phasing
change. This is an **additive** mode -- the default electron-mode content above is unchanged.

- **Input rewire (`config/linac4-8.yaml`).** The `io.sec3_particles` / `io.sec3_summary` paths are
  pointed at the converter output (`logs/diags/converter/main/particles` and its
  `injection_summary.json`) instead of the linac3 sec3 dump, so `upstream_exit_lab_z` continues the
  lab-z chain across the converter. `beam.species` is set to `positrons`. `beam.min_ke_mev` is
  **lowered** -- the converter's positron spectrum is soft (far below the relativistic electron core),
  so the electron-mode 12 MeV captured-core guard would discard most of the beam; the positron cut is
  set low enough to keep the physical spectrum.

- **Charge sign (`sim/linac4-8.py`).** The Impact-T deck header `Bcharge` is set to **`+1.0`** for
  positrons (it is `-1.0` for electrons). **`Bmass` is unchanged** -- the positron has the same rest
  mass as the electron, so only the charge sign flips. The openPMD handoff OUT is written as a
  **`positrons`** group with charge **+e** (via `loadparticles.write_openpmd_particles`, same as the
  electron path but with the positron species spelling).

- **The frozen crest phases DO NOT carry over.** The per-section `crest_phase_deg` are **absolute**
  Impact-T `theta0` values (referenced to the deck's t = 0), and they were derived for the **electron**
  beam. For positrons they are wrong on two counts: (1) flipping the charge sign flips which RF phase
  is **accelerating**, so the on-crest base phase moves by ~180 degrees; and (2) the positron
  **injection energy and velocity are much lower** than the electron core, so the bunch arrives at
  each section at a different phase and the chained-deck phase walk (drifts + finite beta) shifts the
  absolute `theta0` by additional large amounts per section. The two effects do not simply add to a
  clean 180 degrees -- **every section's crest must be re-derived**.

- **Re-derive before running.** Run **`sim/autophase_impact.py`** *after* the converter and *before*
  `sim/linac4-8.py` to re-derive the per-section absolute crests for the positron beam and rewrite
  them into `config/linac4-8.yaml`. Until that is done, a positron run with the electron-mode crest
  phases is meaningless (the beam is off-crest, decelerating in places). Even with the re-derived
  crests, the absolute transmission is a lower bound: the converter's high-divergence positrons need
  a capture optic at the target that is not modelled (see `docs/converter.md`). The robust positron
  deliverable, as in electron mode, is the **longitudinal** physics (per-section delta-E, exit
  `<KE>`), not the unfocused transmission.
