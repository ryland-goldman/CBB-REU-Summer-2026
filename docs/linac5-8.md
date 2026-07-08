# Linac sections 5-8 (Impact-T)

The rest of the line to CHESS, after the WarpX linac sections 1-4 and the e+/e- converter target:
four S-band (2856 MHz) traveling-wave (TW) accelerating sections (CEA 4/5 + CU 3/4), chained into
**one** Impact-T deck and integrated as one time-ordered beam. Reads the **converter positron beam**
(produced after the WarpX section 4) and accelerates the captured positron core on-crest through the
four sections.

```
... linac1-4 (sec1/2/3/4, WarpX RZ) --[4->5 boundary]--> converter (G4beamline) --> linac5-8 (this, Impact-T, 4 TW sections, e+ mode)
```

The converter sits at the **4->5 boundary** (after the WarpX section 4): it drives the section-4
exit electrons into a tungsten target and hands the resulting positron beam to this stage.

## Why Impact-T, not WarpX

The upstream stages are WarpX/pywarpx runs. This stage is an external serial Impact-T run
(`ImpactTexe`) driven through **lume-impact**, in-process (no pywarpx global-geometry binding, so no
per-stage subprocess isolation is needed). Impact-T integrates one beam through one time-ordered
lattice -- the natural fit for the remaining sections that BMAD/LinacSim also treat as one generic
linac.

The run executes in `logs/diags/linac5-8/` (`use_temp_dir=False`, `workdir=` set), so its `fort.18`
lands at a known `logs/diags/linac5-8/fort.18`. The progress bar (`sim.helpers.tqdmwrapper.
impact_progress`) polls that file's column 1 (reference `z` [m]) and advances 0 -> total lattice
length while `I.run()` executes.

## Field model -- generic constant-gradient TW, no field maps

Sections 5-8 have **no GPT/CST field maps** (none exist; LinacSim/BMAD model them with the generic
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
**derived once** and are **hardcoded** in `config/linac5-8.yaml`, then applied directly to the deck:

- An `rf_field_scale` **ControlGroup** is created per section over its 4 solrf cells with factors
  `[1, 1/sin(b0 d), 1/sin(b0 d), 1]`, `absolute=True`. Its value `S` sets entrance/exit = `S`,
  body = `S/sin(b0 d)`, preserving the template body ratio. The group is set to the frozen
  `field_scale`. (The build-time element scales are also seeded, but the group apply is what makes
  the body ratio exact -- and since the group is `absolute=True` defaulting 0, it MUST be set or the
  deck silently runs zero-field.)
- Each section's absolute crest phase is pinned per solrf sub-element via `theta0_deg`
  (entrance +0, body_1 +30, body_2 +90, exit +0 on top of `crest_phase_deg`).

These setpoints were calibrated on-crest at the Fromowitz operating point (17 MW, the sec 5/6
forward power); they are an operating-point artifact, not first-principles values. There are **no
per-run validation gates**.

> **The crests are autophased in the chain.** `sim/main.py` runs `sim/autophase_impact.py` before
> this stage, re-deriving every section's `crest_phase_deg` for the positron beam on the
> section-5-start deck and rewriting the YAML. It is a **numerical** model (~30 s of pure numpy, no
> Impact-T launched): it RK4-integrates the pencil-ised positron core through the exact on-axis Ez it
> reconstructs from the rfdata4-7 Fourier shapes. The shipped values are only the seed for a
> standalone run; if you run `sim/linac5-8.py` alone on a changed upstream beam, run
> `sim/autophase_impact.py` first (see *Positron mode*).

## Operating point & energy budget

`POWER_MW = 17` (the Fromowitz sec 5/6 forward power) for the whole linac. The per-section ΔE target
(recorded in the frozen-calibration table for the section-gains figure) is `sqrt(P_op/15)`-scaled
from the @15 MW `details.md` column:

```
dE_target,i(P_op) = dE_table,i * sqrt(P_op / 15)
```

| Sec | Type  | L (m) | dE@15MW | bore r (mm) |
|-----|-------|-------|---------|-------------|
| 5   | CEA 4 | 5.15  | 55      | 14.7 -> 11.7 |
| 6   | CEA 5 | 5.15  | 55      | 14.7 -> 11.7 |
| 7   | CU 3  | 4.97  | 51      | 14.7 -> 11.7 |
| 8   | CU 4  | 4.97  | 51      | 14.7 -> 11.7 |

(Section 4 — CU 5 — moved upstream into the WarpX stage `sim/linac1-4.py`; see `docs/linac1-4.md`.)

## Captured-core cut (handoff IN)

The input is the **converter positron beam** (a soft, divergent shower product, NOT a relativistic
captured core). The `MIN_KE_MEV` cut (default 2 MeV for positrons) keeps the usable core; the
rigid-crest no-slip TW model is only valid for that relativistic core. The core is matched to `Np`
(reweighted to preserve core charge), drifted to mean `t`, and `z`-zeroed for Impact-T injection.

**Upsampling for statistics (`beam.upsample`).** The converter positron yield caps the core macro
count (~11 k), and capture is ~1%, so a plain run ends with ~100 survivors -- too few for
downstream emittance/energy statistics. Because **space charge is OFF**, each macroparticle is
dynamically independent, so the core is *smeared-upsampled* to `Np`
(`loadparticles.upsample_smeared`): bootstrap-draw parents, then jitter each clone by
`upsample_smear` x its k-NN distance in the **transverse** phase space (x, y, px, py). Plain
duplication would not help -- coincident clones track identically under the deterministic SC-off
optics; the local smear decorrelates them so they sample the local density. pz/z/t are carried from
the parent, so the energy spectrum is preserved to ~0.5% (not exact -- total energy still depends on
the smeared px/py), and the `MIN_KE_MEV` floor is **re-imposed** after the smear (a low-pz,
large-angle parent kept only by its transverse momentum can otherwise leak a clone below the floor). The survival fraction and the survivor moments are preserved; only the sampling
density rises, so `Np` = 120 k injected -> ~1 k survivors. The genuine resolution is still set by
the real converter-core count, recorded as `n_core_raw` in the summary.

**Honest capture denominator:** `injection_summary.json` records `q_injected_C` = the **full**
converter positron-beam charge at the handoff (NOT the post-cut core), so within-stage capture
(`q_out / q_injected`) counts both the dropped tail and in-run loss. The tracked core charge is
recorded separately (`q_core_injected_C`).

## Space charge, quads & capture solenoids

- **SC OFF** (`space_charge: false` => `Bcurr = 0`): the headline. Transverse SC is negligible at
  >25 MeV (`~ 1/gamma^2`, gamma > 49 at entry).
- **Capture-optics focusing ON -- the real CESR optics.** Both the layout and the strengths come
  from the CLASSE BMAD deck (`section_5_8_layout.bmad`, kd324's 2021 linacsim working copy at
  `../Cornell/reference/classe-docs/nfs_acc_user_kd324_documents_mycesr_linacsim/BMAD/`) -- Dan Fromowitz's
  dissertation lattice updated to the machine operating point:
  - **Capture solenoids on sections 5 & 6** (`solenoid_b_tesla`, 0.243 T flat-top). Sections 7 & 8
    have none (cavity 7 is the first section without a solenoid). The deck superimposes three
    coils per section (the `pos_sol_A/B/C.bmad` rotationally-symmetric r-z field grids at their
    264 A supply setpoint); the extracted on-axis Bz is a **+0.2429 T flat-top over ~5.04 m**
    falling to ~27-33% of peak at the section flanges, plus end-trim coils (0.054 / 0.070 T) -- the C trim coil sits at 5.18-5.27 m from the cavity
    start, i.e. **in the exit line**, and the fringe decays across the inter-section gap (0.117 T
    at the cavity start, 0.187 T at the 5.15 m active exit, ~0.007 T at the next cavity). That
    measured profile is committed as `fieldmaps/rfdata/pos_sol_onaxis_264A.txt` and used
    directly: each solenoid is a **solenoid-only `solrf`** (`rf_field_scale = 0`, static
    `solenoid_field_scale` Bz) **overlapping the cavity AND its exit line** (`L = cavity + gap`,
    not advancing the deck z, so the crest geometry is unchanged). The element windows tile, and
    each window sums every string's profile (the neighbour's backward fringe included), so
    nothing is double-counted. Fourier-decomposed (n=120) for Impact-T's static type-105
    paraxial expansion; the sec5 window's periodic-wrap ends nearly match (0.117 / 0.124 T),
    while the sec6 window wraps with a ~0.12 T mismatch whose Gibbs ringing stays within ~5 cm
    of the window ends. `solenoid_b_tesla` scales the profile normalised to the 264 A flat-top.
    The sections' values **track the converter capture field** (config `solenoid_tracking`): each
    `solenoid_b_tesla` is its value at the converter reference field `conv_b_ref` (the CONV_HV
    3300 A setpoint the 264 A string is paired with), and the driver multiplies it by
    `converter solenoid.b_tesla / conv_b_ref` — so retuning or scanning the converter capture
    field scales the whole capture line together (converter block absent or disabled → ratio 1).
    The applied fields are recorded in the summary (`sec56_solenoid_b_tesla`,
    `conv_field_ratio`).
  - **Real exit-optics lines between sections** (config `exit_optics`): the deck's drift/quad line
    after each cavity, element lengths and order verbatim, with gradients from the deck's
    calibrated CU overlay (`grad = polarity * CU * (T/m per A) * (A per CU)`) at the machine
    setpoint. Section 5's Q5 doublet is ~0 (7 / 21 CU -- the solenoid focuses there); sections
    6 and 7 end in **symmetric QH-QV-QH triplets** (-3.86/+3.71/-3.86 and -4.74/+4.62/-4.74 T/m);
    section 8's trailing Q8 doublet (-1.34/+1.15 T/m) **is placed** and ends the modelled line
    (the positron snout continues from there). A per-section `quad_scale` (optimizer knob, 1.0 =
    machine setpoint) multiplies that section's gradients. The exit lines carry the real total
    length, so the cumulative path length -- and thus the bunch arrival time at sections 6-8 --
    matches the deck the ABSOLUTE crest phases are calibrated on (a changed spacing throws
    sections 6-8 off-crest: re-run `sim/autophase_impact.py` after any geometry edit). The
    deck's inter-section pipe radius is recorded per gap (`pipe_radius_m`, 20 / 25.4 mm) and set
    as the element radius -- Impact-T honours it on quads but not on drifts, and the 20 mm
    `xyrad_m` computational box wall remains the binding gap aperture. The drift lengths are the
    deck's flange-frame values (the BMAD sections are flange-to-flange, slightly longer than the
    active lengths modelled here); the 0.129 m D_CONV drift ahead of section 5 is owned by the
    converter stage, not this deck.
- **Transmission folds the real capture acceptance.** The beam scrapes the real tapered bore
  (`bore_aperture_on`, the solrf `radius`); the sec5/6 solenoids + quads collect what the bore
  accepts. Per the thesis the capture is intrinsically inefficient (a few %), so transmission stays
  low. It is measured from the **macroparticle count** (`n_out / n_in`) BEFORE the openPMD charge
  re-imposition, so it can never be masked to 1.0. The taper is per sub-element, not per section:
  `_section_subelements` (`sim/linac5-8.py`) narrows the bore radius linearly from the section's
  entrance radius to its (smaller) exit radius (config `bore_in` diameters, converted to radii), and
  samples each of the 4 solrf sub-elements (entrance/body_1/body_2/exit) at its own z along that
  taper -- entrance at `r_in`, exit at `r_exit`, the two body cells at the taper's midspan -- so the
  exit sub-element (the smallest radius) is the binding aperture within a section, not just the
  overall `xyrad_m` computational-box wall.

## Output

`logs/diags/linac5-8/main/particles` -- the per-section and final beams as WarpX-layout openPMD
slices (sorted by `<z>`). `logs/diags/linac5-8/main/injection_summary.json` -- `q_injected_C`,
`q_core_injected_C`, `z_inject_lab_m` (the lab z of the converter handoff, chained from the converter
summary), `ke_in_mev`, `ke_out_mev`, `transmission_core`, `total_lattice_length_m`, `power_mw`, the
frozen calibration table, and `stat_vs_z` (`z_m`, `ke_mev`, `sigma_ke_mev`, `sigma_x_m`, `sigma_y_m`,
`norm_emit_x`, `norm_emit_y` from `I.stat`).

## Figures (`sim/plot/linac5-8.py` -> `logs/plots/linac5-8/`)

- `energy_spectrum` -- exit KE spectrum.
- `evolution_vs_z` -- mean KE, norm. emittance eps_n,x, sigma_x, and surviving charge vs z.
- `energy_spread` -- absolute sigma_KE and relative sigma_KE/`<KE>` vs z.
- `phase_space_z_KE` -- exit longitudinal phase space (z - <z> vs KE).
- `transverse_r_pr` -- exit transverse phase space (r, p_r).
- `section_gains` -- per-section achieved dE (from the vs-z KE curve) vs the frozen target dE. The
  section-5 bar reads low by construction: dKE is differenced against the whole injected core's
  mean while capture scraping removes the low-energy tail inside section 5 -- a capture artifact,
  not an off-crest or field-scale deficiency (sections 6-8 land on target).

The `evolution_vs_z`/`energy_spread` curves are read from the summary's `stat_vs_z` table (`I.stat`)
when present; a sparser particle-slice fallback (from the openPMD dumps) covers legacy summaries that
predate that table, including one that predates the surviving-charge column. `section_gains` attributes
per-section achieved dE from the vs-z KE curve using the calibration table's `z_entry_m`/`z_exit_m`
(the real deck geometry `sim/linac5-8.py` writes); a legacy summary lacking those falls back to an
even split of the z-grid, which coarsely mis-attributes gain across the inter-section drifts.

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

This stage runs the **positron** beam from the converter (`docs/converter.md`), which sits after the
WarpX section 4. The converted positrons (low-energy, high-divergence, large energy spread) are the
input; the deck accelerates the captured positron core on-crest. The charge sign and the phasing
follow from the positron species.

- **Input (`config/linac5-8.yaml`).** The `io.conv_particles` / `io.conv_summary` paths point at the
  converter output (`logs/diags/converter/main/particles` and its `injection_summary.json`), so
  `upstream_exit_lab_z` continues the lab-z chain across the converter. `beam.species` is `positrons`.
  `beam.min_ke_mev` is **low** (2 MeV) -- the converter's positron spectrum is soft, so a high
  captured-core guard would discard most of the beam; the cut is set low enough to keep the physical
  spectrum.

- **Charge sign (`sim/linac5-8.py`).** The Impact-T deck header `Bcharge` is set to **`+1.0`** for
  positrons (it is `-1.0` for electrons). **`Bmass` is unchanged** -- the positron has the same rest
  mass as the electron, so only the charge sign flips. The openPMD handoff OUT is written as a
  **`positrons`** group with charge **+e** (via `loadparticles.write_openpmd_particles`, same as the
  electron path but with the positron species spelling).

- **The crests are species- and geometry-dependent — autophased each run.** The per-section
  `crest_phase_deg` are **absolute** Impact-T `theta0` values (referenced to the deck's t = 0), so
  they depend on (1) the charge sign — for positrons the accelerating phase moves ~180° from
  electrons; (2) the positron injection energy/velocity, which sets the chained-deck phase walk
  (drifts + finite beta) across sections; and (3) the deck geometry — this deck now **starts at
  section 5** (section 4 is a WarpX stage upstream of the converter). `sim/main.py` therefore runs
  **`sim/autophase_impact.py`** *after* the converter and *before* `sim/linac5-8.py` to re-derive
  every section's crest for the positron beam and rewrite `config/linac5-8.yaml`.

- **How the autophase works (numerical, no Impact-T).** `sim/autophase_impact.py` is a 1D
  longitudinal model — the analog of the WarpX `sim/autophase.py`. It reconstructs each TW section's
  on-axis Ez from the rfdata4-7 Fourier shapes (matching `impact.fieldmaps.ele_field` to float noise:
  the 4-line +0/+30/+90/+0 superposition, the 1/sin(β₀d) body scale, the absolute `cos(2πft + θ0)`
  phase), places the sections at the real deck z (sections + the real exit-optics line lengths, so the
  cumulative arrival time the absolute θ0 depends on is correct), pencil-ises the positron core
  (momentum onto +z, on-axis — the bare ~600 mrad divergence would otherwise scrape before any
  section ends), and RK4-integrates the whole bunch through a phase scan (coarse → fine → parabolic),
  pinning earlier sections to their found crest. ~30 s of numpy versus the old minutes-to-tens-of-
  minutes Impact-T scan, with no physics lost (the Impact-T scan was already pencil/SC-off/quads-off).

- **Crest vs reported energy.** The numeric crests match Impact-T's (re-derived values agree with the
  previously Impact-T-derived YAML crests to ~1°). The per-section `⟨KE⟩` the tool prints is a **model
  energy used only to locate the crest** — its magnitude is ~1.7× Impact-T's actual gain, a *constant*
  offset between `ele_field` and the Fortran Impact-T field. A constant field-scale offset does not
  move the argmax over phase, so the crest is unaffected; the stage energy is what `sim/linac5-8.py`
  reports, not this model number. `sim/autophase_impact.py --verify` confirms the crest directly with
  a tight 3-point Impact-T phase scan (crest and crest ± `VERIFY_DELTA_DEG`) with the aperture/pipe
  scrape disabled, which is immune to the amplitude offset.

- **Standalone runs.** If you invoke `sim/linac5-8.py` alone on a changed upstream beam, run
  `sim/autophase_impact.py` first — otherwise the shipped crest seeds may be off-crest (decelerating
  in places). The sec5/6 capture solenoids + exit-line quads supply the post-target capture optics,
  but the capture stays intrinsically low-efficiency (a few %, per the thesis), so absolute
  transmission is small; the robust positron deliverable is the **longitudinal** physics
  (per-section delta-E, exit `<KE>`).
