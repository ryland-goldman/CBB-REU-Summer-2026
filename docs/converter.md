# Positron converter target (G4beamline / Geant4)

The e+/e- converter that sits at the **4->5 boundary**, between the WarpX linac (sections 1-4) and
the Impact-T linac (sections 5-8). It reads the WarpX section-4 exit **electron** beam (the captured
relativistic core), drives it into a thin high-Z tungsten target, and produces e+/e- pairs by
bremsstrahlung -> pair production. The converted **positron** beam is the input to linac sections
5-8, which accelerate `q = +e`.

```
... linac1-4 (sec4, WarpX RZ) --[4->5 boundary]--> converter (this, G4beamline) --> linac5-8 (Impact-T, e+ mode)
                                  e- in                7 mm W target              e+ out
```

Unlike every other stage it is **not** a self-field/RF tracking run: it is a single-shot
Monte-Carlo particle-shower simulation through a solid radiator with **no field solve**.

```bash
conda activate CBB
python sim/converter.py          # reads the sec4 exit e-, runs g4bl, writes the e+ handoff
python sim/plot/converter.py     # figures from logs/diags/converter/main (sim must run first)
```

---

## Why G4beamline / Geant4, not WarpX/PIC

Every upstream stage is a particle-in-cell or beam-line tracking run: charged particles drift and
respond to self-consistent (WarpX) or applied (Impact-T) electromagnetic fields. **A converter is a
different physics problem.** The energy loss that matters here is **radiative** — a relativistic
electron traversing a high-Z nucleus emits a bremsstrahlung photon, and that photon converts to an
e+/e- pair in the field of another nucleus. Both processes are QED interactions with the target
nuclei (with a hadronic tail at these energies), governed by cross-sections and a stochastic
electromagnetic shower, not by macroscopic E/B fields on a mesh. A PIC code has no model for them.

Geant4 — driven here through **G4beamline** (`g4bl`), its beam-line front end — carries exactly this
physics: a reference physics list assigns each particle/energy the appropriate processes (multiple
Coulomb scattering, ionisation, bremsstrahlung, pair production, annihilation, photoelectric/Compton)
and tracks the resulting shower through the material. The run therefore has **no Poisson/Maxwell
solve at all**: the target is field-free, space charge is irrelevant during the sub-ns transit, and
the only "lattice" is the geometry of the radiator. This is why the stage is an external `g4bl`
process rather than another WarpX/Impact-T driver.

---

## Target design — a thin tungsten radiator

The converter is a single solid cylinder of **tungsten** (`G4_W`, the NIST material), **7 mm thick**,
**radius 10 mm**, with the incident electron beam on-axis. Tungsten is the canonical converter
material: its high atomic number (Z = 74) gives a short radiation length and a large pair-production
cross-section, so the bremsstrahlung -> pair-production cascade develops in a millimetre-scale depth,
and its high density and melting point survive the deposited power.

The mechanism in the radiator is a two-step electromagnetic shower:

1. **Bremsstrahlung.** The incident e- radiates hard photons in the Coulomb field of the W nuclei.
   The radiated spectrum hardens with Z and with target thickness (more nuclei to scatter off).
2. **Pair production.** Those photons convert to e+/e- pairs in the same nuclear field. The
   positrons are the product of interest; the electrons (primaries that survive plus shower e-) and
   the un-converted photons are also tracked and counted but discarded at the handoff.

**Why ~7 mm.** Target thickness is a trade. Too thin and few photons are produced and fewer convert,
so the positron yield is low. Too thick and the shower over-develops: the positrons that *are*
produced re-scatter, lose energy, and spread in angle before they can leave the back face, and the
absorbed dose climbs. The yield-vs-thickness curve for a thin converter peaks around **~2 radiation
lengths** of the radiator (X0 ~ 3.5 mm in tungsten, so ~7 mm), which is the qualitative basis for the
7 mm choice here. (The actual optimum depends on the downstream capture acceptance; this stage fixes
the geometry and reports the yield it gives rather than optimising it — see the caveats.)

---

## Capture solenoid (the `solenoid:` block)

The positrons leaving the back face are **high-divergence** — a converter sprays them over a large
solid angle — so a real e+ source immerses the target exit in a strong axial magnetic field to
re-collect them into the downstream acceptance (the **capture optic**: a high-field solenoid, or an
adiabatic matching device). This stage models it as a **real current coil** of default **5 T peak
over 1 m**, starting at the target front face so the field surrounds the radiator and the drift to
the sampling plane.

It is realised in the g4bl deck as a `coil` + `solenoid`, so g4bl computes the **full Maxwellian
field**, including the **end fringe** — the radial `B_r = -(r/2)·∂B_z/∂z` where the field ramps down
at each coil end. That fringe is what gives a solenoid its lens action: `B_r` crossed with the
azimuthal velocity is the inward focusing kick. (A hard-edge uniform-`B_z` element would have *no*
`B_r` and so no fringe focusing at all — it would bound the beam *radius* but leave the divergence
untouched.) The conductor **current density is auto-solved** (`coil_current_density`, the exact
thick-solenoid centre formula, matched to g4bl's coil to < 1e-5) so the **central (peak)** field
equals `b_tesla`; the field sags to ~half at each coil end and fringes out beyond. The peak field,
length, coil bore/winding radii, upstream start, and the post-coil drift are all set in
`config/converter.yaml`; `enabled: false` reverts to the field-free run.

The sampling plane sits `exit_drift_mm` **past the coil exit**, in field-free space, so the exit
fringe has fully acted before the handoff. There is an inherent **divergence-vs-size trade**: the
fringe converts the positrons' large divergence into a more parallel beam (emittance is conserved,
not reduced), so as the now-parallel beam drifts past the coil it grows in radius. A longer
`exit_drift_mm` gives a more field-free handoff but a larger beam. The coil bounds the captured
phase space and supplies the fringe focusing; it does **not** reduce the emittance, so it makes no
claim to an absolute capture efficiency.

---

## Precision settings (the `physics:` block)

Three knobs in `config/converter.yaml`'s `physics:` block set the Geant4 simulation accuracy. Each
trades fidelity against runtime, and these three are the ones that matter for a thin EM converter:

- **`physics_list: QGSP_BERT_EMZ`** — the Geant4 reference physics list. The `QGSP_BERT` base is the
  standard hadronic list; the **`EMZ`** suffix selects **EM Option 4**, Geant4's high-accuracy
  electromagnetic models for e-/e+/gamma (more precise multiple-scattering, bremsstrahlung-angular,
  and pair-production sampling than the default `EMV`/`EMZ` opt0/3 lists). For a converter the EM
  shower *is* the physics, so the high-accuracy EM models are the right choice despite being slower;
  the hadronic base costs little here (few hadronic interactions at tens of MeV) but is kept for
  correctness.

- **`min_range_cut`** — the Geant4 **production / range cut**. Geant4 expresses secondary-production
  thresholds as a *range* (a distance in the material), not an energy: a secondary whose range would
  be below the cut is **not created as a track** — its energy is deposited locally instead. A small
  cut produces (and tracks) soft secondaries down to a fine scale: more accurate low-energy shower
  detail, more tracks, slower. A large cut suppresses soft secondaries: faster, but it can bias the
  low-energy tail of the positron spectrum. The setpoint is chosen fine enough to resolve the
  positrons that fall in the downstream capture band without tracking the dose-only soft shower.

- **`max_step`** — the maximum Geant4 step length **inside the target**. Capping the step forces
  Geant4 to re-evaluate the shower (energy loss, scattering, secondary production) at a fine depth
  resolution through the radiator, rather than taking one long step that smears the shower's
  longitudinal development. Smaller `max_step` = better-resolved shower depth and exit phase space,
  more steps, slower. It is set to a small fraction of the 7 mm thickness so the cascade is resolved.

Together these are the standard accuracy levers for a converter: **what EM physics** (`physics_list`),
**how soft a secondary to keep** (`min_range_cut`), and **how finely to step the shower**
(`max_step`). Nothing else in the run needs tuning — there is no field solver, no time step.

---

## Yield bookkeeping — one event per incident electron

The yield is the physical quantity the stage exists to deliver, so the charge bookkeeping is kept
exact rather than renormalised:

- Each incident electron is **one g4bl event** (one `EventID`). The driver writes the section-4 exit
  beam as a BLTrackFile with `EventID = 1..N` (one per macroparticle) and runs g4bl with one primary
  per event.
- Every secondary g4bl produces **inherits its primary's `EventID`**. So each output positron knows
  which incident electron produced it.
- The driver maps each output positron back to its incident electron's macroparticle **charge**, so
  the handoff weight is `q_positron = q_incident x yield` — the per-electron positron multiplicity is
  carried through as charge. The handoff therefore preserves the **physical positron yield** (how
  much positron charge a given incident electron charge produced) rather than re-normalising the
  output to a fixed total or to one unit weight per macroparticle.

This is why the **original** BLTrackFile is used (see Gotchas): the alternate g4bl outputs do not
carry the exact `EventID` the back-mapping keys on.

---

## Units and the data handoff

g4bl's `BLTrackFile` ASCII format is **mm / MeV-c / ns**, one row per track, with a `PDGid` column
(`e- = 11`, `e+ = -11`, `gamma = 22`) plus `EventID`/`TrackID`/`ParentID`/`Weight`. The repo's
internal beam convention is openPMD / `pmd_beamphysics.ParticleGroup`: **m**, momentum in **eV/c**
(`= gamma*beta*MC2_EV`), time in **s**, and per-macroparticle **weight = charge [C]**.

The conversion both ways lives in **`sim/helpers/g4bl.py`**:

- `read_bltrackfile` / `write_bltrackfile` — pandas ASCII IO (`comment="#"` skips g4bl's header;
  integer id columns stay int so `PDGid` is `11`, not `11.0`).
- `particlegroup_to_bltrack_df` — the electron input: ParticleGroup -> BLTrackFile DataFrame, tagging
  `PDGid = 11` and `EventID = 1..N` (the incident-electron identity the yield map keys on).
- `bltrack_df_to_particlegroup` — the positron output: BLTrackFile DataFrame -> ParticleGroup, with
  `species = "positron"` and the per-row charge `weight_C` the driver computed from the yield map.
  ParticleGroup derives the signed charge/mass (+e) from the species name.
- `bltrack_ke_mev` — kinetic energy [MeV] from the file momenta (electron/positron rest mass).

The driver selects the **positron** rows (`PDGid == -11`) from g4bl's output, converts them, and
writes the openPMD group as **`positrons`** with charge **+e**. (The electron and photon rows are
read only for the yield-bar figure and the summary counts.)

---

## Lab-z chaining

The converter sits at a known lab z (the WarpX section-4 exit, the 4->5 boundary), and the downstream Impact-T
stage continues the lab-frame z chain via `upstream_exit_lab_z`. The converter therefore records
`z_inject_lab_m` / `z_inject_mean_m` in its `injection_summary.json` the same way the linac sections
do: the lab z at which its output beam sits, so `sim/linac5-8.py`'s `upstream_exit_lab_z` chain
continues unbroken across the converter. The g4bl run works in **local mm coordinates** (the target
front face at the local origin); the driver maps the local-frame exit z back into the converter's
place in the lab frame for the summary.

---

## Output

`logs/diags/converter/main/particles` — the converted **positron** beam as an openPMD `positrons`
group (charge +e), in the repo's WarpX-compatible layout, ready for `sim/linac5-8.py`.

`logs/diags/converter/main/injection_summary.json` — the run's bookkeeping. It holds the **target
geometry** (material, thickness, radius), the **precision settings** (physics list, `min_range_cut`,
`max_step`), the **per-species yields** (e+/e-/gamma counts and per-incident-electron yields), the
**positron charge** (the back-mapped `q_positron`), and the positron beam's **KE and divergence
statistics**, plus the `z_inject_lab_m`/`z_inject_mean_m` lab-z chain fields. (Field *names*
documented here; the *values* are run output and are not quoted.)

---

## Figures (`sim/plot/converter.py` -> `logs/plots/converter/`)

- **positron energy spectrum** — the e+ kinetic-energy distribution leaving the target (the soft,
  broad spectrum a converter produces).
- **e+/e-/gamma yield bars** — the per-species output counts (or per-incident-electron yields): how
  much of each the shower produced.
- **positron transverse divergence / r-p_r** — the e+ angular phase space (a converter beam is
  high-divergence; this is what a capture optic would have to accept).
- **longitudinal z-KE** — the positron longitudinal phase space leaving the target.

---

## Gotchas / caveats

- **PDGid sign convention.** `e+ = -11`, `e- = +11` (positron is the *negative* PDG id). The driver
  selects positrons by `PDGid == -11`; getting the sign backwards silently hands the **electrons**
  to linac5-8.
- **Monte-Carlo determinism is conditional.** g4bl seeds its RNG per `EventID`, so the run is
  reproducible **only** with a fixed deck *and* a fixed input BLTrackFile (same events, same order).
  Any change to the upstream electron beam or to the physics block changes the shower realisation.
- **Use the original BLTrackFile.** g4bl can also emit `BLTrackFile2`/ROOT outputs; the original
  ASCII `BLTrackFile` is used because it carries the exact `EventID` the yield back-mapping needs
  (the alternates renumber or reformat it).
- **The positron beam is not linac-ready as-is.** The converted positrons are **low-energy,
  high-divergence, and large-energy-spread** — a soft shower product, nothing like the relativistic
  captured core the electron linac5-8 was tuned for. Two consequences:
  - **linac5-8's crest phases must be re-derived for positrons.** The frozen per-section
    `crest_phase_deg` were calibrated on the electron beam; for positrons they are wrong (the charge
    sign flips the accelerating phase, and the much lower injection energy/velocity shifts the
    arrival phase further). Run **`sim/autophase_impact.py`** *after* the converter and *before*
    `sim/linac5-8.py` to re-derive them. (See the "Positron mode" section of `docs/linac5-8.md`.)
  - **Capture optic.** A real positron source places a **capture optic** (an adiabatic matching
    device / high-field solenoid) right at the target exit to collect the divergent positrons into
    the linac acceptance. This stage models it as a **uniform axial capture solenoid** (the
    `solenoid:` block — default 5 T over 1 m, starting at the target front face) that focuses the
    positrons before the sampling plane; set `enabled: false` to recover the bare field-free run. It
    is a uniform-field idealisation, not a tapered AMD, so the accepted transmission/yield through
    linac5-8 are still qualitative — the deliverable is the converter physics (spectrum, yield,
    divergence), not an absolute accepted-positron number.
- **g4bl is an external binary.** G4beamline 3.08 is invoked as a separate process (`g4bl`), not a
  pip/conda dependency — it must be installed and on `PATH` independently of the `CBB` environment.
