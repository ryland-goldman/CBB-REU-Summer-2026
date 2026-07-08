# Cornell Linac Beam Simulation

Beam-dynamics simulation of the **Cornell High Energy Synchrotron Source (CHESS)** electron source
front end — a from-first-principles rebuild of Adam Bartnik's
[LinacSim](https://cesrwww.lepp.cornell.edu/wiki/CESR/LinacSim) **cathode → gun → injector → linac**
chain. The first seven stages run in [WarpX](https://warpx.readthedocs.io) (the particle-in-cell
code, via `pywarpx` / lume-warpx); the converter runs in
[G4beamline](http://www.muonsinc.com/) (Geant4) and the final stage in
[Impact-T](https://github.com/impact-lbl/IMPACT-T) (via lume-impact). Each stage reads the previous
stage's openPMD beam, so the stages form one self-consistent accelerator chain. The **4→5 boundary**
(the `linac4` exit) is the slot for the e⁺/e⁻ converter target.

```
cathode ─► gun ─► injector ─► linac1 ─► linac2 ─► linac3 ─► linac4 ─►[4→5]─► converter ─► linac5-8
SCL       CESR    2 prebunchers  SLAC 3 m   CEA 2     CEA 3    CU 5    e+/e-    W target   CU/CEA S-band
diode     gun     + 6 solenoids  TW capture                            target              (Impact-T)
(2D)      (RZ)    (RZ)           (RZ)       (RZ)      (RZ)     (RZ)    (G4bl)               (sections 5–8)
```

## Layout

```
config/    one YAML per stage — every tunable option is hardcoded here (edit to retune)
sim/       main.py + one driver per stage (cathode, gun, injector, linac1-4, linac5-8)
  helpers/ stage-agnostic plumbing: tools, buildfields, loadparticles, metrics, tqdmwrapper
  plot/    one plotter per stage + common (shared figures) + chain (cross-stage)
logs/      diags/<stage>/ (openPMD) · plots/<stage>/ (PNG) · pipeline/log_<date>.log
docs/      per-stage physics notes
fieldmaps/ gdf/ (GPT field-map inputs) · h5/ (built openPMD maps) · rfdata/ (Impact-T templates)
```

The four WarpX linac sections share **one** driver (`sim/linac1-4.py`, section chosen by a CLI
argument); the four Impact-T sections share `sim/linac5-8.py`.

## Setup

All simulations run in the **CBB** conda environment:

```bash
conda activate CBB                 # Miniforge at ~/miniforge3
pip install -r requirements.txt    # warpx / impact-t / openpmd-api are best from conda-forge
```

Field maps are read from `fieldmaps/gdf/` (committed). To *regenerate* a `.gdf` map with GPT (not
needed to run the chain), use the `gpt_remote` wrapper, which runs GPT on a remote host.

## Run

```bash
python sim/main.py                 # the whole chain, in order, with a final beam summary
```

Each stage is a self-contained script you can run alone (from the repo root):

```bash
python sim/cathode.py              # one stage's simulation
python sim/plot/cathode.py         # its figures (from existing diagnostics)
python sim/linac1-4.py 2           # linac section 2 (argument selects the section)
python sim/plot/linac1-4.py 2      # its figures (from existing diagnostics)
```

The WarpX stages each run in a fresh subprocess (pywarpx binds one geometry per interpreter);
Impact-T runs the same way. Subprocess output is captured to `logs/pipeline/log_<date>.log` while
progress bars stay on the terminal.

## GUI

```bash
python sim/gui.py                  # desktop control panel, from the repo root in the CBB env
```

`sim/gui.py` is a Tk desktop control panel for the chain, three panes in one window: a left column
of per-stage pipeline cards (cathode … linac 5-8), each with Edit Config / Run / Plot buttons and,
for the linac stages, an Autophase toggle that runs the same `sim/autophase*.py` pre-step
`sim/main.py` drives before that stage's sim; a right notebook with a Beam Explorer tab (Trends /
1D / 2D views over a stage's existing openPMD dumps, ordered by dump ⟨z⟩) and a Plots tab (browses
the PNGs already written under `logs/plots/<stage>/`); and a bottom console mirroring every
launched subprocess's stdout/stderr. It only reads existing diagnostics/figures on disk — it shells
out to the same `sim/plot`/`sim/autophase*.py` scripts `sim/main.py` uses, and does not simulate
anything itself.

## Stages

| Stage | Driver / config | What it does |
|-------|-----------------|--------------|
| **Cathode** | `sim/cathode.py` · `config/cathode.yaml` | Finite, space-charge-limited (Child–Langmuir) thermionic diode in 2D x–z. The electron source. |
| **Gun** | `sim/gun.py` · `config/gun.yaml` | CESR electrostatic gun in RZ from the `CESR_gun.gdf` field map, with the relativistic EMS self-field. Timed beam release; writes a field-free-pad exit handoff. |
| **Injector** | `sim/injector.py` · `config/injector.yaml` | The LinacSim injector in one RZ space-charge run: two 214 MHz prebunchers (velocity bunching) + six solenoid lenses (focusing) + the 9.547 mm iris, handing off near z ≈ 2.03 m. |
| **Linac 1–4** | `sim/linac1-4.py` · `config/linac{1,2,3,4}.yaml` | Four SLAC-design 3 m, 2π/3 traveling-wave sections (RZ, WarpX) reusing the SLAC field maps. Section 1 captures the injector beam through the iris; sections 2–4 accelerate the captured core. Section 4's exit is the 4→5 boundary. |
| **Converter** | `sim/converter.py` · `config/converter.yaml` | e⁺/e⁻ converter target (G4beamline/Geant4): drives the section-4 exit electrons into a 6.35 mm tungsten radiator (brems → pair production) and hands the resulting positron beam to the Impact-T linac, focused by the real CONV_HV capture-solenoid fieldmap from the CLASSE BMAD deck (0.70 T peak at the target back face). |
| **Linac 5–8** | `sim/linac5-8.py` · `config/linac5-8.yaml` | Four S-band traveling-wave sections (CEA 4/5 + CU 3/4) in one Impact-T deck, using the generic `rfdata4–7` field shape, accelerating the converter positrons. Space charge off; the real CESR capture optics are modelled (sec 5/6 264 A capture solenoids, 0.243 T flat-top, + the inter-section drift/quad lines from the CLASSE BMAD deck). |

Each stage's physics, field model, configuration knobs, and figures are documented in
[`docs/`](docs/).

## Shared helpers

`sim/helpers/` is the stage-agnostic plumbing every driver imports:

- **`buildfields.py`** builds the GDF→openPMD field maps read by the stage drivers. The SLAC
  traveling-wave maps are built once and reused by linac sections 1–4 rather than rebuilt per
  section. Map-geometry constants (gap centres, bore radius, 1-J/1-kW reference voltages) are fixed
  facts of the committed GDF inputs, not tunable operating points — tunables (gun voltage, RF
  power, frequency, Q) live in `config/*.yaml` and are passed in or read by the stage drivers. The
  openPMD files it writes use axis order `["r","z"]`, `m=0`, nodal component centering, and V/m & T
  unit dimensions — a deliberate, reader-validated deviation from WarpX's native RZ diag schema,
  kept because it matches how the stage drivers load these maps. Its prebuncher z-grid uniformity
  check exists because the GDF export's per-interval z spacing wobbles (worst case ~1.6% on the
  prebuncher map); using the mean step (rather than the first interval) avoids a drift of the far
  end of the 305 mm prebuncher map that otherwise mis-sizes the gap voltage.
- **`loadparticles.py`** handles beam handoff between stages. Gun's `anode_beam_mask` takes the
  anode-crossing electron flux as a single crest-time slab (z in the top fraction of the gap,
  moving forward) rather than an id-tracked transit count, since cathode particle dumps are written
  far sparser in time than the ~62 ps gap transit — this also excludes the dense near-cathode
  space-charge pileup and the reflected half of over-injected charge that never exits. Linac 5–8
  upsamples the converter/positron beam with `upsample_smeared` rather than plain bootstrap
  resampling: it bootstrap-draws parents with replacement, then jitters each clone's transverse
  phase space (x, y, px, py by default) by an amount proportional to its distance to its k-th
  nearest neighbour in standardized phase space (a local, KDE-style bandwidth), so clones decorrelate
  instead of landing on exact duplicates that would track identically under SC-off deterministic
  optics. Only the transverse columns are smeared — pz/z/t carry from the parent unchanged — so
  total energy is only approximately preserved (~0.5%), and a clone of a low-pz or large-angle
  parent can drop below an upstream KE floor after smearing (any KE cut that must hold exactly needs
  to be re-imposed downstream of the upsample). This upsampling is safe only with self-fields off,
  since the synthetic clones would otherwise inject a spurious self-field.
- **`metrics.py`**'s `screen_profile` treats each macroparticle's pooled `(id, z, quantity)` rows
  across all volumetric dumps as samples of a single continuous id-trajectory, valid only when
  motion is forward and monotonic in z; it interpolates each id's trajectory onto fixed z-screens
  and charge-weights the result, giving a true local-in-z phase-space diagnostic (no z-binning). It
  backs the gun `energy_gain` and `beam_envelope` figures.
- **`sim/plot/common.py`**'s `evolution_screens` assumes each `ParticleGroup` weight is real Coulomb
  charge for its 4th (charge) panel. The 2D cathode-stage diagnostic is a slab (not RZ), so its
  particle weight is charge-per-unit-out-of-plane-length rather than true Coulombs — the
  charge-vs-z panel from `evolution_screens`/`evolution_vs_z` is omitted for the cathode stage as
  not physically meaningful there.

## Configuration & frozen RF setpoints

Every option lives in `config/*.yaml` — edit a value there to retune a stage. The linac RF
setpoints (each section's crest phase and field scale) were derived once and **hardcoded** into the
linac YAMLs, so the drivers simply read and apply them (no runtime crest search or calibration). If
you change an upstream knob that shifts the beam, re-derive the affected setpoint.

## Outputs

- `logs/diags/<stage>/` — openPMD particle/field diagnostics + `injection_summary.json` sidecars.
- `logs/plots/<stage>/` — per-stage figures.
- `logs/pipeline/log_<date>.log` — the full run transcript.

Diagnostics and figures are regenerated by re-running; commit any PNGs you want to keep with
`git add -f`.

## Reference materials

`reference/` (in the original repo) holds documentation for the accelerator-physics tools
considered — WarpX, IMPACT-T, GPT, BMAD, G4beamline, LinacSim, lume-impact, lume-gpt, distgen,
openPMD-beamphysics / openPMD-viewer, easygdf — plus papers.
