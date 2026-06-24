# Injector (WarpX RZ)

Third stage of the Cornell Linac chain:

```
cathode  ->  gun  ->  injector (this)  ->  linac_sec1
```

The injector is the full injector subsection modelled in **one self-consistent RZ
space-charge run**: two 214 MHz prebuncher cavities (Preb 2 reversed) velocity-bunch the gun
exit beam while six solenoid lenses focus it, handing a focused, collimated beam to the linac
at the true structure entrance **z ~= 2.03 m**. Modelling every element in one drift is
essential -- the bunching, two-cavity phasing, and transverse focusing are coupled through the
self-field.

```
Lens 0A  ->  Prebuncher 1  ->  Prebuncher 2 (reversed)  ->  Sol 0 / Lens 0E  ->  9.547 mm iris
 6 A          8 kW                10 kW                       40 A / 10 A          @ 1.922 m
 @0.225 m     @0.534 m            @1.318 m                    @1.897 / 1.914 m
```

Driven through **lume-warpx**: every constant lives in `config/injector.yaml` (operating point +
the solenoid/prebuncher knobs in its `params:` block) and `sim/injector.py` reads them back,
imports the gun handoff via `WarpX(initial_particles=...)`, and overrides only the
runtime-computed values (the per-field RF/solenoid time functions, step count, `dt`, diagnostic
period). Edit `config/injector.yaml` to retune.

## Running

```
python -m sim.injector          # build fields + run the sim (writes logs/diags/injector/main)
python -m sim.plot.injector     # regenerate figures from the existing diags
```

`sim/injector.py` calls `build_injector_fields()` (idempotent) to (re)build the openPMD maps,
then reads the gun exit beam from **`logs/diags/gun/handoff`** (the full reconstructed
time-release exit beam) when present, else from `logs/diags/gun/particles` (the legacy
snapshot). The beam is already RZ; the run shifts its tail to the entrance and imports it.

## Field maps

`sim/helpers/buildfields.build_injector_fields()` writes eight openPMD maps from the `.gdf`
sources in `fieldmaps/gdf/`:

- `preb1_EB.h5` -- Prebuncher 1, forward 1-J cavity field, gap at `Z_GAP_CENTER_1 = 0.534 m`.
- `preb2_EB.h5` -- the **same forward field**, gap at `Z_GAP_CENTER_2 = 1.318 m`. The two
  cavities differ only in lab-z placement (`grid_global_offset`) and the run-time reversal
  phase; there is **no mirrored/negated map**.
- `lens0a.h5` / `lens0b.h5` / `lens0c.h5` / `lens0d.h5` / `sol0.h5` / `lens0e.h5` -- the six
  static, per-Ampere B-only solenoid maps.

> **Descending-z GDF -> Er sign.** `prebuncher_25D.gdf` stores its z column **descending**
> (+152.4 -> -152.4 mm), unlike the solenoid maps (ascending). The flat GDF rows are reversed to
> the ascending axis *before* the `(nz,nr).T` transpose. Skipping the reversal z-flips the map
> relative to its own axis, which negates the **odd** `Er` (Ez/Bphi are even and unaffected) -- a
> silent bug: on-axis bunching is unchanged, but every off-axis transverse RF force is
> wrong-signed. The gap-parity asserts do not catch it (parity is invariant under the flip), so
> the build adds a raw-GDF orientation check at a fixed off-axis +z point.

## RF drive: zero-crossing, centroid-referenced velocity bunching

Each cavity drives its 1-J map as a standing-wave TM mode: `Er,Ez(t) = map*scale*cos(wt+phi)`,
`Bphi(t) = map*scale*sin(wt+phi)`, with field scale `sqrt(1e3*Q*P/(2pi f_RF))`. The drive phase is

```
phi = -w*t_gap + base + radians(phi_off_deg) + rev_phase
```

where `t_gap` is the arrival time of the **bunch centroid** (`z_ref = z_centroid`, not the tail)
at the gap. The faithful default is **`PHASE="zc"`** (`base = pi/2`) with `phi_off = 0`: this
lands the **centroid on the RF zero-crossing**, so the net mean-energy kick is zero and each
cavity acts as a **pure velocity buncher** -- bunching does not change the mean energy.
`base = pi` (`"crest"`) is the legacy net-accelerating reference, kept for comparison.

**Why centroid-referenced.** The physical beam is the ~2 ns time-release stream spanning a wide
slice of RF phase (the centroid sits well behind the tail). Phasing the *tail* to the
zero-crossing would leave the bulk on the decelerating slope, so the reference plane is the
centroid. The cavity is energy-neutral by the symmetry of -cos(phi) about the zero-crossing; the
head is decelerated and the tail accelerated => a compressive chirp. A single zero-crossing
imparts a *sinusoidal* (not linear) chirp across the wide bunch, so the longitudinal phase space
folds into a checkmark at the waist rather than collapsing to a point (fold-limited bunching).

### Reversed install (`PREB2_REV_PHASE = pi`)

GPT installs Preb 2 with a 180-degree rotation (`-1,0,0`). For this map's definite parity --
Ez EVEN, Er ODD, Bphi EVEN about the gap (asserted in the build; also forced by Maxwell for a
TM0 mode) -- the rotation flips all three lab components, i.e. a global E,B sign flip **= +pi in
absolute drive phase**. In the `zc + phi_off = 0` parametrization, `phi_off` carries no reversal
information, so Preb 2 must carry the genuine geometric +pi itself: **`PREB2_REV_PHASE = pi`**.

Both `rev = 0` and `rev = pi` are energy-flat (the centroid sits on a zero-crossing either way);
they differ only in the chirp **slope**. `rev = 0` bunches harder and over-compresses to an
earlier waist that re-expands by the handoff; the faithful `rev = pi` is the less-aggressive
slope whose waist lands at the 2.03 m handoff. `rev = pi` is therefore both geometrically and
operationally correct -- do not "fix" it to `rev = 0`. The value is convention-dependent: under
the legacy `crest`+GUI convention the reversal is absorbed into the crest reference, so set
`rev = 0` if you switch `PHASE` back to `"crest"` with GUI `phi_off`.

### Preb-2 timing caveat

Preb-2's arrival is timed in two segments -- `v_beam` to gap 1, then the post-Preb-1 speed over
gap 1 -> gap 2 -- using an analytic estimate of Preb-1's net kick. At the `zc`/centroid default
Preb-1's net kick is ~0 (energy-flat), so `v_after_preb1 ~= v_beam` and the residual is small.
The exact fix is a two-pass run (read the post-Preb-1 beta from a diagnostic, rebuild the Preb-2
timing); only needed if a future study gives Preb-1 a non-zero net kick (a Preb-1 power scan, or
`PHASE="crest"` which net-accelerates).

## Solenoid lenses (transverse focusing, native placement)

All six injector solenoids are built and wired (in z-order). LENS_0A / SOL_0 / LENS_0E carry the
default currents (6 / 40 / 10 A); **LENS_0B / 0C / 0D default to 0 A** -- inert at the default
operating point (a 0-A lens is skipped) but real magnets, `config`-overridable for transport
studies. Each lens is a separate per-Ampere B-only openPMD map; the 1-A maps scale linearly with
current.

- **Placement is NATIVE absolute machine-z**, matching `gpt_master.in`, which installs every
  solenoid with `Map2D_B("wcs", "z", 0.0, ...)` -- the GDF's stored Z column *is* absolute
  machine z, so `grid_global_offset` is simply the native origin `z[0]`. **Do NOT align
  argmax -> GUI z:** that is correct only for the narrow lenses (argmax ~= GUI) but **wrong for
  the flat-top SOL_0** -- its argmax is an arbitrary point on its broad plateau and the GUI value
  is a center/edge label, not the peak. Forcing argmax -> GUI shifts the whole channel ~+1 m
  (mostly past the 2.03 m handoff) and leaves the prebunchers unfocused.
- **SOL_0 is a long channel solenoid**: its strong-field region spans the entire prebuncher
  section, focusing both cavities. LENS_0B/0C/0D cluster near Lens 0E around 1.9 m. Lens 0A sets
  the early envelope and the Lens 0E telescope adds the final squeeze into the iris.
- **Ordering gotcha:** picmi forces the global `E_ext_particle_init_style` to "none" if the
  last-added applied field has `load_E=False`, so the B-only solenoids are added **before** the
  RF cavities; an `assert fields[-1].load_E` guards it.
- The build asserts each solenoid's in-domain lab-z peak is in [0, 2.10 m] and **upstream of the
  2.03 m handoff** (so the linac never inherits a beam still inside a lens). A loose
  peak-within-half-FWHM-of-GUI-z cross-check runs only for the narrow lenses.

> **Handoff-seam residual field (approximation).** Lens 0E peaks just upstream of the handoff and
> its field tail is still substantial at the plane; the native-placed Sol 0 peaks far upstream so
> its tail there is small. The linac stage models no solenoid, so this continuing focus is dropped
> at the seam -- an unphysical discontinuity that makes capture *more* conservative (the real
> continuous field would hold the envelope tighter into the structure). A fully faithful treatment
> would carry the tails into the linac stage or move the handoff downstream -- a follow-up.

## The 9.547 mm collimator

The prebuncher subsection carries a `scatteriris` of radius **9.547 mm at z = 1.922 m** followed
by a 9.547 mm pipe to 2.1 m -- the injector->linac aperture. Past 1.922 m the restriction is the
SLAC bore, and Sol 0 / Lens 0E peak just upstream precisely to squeeze the beam through. The
faithful success metric is **transmission through the 9.547 mm iris**, not "contained within the
36 mm domain."

It is applied as a **multi-plane id scrape** (`pipe_violator_ids` / `survivor_mask`), not an
in-run scrape: this pywarpx RZ build's particle-position SoA accessors do not expose the radial
position as a named component, so an afterstep weight-zeroing callback is unavailable. A single
radial cut at the 2.03 m handoff would be wrong: the Sol 0 / Lens 0E telescope focuses the beam
*hard* across the 1.922 -> 2.03 m tail -- it converges, not diverges, so a particle outside the
iris at 1.922 m (scraped by the real machine) can converge back inside it by 2.03 m and be
wrongly kept. Instead the continuous pipe is emulated by tracking particle `id` across every dump
from z = 1.922 m on: a particle outside the aperture at *any* plane in the pipe is removed; only
the survivors are injected into the linac. This is exact in the dense-dump limit and reduces to
the entrance-plane cut when the envelope happens to be monotone. The sim writes the collimated
handoff charge / transmission to `logs/diags/injector/main/injection_summary.json`.

## Domain / grid

- **z:** `ZMIN=0`, `ZMAX=2.10 m`, with a field-free exit drift so the handoff beam coasts.
  Handoff snapshot at **z ~= 2.03 m**.
- **r:** `RMAX=0.036 m`, `NR=80` (dr ~= 0.45 mm) -- the RF map reaches 36 mm and needs the radial
  resolution; do not copy the linac's NR=16.
- **Cell aspect (binding):** `NZ=1664` gives dz ~= 1.262 mm => dz/dr ~= 2.80:1 and is divisible
  by the blocking factor (8). Do not coarsen NZ -- this long-thin box is convergence-bound, so
  coarsening slows the per-step MLMG solve faster than it removes cells and under-resolves the
  ~1 mm bunch. Speed it via `CFL`, `maximum_iterations`, `required_precision`.
- **Handoff diagnostic:** the dump cadence (`period`) is sized from the post-Preb-2 speed (the
  speed *at* the handoff plane) so the spacing near 2.03 m is <= `HANDOFF_DZ`, landing a snapshot
  within ~1 mm of the plane. picmi exposes only a uniform `period` (a true z-station diagnostic
  isn't available in this build).

## Self-field solver (relativistic electromagnetostatic)

The beam self-field uses WarpX's **electromagnetostatic** solver (`ES_MLMG_EMS`), not the plain
lab-frame electrostatic Poisson solve. In addition to `div^2 phi = -rho/eps0` it solves the
Coulomb-gauge vector potential from the beam current (`div^2 A = -mu0 j`, `B = curl A`), so the
self magnetic field is included and the relativistic magnetic-pinch term `q*beta x B` partially
cancels the radial space-charge repulsion: the net transverse self-force is `q E_r / gamma^2`
rather than the pure-electrostatic `q E_r`. This removes the ~gamma^2 transverse-SC over-repulsion
the lab-frame solver incurs over this long line. Matches the gun's solver.

> **Outer radial wall BC is `dirichlet`, not `neumann`.** The magnetostatic vector-Poisson's
> dominant A_z component (driven by `j_z`) has an all-Neumann, *singular* operator unless the
> outer wall is grounded. At RMAX=36 mm (well outside the beam) the self-field has decayed, so the
> beam dynamics are unaffected; this only makes A_z well-posed.

## Config knobs (`params:` in `config/injector.yaml`)

| Key | Meaning |
|-----|---------|
| `CFL` | `dt = CFL * dz / v_beam` |
| `MAX_PART` | downsample the gun snapshot (reweighted) |
| `MAX_STEPS` | 0 -> auto-derive from transit; >0 -> fixed |
| `N_DIAGS` / `HANDOFF_DZ` | diagnostic cadence (handoff-spacing capped) |
| `TRANSIT_MARGIN` | stop just before the bunch centre reaches the exit |
| `Z_INJECT` | lab z of the bunch tail |
| `PHASE` | `zc` (centroid on zero-crossing) or `crest` (legacy net-accelerating) |
| `F_RF` / `Q_L_1` / `Q_L_2` | RF operating point (frequency, loaded Q per cavity) |
| `PREB1_KW` / `PREB2_KW` | prebuncher powers (0 -> cavity off) |
| `PREB1_PHI_OFF` / `PREB2_PHI_OFF` | phase offset [deg] from the base |
| `PREB2_REVERSED` / `PREB2_REV_PHASE` | reversed-install +pi handling |
| `I_LENS0A`..`I_LENS0E` / `I_SOL0` | solenoid currents [A] (0 -> lens off) |
| `COLLIM_R` / `COLLIM_Z` / `COLLIMATE` | iris radius / start / enable the scrape report |

## Figures

`python -m sim.plot.injector` reads `logs/diags/injector/main/particles` and writes to
`logs/plots/injector/` (all filenames keep the `injector_` prefix). Three layers: lume-warpx's
plotting helpers, the shared custom figures in `sim/plot/common.py`, and the stage-specific rich
figures (raw openPMD over the whole run; the prebuncher field lobes read the `preb1/2_EB.h5` maps):

- `injector_phase_space_z_KE.png` -- `plot2D("z","kinetic_energy")`: the velocity-bunched exit beam.
- `injector_transverse_x_px.png` -- `plot2D("x","px")`: the solenoid-focused exit transverse phase space.
- `injector_centroid_vs_t.png` -- `plot1D("t","mean_z")`: the bunch traversing the line.
- `injector_bunch_length_vs_t.png` -- `plot1D("t","sigma_z")`: sigma_z compressing to its waist.
- `injector_emittance_vs_t.png` -- `plot1D("t","norm_emit_x")`: transverse emittance over the run.
- `injector_beamsize_vs_t.png` -- `plot1D("t","sigma_x")`: the transverse envelope (solenoid focusing).
- `injector_energy_spectrum.png` -- charge-weighted handoff KE histogram.
- `injector_current_profile.png` -- longitudinal current I(z), the velocity-bunched current peak.
- `injector_energy_chirp.png` -- slice <KE> vs z, the chirp that drives the bunching.
- `injector_beam_spot_xy.png` -- the transverse x-y spot (RZ-reconstructed y) with marginals.
- `injector_cavity.png` -- on-axis E_z lobes of both prebunchers (scaled by drive amplitude), lab z, handoff marked.
- `injector_line.png` -- sigma_z, peak current, and mean KE vs <z> along the line.
- `injector_bunch_profile.png` -- the line-charge density lambda(z) at four stations.
- `injector_phasespace.png` -- the longitudinal phase space (z-<z>, KE-<KE>) at the same four stations.
