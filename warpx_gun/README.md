# CESR gun in WarpX (RZ)

The second stage of the Cornell Linac electron source, simulated in WarpX. Stage 1
(`../warpx_cathode/`) is the thermionic cathode at the Child–Langmuir limit; here we take
its emitted electrons and accelerate them through the **CESR gun** — the electrostatic
accelerating structure modelled in Adam Bartnik's Linac GUI with the Poisson–Superfish field
map `CESR_gun.gdf` (the "Chili Gun Mk II", ~150 kV).

The gun field is applied as an external **electrode field** on the particles; WarpX's
electrostatic solver supplies the self-consistent beam **space charge** on top. Geometry is
**RZ (cylindrical)**, matching the field map's native symmetry.

## Pipeline

```bash
conda activate CBB
python warpx_gun/build_gun_field.py   # CESR_gun.gdf  ->  gun_field/gun_E.h5 (openPMD)
python warpx_gun/gun_sim.py           # RZ WarpX run  ->  diags/{fields,particles}/
python warpx_gun/plot_gun.py          # figures       ->  results/*.png
```

`build_gun_field.py` reads `../fieldmaps/CESR_gun.gdf`; `gun_sim.py` reads the cathode output
from `../warpx_cathode/diags/particles/`. Both paths are set near the top of each script.

## The gun field map

`CESR_gun.gdf` is a 2D cylindrical `(R, Z)` map of the gun's electrostatic field from
Poisson–Superfish, read with `easygdf`:

| quantity | value |
|----------|-------|
| R grid   | 151 points, 0 → 15 mm (ΔR = 0.1 mm) |
| Z grid   | 521 points, 0 → 51.77 mm (ΔZ ≈ 99.5 µm) |
| fields   | `Er`, `Ez` (V/m), normalized to a **1 kV** cathode→exit drop |
| magnetic | none — purely electrostatic gun |

**Voltage scaling and sign.** The map is normalized to a *+1 kV* cathode (V = +1000 at the
cathode, 0 at the exit), so its on-axis `Ez = -dV/dz` is *positive* — which would push
electrons back into the cathode. A real gun holds the cathode at *negative* high voltage with
the anode grounded, so we scale by a **negative** factor, `SCALE = -150` → a **−150 kV**
cathode. After scaling the on-axis field is `Ez(cathode) ≈ -1.94 MV/m`, peaking at
`-4.88 MV/m` near z ≈ 28 mm, and the 150 kV potential drop accelerates electrons in +z.

`build_gun_field.py` writes the scaled field as an openPMD file in the layout WarpX's
`read_from_file` external-field reader requires for RZ: geometry `thetaMode` with a single
`m = 0` mode, mesh record `E` with components `r`,`t`,`z`, axis labels `["r","z"]`, dataset
shape `(1, nr, nz)`. `gun_sim.py` then loads it via the raw WarpX inputs

```
particles.E_ext_particle_init_style = read_from_file
particles.read_fields_from_path     = warpx_gun/gun_field/gun_E.h5
particles.B_ext_particle_init_style = none
```

(PICMI has no class for a tabulated particle-applied field; `LoadInitialField` only sets a
one-time grid initial condition, which the Poisson solve overwrites — wrong for a static
electrode field.)

## Beam source — chaining the cathode output

The cathode run is a **continuous (DC) emitter**, so the weights in its last particle
snapshot encode the steady-state population *in transit through the diode* (~102 nC), not a
bunch charge. We:

1. Import the emitted **phase-space distribution** (positions + momenta) from the last
   cathode snapshot.
2. Remap the 2D `(x, z)` slab into RZ: treat `|x|` as the radius `r` and smear the particles
   uniformly in azimuth (`x = r cosθ, y = r sinθ`), rotating the transverse momentum
   accordingly.
3. **Renormalize** the total weight to a physical gun bunch charge `BUNCH_CHARGE = 0.1 nC`
   (the CESR gun is grid-pulse gated; 0.1 nC matches the prior IMPACT-T gun model).

**Why renormalize:** injecting the full 102 nC as one instantaneous bunch is unphysical — its
radial space-charge field (~50 MV/m) dwarfs the gun field and blows the beam apart before it
accelerates (observed directly: the beam is absorbed within ~50 steps). At 0.1 nC the beam
transports cleanly and accelerates. Set `BUNCH_CHARGE` at the top of `gun_sim.py` to explore
the space-charge regime.

**Approximations.** The cathode model is a 2D Cartesian slab, not RZ, so the slab→radius
remap is an approximation (it preserves the radial profile and momentum spread but not the
true cylindrical emission). The DC beam is treated as a single injected bunch.

## Simulation parameters (`gun_sim.py`)

| parameter | value |
|-----------|-------|
| geometry | RZ, `n_azimuthal_modes = 1` |
| grid | 96 (r) × 384 (z), r ∈ [0, 15 mm], z ∈ [0, 51.77 mm] |
| solver | electrostatic, lab frame, Multigrid (self-field only) |
| applied field | scaled `CESR_gun.gdf`, −150 kV, read from file |
| bunch | 0.1 nC, imported cathode phase space, ~222k macroparticles |
| time step | `dt = 0.4·Δz/v_exit` (v_exit ≈ 0.63 c at 150 keV) |
| duration | ~2× gun-transit time |

The lab-frame electrostatic solver is non-relativistic in its self-field treatment; at the
gun exit (β ≈ 0.63) this mildly overestimates the space-charge field — acceptable for a
single-pass gun demo, but note it if pushing to higher voltage.

## Figures (`results/`)

1. **`gun_field.png`** — on-axis `Ez(z)` and implied potential of the scaled field map: the
   accelerating field the beam sees (Ez < 0, 150 kV total drop).
2. **`beam_rz.png`** — `r–z` beam distribution at launch / mid-gun / exit: transport through
   the gun, including the near-cathode radial focusing.
3. **`energy_gain.png`** — mean and max kinetic energy vs. ⟨z⟩, climbing toward ~150 keV.
4. **`exit_phase_space.png`** — longitudinal `z–KE` phase space and the energy spectrum at the
   last dump.

## Notes / extensions

- The beam energy gain tracks `∫ e·|Ez| dz` (≈ 7.5 keV by z ≈ 4 mm), reaching the full 150 keV
  set by the cathode→exit potential drop.
- To approach the continuous-emission picture, inject a train of bunches or feed the cathode
  current directly rather than a single snapshot.
- A solenoid (magnetic focusing) could be added via a second `read_from_file` B map if the
  downstream Linac optics are included.
