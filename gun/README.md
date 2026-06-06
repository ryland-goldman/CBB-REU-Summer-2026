# CESR gun in WarpX (RZ)

The second stage of the Cornell Linac electron source, simulated in WarpX. Stage 1
(`../cathode/`) is the thermionic cathode at the Child–Langmuir limit; here we take
its emitted electrons and accelerate them through the **CESR gun** — the electrostatic
accelerating structure modelled in Adam Bartnik's Linac GUI with the Poisson–Superfish field
map `CESR_gun.gdf` (the "Chili Gun Mk II", ~150 kV).

The gun field is applied as an external **electrode field** on the particles; WarpX's
electrostatic solver supplies the self-consistent beam **space charge** on top. Geometry is
**RZ (cylindrical)**, matching the field map's native symmetry.

## Pipeline

```bash
conda activate CBB
python -c "import gun; gun.run()"   # build field map + sim + plots in one call
```

or, equivalently, the individual scripts:
```bash
python gun/build_gun_field.py   # CESR_gun.gdf  ->  gun_field/gun_E.h5 (openPMD)
python gun/gun_sim.py           # RZ WarpX run  ->  diags/{fields,particles}/
python gun/plot_gun.py          # figures       ->  results/*.png
```

To override the gun voltage or bunch charge: `gun.config(GUN_VOLTAGE=150e3,
BUNCH_CHARGE=0.1e-9)` before `gun.run()`. Keys must match the module-level constants in
`gun/build_gun_field.py` and `gun/gun_sim.py`. `build_gun_field.py` reads
`fieldmaps/CESR_gun.gdf`; `gun_sim.py` reads the cathode output from
`cathode/diags/particles/`. All paths are repo-root-relative.

**Performance knobs** (`config()`-overridable module constants; defaults reproduce the
original run): `REQUIRED_PRECISION` (1e-5) and `MAX_ITERS` (None) for the MLMG solve;
`CFL` (0.4, `dt = CFL·dz/v_exit`), `TRANSIT_MARGIN` (1.15) and `AVG_SPEED_FRAC` (0.6) for
the auto-derived run length, or `MAX_STEPS` (>0) to fix it; `N_DIAGS` (40) for the openPMD
dump count; `MAX_PART` (0 = no cap) to downsample the imported cathode bunch (reweighted,
charge-preserving); and the grid `nr, nz`. Runtime ≈ `nz²` (per-step cost ∝ cells, and
`dz = ZMAX/nz` ⇒ fewer steps as `nz` drops), so halving `nz` ≈ 4× faster. This holds because the
gun's cells are near-isotropic (`dz/dr ≈ 1.3`) so the MLMG solve stays well-conditioned as `nz`
drops — **unlike the prebuncher's long-thin box**, where coarsening `NZ` slows the solve instead
(see `prebuncher/README.md`). Keep `N_DIAGS ≥ 20` so `space_charge.png` still finds its
near-launch field snapshot (it self-skips otherwise).

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
particles.read_fields_from_path     = gun/gun_field/gun_E.h5
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
   accordingly. Crucially, the revolution carries a **2πr Jacobian**: a slab uniform in `x`
   has a flat `dN/dr`, which—revolved naively with `r = |x|` and unchanged weight—would give
   areal density `n(r) ∝ 1/r`, a spurious on-axis charge cusp. We therefore
   **importance-resample by `r`** (draw particles with probability ∝ `r`, with replacement),
   so `dN/dr → r·dN/dr` and `n(r)` matches the cathode's true radial profile (a flat-top
   emitting strip → a uniform-density disc). This keeps the macroparticle weights uniform.
3. **Renormalize** the total weight to a physical gun bunch charge `BUNCH_CHARGE = 0.1 nC`
   (the CESR gun is grid-pulse gated; 0.1 nC matches the prior IMPACT-T gun model).

**Why renormalize:** injecting the full 102 nC as one instantaneous bunch is unphysical — its
radial space-charge field (~50 MV/m) dwarfs the gun field and blows the beam apart before it
accelerates (observed directly: the beam is absorbed within ~50 steps). At 0.1 nC the beam
transports cleanly and accelerates. Set `BUNCH_CHARGE` at the top of `gun_sim.py` to explore
the space-charge regime.

**Approximations.** The cathode model is a 2D Cartesian slab, not RZ, so the slab→radius
remap is an approximation: the `r`-importance resample (step 2) makes the **areal density**
match the cathode's radial profile, but the reconstructed azimuthal distribution is assumed
uniform (it cannot recover the true cylindrical emission, which the 2D slab never had). The DC
beam is treated as a single injected bunch.

## Simulation parameters (`gun_sim.py`)

| parameter | value |
|-----------|-------|
| geometry | RZ, `n_azimuthal_modes = 1` |
| grid | 96 (r) × 384 (z), r ∈ [0, 15 mm], z ∈ [0, 51.77 mm] |
| solver | electrostatic, lab frame, Multigrid (self-field only) |
| applied field | scaled `CESR_gun.gdf`, −150 kV, read from file |
| bunch | 0.1 nC, imported cathode phase space, ~222k macroparticles (optionally capped by `MAX_PART`, reweighted) |
| time step | `dt = CFL·Δz/v_exit` (`CFL`=0.4; v_exit ≈ 0.63 c at 150 keV) |
| duration | `TRANSIT_MARGIN`×gun-transit time (=1.15; bunch average speed ≈ `AVG_SPEED_FRAC`·v_exit, =0.6), or fixed via `MAX_STEPS`; stops as the beam reaches the exit — running longer empties the domain and aborts the Multigrid self-field solve |

The lab-frame electrostatic solver is non-relativistic in its self-field treatment: it applies
the rest-frame Coulomb field `qE_r` with no magnetic-pinch cancellation, whereas the true net
transverse force is `qE_r/γ²`. At the gun exit (148 keV, γ ≈ 1.29) this **overestimates the
transverse space-charge force by ≈ γ² = 1.66×, i.e. ~66 %** — ramping from a few % near the
cathode (10 keV) to ~66 % at exit. (The genuine fix is WarpX's relativistic ES mode, out of
scope for this single-pass demo; the same caveat applies, but shrinks, in the more relativistic
prebuncher/linac stages — see those READMEs.) Acceptable for the demo, but note it if pushing
to higher voltage or interpreting the absolute σ_r / emittance.

## Figures (`results/`)

1. **`gun_field.png`** — on-axis `Ez(z)` and implied potential of the scaled field map: the
   accelerating field the beam sees (Ez < 0, 150 kV total drop).
2. **`beam_rz.png`** — `r–z` beam distribution at launch / mid-gun / exit: transport through
   the gun, including the near-cathode radial focusing.
3. **`energy_gain.png`** — mean and max kinetic energy vs. ⟨z⟩, climbing toward ~150 keV.
4. **`exit_phase_space.png`** — longitudinal `z–KE` phase space and the energy spectrum at the
   last dump.
5. **`beam_envelope.png`** — RMS radial size `σ_r = √⟨x²⟩` and normalized transverse emittance
   `εn,x` vs. `⟨z⟩`: the near-cathode radial focusing of `beam_rz.png` made quantitative, plus the
   space-charge / aberration emittance growth along the gun.
6. **`space_charge.png`** — `r–z` maps of the beam **self-field** (`ρ` and the space-charge
   potential well `φ`, ≈ −250 V) at a near-launch snapshot — the dumped self-field nothing else
   plots, and the well that motivates renormalizing the bunch to 0.1 nC.

## Notes / extensions

- The beam energy gain tracks `∫ e·|Ez| dz` (≈ 7.5 keV by z ≈ 4 mm), approaching the ~150 keV
  set by the cathode→exit potential drop (the space-charge-loaded beam lands at ~148 keV mean).
- To approach the continuous-emission picture, inject a train of bunches or feed the cathode
  current directly rather than a single snapshot.
- A solenoid (magnetic focusing) could be added via a second `read_from_file` B map if the
  downstream Linac optics are included.
