# CESR Gun — Electrostatic Accelerating Structure (WarpX RZ)

The second stage of the Cornell linac electron source. Stage 1 (`docs/cathode.md`) is the
thermionic cathode at the Child–Langmuir limit; here we take its emitted electrons and
accelerate them through the **CESR gun** — the electrostatic accelerating structure modelled
in Adam Bartnik's Linac GUI with the Poisson–Superfish field map `CESR_gun.gdf` (the "Chili
Gun Mk II", ~150 kV). Built on `pywarpx`, driven through lume-warpx: every constant lives in
`config/gun.yaml`, and `sim/gun.py` reads them back, overriding only the runtime-computed
values (`dt`, step count, the seed arrays, diagnostic periods).

Geometry is **RZ (cylindrical)**, matching the field map's native symmetry. The gun field is
applied as an external **electrode field** on the particles; WarpX's **electromagnetostatic**
solver supplies the self-consistent beam **space charge** on top.

Run:
```bash
conda activate CBB
python sim/gun.py        # builds the field map, then writes openPMD to logs/diags/gun/{fields,particles}
python sim/plot/gun.py   # writes the figures to logs/plots/gun/
```

`sim/gun.py main()` runs only the simulation (it builds `fieldmaps/h5/gun_E.h5` first, so a
fresh checkout works); `sim/plot/gun.py main()` runs only the plotting (the sim must have been
run first).

---

## The applied gun field

`CESR_gun.gdf` is a 2D cylindrical `(R, Z)` map of the gun's electrostatic field from
Poisson–Superfish, read with `easygdf`. Its `Er`, `Ez` columns are normalized to a **+1 kV**
cathode→exit drop (V = +1000 at the cathode, 0 at the exit), so the native on-axis
`Ez = -dV/dz` is *positive* — which would push electrons back into the cathode. A real gun
holds the cathode at *negative* high voltage with the anode grounded, so the build scales by a
**negative** factor, `SCALE = -GUN_VOLTAGE/1 kV` (= −150 for the shipped 150 kV point), giving
an on-axis field that *accelerates* electrons in +z and a 150 kV potential drop. The map is
purely electrostatic — no magnetic field.

`sim.helpers.buildfields.build_gun_field(GUN_VOLTAGE)` reads `fieldmaps/gdf/CESR_gun.gdf` and
writes the scaled field to `fieldmaps/h5/gun_E.h5` as an openPMD `thetaMode` (RZ, `m = 0`)
mesh record `E` with components `r`, `t`, `z`. `sim/gun.py` calls it at the start of `main()`
(idempotent, so a fresh checkout rebuilds the map), then loads it via the
`AppliedFromFile` field in the YAML — WarpX applies it directly to the particles as an
external electrode field (PICMI has no class for a tabulated particle-applied field, and a
grid initial condition would be overwritten by the Poisson solve, so the raw `read_from_file`
applied-field path is used).

`GUN_VOLTAGE` is a single config value (`params.GUN_VOLTAGE` in `config/gun.yaml`): it sets
*both* the field-map scale (passed to `build_gun_field`) and the exit kinematics (γ, `v_exit`,
the time step and run length), so the two cannot drift apart.

---

## Beam source — chaining the cathode output

The cathode run is a **continuous (DC) emitter**, so the weights in its last particle snapshot
encode the steady-state population *in transit through the diode*, not a bunch charge. The gun
imports that phase space and reshapes it (`load_cathode_bunch`):

1. **Import** the emitted phase-space distribution (positions + momenta) from the last cathode
   snapshot at `logs/diags/cathode/particles`, restricted (`gap_d`/`anode_frac`) to the
   forward-moving beam crossing the anode plane — the delivered flux only. This deliberately
   drops the near-cathode charge pileup and any reflected/over-injected population, which would
   otherwise inflate the seeded charge.
2. **2D-slab → RZ remap.** The cathode is a 2D Cartesian (x, z) slab; treat `|x|` as the radius
   `r` and smear the particles uniformly in azimuth (`x = r cosθ, y = r sinθ`), rotating the
   transverse momentum accordingly (radial component `ux·sign(x)`, azimuthal `uy`). The
   revolution carries a **2πr Jacobian**: a slab uniform in `x` has flat `dN/dr`, which revolved
   naively (`r = |x|`, unchanged weight) would give areal density `n(r) ∝ 1/r`, a spurious
   on-axis charge cusp. We therefore **importance-resample by `r·w`** (draw with replacement
   with probability ∝ `r·w`), so `dN/dr → r·dN/dr` and `n(r)` matches the cathode's true radial
   profile (a flat-top emitting strip → a uniform-density disc), keeping macroparticle weights
   uniform. Because the draw is with replacement, the effective independent sample count is
   below the drawn count (relevant if `MAX_PART` is set small); each resampled copy still draws
   its own independent azimuthal angle `θ`, so duplicate copies from the resampling spread around
   the ring rather than overlapping at the same azimuthal position.
3. **Renormalize** the total weight to a physical gun bunch charge `BUNCH_CHARGE` (1 nC). The
   CESR gun is grid-pulse gated; injecting the full DC-transit charge as one instantaneous
   bunch is unphysical — its radial space-charge field dwarfs the gun field and blows the beam
   apart before it accelerates. At 1 nC the beam transports and accelerates cleanly.

This slab→radius remap is an approximation: the `r`-importance resample makes the **areal
density** match the cathode's radial profile, but the reconstructed azimuthal distribution is
assumed uniform (it cannot recover a true cylindrical emission the 2D slab never had).

### Timed release

Even at 1 nC, *how* the bunch is fed to the gun matters. The CESR gun is gated by a **2 ns grid
pulse** (`PULSE_WIDTH`), which is several gun-transit-times long, so the physical beam is a
**long, low-density, quasi-DC stream**. The imported macroparticles are released over
`PULSE_WIDTH` by a per-step `beforestep` injection callback
(`ParticleContainerWrapper.add_particles`), so only the fraction of charge already emitted is
present at any instant — the physical low line-density. PICMI needs a non-empty initial
distribution, so the YAML seeds a single macroparticle; the callback injects the rest, walking
a t-sorted emission-time array in one pass.

This reproduces the *mechanism* of the gun's time-release, not the pulse *profile*: emission
times are drawn uniformly (flat-top), whereas a real pulse has edge ramps that emit less at the
edges. The cathode run is a *stationary* DC emitter with no per-particle birth time, so every
release time carries the **same** steady-state cathode phase space — a DC line-density
reconstruction, not a resolved emission history. (The older "snapshot" representation, which
crammed all 1 nC in at t=0 and over-stated the space-charge force, is not modelled here — only
the physical timed release.)

### Exit-beam handoff

The released beam's ballistic z-extent (~`v_exit·PULSE_WIDTH`, ≈ 0.4 m) is many times the gun
domain, so no single volumetric snapshot can hold it. `build_exit_handoff()` reconstructs the
full exit beam by particle **id** across the volumetric dumps: it samples each particle in the
**field-free pad past the field map** (z ≥ `ZMAX_FIELD`) — its first appearance there is its
exit-plane phase space — and ballistically drifts all samples to a common reference time to
rebuild one consistent snapshot (head downstream, tail at the entrance).

Sampling in the **field-free** pad (not at the particle's last in-field dump) is essential: a
field-free drift preserves εn,x, but drifting a still-in-field particle as if field-free
manufactures a spurious x–u correlation and inflates the emittance — the reason the domain
carries the `ZPAD` field-free drift pad. Particles that never reach the pad are classified as a
radial (`r = RMAX`) loss or, for the last sliver of the pulse, an un-flushed tail (the run ends
first), and counted. The result is written to `logs/diags/gun/handoff` (openPMD), which the
injector stage reads.

Note the εn,x of the reconstructed *instantaneous* beam is the **projected** emittance of the
~0.4 m drifting bunch (the head drifted longer than the tail), physically several times the
**per-slice** beam quality — that head–tail projection is real (an instantaneous diagnostic of
the 2 ns beam would measure it), not an artifact, and the injector's prebunchers compress the
beam longitudinally and recover it.

---

## Space-charge model — electromagnetostatic self-field

The self-field uses WarpX's **lab-frame electromagnetostatic** solver
(`ES_MLMG_EMS`): on top of the electrostatic Poisson solve (∇²φ = −ρ/ε₀, **E** = −∇φ) it also
solves the Coulomb-gauge vector potential from the beam current (∇²**A** = −μ₀**j**,
**B** = ∇×**A**), so the beam's **self magnetic field** is included. The resulting
magnetic-pinch force `qβ×B` partially cancels the radial electric repulsion, giving the correct
relativistic net transverse self-force `qE_r/γ²` rather than the pure-electrostatic `qE_r`. At
the gun exit (γ ≈ 1.3) a plain lab-frame electrostatic solver would **overestimate the
transverse space-charge force by ≈ γ²** (ramping from a few % near the cathode to tens of % at
exit); the electromagnetostatic solver removes that error self-consistently, and generalizes to
a multi-velocity beam (unlike WarpX's per-species relativistic ES mode).

The magnetostatic vector-Poisson solve adds a 3-component MLMG solve per step (≈2× slower per
step than an electrostatic-only gun). **Boundary requirement:** the outer radial wall is
`dirichlet` (not `neumann`) — the `A_z` component (driven by the dominant beam current `j_z`)
would otherwise have an all-Neumann singular operator and the MLMG bottom solve **diverges**.
Grounding the pipe at r = 15 mm — well outside the beam — makes it well-posed (and models the
real beampipe as a grounded conductor). The z = 0 cathode plate is also `dirichlet`: the metal
cathode is a Dirichlet plane whose opposite-sign image partly cancels the bunch self-field (a
Neumann plane would give a same-sign image that wrongly adds to it).

---

## Time step and run length

The exit kinematics follow from `GUN_VOLTAGE`: γ = 1 + e·V/(m_e c²), `v_exit = c√(1 − γ⁻²)`.
The time step is `dt = CFL·Δz/v_exit` (`CFL = 0.7`, the shipped Fast value — Balanced is 0.4;
`Δz = ZMAX/nz`). The run length is sized on
the **field transit** (`ZMAX_FIELD`), not the padded domain: `PULSE_WIDTH +
TRANSIT_MARGIN × (ZMAX_FIELD / (AVG_SPEED_FRAC·v_exit))`, so the last-released particle clears
the gun while the run stops with the beam still in the pad (over-running drains the padded
domain and aborts the MLMG solve). `MAX_STEPS > 0` fixes the count instead.

`AVG_SPEED_FRAC = 0.6` is hand-tuned for the 150 kV point — `v_exit` is recomputed from
`GUN_VOLTAGE` but the average-speed fraction is not, so re-check it if `GUN_VOLTAGE` changes
substantially.

---

## Configuration knobs (`config/gun.yaml`)

Physics / numerics:
- `grid.number_of_cells` `[nr, nz]` — shipped **Balanced** default `[128, 384]` (dz ≈ 0.19 mm).
  The cross-code-converged grid is `[128, 712]` (dz ≈ 0.10 mm; runtime ≈ `nz²`, so the coarser
  grid is ~4× faster). The gun's cells are near-isotropic, so the MLMG solve stays
  well-conditioned at either `nz`.
- `solver.required_precision` / `solver.warpx_magnetostatic_required_precision` — MLMG
  tolerances; ship **1e-3** (the Fast value; Balanced is 1e-4, Conservative/benchmark 1e-5).
- `species[0].warpx_do_not_deposit` — `false` (space charge ON). `true` turns the self-field
  off (a large physics change, not a mild diagnostic — the self-field is dominant here).

Operating point (`params:`):
- `GUN_VOLTAGE` (150 kV) — cathode high voltage; sets the field-map scale and exit kinematics.
- `BUNCH_CHARGE` (1 nC) — renormalized gun bunch charge.
- `PULSE_WIDTH` (2 ns) — the grid-pulse emission window the beam is released over.
- `CFL` (0.7), `TRANSIT_MARGIN` (1.15), `AVG_SPEED_FRAC` (0.6), `MAX_STEPS` (0 = auto) — time
  step and run-length controls.
- `ZMAX_FIELD` (51.765 mm), `ZPAD` (20 mm) — the field-map z-extent (exit plane) and the
  field-free drift pad past it (must equal the grid `upper_bound[1] − ZMAX_FIELD`).
- `N_DIAGS` (40) — openPMD dump count.
- `MAX_PART` (50000; 0 = no cap) — downsample the imported cathode bunch (reweighted,
  charge-preserving).
- `RNG_SEED` (0) — seeds the resample / emission-time draws.

---

## Figures (`logs/plots/gun/`)

Generated by `sim/plot/gun.py` in three layers — lume-warpx's plotting helpers, the shared
custom figures in `sim/plot/common.py`, and the stage-specific rich figures (raw openPMD):

1. **`phase_space_z_KE.png`** — `plot2D("z","kinetic_energy")`: exit longitudinal phase space.
2. **`transverse_x_px.png`** — `plot2D("x","px")`: exit transverse phase space.
3. **`Ez_rz.png`** — `plot_fields("E","z","r")`: field in the gun gap (applied + self-field).
4. **`self_charge_rz.png`** — `plot_fields("rho","z","r")`: the beam self charge density.
5. **`centroid_vs_t.png`** — `plot1D("t","mean_z")`: the bunch marching down the gun.
6. **`emittance_vs_t.png`** — `plot1D("t","norm_emit_x")`: normalized transverse emittance.
7. **`beamsize_vs_t.png`** — `plot1D("t","sigma_x")`: the transverse envelope σ_x.
8. **`energy_spectrum.png`** — `common.energy_spectrum`: charge-weighted exit KE histogram.
9. **`current_profile.png`** — `common.current_profile`: longitudinal current I(z) of the pulse.
10. **`beam_spot_xy.png`** — `common.beam_spot`: the transverse x–y spot (RZ-reconstructed y).
11. **`gun_field.png`** — the on-axis applied E_z(z) and implied potential, from `gun_E.h5`.
12. **`beam_rz.png`** — the beam shape in r–z at launch / mid-gun / exit.
13. **`energy_gain.png`** — mean and max KE vs z on **fixed-z virtual screens**
    (`metrics.screen_profile`): local-in-z energy gain saturating at the gun voltage. The screen
    reconstruction (id-track each particle and interpolate to every z-plane it crosses) avoids
    the bin-to-bin εn,x jitter a z-histogram of the timed quasi-DC stream would produce.
14. **`beam_envelope.png`** — σ_x and εn,x vs z on the same virtual screens.
15. **`space_charge.png`** — the beam **self-field** ρ(r,z) and φ(r,z) near launch, separate
    from the applied gun field.

---

## Notes

- **Fresh diags on rerun:** WarpX *appends* one openPMD file per dump, so `sim/gun.py` removes
  the existing `logs/diags/gun/{fields,particles}` at the start of each run. Without this,
  re-running with a different grid/step count leaves stale files that interleave with the new
  ones and the plots read both runs as a single series.
- A solenoid (magnetic focusing) could be added via a second applied-field B map if the
  downstream linac optics are included.
