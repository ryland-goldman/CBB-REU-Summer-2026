# CESR gun in WarpX (RZ)

The second stage of the Cornell Linac electron source, simulated in WarpX. Stage 1
(`../cathode/`) is the thermionic cathode at the Child–Langmuir limit; here we take
its emitted electrons and accelerate them through the **CESR gun** — the electrostatic
accelerating structure modelled in Adam Bartnik's Linac GUI with the Poisson–Superfish field
map `CESR_gun.gdf` (the "Chili Gun Mk II", ~150 kV).

Driven through **lume-warpx**: every constant lives in `gun.yaml` and `gun_sim.py` reads them
back, overriding only the runtime-computed values (`dt`, step count, the seed arrays, diagnostic
periods) and building the time-release beam with a `beforestep` callback. Edit `gun.yaml` to
retune (the `config()` knob API is bypassed for this stage).

The gun field is applied as an external **electrode field** on the particles; WarpX's
**electromagnetostatic** solver supplies the self-consistent beam **space charge** on top —
both the electrostatic self-field (∇²φ = −ρ/ε₀) and, from the beam current, the self
**magnetic** field (Coulomb-gauge ∇²A = −μ₀j, B = ∇×A), so the relativistic magnetic pinch is
included (see *Space-charge model* below). Geometry is **RZ (cylindrical)**, matching the field
map's native symmetry.

## Pipeline

```bash
conda activate CBB
python -c "import gun; gun.run()"   # build field map + sim + plots in one call
```

or, equivalently, the individual scripts (module form — `gun_sim.py` and `plot_gun.py`
import `pipeline.*`, which is only on `sys.path` when launched from the repo root):
```bash
python -m gun.build_gun_field   # CESR_gun.gdf  ->  gun_field/gun_E.h5 (openPMD)
python -m gun.gun_sim           # RZ WarpX run  ->  diags/{fields,particles}/
python -m gun.plot_gun          # figures       ->  results/*.png
```

`GUN_VOLTAGE` is a `build_gun_field.py` module constant (it sets the field-map scale), so
`gun.config(GUN_VOLTAGE=150e3)` before `gun.run()` still works. **Every other knob lives in
`gun/gun.yaml`** — the `config()` API is bypassed for this WarpX stage; to retune, edit the YAML.
`build_gun_field.py` reads `fieldmaps/CESR_gun.gdf`; `gun_sim.py` reads the cathode output from
`cathode/diags/particles/`. All paths are repo-root-relative.

**Performance knobs** (in `gun.yaml`): `required_precision` / `warpx_magnetostatic_required_precision`
(`solver:`, ship **1e-4** — the Balanced profile; the Conservative/benchmark value is 1e-5) for the
MLMG solve; in `params:`, `CFL` (0.4, `dt = CFL·dz/v_exit`), `TRANSIT_MARGIN` (1.15) and
`AVG_SPEED_FRAC` (0.6) for the auto-derived run length, or `MAX_STEPS` (>0) to fix it
(`AVG_SPEED_FRAC=0.6` is hand-tuned for the 150 kV point — `v_exit` is recomputed from
`GUN_VOLTAGE`, but the average-speed fraction is not, so re-check it if `GUN_VOLTAGE` changes
substantially); `N_DIAGS` (40) for the openPMD dump count; `MAX_PART` (0 = no cap) to downsample
the imported cathode bunch (reweighted, charge-preserving); `BEAM_RELEASE` (`"timed"`/`"snapshot"`)
and `PULSE_WIDTH` (2 ns) for the beam representation (see *Beam source*; `"timed"` runs ~5× longer
than `"snapshot"`); and the grid `number_of_cells` `[nr, nz]`. The species
`warpx_do_not_deposit` flag (default `false` = space charge ON): `true` turns the self-field off
(only the applied gun field acts) — but the self-field is *dominant* here at 149 keV (it "dwarfs
the gun field," ~17% of charge is already lost to it), so it is a large physics change, not a mild
diagnostic. Runtime ≈ `nz²` (per-step cost ∝ cells, and `dz = ZMAX/nz` ⇒ fewer steps as `nz`
drops), so halving `nz` ≈ 4× faster — the reason the **shipped Balanced default is `nz = 384`**
(dz ≈ 0.19 mm). The cross-code-converged grid is `nz = 712` (dz ≈ 0.10 mm; see *Cross-code
validation* and *Simulation parameters*) — set it for a Conservative near-cathode-εn,x run. The
gun's cells are near-isotropic so the MLMG solve stays well-conditioned at either `nz` — **unlike
the injector's long-thin box**, where coarsening `NZ` slows the solve (see `injector/README.md`).
The `fields` diagnostic (`Ez_rz.png`, `self_charge_rz.png`) needs a few dumps to be meaningful, so
keep `N_DIAGS` reasonable.

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

`build_gun_field.py` writes the scaled field as an openPMD file (via the shared
`pipeline/fieldio.py` writer) in this RZ layout: geometry `thetaMode` with a single
`m = 0` mode, mesh record `E` with components `r`,`t`,`z`, axis labels `["r","z"]`, dataset
shape `(1, nr, nz)`. This `["r","z"]`/`m=0` order is a **deliberate, reader-validated
deviation** from WarpX's own RZ field *diagnostic*, which emits the opposite (`["z","r"]`,
`m=1;imag=+`, shape `(modes, nz, nr)`); both load correctly because the `read_from_file`
reader is axis-labels-aware, so this convention is a choice (one shared across all three
build scripts), not something the reader strictly requires. `gun_sim.py` then loads it via
the raw WarpX inputs

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
snapshot encode the steady-state population *in transit through the diode* (~82 nC), not a
bunch charge. We:

1. Import the emitted **phase-space distribution** (positions + momenta) from the last
   cathode snapshot.
2. Remap the 2D `(x, z)` slab into RZ: treat `|x|` as the radius `r` and smear the particles
   uniformly in azimuth (`x = r cosθ, y = r sinθ`), rotating the transverse momentum
   accordingly. Crucially, the revolution carries a **2πr Jacobian**: a slab uniform in `x`
   has a flat `dN/dr`, which—revolved naively with `r = |x|` and unchanged weight—would give
   areal density `n(r) ∝ 1/r`, a spurious on-axis charge cusp. We therefore
   **importance-resample by `r`** (draw particles with probability ∝ `r·w`, with
   replacement — ≡ ∝ `r` for the cathode's uniform weights), so `dN/dr → r·dN/dr` and
   `n(r)` matches the cathode's true radial profile (a flat-top emitting strip → a
   uniform-density disc). This keeps the macroparticle weights uniform. Because the draw
   is with replacement, macroparticles are duplicated and the effective independent
   sample count is below the drawn count — relevant if `MAX_PART` is set small.
3. **Renormalize** the total weight to a physical gun bunch charge `BUNCH_CHARGE = 1 nC`
   (the CESR gun is grid-pulse gated; 1 nC matches the original LinacSim `gpt_master.in`
   `total_charge = -1e-9`).

**Why renormalize:** injecting the full 82 nC as one instantaneous bunch is unphysical — its
radial space-charge field dwarfs the gun field and blows the beam apart before it
accelerates (observed directly: the beam is absorbed within ~50 steps). At 1 nC the beam still
transports and accelerates to ~149 keV. Set `BUNCH_CHARGE` at the top of `gun_sim.py` to
explore the space-charge regime.

### Beam representation — time-release vs snapshot (`BEAM_RELEASE`)

Even at 1 nC, *how* the bunch is fed to the gun matters. The CESR gun is gated by a **2 ns
grid pulse** (`cathode_master.in` `twidth=2`), which is ~four gun-transit-times long
(transit ≈ 0.47 ns), so the physical beam is a **long, low-density, quasi-DC stream** — the
original LinacSim GPT deck emits it that way via `settdist("beam","F",…,"t",…)`. Two modes:

- **`"timed"` (default — the realistic representation).** The imported macroparticles are
  released over `PULSE_WIDTH` (=2 ns) by a per-step `installbeforestep` injection callback
  (`ParticleContainerWrapper.add_particles`), so only the fraction of charge already emitted
  is present at any instant — the physical low line-density. This reproduces the *mechanism* of
  GPT's native time-release (`settdist`) and the cross-code benchmark's `warpx_tr.py`. It does
  **not** reproduce the pulse *profile*: emission times are drawn uniformly (flat-top), whereas
  the real 2 ns pulse has ~30 V/ns edge ramps (`cathode_master.in` `Vp=30`) that emit less at
  the edges; the dominant correction is the line-density drop, not the edge shape, but the
  trapezoidal profile is not modelled. The cathode run is a *stationary* DC emitter with no
  per-particle birth time, so every release time carries the **same** steady-state cathode phase
  space — a DC line-density reconstruction, not a resolved emission history.
- **`"snapshot"`.** All 1 nC is present at t=0, initially crammed into the ~0.2 mm cathode-exit
  z-extent — a line-charge density orders of magnitude above the real 2 ns beam, which
  **over-states the space-charge force.** Two consequences are *confounded* here and should not
  be read as one clean number: the over-dense snapshot both (a) inflates emittance (the WarpX–GPT
  benchmark's *controlled* beam-representation term is ~28 % on a matched core) **and** (b) blows
  the high-amplitude halo to the pipe wall, so the snapshot loses ~19 % of the beam that the
  time-released run keeps (transmission ≈ 81 % snapshot vs ≈ 100 % timed, same grid/charge).
  Because the two modes' εn,x are then computed on *different* particle sets (timed keeps the
  halo snapshot scrapes), their raw εn,x difference is **not** a clean measure of the
  space-charge over-statement — the controlled benchmark number is. Snapshot is kept for speed
  and back-compat only — it is **not** a realistic operating point.

**Exit-beam handoff (timed mode).** The released beam's ballistic z-extent
(~v_exit·`PULSE_WIDTH` ≈ 0.4 m) is many times the gun domain, so no single volumetric snapshot
can hold it. `build_exit_handoff()` reconstructs the full exit beam by particle **id** across
the volumetric dumps (the `pipeline/collimator.py` idiom), sampling each particle in the
**field-free pad past the field map** (z ≥ `ZMAX_FIELD`): its first appearance there is its
exit-plane phase space, and a ballistic drift to a common reference time then rebuilds one
consistent snapshot (head downstream / tail at the entrance, matching the injector's
`z − z.min() + Z_INJECT`). Sampling in the **field-free** pad — not at the particle's last
in-field dump — is essential: a field-free drift preserves εn,x, but drifting a still-in-field
particle as if field-free manufactures a spurious x–u correlation and inflates εn,x ~8× (the
reason the domain carries the `ZPAD` drift pad). Note εn,x of the reconstructed *instantaneous*
beam is the **projected** emittance of the ~0.4 m drifting bunch (head drifted longer than tail):
≈ 45 mm·mrad, physically ~8× the **per-slice** beam quality (≈ 5.7 mm·mrad steady-state — the
drift preserves per-slice εn,x). That head–tail projection is real (an instantaneous diagnostic of the
2 ns beam would measure it), not an artifact; the injector's prebunchers compress the beam
longitudinally and recover it. Do **not** read the 45 mm·mrad handoff number as the beam quality. Particles that never reach the pad are classified as a
radial (`r=RMAX`) loss or, for the last sliver of the pulse, an un-flushed tail (run ends first,
~2 %), and counted. The result is written to `gun/diags/handoff` (openPMD via
`pipeline.impact_io.write_openpmd_particles`); the injector's `load_gun_bunch` reads it when
present, else the volumetric `gun/diags/particles` (legacy snapshot). **Caveat:** the injector's
prebuncher phases were tuned against the compact snapshot handoff, so the longer timed beam
shifts its operating point and should be re-validated (the LinacSim input-reconciliation
backlog) — the handoff is wired, the re-tuning is a follow-up.

**Approximations.** The cathode model is a 2D Cartesian slab, not RZ, so the slab→radius
remap is an approximation: the `r`-importance resample (step 2) makes the **areal density**
match the cathode's radial profile, but the reconstructed azimuthal distribution is assumed
uniform (it cannot recover the true cylindrical emission, which the 2D slab never had). The DC
beam is treated as a single injected bunch.

## Simulation parameters (`gun_sim.py`)

| parameter | value |
|-----------|-------|
| geometry | RZ, `n_azimuthal_modes = 1` |
| grid | 128 (r) × **384 (z)** shipped Balanced default (dz≈0.19 mm); **712** is the cross-code-converged grid (dz≈0.10 mm), set it in `gun.yaml` `number_of_cells` for a Conservative near-cathode-εn,x run. r ∈ [0, 15 mm], z ∈ [0, 71.77 mm] = field map (`ZMAX_FIELD`=51.77 mm) + `ZPAD`=20 mm field-free drift pad, padded so the exit beam is sampled in field-free space (see *Cross-code validation* and *Beam source*) |
| boundaries (fields) | axis at r=0; **dirichlet** (grounded) outer radial wall + both z plates |
| solver | electromagnetostatic, lab frame, Multigrid (self E **and** B from beam current) |
| applied field | scaled `CESR_gun.gdf`, −150 kV, read from file |
| beam representation | `BEAM_RELEASE="timed"` (default): release the 1 nC over the 2 ns grid pulse (`PULSE_WIDTH`); `"snapshot"`: all at t=0 (over-states space charge) — see *Beam source* |
| bunch | 1 nC, imported cathode phase space, ~133k macroparticles (optionally capped by `MAX_PART`, reweighted) |
| time step | `dt = CFL·Δz/v_exit` (`CFL`=0.4; v_exit ≈ 0.63 c at 150 keV) |
| duration | snapshot: `TRANSIT_MARGIN`×gun-transit (=1.15; bunch average speed ≈ `AVG_SPEED_FRAC`·v_exit, =0.6); timed: `PULSE_WIDTH + TRANSIT_MARGIN`×transit so the last-released particle clears the gun. Or fixed via `MAX_STEPS`. Snapshot stops as the beam reaches the exit (running longer empties the domain and aborts the Multigrid solve); timed stays populated until the pulse flushes, so the solve is never charge-starved mid-run |

### Space-charge model — electromagnetostatic self-field

The self-field uses WarpX's **lab-frame electromagnetostatic** solver (`warpx_magnetostatic=True`
on the PICMI `ElectrostaticSolver`): on top of the electrostatic Poisson solve (∇²φ = −ρ/ε₀,
**E** = −∇φ) it also solves the Coulomb-gauge vector potential from the beam current
(∇²**A** = −μ₀**j**, **B** = ∇×**A**), so the beam's **self magnetic field** is included. The
resulting magnetic-pinch force `qβ×B` partially cancels the radial electric repulsion, giving the
correct relativistic net transverse self-force `qE_r/γ²` rather than the pure-electrostatic `qE_r`.
At the gun exit (149 keV, γ ≈ 1.29) the plain labframe-electrostatic solver would **overestimate
the transverse space-charge force by ≈ γ² = 1.66×, i.e. ~66 %** (ramping from a few % near the
cathode at 10 keV to ~66 % at exit); the electromagnetostatic solver removes that error
self-consistently. (WarpX's per-species *relativistic* ES mode is an alternative for a single
drifting species; the electromagnetostatic mode generalizes to multi-velocity beams.)

The magnetostatic vector-Poisson solve reuses the `REQUIRED_PRECISION` / `MAX_ITERS` knobs (passed
explicitly as `warpx_magnetostatic_required_precision` / `warpx_magnetostatic_max_iters`). It adds
a 3-component MLMG solve per step, so it is **≈2× slower per step** than an electrostatic-only
gun. (Absolute runtime is dominated by the `timed` default — full pulse + the finer 128×712
padded grid run in minutes, not the old ~80 s; use `BEAM_RELEASE="snapshot"` and a coarser grid
for a quick check.) **Boundary requirement:** the outer radial wall is
`dirichlet` (not `neumann`) — the `A_z` component (driven by the dominant beam current `j_z`)
would otherwise have an all-Neumann singular operator and the MLMG bottom solve **diverges**
(`MLMG failed`). Grounding the pipe at r = 15 mm — well outside the r ≲ 8 mm beam — makes it
well-posed (it also models φ as a grounded conductor rather than a Neumann mirror, physically
the real beampipe); the headline exit energy is unchanged (149 keV). Transmission depends on the
beam representation: ≈ 81 % for the over-dense `snapshot` beam (halo blown to the wall), ≈ 100 %
for the realistic `timed` release (the low line-density beam barely diverges) — see *Beam
representation*.

### Cross-code validation (WarpX vs GPT)

A from-scratch WarpX–GPT cross-code benchmark of this exact 150 kV gun
(`CornellMisc/.../bench/writeup`, "Reconciling WarpX and GPT on the 150 kV gun emittance")
reconciled a 30–40 % emittance disagreement between the two codes and, in doing so, validated
the physics choices this stage makes:

- **Relativistic / magnetic-pinch space charge.** A lab-frame ES solver over-states the
  transverse space-charge force by ≈γ² (the dropped magnetic cancellation). The benchmark's fix
  is WarpX's *relativistic* ES; this stage's **electromagnetostatic** solver achieves the same
  by solving the self-**B** explicitly (qβ×B pinch ⇒ net qE_r/γ²), and generalizes to a
  multi-velocity beam. (The benchmark integrates the γ² error to ~+2 % on its short 10 ps bunch;
  it grows with γ and line length.)
- **Dirichlet cathode image, *not* Neumann.** The metal cathode is a Dirichlet plane whose
  opposite-sign image partly cancels the bunch self-field; a Neumann plane gives a *same*-sign
  image that adds to it (a +12 % error in the benchmark). This stage's z = 0 plate is dirichlet.
- **Grid must resolve the near-cathode dynamics.** εn,x falls and only converges once the grid
  is fine enough near the cathode (the benchmark converged by nz ≈ 720 over 55 mm); nz = 384 over
  51.77 mm sat on the unconverged side, so the **converged grid is 128 × 712** (over the 71.77 mm
  field-free-padded domain; dz≈0.10 mm). The shipped **Balanced** default coarsens this to
  `nz = 384` (dz≈0.19 mm) for ~4× speed — the same effective default the pre-YAML pipeline ran via
  its Balanced `config()` block; set `nz = 712` in `gun.yaml` for a converged-εn,x Conservative run.
- **Match the beam representation.** The benchmark's single largest term (~28 %) is snapshot vs
  time-release — implemented here as `BEAM_RELEASE` (see *Beam source*).

With matched physics and converged numerics the benchmark found WarpX and GPT consistent to
within statistics (~3 %) on a controlled snapshot, so these settings are the cross-validated
recipe, not just a WarpX convention.

## Figures (`results/`)

Generated by `plot_gun.py` entirely with lume-warpx's plotting helpers:

1. **`phase_space_z_KE.png`** — `plot2D("z","kinetic_energy")`: exit longitudinal phase space, the beam at ~150 keV.
2. **`transverse_x_px.png`** — `plot2D("x","px")`: exit transverse phase space.
3. **`Ez_rz.png`** — `plot_fields("E","z","r")`: field in the gun gap (applied electrode + self-field).
4. **`self_charge_rz.png`** — `plot_fields("rho","z","r")`: the beam self charge density.
5. **`centroid_vs_t.png`** — `plot1D("t","mean_z")`: the bunch marching down the gun.
6. **`emittance_vs_t.png`** — `plot1D("t","norm_emit_x")`: normalized transverse emittance over the run.

## Notes / extensions

- The beam energy gain tracks `∫ e·|Ez| dz` (≈ 7.5 keV by z ≈ 4 mm), approaching the ~150 keV
  set by the cathode→exit potential drop (the space-charge-loaded beam lands at ~149 keV mean,
  cross-validated against GPT at 148.9 keV — see *Cross-code validation*).
- To approach the continuous-emission picture, inject a train of bunches or feed the cathode
  current directly rather than a single snapshot.
- A solenoid (magnetic focusing) could be added via a second `read_from_file` B map if the
  downstream Linac optics are included.
- **Fresh diags on rerun:** WarpX *appends* one openPMD file per dump, so `gun_sim.py`
  `shutil.rmtree`s `gun/diags/` at the start of each run. Without this, re-running with a
  different grid/step count (hence different diag step numbers) leaves stale files that
  interleave with the new ones; the plots then read both runs as a single series and show a
  fan of overlapping curves. (Mirrors `injector_sim.py` / `linac_sec1_sim.py`.)
