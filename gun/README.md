# CESR gun in WarpX (RZ)

The second stage of the Cornell Linac electron source, simulated in WarpX. Stage 1
(`../cathode/`) is the thermionic cathode at the Child–Langmuir limit; here we take
its emitted electrons and accelerate them through the **CESR gun** — the electrostatic
accelerating structure modelled in Adam Bartnik's Linac GUI with the Poisson–Superfish field
map `CESR_gun.gdf` (the "Chili Gun Mk II", ~150 kV).

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

To override the gun voltage or bunch charge: `gun.config(GUN_VOLTAGE=150e3,
BUNCH_CHARGE=1.0e-9)` before `gun.run()`. Keys must match the module-level constants in
`gun/build_gun_field.py` and `gun/gun_sim.py`. `build_gun_field.py` reads
`fieldmaps/CESR_gun.gdf`; `gun_sim.py` reads the cathode output from
`cathode/diags/particles/`. All paths are repo-root-relative.

**Performance knobs** (`config()`-overridable module constants; defaults reproduce the
original run): `REQUIRED_PRECISION` (1e-5) and `MAX_ITERS` (None) for the MLMG solve;
`CFL` (0.4, `dt = CFL·dz/v_exit`), `TRANSIT_MARGIN` (1.15) and `AVG_SPEED_FRAC` (0.6) for
the auto-derived run length, or `MAX_STEPS` (>0) to fix it (`AVG_SPEED_FRAC=0.6` is
hand-tuned for the 150 kV point — `v_exit` is recomputed from `GUN_VOLTAGE`, but the
average-speed fraction is not, so re-check it if `GUN_VOLTAGE` is changed substantially
via `config()`); `N_DIAGS` (40) for the openPMD
dump count; `MAX_PART` (0 = no cap) to downsample the imported cathode bunch (reweighted,
charge-preserving); `BEAM_RELEASE` (`"timed"`/`"snapshot"`) and `PULSE_WIDTH` (2 ns) for the
beam representation (see *Beam source*; `"timed"` runs ~5× longer than `"snapshot"` because the
run spans the full pulse + a transit instead of one transit); the grid `nr, nz`; and
`SPACE_CHARGE` (default `True`). Setting
`SPACE_CHARGE=False` passes `warpx_do_not_deposit` (beam self-field off, only the applied gun field
acts) — but note the self-field is *dominant* here at 146 keV (it "dwarfs the gun field," ~17% of
charge is already lost to it), so SC-off is a large physics change, not a mild diagnostic. Runtime ≈
`nz²` (per-step cost ∝ cells, and
`dz = ZMAX/nz` ⇒ fewer steps as `nz` drops), so halving `nz` ≈ 4× faster. This holds because the
gun's cells are near-isotropic (`dz/dr ≈ 1.3`) so the MLMG solve stays well-conditioned as `nz`
drops — **unlike the injector's long-thin box**, where coarsening `NZ` slows the solve instead
(see `injector/README.md`). Keep `N_DIAGS ≥ 20` so `space_charge.png` still finds its
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
transports and accelerates to ~146 keV. Set `BUNCH_CHARGE` at the top of `gun_sim.py` to
explore the space-charge regime.

### Beam representation — time-release vs snapshot (`BEAM_RELEASE`)

Even at 1 nC, *how* the bunch is fed to the gun matters. The CESR gun is gated by a **2 ns
grid pulse** (`cathode_master.in` `twidth=2`), which is ~four gun-transit-times long
(transit ≈ 0.47 ns), so the physical beam is a **long, low-density, quasi-DC stream** — the
original LinacSim GPT deck emits it that way via `settdist("beam","F",…,"t",…)`. Two modes:

- **`"timed"` (default — the realistic representation).** The imported macroparticles are
  released over `PULSE_WIDTH` (=2 ns) by a per-step `installbeforestep` injection callback
  (`ParticleContainerWrapper.add_particles`), so only the fraction of charge already emitted
  is present at any instant — the physical low line-density. This is the WarpX counterpart of
  GPT's native time-release and of the cross-code benchmark's `warpx_tr.py`.
- **`"snapshot"`.** All 1 nC is present at t=0, initially crammed into the ~0.2 mm cathode-exit
  z-extent — a line-charge density orders of magnitude above the real 2 ns beam. It
  **over-states the space-charge force**: the WarpX–GPT benchmark put the beam-representation
  effect at ~28 % on emittance, and here the over-dense snapshot additionally blows ~19 % of
  the beam off as halo (transmission ≈ 81 % vs ≈ 100 % for the same grid/charge time-released).
  Kept for speed and back-compat only — it is **not** a realistic operating point.

**Exit-beam handoff (timed mode).** The released beam's ballistic z-extent
(~v_exit·`PULSE_WIDTH` ≈ 0.38 m) is many times the 51.77 mm gun domain, so no single
volumetric snapshot can hold it. `build_exit_handoff()` reconstructs the full exit beam the
same way `pipeline/collimator.py` tracks the iris scrape — by particle **id** across the
volumetric dumps: each macroparticle's last appearance is found, particles last seen at the
pipe wall (`r ≥ RMAX`) are dropped as radial loss, and the rest are drifted ballistically to a
common reference time (rebuilding one consistent snapshot, head downstream / tail at the
entrance). The result is written to `gun/diags/handoff` as openPMD (via
`pipeline.impact_io.write_openpmd_particles`); the injector's `load_gun_bunch` reads it when
present, else falls back to the volumetric `gun/diags/particles` (legacy snapshot). **Caveat:**
the injector's prebuncher phases were tuned against the compact snapshot handoff, so the longer
timed beam shifts its operating point and should be re-validated (the LinacSim
input-reconciliation backlog) — the handoff is wired, the re-tuning is a follow-up.

**Approximations.** The cathode model is a 2D Cartesian slab, not RZ, so the slab→radius
remap is an approximation: the `r`-importance resample (step 2) makes the **areal density**
match the cathode's radial profile, but the reconstructed azimuthal distribution is assumed
uniform (it cannot recover the true cylindrical emission, which the 2D slab never had). The DC
beam is treated as a single injected bunch.

## Simulation parameters (`gun_sim.py`)

| parameter | value |
|-----------|-------|
| geometry | RZ, `n_azimuthal_modes = 1` |
| grid | 128 (r) × 512 (z), r ∈ [0, 15 mm], z ∈ [0, 51.77 mm] (finer than the old 96×384 to resolve the near-cathode dynamics — see *Cross-code validation*) |
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
At the gun exit (146 keV, γ ≈ 1.29) the plain labframe-electrostatic solver would **overestimate
the transverse space-charge force by ≈ γ² = 1.66×, i.e. ~66 %** (ramping from a few % near the
cathode at 10 keV to ~66 % at exit); the electromagnetostatic solver removes that error
self-consistently. (WarpX's per-species *relativistic* ES mode is an alternative for a single
drifting species; the electromagnetostatic mode generalizes to multi-velocity beams.)

The magnetostatic vector-Poisson solve reuses the `REQUIRED_PRECISION` / `MAX_ITERS` knobs (passed
explicitly as `warpx_magnetostatic_required_precision` / `warpx_magnetostatic_max_iters`). It adds
a 3-component MLMG solve per step, so the run is **≈2× slower** than the electrostatic-only gun
(~80 s vs ~40 s at the default grid). **Boundary requirement:** the outer radial wall is
`dirichlet` (not `neumann`) — the `A_z` component (driven by the dominant beam current `j_z`)
would otherwise have an all-Neumann singular operator and the MLMG bottom solve **diverges**
(`MLMG failed`). Grounding the pipe at r = 15 mm — well outside the r ≲ 8 mm beam — makes it
well-posed; the headline exit energy is unchanged (146 keV) and transmission is ~81 %.

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
  is fine enough near the cathode (the benchmark converged by nz ≈ 720 over 55 mm); the old
  nz = 384 over 51.77 mm sat on the unconverged side, so the grid here is **128 × 512**.
- **Match the beam representation.** The benchmark's single largest term (~28 %) is snapshot vs
  time-release — implemented here as `BEAM_RELEASE` (see *Beam source*).

With matched physics and converged numerics the benchmark found WarpX and GPT consistent to
within statistics (~3 %) on a controlled snapshot, so these settings are the cross-validated
recipe, not just a WarpX convention.

## Figures (`results/`)

1. **`gun_field.png`** — on-axis `Ez(z)` and implied potential of the scaled field map: the
   accelerating field the beam sees (Ez < 0, 150 kV total drop).
2. **`beam_rz.png`** — `r–z` beam distribution at launch / mid-gun / exit: transport through
   the gun, including the near-cathode radial focusing.
3. **`energy_gain.png`** — mean and max kinetic energy vs. ⟨z⟩, climbing toward ~150 keV.
4. **`exit_phase_space.png`** — longitudinal `z–KE` phase space and the energy spectrum at the
   last dump.
5. **`beam_envelope.png`** — per-plane RMS size `σ_x = √⟨x²⟩` and normalized transverse emittance
   `εn,x` vs. `⟨z⟩`: the near-cathode radial focusing of `beam_rz.png` made quantitative, plus the
   space-charge / aberration emittance growth along the gun. (`σ_x` is the single-plane RMS that
   pairs with `εn,x`; the radial RMS is `√⟨r²⟩ = √2·σ_x`.)
6. **`space_charge.png`** — `r–z` maps of the beam **self-field** (`ρ` and the space-charge
   potential well `φ`, ≈ −250 V) at a near-launch snapshot — the dumped self-field nothing else
   plots, and the well that motivates renormalizing the bunch to 1 nC.

## Notes / extensions

- The beam energy gain tracks `∫ e·|Ez| dz` (≈ 7.5 keV by z ≈ 4 mm), approaching the ~150 keV
  set by the cathode→exit potential drop (the space-charge-loaded beam lands at ~146 keV mean).
- To approach the continuous-emission picture, inject a train of bunches or feed the cathode
  current directly rather than a single snapshot.
- A solenoid (magnetic focusing) could be added via a second `read_from_file` B map if the
  downstream Linac optics are included.
- **Fresh diags on rerun:** WarpX *appends* one openPMD file per dump, so `gun_sim.py`
  `shutil.rmtree`s `gun/diags/` at the start of each run. Without this, re-running with a
  different grid/step count (hence different diag step numbers) leaves stale files that
  interleave with the new ones; the plots then read both runs as a single series and show a
  fan of overlapping curves. (Mirrors `injector_sim.py` / `linac_sec1_sim.py`.)
