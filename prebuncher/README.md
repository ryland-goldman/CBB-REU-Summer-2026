# WarpX CESR Prebuncher (RZ)

Third stage of the Cornell Linac chain modelled in WarpX:

```
cathode (cathode/)  ->  gun (gun/)  ->  prebuncher (this)
```

The gun's exit beam (~148 keV, β ≈ 0.63, 0.1 nC, already RZ) is driven through the CESR
prebuncher — a standing-wave TM RF cavity — and bunches in the downstream drift.

## Running

```python
# from repo root, in the CBB env:
import prebuncher
prebuncher.config(POWER_KW=800, PHASE="zc",
                  OUTDIR="prebuncher/diags/P800_zc")   # optional
prebuncher.run()        # build field + sim + plots
# prebuncher.plot()     # re-generate figures from existing diags/
```

`prebuncher.run()` runs **a single case**. Defaults live at the top of
`prebuncher/prebuncher_sim.py` (`POWER_KW`, `PHASE`, `OUTDIR`); `config()` overrides them
before the run.

`build_prebuncher_field` reads `fieldmaps/prebuncher_25D.gdf`; the sim reads the gun
output from `gun/diags/particles/` (both paths repo-root-relative). To run the whole accelerator chain
(cathode → gun → prebuncher), use **`pipeline/run_pipeline.py`** in the repo root. The
results table below was produced by running several powers (`POWER_KW=160/300/500/800`) and the
`P=0` drift baseline manually. The baseline must be given an explicit
`OUTDIR="prebuncher/diags/P0_drift"`: auto-derive turns `POWER_KW=0` into `P0_zc`/`P0_crest`,
which the plotter treats as a powered case, not the baseline (see **Outputs**).

## Field map

`prebuncher_25D.gdf` (read with easygdf) is a 2.5-D axisymmetric (R, Z) map of the cavity,
**normalised to 1 J of stored energy**:

| R grid | 72 points, 0 → 36.07 mm |
| Z grid | 601 points, −152.4 → +152.4 mm (gap-centred) |
| fields | `Er`, `Ez` (V/m), `H` = Bφ (T) |

**H is the azimuthal magnetic field Bφ in Tesla** (not A/m): with H-as-Tesla the peak
E/cB ≈ 15, sensible for a resonant cavity; the A/m reading gives an unphysically negligible
B. This matches GPT's `Map25D_TM(... "Er","Ez","H" ...)`, which uses the column directly as
Bφ. (The map's energy integral over the r < 36 mm bore is only ~0.05 J — most of the cavity's
magnetic energy lives outside the beam aperture; irrelevant to the 1-J scale formula below.)

`build_prebuncher_field.py` writes the **raw 1-J map** as one openPMD file with two meshes
(`E` and `B`, both `thetaMode` m = 0, components r/t/z, shape `(1, nr, nz)`), in the layout
WarpX's `read_from_file` external-field reader expects. `grid_global_offset` shifts the
gap-centred map to lab z = `Z_GAP_CENTER` (0.20 m).

## RF drive — reproducing GPT's `Map25D_TM`

The reference GPT model drives the map as a standing-wave TM mode:

```
Er,Ez(t) = map · scale · cos(ω t + φ)
Bφ(t)    = H   · scale · sin(ω t + φ)        (E and B 90° out of phase)
```

`prebuncher_sim.py` reproduces this with `picmi.LoadAppliedField` (one openPMD path for E+B;
independent `warpx_E_time_function` = `cos`, `warpx_B_time_function` = `sin`). Constants are
baked into the parser strings because the `LoadAppliedField` wrapper forwards extra kwargs to
the picmistandard base class, which rejects them.

Parameters fixed by `reference/Linac Simulation Documentation/details.md`:

- **Frequency** `f_RF = 18 × master RF = 18 × (499.7645 MHz / 42) = 214.18 MHz`.
- **Scale** `scale = sqrt(1e3 · Q · P / (2π f_RF) / 1 J)`, loaded **Q = 3000** (prebuncher 1).
  The map's 1-J transit-weighted effective gap voltage is **V1J ≈ 439 kV**, so the physical
  gap voltage of a case is `V_gap = scale · V1J`.

Two inputs are left to the operator by design ("the user specifies dissipated power and
relative phase"):

- **Power P** — set per run (default `800 kW`; characterised over `[160, 300, 500, 800] kW`,
  ⇒ V_gap ≈ 262–586 kV; see the intrinsic-chirp threshold below for why these are higher than a
  textbook prebuncher).
- **Phase** — `zc` or `crest` (both characterised below):
  - `zc` (zero-crossing): bunch centre crosses the gap at the field zero-crossing → head
    decelerated, tail accelerated → velocity (ballistic) bunching downstream.
  - `crest`: maximum energy gain at the gap → mostly accelerates, little bunching.
  The phase is computed from the bunch-centre arrival time `t_gap = (Z_GAP_CENTER − z_inject)/v`
  so the cavity is correctly phased relative to the launched beam.

## Ballistic bunching, and the intrinsic-chirp threshold

The gun bunch is already short — σ_z ≈ 1 mm ≈ 0.1 % of the 214 MHz RF wavelength (βλ ≈
0.88 m). A 214 MHz cavity therefore cannot micro-bunch it at the RF fundamental (the bunch
sits at a single RF phase). Instead the prebuncher acts as a **ballistic buncher**: at the
zero-crossing it imparts a head→tail energy chirp across the 1 mm bunch, which compresses it
in the drift.

**Crucially, the gun beam arrives with an intrinsic +1.40 keV/mm energy chirp** (head — larger
z — already higher energy), which *debunches* on its own in a drift. The zero-crossing cavity
adds **−3.05 keV/mm per unit field scale**, so the net chirp is

```
c_net = 1.40 − 3.05 · scale   [keV/mm]   (scale = sqrt(stored_energy / 1 J))
```

A single-particle estimate gives a bunching threshold of `scale ≳ 0.46` (≈95 kW, V_gap ≈
202 kV) — already higher than a textbook low-voltage prebuncher, because the cavity must
first cancel the gun's intrinsic chirp — which is why the scanned powers are 160–800 kW.

**Space charge dominates the real result.** At 0.1 nC the ~1 mm bunch is dense: in free drift
it expands to ~19 mm over the 1.3 m line, and even the strongest cavity (800 kW) cannot
compress it *below* the injected 0.985 mm at these powers. So the meaningful, space-charge-
honest metric is **bunching relative to a drift-only baseline (P = 0)**: `σ_z,drift(z) /
σ_z,cavity(z)` at matched ⟨z⟩. The diagnostics track **σ_z(z)** (vs. the drift baseline), that
**bunching ratio**, **peak current**, and the **z–KE phase space** (the chirp flipping
negative through the cavity), not the (≈1, flat) RF-fundamental bunching factor.

## Simulation parameters

| parameter | value |
|-----------|-------|
| geometry | RZ, `n_azimuthal_modes = 1` |
| grid | `NR`=80 (r) × `NZ`=1024 (z), r ∈ [0, 36 mm], z ∈ [0, 1.30 m] |
| solver | electrostatic, lab frame, Multigrid (self-field only), `REQUIRED_PRECISION`=1e-4, `MAX_ITERS`≤500 |
| applied field | 1-J `prebuncher_25D.gdf` map × `scale` × cos/sin(ω t + φ), E + B, read from file |
| beam | gun-exit snapshot, downsampled to `MAX_PART`=50k macroparticles (reweighted), 0.1 nC, z ≈ 5 mm |
| time step | `dt = CFL · Δz / v_beam` (`CFL`=0.8; ≈ 5 ps; RF period 4.67 ns) |
| duration | `TRANSIT_MARGIN`×bunch transit of the 1.30 m domain (=0.97; crest uses the accelerated v), or fixed via `MAX_STEPS`; `N_DIAGS`=60 openPMD dumps |

**Performance knobs** (`config()`-overridable module constants; defaults reproduce the
original run). This stage is ≈75% of pipeline runtime — its self-field MLMG solve dominates
— so it is the main optimization target. Tunables: `NR, NZ` (grid), `CFL` (dt / step count),
`MAX_STEPS`/`TRANSIT_MARGIN` (run length), `REQUIRED_PRECISION` and `MAX_ITERS` (MLMG solve),
`MAX_PART` (gun-snapshot downsample), `N_DIAGS` (dump count; the plotter is robust down to
~6–8).

**Gotcha — do NOT coarsen `NZ` to go faster (measured).** Unlike the gun (near-isotropic
cells, where `nz` scales ≈ `nz²`), this box is long and thin (1.3 m × 36 mm), so its cells are
anisotropic and the MLMG Poisson solve is **convergence-bound, not cell-bound**. Lowering `NZ`
*raises* `dz/dr` (1024→512: 2.8→5.6) and slows MLMG convergence faster than it removes cells —
measured **1.37× slower per step** at `NZ=512` than at `NZ=1024`, so the fewer-steps win is more
than cancelled. And `dz` at `NZ=512` is 2.5 mm > the ~1 mm bunch, so it can't resolve the
bunching being measured. **Keep `NZ=1024`.** The effective levers are `CFL` (0.8→0.95 ⇒ ~16%
fewer steps, no per-step penalty) and `MAX_ITERS`/`REQUIRED_PRECISION` (500/1e-4 → 150/1e-3 ⇒
~20% cheaper per-step solve). `MAX_PART` barely affects runtime here (the run is solve-bound,
not particle-bound) — keep it generous for accuracy. The default Balanced profile in
`run_pipeline.py` uses `CFL=0.95, MAX_ITERS=150, REQUIRED_PRECISION=1e-3` (the space-charge field
is a small perturbation on the 148 keV beam, so the looser solve shifts the bunching only
slightly); since space charge *drives* the bunching, check the σ_z(z) figure if you loosen it
further.

WarpX's MLMG Poisson solve is memory-bandwidth bound, so these small grids gain nothing from
OpenMP threads — threads contend for the same memory bus and add fork/join overhead with no
speedup (full Balanced chain ~1.1 min at `OMP_THREADS=1`; `OMP_THREADS=14` showed no gain). The
pipeline therefore runs **single-threaded by default** (`OMP_THREADS=1`); raise `OMP_THREADS`
only for the much larger original-config grids, where per-thread work outgrows the overhead.
To compare several powers, call `prebuncher.run(plots=False)` once per power (e.g. in
a Python loop), each with its own `OUTDIR`; a final `prebuncher.plot()` then picks up every `diags/P*`
directory.

## Results

Bunching ratio = max over z of σ_z,drift / σ_z,cavity (>1 ⇒ the cavity is bunching):

| P [kW] | V_gap [kV] | zero-crossing bunching | focus z | on-crest |
|-------|-----------|------------------------|---------|----------|
| 160 | 262 | 2.05× | 1.09 m | accelerates to 402 keV |
| 300 | 359 | 2.84× | 1.09 m | → 499 keV |
| 500 | 463 | 3.89× | 0.62 m | → 604 keV |
| 800 | 586 | 5.40× | 0.49 m | → 721 keV |

The drift beam expands 0.985 → ~19 mm; the zero-crossing cavity suppresses this increasingly
with power, reaching a transient focus that moves upstream (1.09 → 0.49 m) and nearly recovers
the injected length at 800 kW. **On-crest** mainly *accelerates* the beam (KE up to 721 keV) with
no chirp-flip focus — its apparent ratio is relativistic suppression of debunching, not true
bunching (see the z–KE phase-space panels: zc flips the chirp negative; crest just shifts up in
energy). The phase-space and σ_z curves show some space-charge filamentation near the focus.

## Outputs

Each `prebuncher.run(...)` call writes `diags/P{P}_{phase}/particles/` when `OUTDIR`
is left unset — the dir is auto-derived as `P{POWER_KW:g}_{PHASE}` (so `POWER_KW=0` yields
`P0_zc`/`P0_crest`, *not* `P0_drift`). The plotter recognises the drift baseline **only** by the
exact name `P0_drift`, so to get a baseline that the comparison treats as the P=0 reference you
must pass it explicitly: `prebuncher.config(POWER_KW=0, OUTDIR="prebuncher/diags/P0_drift")`.
`prebuncher.plot()` reads every `diags/P*` directory present and writes to `results/`:

The per-case figures use **config-independent filenames** (the power/phase lives in the figure
titles and the `diags/<case>` input dir, not the filename), so changing the operating point
overwrites them in place. With several `diags/P*` cases present they are overwritten (last case
wins) — use `compare_power_phase.png` for the cross-case scan.

- `prebuncher_line.png` — σ_z(z) (vs. the drift baseline, when present, with the max-bunching point
  marked) and peak current / mean energy.
- `prebuncher_phasespace.png` — z–KE at injection / cavity exit / best focus (the
  σ_drift/σ_cavity maximum with a drift baseline, else the post-cavity σ_z minimum).
- `prebuncher_cavity.png` — the RF drive: on-axis Ez(z) of the scaled 1-J map placed at the lab
  gap, and the cos/sin RF waveform vs. time around the gap arrival (bunch centre on the field
  zero-crossing for `zc`, crest for `crest`).
- `prebuncher_bunch_profile.png` — the real longitudinal line-charge density λ(z) at the same
  three snapshots as the phase-space figure (compression and space-charge filamentation/spikes the
  scalar σ_z curve cannot show).
- `compare_power_phase.png` — σ_z(z) for the drift baseline vs. all zc powers, and the bunching
  ratio vs. power (zc vs. crest).
- a printed summary table (P, phase, σ_z0, σ_z,min, bunching ratio, focus z, I_peak, final KE).

## Notes / caveats

- We model **one** prebuncher. `details.md` describes two sharing this map (the second
  installed reversed, direction `−1,0,0`).
- **P and φ are undocumented operating points** (GUI inputs); we scan/compare rather than fix
  them.
- The **gun→prebuncher drift distance** (`Z_prebuncher1`) is undocumented; we use a short
  entrance drift with the gap at `Z_GAP_CENTER = 0.20 m`. Adjust if a real spacing is known.
- The lab-frame electrostatic self-field is mildly non-relativistic at β ≈ 0.63 (same caveat
  as the gun); acceptable for this stage.
- **Absolute compression below the injected σ_z would need lower bunch charge, a shorter
  gun→cavity drift, or higher V_gap.** Here the bunch is already short and space-charge dense,
  so we report bunching relative to the drift baseline. Set `MAX_PART` or `Z_GAP_CENTER`
  in `prebuncher_sim.py` to explore.
- If a lower-power focus falls outside the 1.30 m domain, raise `ZMAX` in `prebuncher_sim.py`.
- The drift baseline and crest curves drop sharply at the last 1–2 points as the beam clears the
  absorbing exit (few particles left); these near-empty dumps are ignored in the analysis.
