# Cathode — Space-Charge-Limited Thermionic Diode (WarpX 2D)

A WarpX model of the **electron source** at the front of the Cornell linac — Adam Bartnik's
"Region 1": a hot thermionic cathode a short distance from a positively biased grid/anode,
operating in the **space-charge-limited (SCL)** regime. Built on `pywarpx`, driven through
lume-warpx: every constant lives in `config/cathode.yaml`, and `sim/cathode.py` reads them
back, overriding only the runtime-computed values (flux, thermal velocity, `dt`, diagnostic
periods).

The cathode has a **finite transverse extent** and is simulated in 2D (x–z). The emitting
strip is much wider than the gap, so on axis we recover the 1D Child–Langmuir physics cleanly,
while the 2D run still resolves the finite-cathode edges.

Run:
```bash
conda activate CBB
python sim/cathode.py        # writes openPMD to logs/diags/cathode/{fields,particles}
python sim/plot/cathode.py   # writes the figures to logs/plots/cathode/
```

`sim/cathode.py main()` runs only the simulation; `sim/plot/cathode.py main()` runs only the
plotting (the sim must have been run first).

---

## The physics: Child–Langmuir / space-charge-limited emission

A hot cathode can supply far more current than a diode can actually transport. As electrons
leave the cathode they pile up just in front of it, and their own negative space charge
**drives the electric field at the cathode surface to zero**. This forms a *virtual cathode*
that reflects any excess emission, so the transmitted current self-regulates to the
**Child–Langmuir limit**:

$$ J_{CL} = \frac{4}{9}\,\varepsilon_0\sqrt{\frac{2e}{m_e}}\,\frac{V^{3/2}}{d^2} $$

In steady state the 1D solution has the characteristic shapes

$$ \phi(z) = V\left(\frac{z}{d}\right)^{4/3}, \qquad
   E_z(z) = -\frac{4V}{3d}\left(\frac{z}{d}\right)^{1/3} $$

— the potential is **depressed below the vacuum (linear) ramp**, and the field is **zero at
the cathode** instead of uniform.

This demo verifies WarpX reproduces this from first principles: we deliberately **over-inject
at 2 × J_CL** and let the self-consistent fields do the limiting — we do not impose the answer.
`sim/cathode.py` computes `J_CL = child_langmuir_current_density(V_anode, gap_d)` and injects a
flux of `over_inject · J_CL / q_e`; the SCL mechanism is what limits the transmitted current.

---

## What the simulation does (`sim/cathode.py`)

- **Geometry**: 2D x–z, cathode plane at `z = 0` held at 0 V, anode at `z = d` held at the
  configured anode bias. Electrons are emitted only from the finite cathode patch (the
  `lower_bound`/`upper_bound` of the flux distribution).
- **Emission**: continuous flux injection (PICMI `UniformFluxDistribution`) at `over_inject ×
  J_CL`, with a small thermal velocity spread set by the cathode temperature and a
  half-Maxwellian normal-momentum distribution (`gaussian_flux_momentum_distribution`). The
  per-component RMS thermal velocity is `thermal_velocity_sigma(T_cathode)`.
- **Solver**: electrostatic lab frame, **Multigrid** Poisson solver (`ES_MLMG_LF`) with
  Dirichlet plate potentials (`warpx_potential_lo_z` / `warpx_potential_hi_z`) and Neumann
  transverse walls. (This differs from the FFT/IGF solver used for open-boundary relativistic
  beams; here we have fixed-potential plates and non-relativistic electrons.)
- **Timestep**: `dt = CFL · (gap_d / nz) / v_final`, where `v_final = sqrt(2 q_e V / m_e)` is
  the cold final velocity through the full bias — i.e. a CFL condition on the fastest electron
  crossing one longitudinal cell.
- **Diagnostics**: field snapshots (`phi`, `rho`, `E`, `J`) and electron particle data. With
  `DIAG_PERIOD: null` the field diagnostic uses a **dense-early union slice**
  (`0:470:5, 470:MAX_STEPS:80`) so the turn-on transient is well sampled; an integer
  `DIAG_PERIOD` instead applies one uniform period to both diagnostics.
- **Output**: openPMD into `logs/diags/cathode/{fields,particles}/` (the diagnostics'
  `write_dir` is set to `logs/diags/cathode` in the YAML). Stale diagnostics are removed before
  each run because the h5 backend appends one file per dump.

### Configured operating point (`config/cathode.yaml`)

The shipped YAML values are the single Balanced operating point — there is no profile/override
machinery for this stage; **edit the YAML to retune**.

| Knob | YAML location |
|------|---------------|
| Anode / grid bias `V` | `grid: warpx_potential_hi_z` (30 V) |
| Cathode potential | `grid: warpx_potential_lo_z` (0 V) |
| Gap `d` | `grid: upper_bound[1]` (200 µm) |
| Domain half-width `W` | `grid: upper_bound[0]` (±16 mm) |
| Cathode patch `±R_cathode` | species `lower_bound`/`upper_bound[0]` (±8 mm) |
| Grid | `grid: number_of_cells` `[nx, nz]` = 128 × 64 |
| Solver tolerance | `solver: required_precision` (3e-5) |
| Macroparticles/cell | species `n_macroparticles_per_cell` (6) |
| Cathode temperature | `params: T_cathode` (1425 K) |
| Over-injection factor | `params: over_inject` (2×) |
| CFL number | `params: CFL` (0.4) |
| Diagnostic period | `params: DIAG_PERIOD` (null → dense-early union slice) |
| Max steps | `simulation: max_steps` (2000) |

The 16 mm / 0.2 mm geometry sits deep in the 1D limit (80× the gap), so the on-axis result
recovers planar Child–Langmuir while the 2D run still resolves the finite-cathode edges. It is
a **DC** demo: Region 1 actually pulses the grid voltage to chop out a bunch; this configured
anode bias is the peak cathode→grid potential difference and the grid pulsing is not modelled.

There is also a `warpx_do_not_deposit` flag on the species (default `false` = space charge ON):
**keep it `false`** — `true` disables the space-charge-limited (Child–Langmuir) mechanism this
stage exists to demonstrate, so the diode would pass the full over-injected current and the
validation figures become invalid. It is a forces-off sanity check only (the run prints a
warning).

---

## The figures (`sim/plot/cathode.py` → `logs/plots/cathode/`)

Generated in three layers — lume-warpx's plotting helpers, the shared custom figures in
`sim/plot/common.py`, and the stage-specific emission-physics figures (raw openPMD):

- **`phase_space_z_KE.png`** — `plot2D("z","kinetic_energy")`: longitudinal phase space across the gap.
- **`transverse_x_px.png`** — `plot2D("x","px")`: transverse phase space (the source's thermal emittance).
- **`potential_xz.png`** — `plot_fields("phi","x","z")`: gap potential, depressed in the beam column.
- **`charge_density_xz.png`** — `plot_fields("rho","x","z")`: the space-charge / virtual-cathode layer.
- **`centroid_vs_t.png`** — `plot1D("t","mean_z")`: the emitted cloud filling the gap.
- **`charge_vs_t.png`** — `plot1D("t","charge")`: tracked charge as emission self-limits at J_CL.
- **`energy_spectrum.png`** — `common.energy_spectrum`: charge-weighted KE histogram (⟨E⟩/σ_E) — the broad low-energy emitted spectrum.
- **`current_profile.png`** — `common.current_profile`: longitudinal current I(z) = Σ(w·v_z)/dz — the flat continuous-DC emission stream (the SCL diode emits a steady current, not a bunch).
- **`child_langmuir.png`** — on-axis φ(z) and E_z(z) vs the planar Child–Langmuir law and the vacuum reference: the space-charge depression below vacuum and the CL z^(4/3)/z^(1/3) shape (with the near-cathode field reversal).
- **`current_saturation.png`** — transmitted current density at the anode vs time, held toward J_CL despite the over-injected flux.
- **`emission_phase_space.png`** — the intrinsic thermal transverse phase space (x, p_x) + p_x histogram with ε_n,x and the ±√(kT·mₑc²) thermal scale (the source quality the gun inherits; the gun's disc remap receives ×√(3/4)).

---

## Notes & possible extensions

- The cathode is much wider than the gap (2R ≫ d), so on axis it sits in the ideal 1D planar
  limit and the J_CL agreement is tight. Shrink `R_cathode` toward `gap_d` to bring out the
  finite-cathode edge effects instead.
- Region 1 actually *pulses* the grid voltage to chop out a bunch. That can be added with a
  time-dependent potential / `AnalyticFluxDistribution`; this demo uses a DC bias to keep the
  Child–Langmuir validation clean.
