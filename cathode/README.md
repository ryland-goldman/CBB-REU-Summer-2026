# WarpX Thermionic Cathode (Space-Charge-Limited Diode)

A WarpX model of the **electron source** at the front of the Cornell Linac —
Adam Bartnik's "Region 1": a hot thermionic cathode a short distance from a
positively biased grid/anode, operating in the **space-charge-limited (SCL)**
regime. Built on `pywarpx`, driven through **lume-warpx**: every constant lives in
`cathode.yaml` and `cathode_diode.py` reads them back, overriding only the
runtime-computed values (flux, thermal velocity, `dt`, diagnostic periods). Edit
`cathode.yaml` to retune (the `config()` knob API is bypassed for this stage).

Unlike the canonical 1D [Pierce-diode example](../reference/WarpX%20Documentation/usage/examples/pierce_diode/README.md),
the cathode here has a **finite transverse extent** and is simulated in 2D (x–z).
The emitting strip (`|x| < 8 mm`) is much wider than the 0.2 mm gap, so on axis we
recover the 1D Child–Langmuir physics cleanly, while the 2D run still resolves the
finite-cathode edges.

Run with:
```bash
conda activate CBB
python -c "import cathode; cathode.run()"   # sim + plots in one call
```

or, equivalently, the individual scripts:
```bash
python cathode/cathode_diode.py   # ~1 min, writes openPMD to diags/
python cathode/plot_cathode.py    # writes the figures to results/
```

To retune the operating point, **edit `cathode/cathode.yaml`** — the `config()` knob API is
bypassed for this WarpX stage. The anode bias is `grid: warpx_potential_hi_z` (30 V); the
gap is `grid: upper_bound[1]` (200 µm); the cathode patch is the flux `lower/upper_bound`
(±8 mm); and `T_cathode` (1425 K), `over_inject` (2×), `CFL` (0.4), `DIAG_PERIOD` live in the
`params:` block.

**Performance knobs** (in `cathode.yaml`): `required_precision` (`solver:`, ships **3e-5** —
the Balanced profile; the Conservative/benchmark value is 1e-5), `n_macroparticles_per_cell`
(`species:`, ships **6** — Balanced; Conservative 10), `CFL` and `DIAG_PERIOD` (`params:`), and
the grid `number_of_cells` `[nx, nz]`. There is also a `warpx_do_not_deposit` flag on the species
(default `false` = space charge ON): **keep it `false`** — `true` disables the space-charge-limited
(Child–Langmuir) mechanism this stage exists to demonstrate, so the diode passes the full 2×J_CL
over-injection (~double the physical current) and the validation figures become invalid. It is a
forces-off sanity check only, not a meaningful cathode operating point (the run prints a warning).
The cathode is only ~7% of pipeline runtime;
**leave `DIAG_PERIOD: null`** — the field diagnostic iterates every dump over the 0–0.15 ns
turn-on window and needs the default dense-early union slice (`0:470:5, 470:MAX_STEPS:80`). An
integer `DIAG_PERIOD` applies one uniform period to both diagnostics and under-resolves it.

---

## The physics: Child–Langmuir / space-charge-limited emission

A hot cathode can supply far more current than a diode can actually transport. As
electrons leave the cathode they pile up just in front of it, and their own
negative space charge **drives the electric field at the cathode surface to zero**.
This forms a *virtual cathode* that reflects any excess emission, so the
transmitted current self-regulates to the **Child–Langmuir limit**:

$$ J_{CL} = \frac{4}{9}\,\varepsilon_0\sqrt{\frac{2e}{m_e}}\,\frac{V^{3/2}}{d^2} $$

In steady state the 1D solution has the characteristic shapes

$$ \phi(z) = V\left(\frac{z}{d}\right)^{4/3}, \qquad
   E_z(z) = -\frac{4V}{3d}\left(\frac{z}{d}\right)^{1/3} $$

— the potential is **depressed below the vacuum (linear) ramp**, and the field is
**zero at the cathode** instead of uniform.

This demo verifies WarpX reproduces this from first principles: we deliberately
**over-inject at 2 × J_CL** and let the self-consistent fields do the limiting —
we do not impose the answer.

---

## What the simulation does (`cathode_diode.py`)

- **Geometry**: 2D x–z, cathode plane at `z = 0` held at 0 V, anode at `z = d = 0.2 mm`
  (200 µm) held at `+30 V`. Electrons are emitted only from the finite cathode patch
  `|x| < 8 mm` (the `lower_bound`/`upper_bound` of the flux distribution).
- **Emission**: continuous flux injection (PICMI `UniformFluxDistribution`) at `2 × J_CL`,
  with a small thermal velocity spread set by a 1425 K cathode and a half-Maxwellian
  normal-momentum distribution (`gaussian_flux_momentum_distribution`).
- **Solver**: electrostatic lab frame, **Multigrid** Poisson solver with Dirichlet
  plate potentials (`warpx_potential_lo_z` / `warpx_potential_hi_z`) and Neumann
  transverse walls. (This differs from the FFT/IGF solver used for open-boundary
  relativistic beams; here we have fixed-potential plates and non-relativistic
  electrons.)
- **Output**: openPMD field snapshots (`phi`, `rho`, `E`, `j`) and electron
  particle data every 200 steps, into `diags/`.

| Parameter | Value |
|-----------|-------|
| Anode bias `V` | 30 V (peak grid potential = `Voff` + `Vpulse`) |
| Gap `d` | 0.2 mm (200 µm) |
| Cathode width `2R` | 16 mm (80× the gap → 1D limit on axis) |
| Cathode temperature | 1425 K |
| Injected current | 2 × J_CL ≈ 1.92 × 10⁴ A/m² |
| Child–Langmuir J_CL | ≈ 9.59 × 10³ A/m² |
| Grid | 128 × 64 cells (x, z), domain ±16 mm × 0.2 mm |
| Steps | 2000 (gap-fill ≈ 480 steps) |

These parameters now match Adam's Region-1 cathode geometry from the original LinacSim inputs
(`reference/Linac Simulation Documentation/input_files/cathode_master.in`): cathode diameter 16 mm,
cathode–grid distance 0.2 mm, and cathode temperature 1425 K. The 16 mm / 0.2 mm geometry sits deep
in the 1D limit (80× the gap) so the on-axis result recovers planar Child–Langmuir. It remains a
**DC** 2D demo — the grid voltage pulsing is not modelled.

**Gap voltage = 30 V, not `Vpulse` = 60 V.** LinacSim drives the grid with a *pulse*
`V(t) = Voff + Vpulse·f(t)` (`details.md` voltage-pulse model), where `Voff = −30 V` is the
off-bias and `Vpulse = 60 V` is the peak *swing*. The pulse shape `f` peaks at 1, so the actual
peak cathode→grid potential difference is `Voff + Vpulse = −30 + 60 = +30 V`. This DC demo uses
that peak (`V_anode = 30 V`); `Vpulse` alone (60 V) is the swing amplitude, not the absolute bias.

---

## The figures (`plot_cathode.py` → `results/`)

Generated in three layers — lume-warpx's plotting helpers, the shared custom figures in
`pipeline/plot_extras.py`, and the stage-specific emission-physics figures (raw openPMD):

- **`phase_space_z_KE.png`** — `plot2D("z","kinetic_energy")`: longitudinal phase space across the gap.
- **`transverse_x_px.png`** — `plot2D("x","px")`: transverse phase space (the source's thermal emittance).
- **`potential_xz.png`** — `plot_fields("phi","x","z")`: gap potential, depressed in the beam column.
- **`charge_density_xz.png`** — `plot_fields("rho","x","z")`: the space-charge / virtual-cathode layer.
- **`centroid_vs_t.png`** — `plot1D("t","mean_z")`: the emitted cloud filling the gap.
- **`charge_vs_t.png`** — `plot1D("t","charge")`: tracked charge as emission self-limits at J_CL.
- **`energy_spectrum.png`** — `plot_extras.energy_spectrum`: charge-weighted KE histogram (⟨E⟩/σ_E) — the broad low-energy emitted spectrum.
- **`current_profile.png`** — `plot_extras.current_profile`: longitudinal current I(z) = Σ(w·v_z)/dz — the flat continuous-DC emission stream (the SCL diode emits a steady current, not a bunch).
- **`child_langmuir.png`** — on-axis φ(z) and E_z(z) vs the planar Child–Langmuir law and the vacuum reference: the space-charge depression below vacuum and the CL z^(4/3)/z^(1/3) shape (with the near-cathode field reversal).
- **`current_saturation.png`** — transmitted current density at the anode vs time, held toward J_CL despite 2× J_CL injected.
- **`emission_phase_space.png`** — the intrinsic thermal transverse phase space (x, p_x) + p_x histogram with ε_n,x and the ±√(kT·mₑc²) thermal scale (the source quality the gun inherits; the gun's disc remap receives ×√(3/4)).

---

## Notes & possible extensions

- The cathode is much wider than the gap (2R = 16 mm ≫ d = 0.2 mm), so on axis it
  sits in the ideal 1D planar limit and the J_CL agreement is tight. Shrink
  `R_cathode` toward `gap_d` to bring out the finite-cathode edge effects instead.
- Adam's Region 1 actually *pulses* the grid voltage to chop out a bunch. That can
  be added with a time-dependent potential / `AnalyticFluxDistribution`; this demo
  uses a DC bias to keep the Child–Langmuir validation clean.

## References
- WarpX Pierce-diode example: `reference/WarpX Documentation/usage/examples/pierce_diode/README.md`
- Linac cathode model: `reference/Linac Simulation Documentation/details.md`
- Flux-injection PICMI API: `reference/WarpX Documentation/usage/python.md` (`UniformFluxDistribution`)
