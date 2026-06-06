# WarpX Thermionic Cathode (Space-Charge-Limited Diode)

A WarpX model of the **electron source** at the front of the Cornell Linac —
Adam Bartnik's "Region 1": a hot thermionic cathode a short distance from a
positively biased grid/anode, operating in the **space-charge-limited (SCL)**
regime. Built with the Python/PICMI interface (`pywarpx`).

Unlike the canonical 1D [Pierce-diode example](../reference/WarpX%20Documentation/usage/examples/pierce_diode/README.md),
the cathode here has a **finite transverse extent** and is simulated in 2D (x–z).
The emitting strip (`|x| < 6 mm`) is much wider than the 0.1 mm gap, so on axis we
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

To override the operating point without editing the source, call
`cathode.config(V_anode=..., gap_d=..., R_cathode=..., T_cathode=..., MAX_STEPS=...)`
before `cathode.run()`. Keys must match the module-level constants at the top of
`cathode/cathode_diode.py`.

**Performance knobs** (also `config()`-overridable module constants; defaults reproduce
the original run): `REQUIRED_PRECISION` (1e-5, MLMG tolerance), `MAX_ITERS` (None →
PICMI default), `PPC` (10, macroparticles/cell), `CFL` (0.4, `dt = CFL·dz/v_final`),
`DIAG_PERIOD` (None), and the grid `nx, nz`. The cathode is only ~7% of pipeline runtime;
**leave `DIAG_PERIOD=None`** — `current_saturation.png` and `rho_z_time.png` iterate every
field dump over the 0–0.15 ns turn-on window and need the default dense-early union slice
(`0:470:5, 470:MAX_STEPS:80`). An integer `DIAG_PERIOD` applies one uniform period to both
diagnostics and under-resolves those two figures.

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

- **Geometry**: 2D x–z, cathode plane at `z = 0` held at 0 V, anode at `z = d = 0.1 mm`
  (100 µm) held at `+50 V`. Electrons are emitted only from the finite cathode patch
  `|x| < 6 mm` (the `lower_bound`/`upper_bound` of the flux distribution).
- **Emission**: continuous flux injection (PICMI `UniformFluxDistribution`) at `2 × J_CL`,
  with a small thermal velocity spread set by a 1200 K cathode and a half-Maxwellian
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
| Anode bias `V` | 50 V |
| Gap `d` | 0.1 mm (100 µm) |
| Cathode width `2R` | 12 mm (120× the gap → 1D limit on axis) |
| Cathode temperature | 1200 K |
| Injected current | 2 × J_CL ≈ 1.65 × 10⁵ A/m² |
| Child–Langmuir J_CL | ≈ 8.25 × 10⁴ A/m² |
| Grid | 128 × 64 cells (x, z), domain ±12 mm × 0.1 mm |
| Steps | 2000 (gap-fill ≈ 480 steps) |

These parameters are a deliberately scaled-down 2D demo of Adam's Region-1 cathode (whose
diameter is 16 mm): the 50 V / 100 µm operating point is chosen to sit deep in the 1D limit so
the on-axis result recovers planar Child–Langmuir, and is not his actual operating geometry.

---

## The figures (`plot_cathode.py` → `results/`)

### `child_langmuir.png` — the validation
On-axis (center of the cathode) `φ(z)` and `Ez(z)` overlaid with the
Child–Langmuir laws and the vacuum reference. The WarpX curve sits right on the
4/3-power potential, and the field is **driven to ~0 at the cathode** — the
defining signature of space-charge-limited emission.

### `cathode_2d.png` — the 2D structure
Maps of charge density, potential, and `|E|`. You can see (1) the dense
space-charge / virtual-cathode layer hugging the emitting strip, (2) the potential
depression in the beam column, and (3) the **field transition at the cathode edges**
`x = ±6 mm`, where the field-suppressed emitting strip meets the full vacuum field
outside — the finite-cathode signature absent from planar theory.

### `current_saturation.png` — self-limiting
Transmitted current (integrated across the beam, referenced to the cathode width)
vs. time. Despite injecting **2× J_CL**, the transmitted current self-limits to the
Child–Langmuir scale — it settles near J_CL (slightly above the cold-emission value,
≈ 110% in this run, with the finite cathode temperature and near-1D geometry). The
cathode does **not** pass the 2× current it is fed — space charge regulates it.

### `rho_z_time.png` — space-charge cloud build-up
On-axis charge density `|ρ|(z, t)` (√ scale) over the turn-on transient: the
space-charge cloud building up and filling the gap (gap-fill ≈ 480 steps), drawn
with `pcolormesh` on the true (non-uniform) time coordinates.

### `field_lines.png` — the 2D cathode-edge field transition
φ equipotential contours + E-field streamlines across the gap, with a zoom on the `+x` edge. At the
cathode edges `x = ±6 mm` the equipotentials **crowd together** and the streamlines **splay** as
`|E|` climbs from its space-charge-suppressed value on the emitting surface to the full vacuum field
outside — the field **transition** at the emission edge (monotonic, no overshoot above `V/d`), the
finite-cathode effect the 1D Child–Langmuir picture omits. (Contour companion to the `φ` panel of
`cathode_2d.png`.)

### `emission_phase_space.png` — intrinsic thermal emittance
Transverse phase space `x` vs. `ux = γβ_x` and the histogram of `ux`, from the last particle
snapshot. The RMS normalized emittance `εn,x ≈ 1.57 mm·mrad` (annotated) is the source's intrinsic
thermal emittance, set by the 1200 K cathode — the beam quality handed to the gun. The run
reproduces the expected thermal momentum spread `√(kT/mₑc²)`.

---

## Notes & possible extensions

- The cathode is much wider than the gap (2R = 12 mm ≫ d = 0.1 mm), so on axis it
  sits in the ideal 1D planar limit and the J_CL agreement is tight. Shrink
  `R_cathode` toward `gap_d` to bring out the finite-cathode edge effects instead.
- Adam's Region 1 actually *pulses* the grid voltage to chop out a bunch. That can
  be added with a time-dependent potential / `AnalyticFluxDistribution`; this demo
  uses a DC bias to keep the Child–Langmuir validation clean.

## References
- WarpX Pierce-diode example: `reference/WarpX Documentation/usage/examples/pierce_diode/README.md`
- Linac cathode model: `reference/Linac Simulation Documentation/details.md`
- Flux-injection PICMI API: `reference/WarpX Documentation/usage/python.md` (`UniformFluxDistribution`)
