# WarpX Thermionic Cathode (Space-Charge-Limited Diode)

A WarpX model of the **electron source** at the front of the Cornell Linac —
Adam Bartnik's "Region 1": a hot thermionic cathode a short distance from a
positively biased grid/anode, operating in the **space-charge-limited (SCL)**
regime. Built with the Python/PICMI interface (`pywarpx`).

Unlike the canonical 1D [Pierce-diode example](../reference/WarpX%20Documentation/usage/examples/pierce_diode/README.md),
the cathode here has a **finite transverse extent** and is simulated in 2D (x–z),
so we can see the on-axis Child–Langmuir physics *and* the 2D field enhancement at
the cathode edges in the same run.

Run with:
```bash
conda activate CBB
python warpx_cathode/cathode_diode.py   # ~1 min, writes openPMD to diags/
python warpx_cathode/plot_cathode.py    # writes the 3 figures to results/
```

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

- **Geometry**: 2D x–z, cathode plane at `z = 0` held at 0 V, anode at `z = d = 4 mm`
  held at `+500 V`. Electrons are emitted only from the finite cathode patch
  `|x| < 6 mm` (the `lower_bound`/`upper_bound` of the flux distribution).
- **Emission**: continuous flux injection (`NFluxPerCell`) at `2 × J_CL`, with a
  small thermal velocity spread set by a 1500 K cathode and a half-Maxwellian
  normal-momentum distribution (`gaussian_flux_momentum_distribution`).
- **Solver**: electrostatic lab frame, **Multigrid** Poisson solver with Dirichlet
  plate potentials (`warpx_potential_lo_z` / `warpx_potential_hi_z`) and Neumann
  transverse walls. (This differs from the FFT/IGF solver used in `warpx_test/` —
  that one is for open-boundary relativistic beams; here we have fixed-potential
  plates and non-relativistic electrons.)
- **Output**: openPMD field snapshots (`phi`, `rho`, `E`, `j`) and electron
  particle data every 200 steps, into `diags/`.

| Parameter | Value |
|-----------|-------|
| Anode bias `V` | 500 V |
| Gap `d` | 4 mm |
| Cathode width `2R` | 12 mm (3× the gap) |
| Cathode temperature | 1500 K |
| Injected current | 2 × J_CL ≈ 3260 A/m² |
| Child–Langmuir J_CL | ≈ 1630 A/m² |
| Grid | 128 × 64 cells (x, z), domain ±12 mm × 4 mm |
| Steps | 2000 (gap-fill ≈ 480 steps) |

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
depression in the beam column, and (3) **bright field enhancement at the cathode
edges** `x = ±6 mm`, where the equipotentials crowd — a genuinely 2D effect with no
counterpart in the planar theory.

### `current_saturation.png` — self-limiting
Transmitted current (integrated across the beam, referenced to the cathode width)
vs. time. Despite injecting **2× J_CL**, the transmitted current saturates near
J_CL (≈ 87% in this run; the small deficit is from transverse beam spreading and
edge losses in the finite 2D geometry). The cathode does **not** pass the current
it is fed — the space charge regulates it.

---

## Notes & possible extensions

- The cathode is intentionally only 3× wider than the gap so the edge effects are
  visible. Make `R_cathode` much larger than `gap_d` to approach the ideal 1D
  planar limit and tighten the J_CL agreement.
- Adam's Region 1 actually *pulses* the grid voltage to chop out a bunch. That can
  be added with a time-dependent potential / `AnalyticFluxDistribution`; this demo
  uses a DC bias to keep the Child–Langmuir validation clean.

## References
- WarpX Pierce-diode example: `reference/WarpX Documentation/usage/examples/pierce_diode/README.md`
- Linac cathode model: `reference/Linac Simulation Documentation/details.md`
- Flux-injection PICMI API: `reference/WarpX Documentation/usage/python.md` (`UniformFluxDistribution`)
