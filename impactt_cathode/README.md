# SPEBS Gun in IMPACT-T (Bunched Cathode Emission)

An IMPACT-T model of the **same** Cornell SPEBS Pierce-gun cathode that
`warpx_cathode/spebs_gun.py` modeled in WarpX — V = 50 kV, cathode→anode gap
d = 23.1 mm, 4 mm-radius thermionic cathode — but in IMPACT-T's native regime:
a finite electron **bunch** emitted from the cathode and tracked across the gap with
3D space charge **and the cathode image charge**. Driven entirely from Python via
[lume-impact](../reference/lume-impact%20Documentation/README.md).

Run with:
```bash
conda activate CBB
python impactt_cathode/spebs_gun_impactt.py img     # image charge ON  (~1-2 min)
python impactt_cathode/spebs_gun_impactt.py noimg   # image charge OFF
python impactt_cathode/plot_impactt.py              # writes 3 figures to results/
```

---

## Why IMPACT-T as well as WarpX?

The WarpX run solved this gun as a **steady-state space-charge-limited diode** and
validated the Child–Langmuir *current limit* — how much current the gap can carry.
IMPACT-T is a **photoinjector / bunched-beam** tracking code: it generates a finite
particle bunch behind the cathode and emits it over a short pulse (`Nemission` steps
within `Temission`), tracking it through the gap. So instead of a current limit, the
IMPACT-T deliverables are **single-bunch dynamics**:

- energy gain across the 50 kV gap,
- transverse beam-size evolution,
- **thermal-emittance** growth from the 1500 K cathode,
- and the **cathode image charge** — the attractive pull of the conductor on the
  emitted charge near `z = 0`, which IMPACT-T includes when `Flagimg = 1` (the
  cathode is fixed at `z = 0`). This effect has no analogue in the steady-state
  WarpX diode picture.

The two models are complementary: WarpX answers *"how much current?"*, IMPACT-T
answers *"what does one bunch look like coming off the cathode?"*.

---

## What the driver does (`spebs_gun_impactt.py`)

- **Accelerating field**: a uniform DC field `E0 = V/d = 2.165 MV/m` is built as an
  on-axis profile with `pmd_beamphysics.FieldMesh.from_onaxis` (flat across the gap,
  with smooth cosine ramps placed *inside* the cathode at `z<0` and *past* the anode,
  so the tracked beam only ever sees the flat region). It is attached as a standard
  `solrf` (type 105) gun element with **RF frequency 0** (i.e. static). `Ez` is
  negative so the force `−eE` accelerates electrons toward `+z`.
- **Emission**: cathode model on (`Nemission = 100`, `Temission = 120 ps`),
  Gaussian transverse spot on a 4 mm cathode, thermal momentum spread
  `σ_p = √(k_B·1500 K / m_e c²)`, launched essentially at rest (≈ 0.4 eV).
- **Charge**: SPEBS runs ~1 A continuously; we emit a **0.1 nC** slice
  (`Bcurr/Bfreq = 1 A / 10 GHz`) so the bunch carries a representative current.
- **Image charge**: `Flagimg = 1` (`img`) vs `0` (`noimg`), toggled from the CLI.
- **Output**: lume-impact parses `fort.18`/`fort.24`–`26` into `I.output['stats']`
  and the final phase space into `I.particles['final_particles']`; the driver dumps
  `results/stats_<tag>.json` plus the exact `ImpactT.in` + `rfdata1001` it ran.

### One non-obvious IMPACT-T detail

The uniform field is written in IMPACT-T's **discrete derivative** field-map format
(`Ez, Ez', Ez'', Ez'''`), which avoids the Gibbs ringing a Fourier representation
would give for a flat-top. IMPACT-T only reads that format when the `solrf` **file ID
is > 1000** — otherwise it interprets the file as Fourier coefficients and the field
collapses to ≈ 0 (the beam then never accelerates). The driver therefore uses
`file_id = 1001` → `rfdata1001`.

---

## Figures (`plot_impactt.py` → `results/`)

- **`gun_dynamics.png`** — mean kinetic energy, RMS beam size, and normalized
  emittance vs z. The energy rises **linearly** to ≈ 50 keV, overlaid with the
  uniform-field prediction `eE0·z`: a direct check that the DC gun field is correct.
- **`final_phasespace.png`** — `x–px` and `z–pz` of the final bunch at the anode.
- **`image_charge.png`** — normalized emittance (zoomed to the first few mm) and
  energy gain for image charge **ON vs OFF**, isolating the cathode image-charge
  effect near `z = 0`.

---

## Validation

| Check | Expectation | Result |
|-------|-------------|--------|
| Exit kinetic energy | ≈ 50 keV (= eV across the gap) | **49.3 keV** ✓ |
| Energy vs z | linear (uniform field) | KE(d/2) ≈ 25 keV, overlays `eE0·z` ✓ |
| Initial energy | thermal launch, ≪ 1 eV | ≈ 0.4 eV ✓ |
| Transmission | most of the bunch reaches the anode | **75%** (rest lost to image charge / space charge in the first emission step at the cathode) |
| Emittance growth | smooth rise from thermal floor | ~1.0 µm → ~4.6 µm, monotonic ✓ |
| Image charge | resolvable ON-vs-OFF difference near cathode | ε_n,x = **4.6 µm (ON)** vs **5.3 µm (OFF)** ✓ |

The image charge *lowers* the emittance growth in the first few mm (it pulls the
expanding bunch back toward the axis near the cathode) — visible as the blue curve
sitting below the orange one in `image_charge.png`.

> **Field-map note:** the uniform field is kept flat out to 40 mm — well past the
> 23.1 mm anode — because the bunch is long (early-emitted electrons accelerate
> longer, so the head runs ~5 mm ahead of the centroid). If the field ramps down
> where the bunch head still is, the ramp's `dEz/dz` gives off-axis radial kicks that
> artificially inflate the *normalized* emittance near the anode. Extending the flat
> region past the tracked beam removes that artifact.

---

## Files

| File | Purpose |
|------|---------|
| `spebs_gun_impactt.py` | lume-impact driver — builds the DC field, runs IMPACT-T, saves stats |
| `plot_impactt.py` | renders the three figures from `results/stats_*.json` |
| `results/deck_<tag>/ImpactT.in`, `rfdata1001` | the exact deck + field map that ran (inspectable) |
| `results/stats_*.json` | parsed beam statistics + final phase space |

Sibling demo: [`warpx_cathode/`](../warpx_cathode/README.md) — the same gun as a
self-consistent Child–Langmuir diode in WarpX.
