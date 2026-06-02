# 7 Output Data

---

## fort.18

| Column | Description |
|--------|-------------|
| 1st | time (secs) |
| 2nd | distance (m) |
| 3rd | gamma |
| 4th | kinetic energy (MeV) |
| 5th | beta |
| 6th | Rmax (m) — R is measured from the axis of pipe |
| 7th | rms energy deviation normalized by MC² |

---

## fort.24, fort.25 — X and Y RMS size information

| Column | Description |
|--------|-------------|
| 1st | time (secs) |
| 2nd | z distance (m) |
| 3rd | centroid location (m) |
| 4th | RMS size (m) |
| 5th | Centroid momentum normalized by MC |
| 6th | RMS momentum normalized by MC |
| 7th | Twiss parameter |
| 8th | normalized RMS emittance (m-rad) |

---

## fort.26 — Z RMS size information

| Column | Description |
|--------|-------------|
| 1st | time (secs) |
| 2nd | centroid location (m) |
| 3rd | RMS size (m) |
| 4th | Centroid momentum normalized by MC |
| 5th | RMS momentum normalized by MC |
| 6th | Twiss parameter |
| 7th | normalized RMS emittance (m-rad) |

---

## fort.27 — Maximum amplitude information

| Column | Description |
|--------|-------------|
| 1st | time (secs) |
| 2nd | z distance (m) |
| 3rd | Max. X (m) |
| 4th | Max. Px (MC) |
| 5th | Max. Y (m) |
| 6th | Max. Py (MC) |
| 7th | Max. Z (m) (with respect to centroid) |
| 8th | Max. Pz (MC) |

---

## fort.28 — Load balance and loss diagnostic

| Column | Description |
|--------|-------------|
| 1st | time (secs) |
| 2nd | z distance (m) |
| 3rd | min # of particles on a PE |
| 4th | max # of particles on a PE |
| 5th | total # of particles in the bunch |

---

## fort.29 — Cubic root of 3rd moments of the beam distribution

| Column | Description |
|--------|-------------|
| 1st | time (secs) |
| 2nd | z distance (m) |
| 3rd | X (m) |
| 4th | Px (MC) |
| 5th | Y (m) |
| 6th | Py (MC) |
| 7th | Z (m) |
| 8th | Pz (MC) |

---

## fort.30 — Square root of square root of 4th moments of the beam distribution

| Column | Description |
|--------|-------------|
| 1st | time (secs) |
| 2nd | z distance (m) |
| 3rd | X (m) |
| 4th | Px (MC) |
| 5th | Y (m) |
| 6th | Py (MC) |
| 7th | Z (m) |
| 8th | Pz (MC) |

---

## fort.31 — 4×4 sigma matrix

| Column | Description |
|--------|-------------|
| 1st | time (secs) |
| 2nd | z distance (m) |
| 3rd | length scale (m) |
| 4th | `<x**2>` |
| 5th | `<xpx>` |
| 6th | `<xy>` |
| 7th | `<xpy>` |
| 8th | `<px**2>` |
| 9th | `<pxy>` |
| 10th | `<pxpy>` |
| 11th | `<y**2>` |
| 12th | `<ypy>` |
| 13th | `<py**2>` |

---

## fort.32 — 6×6 sigma matrix

| Column | Description |
|--------|-------------|
| 1st | time (secs) |
| 2nd | z distance (m) |
| 3rd | length scale (m) |
| 4th | `<x**2>` |
| 5th | `<xpx>` |
| 6th | `<xy>` |
| 7th | `<xpy>` |
| 8th | `<xz>` |
| 9th | `<xpz>` |
| 10th | `<px**2>` |
| 11th | `<pxy>` |
| 12th | `<pxpy>` |
| 13th | `<pxz>` |
| 14th | `<pxpz>` |
| ... | ... |
| 23rd | `<zpz>` |
| 24th | `<pz**2>` |

---

## fort.34, fort.35 — X and Y output information in dipole reference coordinate system (inside dipole ONLY)

| Column | Description |
|--------|-------------|
| 1st | time (secs) |
| 2nd | z (m) |
| 3rd | x or y (m) |
| 4th | xrms or yrms (m) |
| 5th | Centroid momentum normalized by MC |
| 6th | RMS momentum normalized by MC |
| 7th | correlation parameter |
| 8th | normalized RMS emittance (m-rad) |

---

## fort.36 — Z output information in dipole reference coordinate system (inside dipole ONLY)

| Column | Description |
|--------|-------------|
| 1st | time (secs) |
| 2nd | z (m) |
| 3rd | rms (m) |
| 4th | Centroid momentum normalized by MC |
| 5th | RMS momentum normalized by MC |
| 6th | correlation parameter |
| 7th | normalized RMS emittance (m-rad) |

---

## fort.37 — Maximum amplitude information in dipole reference coordinate system (inside dipole ONLY)

| Column | Description |
|--------|-------------|
| 1st | time (secs) |
| 2nd | z distance (m) |
| 3rd | Max. X (m) |
| 4th | Max. Px (MC) |
| 5th | Max. Y (m) |
| 6th | Max. Py (MC) |
| 7th | Max. Z (m) (with respect to centroid) |
| 8th | Max. Pz (MC) |

---

## fort.38 — Reference particle information in dipole reference coordinate system (inside dipole ONLY)

| Column | Description |
|--------|-------------|
| 1st | time (secs) |
| 2nd | x distance (m) |
| 3rd | Px/MC |
| 4th | y (m) |
| 5th | Py/MC |
| 6th | z (m) |
| 7th | Pz/MC |

---

## fort.40 — Initial particle distribution at t = 0

Particle coordinates are x (m), Px/MC, y (m), Py/MC, z (m), Pz/MC.

---

## fort.50 — Final particle distribution

Final particle distribution projected to the centroid location of the bunch.

---

## fort.60 — Slice information of the initial distribution

| Column | Description |
|--------|-------------|
| 1st | bunch length (m) |
| 2nd | number of macroparticles per cell |
| 3rd | current profile |
| 4th | x slice emittance (m-rad) |
| 5th | y slice emittance (m-rad) |
| 6th | energy spread per cell without taking out correlation (eV) |
| 7th | uncorrelated energy spread per cell (eV) |

---

## fort.70 — Slice information of the final distribution

Output file for slice information of the final distribution (same format as fort.60).
