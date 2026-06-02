# 8 Physical Models

The general equations of motion used in the IMPACT-T code are:

```
ṙ = p / (mγ)

ṗ = q(E + p/(mγ) × B)
```

where `γ = 1/sqrt(1 - β²)`, `βᵢ = vᵢ/c` with `i = x, y, z`, `c` is the speed of light, `m` is the rest mass of particle, `q` is the charge of particle. The electric field, **E**, and the magnetic field, **B**, include the contributions from the external focusing and accelerating fields and the space-charge fields of intra-particle Coulomb interactions.

Given electric and magnetic fields, the equations of motion are solved using a second-order leap-frog algorithm: the particles are drifted half time step; the particles are collected and deposited onto a three-dimensional grid; the Poisson equation is solved in the beam frame; the electric and magnetic fields are obtained in the laboratory frame through the Lorentz transformation; the particle momenta are updated using both the space-charge fields and external fields for one time step according to Eq. 2; the particles are drifted another half time step. This procedure is repeated for many time steps until the beam is out of the computational domain of beam line elements.

---

| Section | File |
|---------|------|
| 8.1 Internal Coordinates of Particle | [8.1_internal_coordinates.md](8.1_internal_coordinates.md) |
| 8.2 Particle Emission from Cathode | [8.2_particle_emission.md](8.2_particle_emission.md) |
| 8.3 Space-Charge Effects | [8.3_space_charge.md](8.3_space_charge.md) |
| 8.4 Short Range Longitudinal and Transverse Wakefields | [8.4_short_range_wakefields.md](8.4_short_range_wakefields.md) |
| 8.5 Longitudinal CSR Wakefield | [8.5_longitudinal_csr_wakefield.md](8.5_longitudinal_csr_wakefield.md) |
| 8.6 RF Fields in Standing Wave Structures | [8.6_rf_fields_standing_wave.md](8.6_rf_fields_standing_wave.md) |
| 8.7 Traveling Wave Structures | [8.7_traveling_wave_structures.md](8.7_traveling_wave_structures.md) |
| 8.8 Solenoid | [8.8_solenoid.md](8.8_solenoid.md) |
| 8.9 Bending Magnet | [8.9_bending_magnet.md](8.9_bending_magnet.md) |
| 8.10 Quadrupole | [8.10_quadrupole.md](8.10_quadrupole.md) |
| 8.11 Multipole | [8.11_multipole.md](8.11_multipole.md) |
