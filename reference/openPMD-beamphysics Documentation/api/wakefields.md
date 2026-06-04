# Wakefields

The `beamphysics.wakefields` module. (Not used by this repo; listed for completeness.
Full signatures live in the package source / docs site — these were build-time
mkdocstrings directives upstream.) Applied to particles via
`ParticleGroup.apply_wakefield(...)` and visualized with `ParticleGroup.wakefield_plot`.

## Resistive-wall wakefield classes
- `ResistiveWallWakefield`
- `ResistiveWallPseudomode`

## Base classes
- `WakefieldBase`
- `PseudomodeWakefield`
- `ImpedanceWakefield`
- `TabularWakefield`
- `Pseudomode`

## Low-level functions
- `longitudinal_impedance_round`
- `longitudinal_impedance_flat`
- `wakefield_from_impedance`
- `wakefield_from_impedance_fft`
- `ac_conductivity`
- `surface_impedance`
- `characteristic_length`
