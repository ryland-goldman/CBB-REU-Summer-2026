# `ImpactZOutput` and Output Data Classes

`from impact.z import ImpactZOutput`

## Top-Level Output Classes

| Class | Description |
|-------|-------------|
| `ImpactZOutput` | Top-level output container |
| `RunInfo` | Run metadata (wall time, MPI rank count, etc.) |
| `OutputStats` | Beam statistics vs. z position (emittance, sizes, centroids) |
| `ImpactZSlices` | Slice-by-slice longitudinal diagnostics |

## Raw Fortran Output File Classes

These map directly to IMPACT-Z's `fort.*` output files. Parsed automatically by `ImpactZOutput`.

| Class | Fort file | Contents |
|-------|-----------|----------|
| `ReferenceParticles` | fort.18 | Reference particle coordinates |
| `RmsX` | fort.24 | RMS x phase-space moments |
| `RmsY` | fort.25 | RMS y phase-space moments |
| `RmsZ` | fort.26 | RMS z phase-space moments |
| `MaxAmplitudeStandard` | fort.27 | Maximum amplitude (standard) |
| `MaxAmplitudeExtended` | fort.28 | Maximum amplitude (extended) |
| `LoadBalanceLossDiagnostic` | fort.29 | Load balance and particle loss |
| `BeamDistribution3rdStandard` | fort.30 | 3rd-order moments (standard) |
| `BeamDistribution3rdExtended` | fort.32 | 3rd-order moments (extended) |
| `BeamDistribution4th` | fort.34 | 4th-order moments |
| `ParticlesAtChargedState` | fort.40+ | Per-charge-state particle data |
| `FortranOutputFileData` | (base) | Base class for all fort.* parsing |

For column-by-column descriptions of the fort.* files, see the [IMPACT-Z Documentation](../../Impact-Z%20Documentation/README.md).
