# `ImpactZInput` and Beamline Elements

`from impact.z import ImpactZInput`

`ImpactZInput` holds the beamline lattice as a typed list of element objects.

## Beamline Elements

| Class | Description |
|-------|-------------|
| `Drift` | Field-free drift section |
| `Quadrupole` | Magnetic quadrupole |
| `Solenoid` | Solenoid magnet |
| `Dipole` | Dipole bending magnet |
| `Multipole` | General multipole kick |
| `DTL` | Drift-tube linac cavity |
| `CCL` | Coupled-cavity linac |
| `CCDTL` | Coupled-cavity drift-tube linac |
| `ConstantFocusing` | Uniform focusing channel |
| `SolenoidWithRFCavity` | Combined solenoid + RF cavity |
| `SuperconductingCavity` | SRF cavity |
| `TravelingWaveRFCavity` | Traveling-wave RF cavity |
| `UserDefinedRFCavity` | Custom RF fieldmap |

## Control / Diagnostic Elements

| Class | Description |
|-------|-------------|
| `WriteFull` | Write full particle distribution to file |
| `WritePhaseSpaceInfo` | Write 6D phase-space statistics |
| `WriteSliceInfo` | Write slice statistics |
| `CollimateBeam` | Apply aperture collimation |
| `ToggleSpaceCharge` | Enable or disable space charge |
| `RotateBeam` | Rotate beam in phase space |
| `ShiftBeamCentroid` / `ShiftCentroid` | Shift beam centroid coordinates |
| `BeamShift` | General beam coordinate shift |
| `BeamEnergySpread` | Apply an energy spread |
| `EnergyModulation` | Modulate beam energy |
| `BeamKickerByRFNonlinearity` | RF nonlinearity kick |
| `KickBeamUsingMultipole` | Multipole kick |
| `Density3D` / `DensityProfile` | Density diagnostic output |
| `DensityProfileInput` | Input density profile |
| `Projection2D` | 2D phase-space projection output |
| `ScaleMismatchParticle6DCoordinates` | Rescale 6D coordinates |
| `IntegratorTypeSwitch` | Switch integration method mid-lattice |
| `HaltExecution` | Halt simulation at this point |
| `RfcavityStructureWakefield` | RF structure wakefield |

## Internal / Helper Classes

| Class | Description |
|-------|-------------|
| `ElementListProxy` | List-like proxy for the lattice element list |
| `ZElement` | Base class for all IMPACT-Z elements |
