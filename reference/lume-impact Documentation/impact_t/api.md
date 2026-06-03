# `Impact` Class — IMPACT-T Interface

`from impact import Impact`

The `Impact` class is the primary entry point for running IMPACT-T simulations. Loading an `ImpactT.in` file calls `configure()` automatically.

## Key Attributes

| Attribute | Description |
|-----------|-------------|
| `I.header` | Dict-like access to `ImpactT.in` header parameters (e.g., `Np`, `Nx`, `Ny`, `Nz`) |
| `I.lattice` | List of beamline element objects |
| `I.particles` | Dict of `ParticleGroup` objects (openPMD-beamphysics) keyed by name (e.g., `"final_particles"`) |
| `I.output` | Parsed output data (statistics, particle distributions) |
| `I.verbose` | Toggle console output |
| `I.workdir` | Temporary working directory for the run |

## Key Methods

| Method | Description |
|--------|-------------|
| `I.configure()` | Prepare working directory and write input files |
| `I.run()` | Execute the `ImpactTexe` (or MPI) binary |
| `I.plot()` | Generate summary plots of beam statistics |
| `I.archive(file)` | Save all input and output data to an HDF5 file |
| `I.load_archive(file)` | Restore a previous run from an HDF5 archive |

## BMAD Interface

LUME-Impact can convert BMAD lattice elements to IMPACT-T input. See `docs/examples/bmad_interface.ipynb` and `docs/examples/bmad_to_impact/` for examples.
