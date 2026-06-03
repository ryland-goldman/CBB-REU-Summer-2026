# `ImpactZ` Class — IMPACT-Z Interface

`from impact.z import ImpactZ`

`ImpactZ` wraps an `ImpactZInput` + `ImpactZOutput` pair with the same run/plot/archive interface as the `Impact` class.

## Key Methods

| Method | Description |
|--------|-------------|
| `I.configure()` | Prepare working directory and write input files |
| `I.run()` | Execute the `ImpactZexe` (or MPI) binary |
| `I.plot()` | Generate summary plots of beam statistics |
| `I.archive(file)` | Save all input and output data to an HDF5 file |
| `I.load_archive(file)` | Restore a previous run from an HDF5 archive |

## Key Attributes

| Attribute | Description |
|-----------|-------------|
| `I.input` | `ImpactZInput` object — see [input.md](input.md) |
| `I.output` | `ImpactZOutput` object — see [output.md](output.md) |
| `I.verbose` | Toggle console output |
