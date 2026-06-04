# Cornell Linac Beam Simulation

Beam-dynamics simulations for the **Cornell High Energy Synchrotron Source (CHESS)** electron
source, built during a Research Experience for Undergraduates (REU) at the Cornell Center for
Bright Beams (CBB) / Cornell Laboratory for Accelerator ScienceS and Education (CLASSE).

The project rebuilds the front end of the **Cornell Linac** — Adam Bartnik's thermionic
source → gun → prebuncher chain — from first principles in [WarpX](https://warpx.readthedocs.io),
the massively-parallel particle-in-cell code, using its Python/PICMI interface (`pywarpx`). Each
stage reads the previous stage's openPMD beam as input, so the simulations form a single
self-consistent accelerator chain.

```
cathode  ─►  gun  ─►  prebuncher
(SCL diode)  (~148 keV)  (RF bunching)
```

## Setup

All simulations run in the **CBB** conda environment:

```bash
conda activate CBB          # Miniforge at ~/miniforge3
# if conda isn't on PATH: source ~/miniforge3/bin/activate
```

## Run the full chain

The end-to-end driver runs every stage in order with live progress bars and a final-beam summary:

```bash
python pipeline/run_pipeline.py
```

Stage toggles and the prebuncher operating point are configured in the `CONFIG` block at the top
of `run_pipeline.py`. See [`pipeline/README.md`](pipeline/README.md) for details.

## Components

| Stage | Directory | What it does |
|-------|-----------|--------------|
| **1. Cathode** | [`warpx_cathode/`](warpx_cathode/README.md) | Thermionic cathode as a finite-extent, space-charge-limited (Child–Langmuir) diode in 2D x–z. The electron source. |
| **2. Gun** | [`warpx_gun/`](warpx_gun/README.md) | CESR electrostatic gun (~150 kV) in RZ, using the `CESR_gun.gdf` Poisson–Superfish field map. Accelerates the cathode beam to ~148 keV. |
| **3. Prebuncher** | [`warpx_prebuncher/`](warpx_prebuncher/README.md) | CESR standing-wave RF prebuncher (RZ) that velocity-bunches the gun's exit beam in the downstream drift. |
| **Pipeline** | [`pipeline/`](pipeline/README.md) | Chains stages 1–3 as subprocesses; each reads the prior stage's openPMD output. |
| **Tutorials** | [`warpx_test/`](warpx_test/README.md) | Introductory WarpX/PICMI demos (single positron, Gaussian bunch, space charge) — the learning warm-up, not part of the chain. |

Each directory's `README.md` documents its physics, field maps, and outputs.

## Reference materials

`reference/` holds documentation for the accelerator-physics tools considered for the project
(WarpX, IMPACT-T/Z, GPT, BMAD, G4beamline, the Linac GUI, LUME-Impact, openPMD tools, easygdf)
plus papers in `reference/Papers/`. See [`CLAUDE.md`](CLAUDE.md) for the full index.

## Notes

- Simulation outputs (`diags/`, `results/`, `*.h5`, `*.gdf`, logs, etc.) are git-ignored — clone
  and re-run to regenerate them.
- Field maps (`CESR_gun.gdf`, `prebuncher_25D.gdf`) live in [`fieldmaps/`](fieldmaps/) and are
  read from there by the `build_*_field.py` scripts; paths are set near the top of each script.
