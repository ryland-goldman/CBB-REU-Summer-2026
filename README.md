# Cornell Linac Beam Simulation

Beam-dynamics simulations for the **Cornell High Energy Synchrotron Source (CHESS)** electron
source, built during a Research Experience for Undergraduates (REU) at the Cornell Center for
Bright Beams (CBB) / Cornell Laboratory for Accelerator ScienceS and Education (CLASSE).

The project rebuilds the front end of the **Cornell Linac** — Adam Bartnik's
[LinacSim](https://cesrwww.lepp.cornell.edu/wiki/CESR/LinacSim) thermionic
source → gun → prebuncher → linac chain — from first principles in [WarpX](https://warpx.readthedocs.io),
the massively-parallel particle-in-cell code, using its Python/PICMI interface (`pywarpx`). Each
stage reads the previous stage's openPMD beam as input, so the simulations form a single
self-consistent accelerator chain.

```
cathode  ─►  gun  ─►  prebuncher  ─►  linac_sec1
(SCL diode)  (~146 keV)  (RF bunching)  (~15 MeV captured)
```

## Setup

All simulations run in the **CBB** conda environment. Activate it and install the Python
dependencies:

```bash
conda activate CBB          # Miniforge at ~/miniforge3
# if conda isn't on PATH: source ~/miniforge3/bin/activate

pip install -r requirements.txt
```

`pywarpx` and `openpmd-api` are best installed via conda/mamba (`mamba install -c conda-forge
warpx openpmd-api`); the rest come from `requirements.txt`. See
[`requirements.txt`](requirements.txt) for the pinned versions.

## Run the full chain

The end-to-end driver runs every stage in order with live progress bars and a final-beam summary:

```bash
python pipeline/run_pipeline.py
```

Each stage is also a top-level Python package — `import cathode; cathode.run()` (likewise
`gun.run()`, `prebuncher.run()`, `linac_sec1.run()`) runs that stage alone. Use
`cathode.config(V_anode=60)` etc. to override the module-level parameters before calling `.run()`
(e.g. `linac_sec1.config(I_SOL=0)` for the unfocused linac case). See
[`pipeline/README.md`](pipeline/README.md) for details.

## Components

| Stage | Directory | What it does |
|-------|-----------|--------------|
| **1. Cathode** | [`cathode/`](cathode/README.md) | Thermionic cathode as a finite-extent, space-charge-limited (Child–Langmuir) diode in 2D x–z. The electron source. |
| **2. Gun** | [`gun/`](gun/README.md) | CESR electrostatic gun (~150 kV) in RZ, using the `CESR_gun.gdf` Poisson–Superfish field map. Accelerates the cathode beam to ~146 keV. |
| **3. Prebuncher** | [`prebuncher/`](prebuncher/README.md) | CESR standing-wave RF prebuncher (RZ) that velocity-bunches the gun's exit beam in the downstream drift. |
| **4. Linac Sec 1** | [`linac_sec1/`](linac_sec1/README.md) | SLAC-design 3 m, 2π/3 traveling-wave accelerating section (RZ) with solenoid focusing. The injected 8 kW prebuncher beam (0.83 nC) has diverged to r_max ≈ 26 mm, so only ~32% enters the 12 mm domain; at the original LinacSim point (40 A, 11 MW) capture is ~0.7% of injected (≈5.7 pC) to ⟨KE⟩ ≈ 15.5 MeV (max ~30 MeV). `I_SOL≈1000 A` raises it to ~7% of injected — capture is injection-limited (bore fit), not focusing-limited. |
| **Pipeline** | [`pipeline/`](pipeline/README.md) | Driver + shared `Stage` runner: orchestrates the four stages in order, spawning a fresh Python subprocess per simulation so pywarpx's per-process geometry binding doesn't trip between stages. |

Each directory's `README.md` documents its physics, field maps, and outputs.

## Reference materials

`reference/` holds documentation for the accelerator-physics tools considered for the project
([WarpX](https://warpx.readthedocs.io),
[IMPACT-T](https://github.com/impact-lbl/IMPACT-T)/[IMPACT-Z](https://github.com/impact-lbl/IMPACT-Z),
[GPT](https://www.pulsar.nl/gpt/), [BMAD](https://www.classe.cornell.edu/bmad/),
[G4beamline](https://www.muonsinternal.com/muons3/G4beamline),
[Adam Bartnik's LinacSim](https://cesrwww.lepp.cornell.edu/wiki/CESR/LinacSim),
[LUME-Impact](https://github.com/ChristopherMayes/lume-impact),
[openPMD-beamphysics](https://github.com/ChristopherMayes/openPMD-beamphysics) /
[openPMD-viewer](https://github.com/openPMD/openPMD-viewer),
[easygdf](https://gitlab.com/chris.pierce/easygdf))
plus papers in `reference/Papers/`. See [`CLAUDE.md`](CLAUDE.md) for the full index.

## Notes

- Simulation outputs (`diags/`, `results/`, `*.h5`, `*.gdf`, logs, etc.) are git-ignored — clone
  and re-run to regenerate them.
- Field maps (`CESR_gun.gdf`, `prebuncher_25D.gdf`, plus the linac_sec1 SLAC traveling-wave and
  solenoid/lens maps `SLAC-3mLinac-field1/field2.gdf`, `SOL_0.gdf`, `LENS_0A..0E.gdf`) live in
  [`fieldmaps/`](fieldmaps/) and are read from there by the `build_*_field.py` scripts; paths are
  set near the top of each script.
