# Cornell Linac Beam Simulation

Beam-dynamics simulations for the **Cornell High Energy Synchrotron Source (CHESS)** electron
source, built during a Research Experience for Undergraduates (REU) at the Cornell Center for
Bright Beams (CBB) / Cornell Laboratory for Accelerator ScienceS and Education (CLASSE).

The project rebuilds the front end of the **Cornell Linac** — Adam Bartnik's
[LinacSim](https://cesrwww.lepp.cornell.edu/wiki/CESR/LinacSim) thermionic
source → gun → injector → linac chain — from first principles. The first four stages are built
in [WarpX](https://warpx.readthedocs.io), the massively-parallel particle-in-cell code, using its
Python/PICMI interface (`pywarpx`); the final stage (linac sections 2–8) is an
[Impact-T](https://github.com/impact-lbl/IMPACT-T) run driven through
[lume-impact](https://github.com/ChristopherMayes/lume-impact). Each stage reads the previous
stage's openPMD beam as input, so the simulations form a single self-consistent accelerator chain.

```
cathode  ─►  gun  ─►  injector  ─►  linac_sec1  ─►  linac_rest
(SCL diode)  (~146 keV)  (2 prebunchers + 3 solenoids)  (~25 MeV captured)  (sections 2–8, ~308 MeV, Impact-T)
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
`gun.run()`, `injector.run()`, `linac_sec1.run()`, `linac_rest.run()`) runs that stage alone. Use
`cathode.config(V_anode=60)` etc. to override the module-level parameters before calling `.run()`
(e.g. `injector.config(I_SOL0=0)` to disable Sol 0). See
[`pipeline/README.md`](pipeline/README.md) for details. After the five stages,
`pipeline.plot_chain()` writes cross-stage figures to the repo-root `results/`.

## Components

| Stage | Directory | What it does |
|-------|-----------|--------------|
| **1. Cathode** | [`cathode/`](cathode/README.md) | Thermionic cathode as a finite-extent, space-charge-limited (Child–Langmuir) diode in 2D x–z. The electron source. |
| **2. Gun** | [`gun/`](gun/README.md) | CESR electrostatic gun (~150 kV) in RZ, using the `CESR_gun.gdf` Poisson–Superfish field map. Accelerates the cathode beam to ~146 keV. |
| **3. Injector** | [`injector/`](injector/README.md) | The full LinacSim injector subsection in one RZ space-charge run (RZ): Lens 0A → Prebuncher 1 (8 kW) → Prebuncher 2 (10 kW, reversed) → Sol 0 / Lens 0E, then the 9.547 mm collimator. Two-cavity velocity bunching + solenoid focusing; hands a focused, collimated beam to the linac at z ≈ 2.03 m. (Replaced the earlier single-cavity `prebuncher/` stage.) |
| **4. Linac Sec 1** | [`linac_sec1/`](linac_sec1/README.md) | SLAC-design 3 m, 2π/3 traveling-wave accelerating section (RZ). Reads the injector's focused beam at the z ≈ 2.03 m handoff and applies the multi-plane 9.547 mm iris scrape at injection; no in-stage solenoid (focusing is upstream now). At the faithful 11 MW point capture is ~7% of true injected to ⟨KE⟩ ≈ 25 MeV — a conservative (γ²) lower bound, tune-sensitive to the upstream lens currents. |
| **5. Linac Rest** | [`linac_rest/`](linac_rest/README.md) | Cornell linac **sections 2–8** (CEA 2/3/4/5 + CU 3/4/5) — seven S-band traveling-wave sections chained into one **Impact-T** deck (via lume-impact, not WarpX). Reads `linac_sec1`'s captured ~25 MeV exit beam, keeps the β > 0.999 core (a model-validity cut), and accelerates on-crest to ≈308 MeV at 11 MW (307.97 survivors through the real bore). No field maps (generic constant-gradient TW, reusing the lume-impact `rfdata4–7` shape); space charge off; quads OFF for the headline. The e+ compressor (CU 2) is out of scope. |
| **Pipeline** | [`pipeline/`](pipeline/README.md) | Driver + shared runner: orchestrates the five stages in order. The four WarpX stages each spawn a fresh Python subprocess per simulation (so pywarpx's per-process geometry binding doesn't trip between stages); the Impact-T `linac_rest` stage runs in-process (external `ImpactTexe`, no pywarpx binding). Then `plot_chain` writes the cross-stage figures to the repo-root `results/`. |
| **Ch-10 comparison** | [`pipeline/plot_chapter10.py`](pipeline/plot_chapter10.py) | Cheng-Yang Tan dissertation Ch.10 (PARMELA) comparison figures + table; cross-stage post-processing (no sim) reading the four location snapshots (gun exit / injector pre-Preb2 entrance / injector handoff / linac_sec1 exit) → Tan Fig 10.2–10.5 + a Table-10.2 comparison. `_repo_default` = **qualitative** comparison (the 8/10 kW default is NOT Tan cond (i)); `_cond_i` = the **true Tan operating point** (50/150 kV, numeric comparison), which needs an injector + linac_sec1 re-run into `*/diags/cond_i` (see [`CLAUDE.md`](CLAUDE.md) Commands). Run separately, not in the default chain. |

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
- Field maps (`CESR_gun.gdf`; the injector's `prebuncher_25D.gdf` cavity + `SOL_0.gdf` /
  `LENS_0A..0E.gdf` solenoid lenses; the linac's `SLAC-3mLinac-field1/field2.gdf` traveling-wave
  pair) live in [`fieldmaps/`](fieldmaps/) and are read from there by the `build_*_field.py`
  scripts; paths are set near the top of each script. The `linac_rest` stage (sections 2–8) has
  **no** field maps — it reuses the lume-impact S-band TW field *shape* `rfdata4–7`, vendored
  into [`linac_rest/rfdata/`](linac_rest/rfdata/) (committed).
- The repo-root [`results/`](results/) holds the cross-stage figures from `pipeline.plot_chain()`
  (also git-ignored; regenerate by re-running). Commit result PNGs with `git add -f`.
