# Example Templates

Example input files and Jupyter notebooks are bundled in the lume-impact source tree under `docs/examples/`.

## IMPACT-T Example Notebooks

| Notebook | Description |
|----------|-------------|
| `basic_impact_examples.ipynb` | Basic IMPACT-T run and output inspection |
| `autophase_example.ipynb` | Automatic RF cavity phasing |
| `distgen_example.ipynb` | Beam distribution generation with distgen |
| `input_parsing_example.ipynb` | Parsing and modifying ImpactT.in files |
| `output_parsing_example.ipynb` | Reading and plotting fort.* output files |
| `plotting_example.ipynb` | Phase-space and statistics plots |
| `parallel_run_example.ipynb` | MPI parallel run example |
| `functional_impact_run.ipynb` | Functional-style run workflow |
| `compare_covariance.ipynb` | Covariance matrix comparison |
| `point_to_point_spacecharge.ipynb` | Point-to-point space charge |
| `fieldmap_reconstruction.ipynb` | RF fieldmap reconstruction |
| `movie_example.ipynb` | Beam evolution animation |
| `bmad_interface.ipynb` | BMAD lattice → IMPACT-T conversion |

### IMPACT-T Element Notebooks (`examples/elements/`)

| Notebook | Element |
|----------|---------|
| `drift.ipynb` | Drift |
| `quadrupole.ipynb` | Quadrupole |
| `solenoid.ipynb` | Solenoid |
| `dipole.ipynb` | Dipole |
| `corrector_coil.ipynb` | Corrector coil |
| `tesla_9cell_cavity.ipynb` | Tesla 9-cell SRF cavity |
| `traveling_wave_cavity.ipynb` | Traveling-wave cavity |
| `apex_gun.ipynb` | APEX DC gun |
| `awa_flatbeam.ipynb` | AWA flat-beam transformer |
| `wakefield.ipynb` | Wakefield element |
| `3d_field.ipynb` | 3D external fieldmap |

### IMPACT-T Input Templates (`examples/templates/`)

| Template | Description |
|----------|-------------|
| `apex_gun/` | APEX DC gun + solenoid photoinjector |
| `lcls_injector/` | LCLS S-band RF gun injector |
| `awa_flatbeam/` | AWA flat-beam transformer (7 solenoid fieldmaps) |
| `fermilab_fast/` | Fermilab FAST injector |
| `chicane/` | Magnetic chicane with CSR fieldmaps |
| `solenoid/` | Solenoid with 1D and 2D fieldmaps |
| `tesla_9cell_cavity/` | Tesla 9-cell SRF cavity (1D and 2D) |
| `traveling_wave_cavity/` | Traveling-wave cavity |
| `dipole/` | Single dipole element |
| `quadrupole/` | Single quadrupole element |
| `drift/` | Simple drift |
| `corrector_coil/` | Corrector coil |
| `rf_gun/` | RF gun |
| `wakefield/` | Wakefield element (SDDS wake file) |
| `3dfield/` | 3D external fieldmap |

---

## IMPACT-Z Example Notebooks

| Notebook | Description |
|----------|-------------|
| `z/basic/basic_impact_z.ipynb` | Basic ImpactZ run |
| `z/example1/example1.ipynb` | Example 1 (FODO-like lattice) |
| `z/example2/example2.ipynb` | Example 2 (RF cavity) |
| `z/example3/example3.ipynb` | Example 3 (space charge) |

### IMPACT-Z Element Notebooks (`examples/z/elements/`)

| Notebook | Description |
|----------|-------------|
| `drift.ipynb` / `drift-bmad.ipynb` | Drift (standalone and vs. BMAD) |
| `dipole.ipynb` / `dipole-bmad.ipynb` | Dipole |
| `quadrupole-bmad.ipynb` | Quadrupole vs. BMAD |
| `sextupole-bmad.ipynb` | Sextupole vs. BMAD |
| `octupole-bmad.ipynb` | Octupole vs. BMAD |
| `decapole-bmad.ipynb` | Decapole vs. BMAD |
| `solenoid-bmad.ipynb` | Solenoid vs. BMAD |
| `lcavity-bmad.ipynb` | Linac cavity vs. BMAD |
| `traveling_wave_rf_cavity.ipynb` | Traveling-wave RF cavity |
| `wiggler-bmad.ipynb` | Wiggler vs. BMAD |
| `collimator-bmad.ipynb` | Collimator vs. BMAD |
| `spacecharge-benchmark.ipynb` | Space charge benchmark |
| `spacecharge-drift-bmad.ipynb` | Space charge in drift vs. BMAD |
| `csr-bench-bmad.ipynb` | CSR benchmark vs. BMAD |
| `csr-zeuthen.ipynb` | Zeuthen CSR benchmark |
| `optics-matching-bmad.ipynb` | Optics matching via BMAD |
| `compare-particle-stats.ipynb` | Particle statistics comparison |
