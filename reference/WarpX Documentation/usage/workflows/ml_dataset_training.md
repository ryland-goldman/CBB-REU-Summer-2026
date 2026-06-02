<a id="ml-dataset-training"></a>

# Training a Surrogate Model from WarpX Data

Suppose we have a WarpX simulation that we wish to replace with a neural network surrogate model.
For example, a simulation determined by the following input script

### Python Input for Training Simulation

> ```python
> #!/usr/bin/env python3
> import math

> import numpy as np

> from pywarpx import picmi

> # Physical constants
> c = picmi.constants.c
> q_e = picmi.constants.q_e
> m_e = picmi.constants.m_e
> m_p = picmi.constants.m_p
> ep0 = picmi.constants.ep0

> # Number of cells
> dim = "3"
> nx = ny = 128
> nz = 35328  # 17664 #8832
> if dim == "rz":
>     nr = nx // 2

> # Computational domain
> rmin = 0.0
> rmax = 128e-6
> zmin = -180e-6
> zmax = 0.0

> # Number of processes for static load balancing
> # Check with your submit script
> num_procs = [1, 1, 64 * 4]
> if dim == "rz":
>     num_procs = [1, 64]

> # Number of time steps
> gamma_boost = 60.0
> beta_boost = np.sqrt(1.0 - gamma_boost**-2)

> # Create grid
> if dim == "rz":
>     grid = picmi.CylindricalGrid(
>         number_of_cells=[nr, nz],
>         guard_cells=[32, 32],
>         n_azimuthal_modes=2,
>         lower_bound=[rmin, zmin],
>         upper_bound=[rmax, zmax],
>         lower_boundary_conditions=["none", "damped"],
>         upper_boundary_conditions=["none", "damped"],
>         lower_boundary_conditions_particles=["absorbing", "absorbing"],
>         upper_boundary_conditions_particles=["absorbing", "absorbing"],
>         moving_window_velocity=[0.0, c],
>         warpx_max_grid_size=256,
>         warpx_blocking_factor=64,
>     )
> else:
>     grid = picmi.Cartesian3DGrid(
>         number_of_cells=[nx, ny, nz],
>         guard_cells=[11, 11, 12],
>         lower_bound=[-rmax, -rmax, zmin],
>         upper_bound=[rmax, rmax, zmax],
>         lower_boundary_conditions=["periodic", "periodic", "damped"],
>         upper_boundary_conditions=["periodic", "periodic", "damped"],
>         lower_boundary_conditions_particles=["periodic", "periodic", "absorbing"],
>         upper_boundary_conditions_particles=["periodic", "periodic", "absorbing"],
>         moving_window_velocity=[0.0, 0.0, c],
>         warpx_max_grid_size=256,
>         warpx_blocking_factor=32,
>     )


> # plasma region
> plasma_rlim = 100.0e-6
> N_stage = 15
> L_plasma_bulk = 0.28
> L_ramp = 1.0e-9
> L_ramp_up = L_ramp
> L_ramp_down = L_ramp
> L_stage = L_plasma_bulk + 2 * L_ramp

> # focusing
> # lens external fields
> beam_gamma1 = 15095
> lens_focal_length = 0.015
> lens_width = 0.003

> stage_spacing = L_plasma_bulk + 2 * lens_focal_length


> def get_species_of_accelerator_stage(
>     stage_idx,
>     stage_zmin,
>     stage_zmax,
>     stage_xmin=-plasma_rlim,
>     stage_xmax=plasma_rlim,
>     stage_ymin=-plasma_rlim,
>     stage_ymax=plasma_rlim,
>     Lplus=L_ramp_up,
>     Lp=L_plasma_bulk,
>     Lminus=L_ramp_down,
> ):
>     # Parabolic density profile
>     n0 = 1.7e23
>     Rc = 40.0e-6
>     Lstage = Lplus + Lp + Lminus
>     if not np.isclose(stage_zmax - stage_zmin, Lstage):
>         print("Warning: zmax disagrees with stage length")
>     parabolic_distribution = picmi.AnalyticDistribution(
>         density_expression=f"n0*(1.+4.*(x**2+y**2)/(kp**2*Rc**4))*(0.5*(1.-cos(pi*(z-{stage_zmin})/Lplus)))*((z-{stage_zmin})<Lplus)"
>         + f"+n0*(1.+4.*(x**2+y**2)/(kp**2*Rc**4))*((z-{stage_zmin})>=Lplus)*((z-{stage_zmin})<(Lplus+Lp))"
>         + f"+n0*(1.+4.*(x**2+y**2)/(kp**2*Rc**4))*(0.5*(1.+cos(pi*((z-{stage_zmin})-Lplus-Lp)/Lminus)))*((z-{stage_zmin})>=(Lplus+Lp))*((z-{stage_zmin})<(Lplus+Lp+Lminus))",
>         n0=n0,
>         kp=q_e / c * math.sqrt(n0 / (m_e * ep0)),
>         Rc=Rc,
>         Lplus=Lplus,
>         Lp=Lp,
>         Lminus=Lminus,
>         lower_bound=[stage_xmin, stage_ymin, stage_zmin],
>         upper_bound=[stage_xmax, stage_ymax, stage_zmax],
>         fill_in=True,
>     )

>     electrons = picmi.Species(
>         particle_type="electron",
>         name=f"electrons{stage_idx}",
>         initial_distribution=parabolic_distribution,
>     )

>     ions = picmi.Species(
>         particle_type="proton",
>         name=f"ions{stage_idx}",
>         initial_distribution=parabolic_distribution,
>     )

>     return electrons, ions


> species_list = []
> for i_stage in range(1):
>     # Add plasma
>     zmin_stage = i_stage * stage_spacing
>     zmax_stage = zmin_stage + L_stage
>     electrons, ions = get_species_of_accelerator_stage(
>         i_stage + 1, zmin_stage, zmax_stage
>     )
>     species_list.append(electrons)
>     species_list.append(ions)

> # add beam to species_list
> beam_charge = -10.0e-15  # in Coulombs
> N_beam_particles = int(1e6)
> beam_centroid_z = -107.0e-6
> beam_rms_z = 2.0e-6
> beam_gammas = [1960 + 13246 * i_stage for i_stage in range(N_stage)]
> # beam_gammas = [1957, 15188, 28432, 41678, 54926, 68174, 81423,94672, 107922,121171] # From 3D run
> beams = []
> for i_stage in range(N_stage):
>     beam_gamma = beam_gammas[i_stage]
>     sigma_gamma = 0.06 * beam_gamma
>     gaussian_distribution = picmi.GaussianBunchDistribution(
>         n_physical_particles=abs(beam_charge) / q_e,
>         rms_bunch_size=[2.0e-6, 2.0e-6, beam_rms_z],
>         rms_velocity=[8 * c, 8 * c, sigma_gamma * c],
>         centroid_position=[0.0, 0.0, beam_centroid_z],
>         centroid_velocity=[0.0, 0.0, beam_gamma * c],
>     )
>     beam = picmi.Species(
>         particle_type="electron",
>         name=f"beam_stage_{i_stage}",
>         initial_distribution=gaussian_distribution,
>     )
>     beams.append(beam)

> # Laser
> antenna_z = -1e-9
> profile_t_peak = 1.46764864e-13


> def get_laser(antenna_z, profile_t_peak, fill_in=True):
>     profile_focal_distance = 0.0
>     laser = picmi.GaussianLaser(
>         wavelength=0.8e-06,
>         waist=36e-06,
>         duration=7.33841e-14,
>         focal_position=[0.0, 0.0, profile_focal_distance + antenna_z],
>         centroid_position=[0.0, 0.0, antenna_z - c * profile_t_peak],
>         propagation_direction=[0.0, 0.0, 1.0],
>         polarization_direction=[0.0, 1.0, 0.0],
>         a0=2.36,
>         fill_in=fill_in,
>     )
>     laser_antenna = picmi.LaserAntenna(
>         position=[0.0, 0.0, antenna_z], normal_vector=[0.0, 0.0, 1.0]
>     )
>     return (laser, laser_antenna)


> lasers = []
> for i_stage in range(1):
>     fill_in = True
>     if i_stage == 0:
>         fill_in = False
>     lasers.append(
>         get_laser(
>             antenna_z + i_stage * stage_spacing,
>             profile_t_peak + i_stage * stage_spacing / c,
>             fill_in,
>         )
>     )

> # Electromagnetic solver

> psatd_algo = "psatd_JRhom"
> if psatd_algo == "galilean":
>     galilean_velocity = [0.0, 0.0] if dim == "3" else [0.0]
>     galilean_velocity += [-c * beta_boost]
>     n_pass_z = 1
>     psatd_JRhom = None
>     current_correction = True
>     divE_cleaning = False
> elif psatd_algo == "psatd_JRhom":
>     n_pass_z = 4
>     galilean_velocity = None
>     psatd_JRhom = "LL2"
>     current_correction = False
>     divE_cleaning = True
> else:
>     raise Exception(
>         f"PSATD algorithm '{psatd_algo}' is not recognized!\n"
>         "Valid options are 'psatd_JRhom' or 'galilean'."
>     )
> if dim == "rz":
>     stencil_order = [8, 16]
>     smoother = picmi.BinomialSmoother(n_pass=[1, n_pass_z])
>     grid_type = "collocated"
> else:
>     stencil_order = [8, 8, 16]
>     smoother = picmi.BinomialSmoother(n_pass=[1, 1, n_pass_z])
>     grid_type = "hybrid"


> solver = picmi.ElectromagneticSolver(
>     grid=grid,
>     method="PSATD",
>     cfl=0.9999,
>     source_smoother=smoother,
>     stencil_order=stencil_order,
>     galilean_velocity=galilean_velocity,
>     warpx_psatd_update_with_rho=True,
>     warpx_psatd_JRhom=psatd_JRhom,
>     warpx_current_correction=current_correction,
>     divE_cleaning=divE_cleaning,
> )

> # Diagnostics
> diag_field_list = ["B", "E", "J", "rho"]
> diag_particle_list = ["weighting", "position", "momentum"]
> coarse_btd_end = int((L_plasma_bulk + 0.001 + stage_spacing * (N_stage - 1)) * 100000)
> stage_end_snapshots = [
>     f"{int((L_plasma_bulk + stage_spacing * ii) * 100000)}:{int((L_plasma_bulk + stage_spacing * ii) * 100000 + 50)}:5"
>     for ii in range(1)
> ]
> btd_particle_diag = picmi.LabFrameParticleDiagnostic(
>     name="lab_particle_diags",
>     species=beams,
>     grid=grid,
>     num_snapshots=25 * N_stage,
>     # warpx_intervals=', '.join([f':{coarse_btd_end}:1000']+stage_end_snapshots),
>     warpx_intervals=", ".join(["0:0"] + stage_end_snapshots),
>     dt_snapshots=0.00001 / c,
>     data_list=diag_particle_list,
>     write_dir="lab_particle_diags",
>     warpx_format="openpmd",
>     warpx_openpmd_backend="bp5",
> )

> btd_field_diag = picmi.LabFrameFieldDiagnostic(
>     name="lab_field_diags",
>     grid=grid,
>     num_snapshots=25 * N_stage,
>     dt_snapshots=stage_spacing / 25 / c,
>     data_list=diag_field_list,
>     warpx_lower_bound=[-128.0e-6, 0.0e-6, -180.0e-6],
>     warpx_upper_bound=[128.0e-6, 0.0e-6, 0.0],
>     write_dir="lab_field_diags",
>     warpx_format="openpmd",
>     warpx_openpmd_backend="bp5",
> )

> field_diag = picmi.FieldDiagnostic(
>     name="field_diags",
>     data_list=diag_field_list,
>     grid=grid,
>     period=100,
>     write_dir="field_diags",
>     lower_bound=[-128.0e-6, 0.0e-6, -180.0e-6],
>     upper_bound=[128.0e-6, 0.0e-6, 0.0],
>     warpx_format="openpmd",
>     warpx_openpmd_backend="h5",
> )

> particle_diag = picmi.ParticleDiagnostic(
>     name="particle_diags",
>     species=beams,
>     period=100,
>     write_dir="particle_diags",
>     warpx_format="openpmd",
>     warpx_openpmd_backend="h5",
> )

> beamrel_red_diag = picmi.ReducedDiagnostic(
>     diag_type="BeamRelevant", name="beamrel", species=beam, period=1
> )

> # Set up simulation
> sim = picmi.Simulation(
>     solver=solver,
>     warpx_numprocs=num_procs,
>     warpx_compute_max_step_from_btd=True,
>     verbose=2,
>     particle_shape="cubic",
>     gamma_boost=gamma_boost,
>     warpx_charge_deposition_algo="standard",
>     warpx_current_deposition_algo="direct",
>     warpx_field_gathering_algo="momentum-conserving",
>     warpx_particle_pusher_algo="vay",
>     warpx_amrex_the_arena_is_managed=False,
>     warpx_amrex_use_gpu_aware_mpi=True,
>     warpx_grid_type=grid_type,
>     # default: 2 for staggered grids, 8 for hybrid grids
>     warpx_field_centering_order=[16, 16, 16],
>     # only for hybrid grids, default: 8
>     warpx_current_centering_order=[16, 16, 16],
> )

> for species in species_list:
>     if dim == "rz":
>         n_macroparticle_per_cell = [2, 4, 2]
>     else:
>         n_macroparticle_per_cell = [2, 2, 2]
>     sim.add_species(
>         species,
>         layout=picmi.GriddedLayout(
>             grid=grid, n_macroparticle_per_cell=n_macroparticle_per_cell
>         ),
>     )

> for i_stage in range(N_stage):
>     sim.add_species_through_plane(
>         species=beams[i_stage],
>         layout=picmi.PseudoRandomLayout(grid=grid, n_macroparticles=N_beam_particles),
>         injection_plane_position=0.0,
>         injection_plane_normal_vector=[0.0, 0.0, 1.0],
>     )

> for i_stage in range(1):
>     # Add laser
>     (laser, laser_antenna) = lasers[i_stage]
>     sim.add_laser(laser, injection_method=laser_antenna)

> # Add diagnostics
> sim.add_diagnostic(btd_particle_diag)
> # sim.add_diagnostic(btd_field_diag)
> # sim.add_diagnostic(field_diag)
> # sim.add_diagnostic(particle_diag)

> # Add reduced diagnostic
> sim.add_diagnostic(beamrel_red_diag)

> sim.write_input_file(f"inputs_training_{N_stage}_stages")

> # Advance simulation until last time step
> sim.step()
> ```

In this section we walk through a workflow for data processing and model training, using data from this input script as an example.
The simulation output is stored in an online [Zenodo archive](https://zenodo.org/records/10368972), in the `lab_particle_diags` directory.
In the example scripts provided here, the data is downloaded from the Zenodo archive, properly formatted, and used to train a neural network.
This workflow was developed and first presented in Sandberg *et al.* [[1](#id24)], Sandberg *et al.* [[2](#id23)].
It assumes you have an up-to-date environment with PyTorch and openPMD.

## Data Cleaning

It is important to inspect the data for artifacts, to
check that input/output data make sense.
If we plot the final phase space of the particle beam,
shown in [Fig. 23](#fig-unclean-phase-space).
we see outlying particles.
Looking closer at the z-pz space, we see that some particles were not trapped in the accelerating region of the wake and have much less energy than the rest of the beam.

<a id="fig-unclean-phase-space"></a>
![Phase space projections showing partially accelerated beam particles](https://gist.githubusercontent.com/RTSandberg/649a81cc0e7926684f103729483eff90/raw/095ac2daccbcf197fa4e18a8f8505711b27e807a/unclean_stage_0.png)

To assist our neural network in learning dynamics of interest, we filter out these particles.
It is sufficient for our purposes to select particles that are not too far back, setting
`particle_selection={'z':[0.280025, None]}`.
After filtering, we can see in [Fig. 24](#fig-clean-phase-space) that the beam phase space projections are much cleaner – this is the beam we want to train on.

<a id="fig-clean-phase-space"></a>
![Phase space projections of filtered beam particles](https://gist.githubusercontent.com/RTSandberg/649a81cc0e7926684f103729483eff90/raw/095ac2daccbcf197fa4e18a8f8505711b27e807a/clean_stage_0.png)

A particle tracker is set up to make sure
we consistently filter out these particles from both the initial and final data.

```python
iteration = ts.iterations[survivor_select_index]
pt = ParticleTracker(
    ts, species=species, iteration=iteration, select=particle_selection
)
```

This data cleaning ensures that the particle data is distributed in a single blob,
as is optimal for training neural networks.

## Create Normalized Dataset

Having chosen training data we are content with, we now need to format the data,
normalize it, and store the normalized data as well as the normalizations.
The script below will take the openPMD data we have selected and
format, normalize, and store it.

### Python dataset creation

> ```python
> #!/usr/bin/env python3
> #
> # Copyright 2023 The WarpX Community
> #
> # This file is part of WarpX.
> #
> # Authors: Ryan Sandberg
> # License: BSD-3-Clause-LBNL
> #

> import os
> import zipfile
> from urllib import request

> import numpy as np
> import torch
> from openpmd_viewer import OpenPMDTimeSeries, ParticleTracker

> c = 2.998e8
> ###############


> def sanitize_dir_strings(*dir_strings):
>     """append '/' to a string for concatenation in building up file tree descriptions"""
>     dir_strings = list(dir_strings)
>     for ii, dir_string in enumerate(dir_strings):
>         if dir_string[-1] != "/":
>             dir_strings[ii] = dir_string + "/"

>     return dir_strings


> def download_and_unzip(url, data_dir):
>     request.urlretrieve(url, data_dir)
>     with zipfile.ZipFile(data_dir, "r") as zip_dataset:
>         zip_dataset.extractall()


> def create_source_target_data(
>     data_dir,
>     species,
>     source_index=0,
>     target_index=-1,
>     survivor_select_index=-1,
>     particle_selection=None,
> ):
>     """Create dataset from openPMD files

>     Parameters
>     ---
>     data_dir : string, location of diagnostic data
>     source_index : int, which index to take source data from
>     target_index : int, which index to take target data from
>     particle_selection: dictionary, optional, selection criterion for dataset

>     Returns
>     ---
>     source_data:  Nx6 array of source particle data
>     source_means: 6 element array of source particle coordinate means
>     source_stds:  6 element array of source particle coordinate standard deviations
>     target_data:  Nx6 array of target particle data
>     target_means: 6 element array of target particle coordinate means
>     target_stds:  6 element array of source particle coordinate standard deviations
>     relevant times: 2 element array of source and target times
>     """
>     (data_dir,) = sanitize_dir_strings(data_dir)
>     data_path = data_dir
>     print("loading openPMD data from", data_path)
>     ts = OpenPMDTimeSeries(data_path)
>     relevant_times = [ts.t[source_index], ts.t[target_index]]

>     # Manual: Particle tracking START
>     iteration = ts.iterations[survivor_select_index]
>     pt = ParticleTracker(
>         ts, species=species, iteration=iteration, select=particle_selection
>     )
>     # Manual: Particle tracking END

>     #### create normalized source, target data sets ####
>     print("creating data sets")

>     # Manual: Load openPMD START
>     iteration = ts.iterations[source_index]
>     source_data = ts.get_particle(
>         species=species,
>         iteration=iteration,
>         var_list=["x", "y", "z", "ux", "uy", "uz"],
>         select=pt,
>     )

>     iteration = ts.iterations[target_index]
>     target_data = ts.get_particle(
>         species=species,
>         iteration=iteration,
>         var_list=["x", "y", "z", "ux", "uy", "uz"],
>         select=pt,
>     )
>     # Manual: Load openPMD END

>     # Manual: Normalization START
>     target_means = np.zeros(6)
>     target_stds = np.zeros(6)
>     source_means = np.zeros(6)
>     source_stds = np.zeros(6)
>     for jj in range(6):
>         source_means[jj] = source_data[jj].mean()
>         source_stds[jj] = source_data[jj].std()
>         source_data[jj] -= source_means[jj]
>         source_data[jj] /= source_stds[jj]

>     for jj in range(6):
>         target_means[jj] = target_data[jj].mean()
>         target_stds[jj] = target_data[jj].std()
>         target_data[jj] -= target_means[jj]
>         target_data[jj] /= target_stds[jj]
>     # Manual: Normalization END

>     # Manual: Format data START
>     source_data = torch.tensor(np.column_stack(source_data))
>     target_data = torch.tensor(np.column_stack(target_data))
>     # Manual: Format data END

>     return (
>         source_data,
>         source_means,
>         source_stds,
>         target_data,
>         target_means,
>         target_stds,
>         relevant_times,
>     )


> def save_warpx_surrogate_data(
>     dataset_fullpath_filename,
>     diag_dir,
>     species,
>     training_frac,
>     batch_size,
>     source_index,
>     target_index,
>     survivor_select_index,
>     particle_selection=None,
> ):
>     source_target_data = create_source_target_data(
>         data_dir=diag_dir,
>         species=species,
>         source_index=source_index,
>         target_index=target_index,
>         survivor_select_index=survivor_select_index,
>         particle_selection=particle_selection,
>     )
>     (
>         source_data,
>         source_means,
>         source_stds,
>         target_data,
>         target_means,
>         target_stds,
>         times,
>     ) = source_target_data

>     # Manual: Save dataset START
>     full_dataset = torch.utils.data.TensorDataset(
>         source_data.float(), target_data.float()
>     )

>     n_samples = full_dataset.tensors[0].size(0)
>     n_train = int(training_frac * n_samples)
>     n_test = n_samples - n_train

>     train_data, test_data = torch.utils.data.random_split(
>         full_dataset, [n_train, n_test]
>     )

>     torch.save(
>         {
>             "dataset": full_dataset,
>             "train_indices": train_data.indices,
>             "test_indices": test_data.indices,
>             "source_means": source_means,
>             "source_stds": source_stds,
>             "target_means": target_means,
>             "target_stds": target_stds,
>             "times": times,
>         },
>         dataset_fullpath_filename,
>     )
>     # Manual: Save dataset END


> ######## end utility functions #############
> ######## start dataset creation ############

> data_url = "https://zenodo.org/records/10810754/files/lab_particle_diags.zip?download=1"
> download_and_unzip(data_url, "lab_particle_diags.zip")
> data_dir = "lab_particle_diags/lab_particle_diags/"

> # create data set

> source_index = 0
> target_index = 1
> survivor_select_index = 1
> batch_size = 1200
> training_frac = 0.7

> os.makedirs("datasets", exist_ok=True)

> # improve stage 0 dataset
> stage_i = 0
> select = {"z": [0.280025, None]}
> species = f"beam_stage_{stage_i}"
> dataset_filename = f"dataset_{species}.pt"
> dataset_file = "datasets/" + dataset_filename
> save_warpx_surrogate_data(
>     dataset_fullpath_filename=dataset_file,
>     diag_dir=data_dir,
>     species=species,
>     training_frac=training_frac,
>     batch_size=batch_size,
>     source_index=source_index,
>     target_index=target_index,
>     survivor_select_index=survivor_select_index,
>     particle_selection=select,
> )

> for stage_i in range(1, 15):
>     species = f"beam_stage_{stage_i}"
>     dataset_filename = f"dataset_{species}.pt"
>     dataset_file = "datasets/" + dataset_filename
>     save_warpx_surrogate_data(
>         dataset_fullpath_filename=dataset_file,
>         diag_dir=data_dir,
>         species=species,
>         training_frac=training_frac,
>         batch_size=batch_size,
>         source_index=source_index,
>         target_index=target_index,
>         survivor_select_index=survivor_select_index,
>     )
> ```

### Load openPMD Data

First the openPMD data is loaded, using the particle selector as chosen above.
The neural network will make predictions from the initial phase space coordinates,
using the final phase space coordinates to measure how well it is making predictions.
Hence we load two sets of particle data, the source and target particle arrays.

```python
iteration = ts.iterations[source_index]
source_data = ts.get_particle(
    species=species,
    iteration=iteration,
    var_list=["x", "y", "z", "ux", "uy", "uz"],
    select=pt,
)

iteration = ts.iterations[target_index]
target_data = ts.get_particle(
    species=species,
    iteration=iteration,
    var_list=["x", "y", "z", "ux", "uy", "uz"],
    select=pt,
)
```

### Normalize Data

Neural networks learn better on appropriately normalized data.
Here we subtract out the mean in each coordinate direction and
divide by the standard deviation in each coordinate direction,
for normalized data that is centered on the origin with unit variance.

```python
target_means = np.zeros(6)
target_stds = np.zeros(6)
source_means = np.zeros(6)
source_stds = np.zeros(6)
for jj in range(6):
    source_means[jj] = source_data[jj].mean()
    source_stds[jj] = source_data[jj].std()
    source_data[jj] -= source_means[jj]
    source_data[jj] /= source_stds[jj]

for jj in range(6):
    target_means[jj] = target_data[jj].mean()
    target_stds[jj] = target_data[jj].std()
    target_data[jj] -= target_means[jj]
    target_data[jj] /= target_stds[jj]
```

### openPMD to PyTorch Data

With the data normalized, it must be stored in a form PyTorch recognizes.
The openPMD data are 6 lists of arrays, for each of the 6 phase space coordinates
$x, y, z, p_x, p_y,$ and $p_z$.
This data are converted to an $N\times 6$ numpy array and then to a PyTorch $N\times 6$ tensor.

```python
source_data = torch.tensor(np.column_stack(source_data))
target_data = torch.tensor(np.column_stack(target_data))
```

### Save Normalizations and Normalized Data

The data is split into training and testing subsets.
We take most of the data (70%) for training, meaning that data is used to update
the neural network parameters.
The testing data is reserved to determine how well the neural network generalizes;
that is, how well the neural network performs on data that wasn’t used to update the neural network parameters.
With the data split and properly normalized, it and the normalizations are saved to file for
use in training and inference.

```python
full_dataset = torch.utils.data.TensorDataset(
    source_data.float(), target_data.float()
)

n_samples = full_dataset.tensors[0].size(0)
n_train = int(training_frac * n_samples)
n_test = n_samples - n_train

train_data, test_data = torch.utils.data.random_split(
    full_dataset, [n_train, n_test]
)

torch.save(
    {
        "dataset": full_dataset,
        "train_indices": train_data.indices,
        "test_indices": test_data.indices,
        "source_means": source_means,
        "source_stds": source_stds,
        "target_means": target_means,
        "target_stds": target_stds,
        "times": times,
    },
    dataset_fullpath_filename,
)
```

## Neural Network Structure

It was found in Sandberg *et al.* [[2](#id23)] that a reasonable surrogate model is obtained with
shallow feedforward neural networks consisting of about 5 hidden layers and 700-900 nodes per layer.
The example shown here uses 3 hidden layers and 20 nodes per layer
and is trained for 10 epochs.

Some utility functions for creating neural networks are provided in the script below.
These are mostly convenience wrappers and utilities for working with [PyTorch](https://pytorch.org/) neural network objects.
This script is imported in the training scripts shown later.

### Python neural network class definitions

> ```python3
> #!/usr/bin/env python3
> #
> # Copyright 2023 The WarpX Community
> #
> # This file is part of WarpX.
> #
> # Authors: Ryan Sandberg
> # License: BSD-3-Clause-LBNL
> #
> from enum import Enum

> from torch import nn


> class ActivationType(Enum):
>     """
>     Activation class provides an enumeration type for the supported activation layers
>     """

>     ReLU = 1
>     Tanh = 2
>     PReLU = 3
>     Sigmoid = 4


> def get_enum_type(type_to_test, EnumClass):
>     """
>     Returns the enumeration type associated to type_to_test in EnumClass

>     Parameters
>     ----------
>     type_to_test: EnumClass, int or str
>         object whose Enum class is to be obtained
>     EnumClass: Enum class
>         Enum class to test
>     """
>     if type(type_to_test) is EnumClass:
>         return type_to_test
>     if type(type_to_test) is int:
>         return EnumClass(type_to_test)
>     if type(type_to_test) is str:
>         return getattr(EnumClass, type_to_test)
>     else:
>         raise Exception("unsupported type entered")


> class ConnectedNN(nn.Module):
>     """
>     ConnectedNN is a class of fully connected neural networks
>     """

>     def __init__(self, layers):
>         super().__init__()
>         self.stack = nn.Sequential(*layers)

>     def forward(self, x):
>         return self.stack(x)


> class OneActNN(ConnectedNN):
>     """
>     OneActNN is class of fully connected neural networks admitting only one activation function
>     """

>     def __init__(self, n_in, n_out, n_hidden_nodes, n_hidden_layers, act):
>         self.n_in = n_in
>         self.n_out = n_out
>         self.n_hidden_layers = n_hidden_layers
>         self.n_hidden_nodes = n_hidden_nodes

>         self.act = get_enum_type(act, ActivationType)

>         layers = [nn.Linear(self.n_in, self.n_hidden_nodes)]

>         for ii in range(self.n_hidden_layers):
>             if self.act is ActivationType.ReLU:
>                 layers += [nn.ReLU()]
>             if self.act is ActivationType.Tanh:
>                 layers += [nn.Tanh()]
>             if self.act is ActivationType.PReLU:
>                 layers += [nn.PReLU()]
>             if self.act is ActivationType.Sigmoid:
>                 layers += [nn.Sigmoid()]

>             if ii < self.n_hidden_layers - 1:
>                 layers += [nn.Linear(self.n_hidden_nodes, self.n_hidden_nodes)]

>         layers += [nn.Linear(self.n_hidden_nodes, self.n_out)]

>         super().__init__(layers)
> ```

## Train and Save Neural Network

The script below trains the neural network on the dataset just created.
In subsequent sections we discuss the various parts of the training process.

### Python neural network training

> ```python3
> #!/usr/bin/env python3
> #
> # Copyright 2023 The WarpX Community
> #
> # This file is part of WarpX.
> #
> # Authors: Ryan Sandberg
> # License: BSD-3-Clause-LBNL
> #
> import os
> import time

> import neural_network_classes as mynn
> import torch
> import torch.nn.functional as F
> import torch.optim as optim

> ############# set model parameters #################

> stage_i = 0
> species = f"beam_stage_{stage_i}"
> source_index = 0
> target_index = 1
> survivor_select_index = 1

> data_dim = 6
> n_in = data_dim
> n_out = data_dim

> learning_rate = 0.0001
> n_epochs = 10
> batch_size = 1200

> loss_fun = F.mse_loss

> n_hidden_nodes = 20
> n_hidden_layers = 3
> activation_type = "ReLU"

> device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
> print(f"device={device}")
> #################### load dataset ################
> dataset_filename = f"dataset_{species}.pt"
> dataset_file = "datasets/" + dataset_filename

> print(f"trying to load dataset+test-train split in {dataset_file}")

> dataset_with_indices = torch.load(dataset_file)
> train_data = torch.utils.data.dataset.Subset(
>     dataset_with_indices["dataset"], dataset_with_indices["train_indices"]
> )
> test_data = torch.utils.data.dataset.Subset(
>     dataset_with_indices["dataset"], dataset_with_indices["test_indices"]
> )
> source_data = dataset_with_indices["dataset"]
> source_means = dataset_with_indices["source_means"]
> source_stds = dataset_with_indices["source_stds"]
> target_means = dataset_with_indices["target_means"]
> target_stds = dataset_with_indices["target_stds"]
> print("able to load data and test/train split")

> ###### move data to device (GPU) if available ########
> source_device = train_data.dataset.tensors[0].to(
>     device
> )  # equivalently, test_data.tensors[0].to(device)
> target_device = train_data.dataset.tensors[1].to(device)
> full_dataset_device = torch.utils.data.TensorDataset(
>     source_device.float(), target_device.float()
> )

> train_data_device = torch.utils.data.dataset.Subset(
>     full_dataset_device, train_data.indices
> )
> test_data_device = torch.utils.data.dataset.Subset(
>     full_dataset_device, test_data.indices
> )

> train_loader_device = torch.utils.data.DataLoader(
>     train_data_device, batch_size=batch_size, shuffle=True
> )
> test_loader_device = torch.utils.data.DataLoader(
>     test_data_device, batch_size=batch_size, shuffle=True
> )

> test_source_device = test_data_device.dataset.tensors[0]
> test_target_device = test_data_device.dataset.tensors[1]

> training_set_size = len(train_data_device.indices)
> testing_set_size = len(test_data_device.indices)

> ###### create model ###########

> model = mynn.OneActNN(
>     n_in=n_in,
>     n_out=n_out,
>     n_hidden_nodes=n_hidden_nodes,
>     n_hidden_layers=n_hidden_layers,
>     act=activation_type,
> )

> training_time = 0
> train_loss_list = []
> test_loss_list = []

> model.to(device=device)


> ########## train and test functions ####
> # Manual: Train function START
> def train(model, optimizer, train_loader, loss_fun):
>     model.train()
>     total_loss = 0.0
>     for batch_idx, (data, target) in enumerate(train_loader):
>         # evaluate network with data
>         output = model(data)
>         # compute loss
>         # sum the differences squared, take mean afterward
>         loss = loss_fun(output, target, reduction="sum")
>         # backpropagation: step optimizer and reset gradients
>         loss.backward()
>         optimizer.step()
>         optimizer.zero_grad()
>         total_loss += loss.item()
>     return total_loss


> # Manual: Train function END


> def test(model, test_loader, loss_fun):
>     model.eval()
>     total_loss = 0.0
>     with torch.no_grad():
>         for batch_idx, (data, target) in enumerate(test_loader):
>             output = model(data)
>             total_loss += loss_fun(output, target, reduction="sum").item()
>     return total_loss


> # Manual: Test function START
> def test_dataset(model, test_source, test_target, loss_fun):
>     model.eval()
>     with torch.no_grad():
>         output = model(test_source)
>         return loss_fun(output, test_target, reduction="sum").item()


> # Manual: Test function END

> ######## training loop ########

> optimizer = optim.Adam(model.parameters(), lr=learning_rate)

> do_print = True

> t3 = time.time()
> # Manual: Training loop START
> for epoch in range(n_epochs):
>     if do_print:
>         t1 = time.time()
>     ave_train_loss = (
>         train(model, optimizer, train_loader_device, loss_fun)
>         / data_dim
>         / training_set_size
>     )
>     ave_test_loss = (
>         test_dataset(model, test_source_device, test_target_device, loss_fun)
>         / data_dim
>         / training_set_size
>     )
>     train_loss_list.append(ave_train_loss)
>     test_loss_list.append(ave_test_loss)

>     if do_print:
>         t2 = time.time()
>         print(
>             "Train Epoch: {:04d} \tTrain Loss: {:.6f} \tTest Loss: {:.6f}, this epoch: {:.3f} s".format(
>                 epoch + 1, ave_train_loss, ave_test_loss, t2 - t1
>             )
>         )
> # Manual: Training loop END
> t4 = time.time()
> print(f"total training time: {t4 - t3:.3f}s")

> ######### save model #########

> os.makedirs("models", exist_ok=True)

> # Manual: Save model START
> model.to(device="cpu")
> torch.save(
>     {
>         "n_hidden_layers": n_hidden_layers,
>         "n_hidden_nodes": n_hidden_nodes,
>         "activation": activation_type,
>         "model_state_dict": model.state_dict(),
>         "optimizer_state_dict": optimizer.state_dict(),
>         "train_loss_list": train_loss_list,
>         "test_loss_list": test_loss_list,
>         "training_time": training_time,
>     },
>     f"models/{species}_model.pt",
> )
> # Manual: Save model END
> ```

### Training Function

In the training function, the model weights are updated.
Iterating through batches, the loss function is evaluated on each batch.
PyTorch provides automatic differentiation, so the direction of steepest descent
is determined when the loss function is evaluated and the `loss.backward()` function
is invoked.
The optimizer uses this information to update the weights in the `optimizer.step()` call.
The training loop then resets the optimizer and updates the summed error for the whole dataset
with the error on the batch and continues iterating through batches.
Note that this function returns the sum of all errors across the entire dataset,
which is later divided by the size of the dataset in the training loop.

```python
def train(model, optimizer, train_loader, loss_fun):
    model.train()
    total_loss = 0.0
    for batch_idx, (data, target) in enumerate(train_loader):
        # evaluate network with data
        output = model(data)
        # compute loss
        # sum the differences squared, take mean afterward
        loss = loss_fun(output, target, reduction="sum")
        # backpropagation: step optimizer and reset gradients
        loss.backward()
        optimizer.step()
        optimizer.zero_grad()
        total_loss += loss.item()
    return total_loss


```

### Testing Function

The testing function just evaluates the neural network on the testing data that has not been used
to update the model parameters.
This testing function requires that the testing dataset is small enough to be loaded all at once.
The PyTorch dataloader can load data in batches if this size assumption is not satisfied.
The error, measured by the loss function, is returned by the testing function to be aggregated and stored.
Note that this function returns the sum of all errors across the entire dataset,
which is later divided by the size of the dataset in the training loop.

```python
def test_dataset(model, test_source, test_target, loss_fun):
    model.eval()
    with torch.no_grad():
        output = model(test_source)
        return loss_fun(output, test_target, reduction="sum").item()


```

### Training Loop

The full training loop performs `n_epochs` number of iterations.
At each iteration the training and testing functions are called,
the respective errors are divided by the size of the dataset and recorded,
and a status update is printed to the console.

```python
for epoch in range(n_epochs):
    if do_print:
        t1 = time.time()
    ave_train_loss = (
        train(model, optimizer, train_loader_device, loss_fun)
        / data_dim
        / training_set_size
    )
    ave_test_loss = (
        test_dataset(model, test_source_device, test_target_device, loss_fun)
        / data_dim
        / training_set_size
    )
    train_loss_list.append(ave_train_loss)
    test_loss_list.append(ave_test_loss)

    if do_print:
        t2 = time.time()
        print(
            "Train Epoch: {:04d} \tTrain Loss: {:.6f} \tTest Loss: {:.6f}, this epoch: {:.3f} s".format(
                epoch + 1, ave_train_loss, ave_test_loss, t2 - t1
            )
        )
```

### Save Neural Network Parameters

The model weights are saved after training to record the updates to the model parameters.
Additionally, we save some model metainformation with the model for convenience,
including the model hyperparameters, the training and testing losses, and how long the training took.

```python
model.to(device="cpu")
torch.save(
    {
        "n_hidden_layers": n_hidden_layers,
        "n_hidden_nodes": n_hidden_nodes,
        "activation": activation_type,
        "model_state_dict": model.state_dict(),
        "optimizer_state_dict": optimizer.state_dict(),
        "train_loss_list": train_loss_list,
        "test_loss_list": test_loss_list,
        "training_time": training_time,
    },
    f"models/{species}_model.pt",
)
```

## Evaluate

In this section we show two ways to diagnose how well the neural network is learning the data.
First we consider the train-test loss curves, shown in [Fig. 25](#fig-train-test-loss).
This figure shows the model error on the training data (in blue) and testing data (in green) as a function of the number of epochs seen.
The training data is used to update the model parameters, so training error should be lower than testing error.
A key feature to look for in the train-test loss curve is the inflection point in the test loss trend.
The testing data is set aside as a sample of data the neural network hasn’t seen before.
The testing error serves as a metric of model generalizability, indicating how well the model performs
on data it hasn’t seen yet.
When the test-loss starts to trend flat or even upward, the neural network is no longer improving its ability to generalize to new data.

<a id="fig-train-test-loss"></a>
![Training and testing loss evolution over epochs](https://gist.githubusercontent.com/RTSandberg/649a81cc0e7926684f103729483eff90/raw/095ac2daccbcf197fa4e18a8f8505711b27e807a/beam_stage_0_training_testing_error.png)

<a id="fig-train-evaluation"></a>
![Model predictions compared with simulation results](https://gist.githubusercontent.com/RTSandberg/649a81cc0e7926684f103729483eff90/raw/095ac2daccbcf197fa4e18a8f8505711b27e807a/beam_stage_0_model_evaluation.png)

A visual inspection of the model prediction can be seen in [Fig. 26](#fig-train-evaluation).
This plot compares the model prediction, with dots colored by mean-square error, on the testing data with the actual simulation output in black.
The model obtained with the hyperparameters chosen here trains quickly but is not very accurate.
A more accurate model is obtained with 5 hidden layers and 900 nodes per layer,
as discussed in Sandberg *et al.* [[2](#id23)].

These figures can be generated with the following Python script.

### Python visualization of progress training neural network

> ```python3
> #!/usr/bin/env python3
> #
> # Copyright 2023 The WarpX Community
> #
> # This file is part of WarpX.
> #
> # Authors: Ryan Sandberg
> # License: BSD-3-Clause-LBNL
> #
> import neural_network_classes as mynn
> import numpy as np
> import torch
> import torch.nn.functional as F
> from matplotlib import pyplot as plt

> c = 2.998e8


> # open model file
> stage_i = 0
> species = f"beam_stage_{stage_i}"
> model_data = torch.load(f"models/{species}_model.pt", map_location=torch.device("cpu"))
> data_dim = 6
> n_in = data_dim
> n_out = data_dim
> n_hidden_layers = model_data["n_hidden_layers"]
> n_hidden_nodes = model_data["n_hidden_nodes"]
> activation_type = model_data["activation"]
> train_loss_list = model_data["train_loss_list"]
> test_loss_list = model_data["test_loss_list"]
> training_time = model_data["training_time"]
> loss_fun = F.mse_loss


> n_epochs = len(train_loss_list)
> train_counter = np.arange(n_epochs) + 1
> test_counter = train_counter

> do_log_plot = False
> fig, ax = plt.subplots()
> if do_log_plot:
>     ax.semilogy(
>         train_counter, train_loss_list, ".-", color="blue", label="training loss"
>     )
>     ax.semilogy(test_counter, test_loss_list, color="green", label="testing loss")
> else:
>     ax.plot(train_counter, train_loss_list, ".-", color="blue", label="training loss")
>     ax.plot(test_counter, test_loss_list, color="green", label="testing loss")
> ax.set_xlabel("number of epochs seen")
> ax.set_ylabel(" loss")
> ax.legend()
> fig_dir = "figures/"
> ax.set_title(f"final test error = {test_loss_list[-1]:.3e} ")
> ax.grid()
> plt.tight_layout()
> plt.savefig(f"{species}_training_testing_error.png")


> ######### plot phase space comparison #######
> device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
> print(f"device={device}")

> model = mynn.OneActNN(
>     n_in=n_in,
>     n_out=n_out,
>     n_hidden_nodes=n_hidden_nodes,
>     n_hidden_layers=n_hidden_layers,
>     act=activation_type,
> )
> model.load_state_dict(model_data["model_state_dict"])
> model.to(device=device)

> ###### load model data ###############
> dataset_filename = f"dataset_{species}.pt"
> dataset_dir = "datasets/"
> model_input_data = torch.load(dataset_dir + dataset_filename)
> dataset = model_input_data["dataset"]
> train_indices = model_input_data["train_indices"]
> test_indices = model_input_data["test_indices"]
> source_means = model_input_data["source_means"]
> source_stds = model_input_data["source_stds"]
> target_means = model_input_data["target_means"]
> target_stds = model_input_data["target_stds"]
> source_time, target_time = model_input_data["times"]


> source = dataset.tensors[0]
> test_source = source[test_indices]
> test_source_device = test_source.to(device)
> with torch.no_grad():
>     evaluation_device = model(test_source_device.float())
> eval_cpu = evaluation_device.to("cpu")

> target = dataset.tensors[1]
> test_target = target[test_indices]

> target_si = test_target * target_stds + target_means
> eval_cpu_si = eval_cpu * target_stds + target_means
> target_mu = np.copy(target_si)
> eval_cpu_mu = np.copy(eval_cpu_si)
> target_mu[:, 2] -= c * target_time
> eval_cpu_mu[:, 2] -= c * target_time
> target_mu[:, :3] *= 1e6
> eval_cpu_mu[:, :3] *= 1e6


> loss_tensor = torch.sum(loss_fun(eval_cpu, test_target, reduction="none"), axis=1) / 6
> loss_array = loss_tensor.detach().numpy()

> tinds = np.nonzero(loss_array > 0.0)[0]
> skip = 10

> plt.figure()
> fig, axT = plt.subplots(3, 3)
> axes_label = {
>     0: r"x [$\mu$m]",
>     1: r"y [$\mu$m]",
>     2: r"z - %.2f cm [$\mu$m]" % (c * target_time),
>     3: r"$p_x$",
>     4: r"$p_y$",
>     5: r"$p_z$",
> }
> xy_inds = [(0, 1), (2, 0), (2, 1)]


> def set_axes(ax, indx, indy):
>     ax.scatter(
>         target_mu[::skip, indx], target_mu[::skip, indy], s=8, c="k", label="simulation"
>     )
>     ax.scatter(
>         eval_cpu_mu[::skip, indx],
>         eval_cpu_mu[::skip, indy],
>         marker="*",
>         c=loss_array[::skip],
>         s=0.02,
>         label="surrogate",
>         cmap="YlOrRd",
>     )
>     ax.set_xlabel(axes_label[indx])
>     ax.set_ylabel(axes_label[indy])
>     # return


> for ii in range(3):
>     ax = axT[0, ii]
>     indx, indy = xy_inds[ii]
>     set_axes(ax, indx, indy)

> for ii in range(2):
>     indx, indy = xy_inds[ii]
>     ax = axT[1, ii]
>     set_axes(ax, indx + 3, indy + 3)

> for ii in range(3):
>     ax = axT[2, ii]
>     indx = ii
>     indy = ii + 3
>     set_axes(ax, indx, indy)


> ax = axT[1, 2]
> indx = 5
> indy = 4
> ax.scatter(
>     target_mu[::skip, indx], target_mu[::skip, indy], s=8, c="k", label="simulation"
> )
> evalplt = ax.scatter(
>     eval_cpu_mu[::skip, indx],
>     eval_cpu_mu[::skip, indy],
>     marker="*",
>     c=loss_array[::skip],
>     s=2,
>     label="surrogate",
>     cmap="YlOrRd",
> )
> ax.set_xlabel(axes_label[indx])
> ax.set_ylabel(axes_label[indy])

> cb = plt.colorbar(evalplt, ax=ax)
> cb.set_label("MSE loss")

> fig.suptitle(f"stage {stage_i} prediction")

> plt.tight_layout()

> plt.savefig(f"{species}_model_evaluation.png")
> ```

### Surrogate Usage in Accelerator Physics

A neural network such as the one we trained here can be incorporated in other BLAST codes.
Consider this [example using neural network surrogates of WarpX simulations in ImpactX](https://impactx.readthedocs.io/en/latest/usage/examples/pytorch_surrogate_model/README.html).
