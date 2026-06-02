<a id="examples-beam-beam-collision"></a>

# Beam-beam collision

This example shows how to simulate the collision between two ultra-relativistic particle beams.
This is representative of what happens at the interaction point of a linear collider.
We consider a right-propagating electron bunch colliding against a left-propagating positron bunch.

We turn on the Quantum Synchrotron QED module for photon emission (also known as beamstrahlung in the collider community) and
the Breit-Wheeler QED module for the generation of electron-positron pairs (also known as coherent pair generation in the collider community).

To solve for the electromagnetic field we use the nodal version of the electrostatic relativistic solver.
This solver computes the average velocity of each species, and solves the corresponding relativistic Poisson equation (see the WarpX documentation for warpx.do_electrostatic = relativistic for more detail).
This solver accurately reproduces the subtle cancellation that occur for some component of `E + v x B`, which are crucial in simulations of relativistic particles.

This example is based on the following paper Yakimenko *et al.* [[9](../../examples.md#id26)].

## Run

The PICMI input file is not available for this example yet.

For [MPI-parallel](https://www.mpi-forum.org) runs, prefix these lines with `mpiexec -n 4 ...` or `srun -n 4 ...`, depending on the system.

```none
#################################
########## MY CONSTANTS #########
#################################
my_constants.mc2 = m_e*clight*clight
my_constants.nano = 1.0e-9
my_constants.GeV = q_e*1.e9

# BEAMS
my_constants.beam_energy = 125.*GeV
my_constants.beam_uz = beam_energy/(mc2)
my_constants.beam_charge = 0.14*nano
my_constants.sigmax = 10*nano
my_constants.sigmay = 10*nano
my_constants.sigmaz = 10*nano
my_constants.beam_uth = 0.1/100.*beam_uz
my_constants.n0 = beam_charge / (q_e * sigmax * sigmay * sigmaz * (2.*pi)**(3./2.))
my_constants.omegab = sqrt(n0 * q_e**2 / (epsilon0*m_e))
my_constants.mux = 0.0
my_constants.muy = 0.0
my_constants.muz = -0.5*Lz+3.2*sigmaz

# BOX
my_constants.Lx = 100.0*clight/omegab
my_constants.Ly = 100.0*clight/omegab
my_constants.Lz = 180.0*clight/omegab

# for a full scale simulation use: nx, ny, nz = 512, 512, 1024
my_constants.nx = 64
my_constants.ny = 64
my_constants.nz = 64

# TIME
my_constants.T = 0.7*Lz/clight
my_constants.dt = sigmaz/clight/10.

# DIAGS
my_constants.every_red = 1.
warpx.used_inputs_file = warpx_used_inputs.txt

#################################
####### GENERAL PARAMETERS ######
#################################
stop_time = T
amr.n_cell = nx ny nz
amr.max_grid_size = 128
amr.blocking_factor = 2
amr.max_level = 0
geometry.dims = 3
geometry.prob_lo = -0.5*Lx -0.5*Ly -0.5*Lz
geometry.prob_hi =  0.5*Lx  0.5*Ly  0.5*Lz

#################################
######## BOUNDARY CONDITION #####
#################################
boundary.field_lo = PEC PEC PEC
boundary.field_hi = PEC PEC PEC
boundary.particle_lo = Absorbing Absorbing Absorbing
boundary.particle_hi = Absorbing Absorbing Absorbing

#################################
############ NUMERICS ###########
#################################
warpx.do_electrostatic = relativistic
warpx.const_dt = dt
warpx.grid_type = collocated
algo.particle_shape = 3
algo.load_balance_intervals=100
algo.particle_pusher = vay

#################################
########### PARTICLES ###########
#################################
particles.species_names = beam1 beam2 pho1 pho2 ele1 pos1 ele2 pos2

beam1.species_type = electron
beam1.injection_style = NUniformPerCell
beam1.num_particles_per_cell_each_dim = 1 1 1
beam1.profile = parse_density_function
beam1.density_function(x,y,z) = "n0 *  exp(-(x-mux)**2/(2*sigmax**2))  * exp(-(y-muy)**2/(2*sigmay**2)) * exp(-(z-muz)**2/(2*sigmaz**2))"
beam1.density_min = n0 / 1e3
beam1.momentum_distribution_type = gaussian
beam1.uz_m = beam_uz
beam1.uy_m = 0.0
beam1.ux_m = 0.0
beam1.ux_th = beam_uth
beam1.uy_th = beam_uth
beam1.uz_th = beam_uth
beam1.initialize_self_fields = 1
beam1.self_fields_required_precision = 5e-10
beam1.self_fields_max_iters = 10000
beam1.do_qed_quantum_sync = 1
beam1.qed_quantum_sync_phot_product_species = pho1
beam1.do_classical_radiation_reaction = 0

beam2.species_type = positron
beam2.injection_style = NUniformPerCell
beam2.num_particles_per_cell_each_dim = 1 1 1
beam2.profile = parse_density_function
beam2.density_function(x,y,z) = "n0 *  exp(-(x-mux)**2/(2*sigmax**2))  * exp(-(y-muy)**2/(2*sigmay**2)) * exp(-(z+muz)**2/(2*sigmaz**2))"
beam2.density_min = n0 / 1e3
beam2.momentum_distribution_type = gaussian
beam2.uz_m = -beam_uz
beam2.uy_m = 0.0
beam2.ux_m = 0.0
beam2.ux_th = beam_uth
beam2.uy_th = beam_uth
beam2.uz_th = beam_uth
beam2.initialize_self_fields = 1
beam2.self_fields_required_precision = 5e-10
beam2.self_fields_max_iters = 10000
beam2.do_qed_quantum_sync = 1
beam2.qed_quantum_sync_phot_product_species = pho2
beam2.do_classical_radiation_reaction = 0

pho1.species_type = photon
pho1.injection_style = none
pho1.do_qed_breit_wheeler = 1
pho1.qed_breit_wheeler_ele_product_species = ele1
pho1.qed_breit_wheeler_pos_product_species = pos1

pho2.species_type = photon
pho2.injection_style = none
pho2.do_qed_breit_wheeler = 1
pho2.qed_breit_wheeler_ele_product_species = ele2
pho2.qed_breit_wheeler_pos_product_species = pos2

ele1.species_type = electron
ele1.injection_style = none
ele1.self_fields_required_precision = 1e-11
ele1.self_fields_max_iters = 10000
ele1.do_qed_quantum_sync = 1
ele1.qed_quantum_sync_phot_product_species = pho1
ele1.do_classical_radiation_reaction = 0

pos1.species_type = positron
pos1.injection_style = none
pos1.self_fields_required_precision = 1e-11
pos1.self_fields_max_iters = 10000
pos1.do_qed_quantum_sync = 1
pos1.qed_quantum_sync_phot_product_species = pho1
pos1.do_classical_radiation_reaction = 0

ele2.species_type = electron
ele2.injection_style = none
ele2.self_fields_required_precision = 1e-11
ele2.self_fields_max_iters = 10000
ele2.do_qed_quantum_sync = 1
ele2.qed_quantum_sync_phot_product_species = pho2
ele2.do_classical_radiation_reaction = 0

pos2.species_type = positron
pos2.injection_style = none
pos2.self_fields_required_precision = 1e-11
pos2.self_fields_max_iters = 10000
pos2.do_qed_quantum_sync = 1
pos2.qed_quantum_sync_phot_product_species = pho2
pos2.do_classical_radiation_reaction = 0

pho1.species_type = photon
pho1.injection_style = none
pho1.do_qed_breit_wheeler = 1
pho1.qed_breit_wheeler_ele_product_species = ele1
pho1.qed_breit_wheeler_pos_product_species = pos1

pho2.species_type = photon
pho2.injection_style = none
pho2.do_qed_breit_wheeler = 1
pho2.qed_breit_wheeler_ele_product_species = ele2
pho2.qed_breit_wheeler_pos_product_species = pos2

#################################
############# QED ###############
#################################
qed_qs.photon_creation_energy_threshold = 0.

qed_qs.lookup_table_mode = builtin
qed_qs.chi_min = 1.e-3

qed_bw.lookup_table_mode = builtin
qed_bw.chi_min = 1.e-2

# for accurate results use the generated tables with
# the following parameters
# note: must compile with -DWarpX_QED_TABLE_GEN=ON
#qed_qs.lookup_table_mode = generate
#qed_bw.lookup_table_mode = generate
#qed_qs.tab_dndt_chi_min=1e-3
#qed_qs.tab_dndt_chi_max=2e3
#qed_qs.tab_dndt_how_many=512
#qed_qs.tab_em_chi_min=1e-3
#qed_qs.tab_em_chi_max=2e3
#qed_qs.tab_em_chi_how_many=512
#qed_qs.tab_em_frac_how_many=512
#qed_qs.tab_em_frac_min=1e-12
#qed_qs.save_table_in=my_qs_table.txt
#qed_bw.tab_dndt_chi_min=1e-2
#qed_bw.tab_dndt_chi_max=2e3
#qed_bw.tab_dndt_how_many=512
#qed_bw.tab_pair_chi_min=1e-2
#qed_bw.tab_pair_chi_max=2e3
#qed_bw.tab_pair_chi_how_many=512
#qed_bw.tab_pair_frac_how_many=512
#qed_bw.save_table_in=my_bw_table.txt

# if you wish to use existing tables:
#qed_qs.lookup_table_mode=load
#qed_qs.load_table_from = /path/to/my_qs_table.txt
#qed_bw.lookup_table_mode=load
#qed_bw.load_table_from = /path/to/my_bw_table.txt

warpx.do_qed_schwinger = 0.

#################################
######### DIAGNOSTICS ###########
#################################
# FULL
diagnostics.diags_names = diag1

diag1.intervals = 15
diag1.diag_type = Full
diag1.write_species = 1
diag1.fields_to_plot = Ex Ey Ez Bx By Bz rho_beam1 rho_beam2 rho_ele1 rho_pos1 rho_ele2 rho_pos2
diag1.format = openpmd
diag1.dump_last_timestep = 1
diag1.species = pho1 pho2 ele1 pos1 ele2 pos2 beam1 beam2

# REDUCED
warpx.reduced_diags_names = ParticleNumber ColliderRelevant_beam1_beam2

ColliderRelevant_beam1_beam2.type = ColliderRelevant
ColliderRelevant_beam1_beam2.intervals = every_red
ColliderRelevant_beam1_beam2.species = beam1 beam2

ParticleNumber.type = ParticleNumber
ParticleNumber.intervals = every_red
```

## QED tables

The quantum synchrotron and nonlinear Breit-Wheeler modules are based on a Monte Carlo algorithm that computes the probabilities of an event from tabulated values.
WarpX comes with builtin tables (see the input file above), however these are low resolution tables that may not provide accurate results.
There are two ways to generate your own lookup table:

* Inside WarpX, at runtime: the tables are generated by WarpX itself at the beginning of the simulation.
  This requires to compile WarpX with `-DWarpX_QED_TABLE_GEN=ON` and to add the desired tables parameters in WarpX’s input file.
  [Here](https://warpx.readthedocs.io/en/latest/usage/parameters.html#lookup-tables-and-other-settings-for-qed-modules)  are more details.
* Outside of WarpX, using an external table generator: the tables are pregenerated, before running the actual simulation.
  This standalone tool can be compiled at the same time as WarpX using `-DWarpX_QED_TOOLS=ON`.
  The table parameters are then passed to the table generator and do not need to be added to WarpX’s input file.
  [Here](https://warpx.readthedocs.io/en/latest/usage/workflows/generate_lookup_tables_with_tools.html) are more details.

Once the tables have been generated, they can be loaded in the input file using
`qed_qs,bw.lookup_table_mode=load` and `qed_qs,bw.load_table_from=/path/to/your/table`.

## Visualize

The figure below shows the number of photons emitted per beam particle (left) and the number of secondary pairs generated per beam particle (right).

We compare different results for the reduced diagnostics with the literature:
\* (red) simplified WarpX simulation as the example stored in the directory `/Examples/Physics_applications/beam-beam_collision`;
\* (blue) large-scale WarpX simulation (high resolution and ad hoc generated tables ;
\* (black) literature results from Yakimenko *et al.* [[9](../../examples.md#id26)].

The small-scale simulation has been performed with a resolution of `nx = 64, ny = 64, nz = 64` grid cells, while the large-scale one has a much higher resolution of `nx = 512, ny = 512, nz = 1024`.
Moreover, the large-scale simulation uses dedicated QED lookup tables instead of the builtin tables.
For the large-scale simulation we have used the following options (added to the input file):

```ini
qed_qs.lookup_table_mode = generate
qed_bw.lookup_table_mode = generate

qed_qs.tab_dndt_chi_min=1e-3
qed_qs.tab_dndt_chi_max=2e3
qed_qs.tab_dndt_how_many=512
qed_qs.tab_em_chi_min=1e-3
qed_qs.tab_em_chi_max=2e3
qed_qs.tab_em_chi_how_many=512
qed_qs.tab_em_frac_how_many=512
qed_qs.tab_em_frac_min=1e-12
qed_qs.save_table_in=my_qs_table.txt

qed_bw.tab_dndt_chi_min=1e-2
qed_bw.tab_dndt_chi_max=2e3
qed_bw.tab_dndt_how_many=512
qed_bw.tab_pair_chi_min=1e-2
qed_bw.tab_pair_chi_max=2e3
qed_bw.tab_pair_chi_how_many=512
qed_bw.tab_pair_frac_how_many=512
qed_bw.save_table_in=my_bw_table.txt
```

The same table can be also obtained using the table generator with the following lines:

```ini
./qed_table_generator --table QS --mode DP -o my_qs_table.txt \
                      --dndt_chi_min 1e-3 --dndt_chi_max 2e3 --dndt_how_many 512 \
                      --em_chi_min 1e-3 --em_chi_max 2e3 --em_frac_min 1e-12 --em_chi_how_many 512 --em_frac_how_many 512


./qed_table_generator --table BW --mode DP -o my_bw_table.txt \
                      --dndt_chi_min 1e-2 --dndt_chi_max 2e3 --dndt_how_many 512 --pair_chi_min 1e-2 --pair_chi_max 2e3 --pair_chi_how_many 512 --pair_frac_how_many 512
```

![Beam-beam collision benchmark against :cite:t:`ex-Yakimenko2019`.](https://gist.github.com/user-attachments/assets/2dd43782-d039-4faa-9d27-e3cf8fb17352)

Below are two visualizations scripts that provide examples to graph the field and reduced diagnostics.
They are available in the `Examples/Physics_applications/beam-beam_collision/` folder and can be run as simply as `python3 plot_fields.py` and `python3 plot_reduced.py`.

### Field Diagnostics

This script visualizes the evolution of the fields ($|E|, |B|, \rho$) during the collision between the two ultra-relativistic lepton beams.
The magnitude of E and B and the charge densities of the primary beams and of the secondary pairs are sliced along either one of the two transverse coordinates ($x$ and $y$).

```python3
#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from openpmd_viewer import OpenPMDTimeSeries

plt.rcParams.update({"font.size": 16})

series = OpenPMDTimeSeries("./diags/diag1")
steps = series.iterations


for slice_axis in ["x", "y"]:  # slice the fields along x and y
    for n in steps:  # loop through the available timesteps
        fig, ax = plt.subplots(
            ncols=2, nrows=2, figsize=(10, 6), dpi=300, sharex=True, sharey=True
        )

        # get E field
        Ex, info = series.get_field(
            field="E", coord="x", iteration=n, plot=False, slice_across=slice_axis
        )
        Ey, info = series.get_field(
            field="E", coord="y", iteration=n, plot=False, slice_across=slice_axis
        )
        Ez, info = series.get_field(
            field="E", coord="z", iteration=n, plot=False, slice_across=slice_axis
        )
        # get B field
        Bx, info = series.get_field(
            field="B", coord="x", iteration=n, plot=False, slice_across=slice_axis
        )
        By, info = series.get_field(
            field="B", coord="y", iteration=n, plot=False, slice_across=slice_axis
        )
        Bz, info = series.get_field(
            field="B", coord="z", iteration=n, plot=False, slice_across=slice_axis
        )
        # get charge densities
        rho_beam1, info = series.get_field(
            field="rho_beam1", iteration=n, plot=False, slice_across=slice_axis
        )
        rho_beam2, info = series.get_field(
            field="rho_beam2", iteration=n, plot=False, slice_across=slice_axis
        )
        rho_ele1, info = series.get_field(
            field="rho_ele1", iteration=n, plot=False, slice_across=slice_axis
        )
        rho_pos1, info = series.get_field(
            field="rho_pos1", iteration=n, plot=False, slice_across=slice_axis
        )
        rho_ele2, info = series.get_field(
            field="rho_ele2", iteration=n, plot=False, slice_across=slice_axis
        )
        rho_pos2, info = series.get_field(
            field="rho_pos2", iteration=n, plot=False, slice_across=slice_axis
        )

        xmin = info.z.min()
        xmax = info.z.max()
        xlabel = "z [m]"

        if slice_axis == "x":
            ymin = info.y.min()
            ymax = info.y.max()
            ylabel = "y [m]"
        elif slice_axis == "y":
            ymin = info.x.min()
            ymax = info.x.max()
            ylabel = "x [m]"

        # plot E magnitude
        Emag = np.sqrt(Ex**2 + Ey**2 + Ez**2)
        im = ax[0, 0].imshow(
            np.transpose(Emag),
            cmap="seismic",
            extent=[xmin, xmax, ymin, ymax],
            vmin=0,
            vmax=np.max(np.abs(Emag)),
        )
        ax[0, 0].set_title("E [V/m]")
        divider = make_axes_locatable(ax[0, 0])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax, orientation="vertical")

        # plot B magnitude
        Bmag = np.sqrt(Bx**2 + By**2 + Bz**2)
        im = ax[1, 0].imshow(
            np.transpose(Bmag),
            cmap="seismic",
            extent=[xmin, xmax, ymin, ymax],
            vmin=0,
            vmax=np.max(np.abs(Bmag)),
        )
        ax[1, 0].set_title("B [T]")
        divider = make_axes_locatable(ax[1, 0])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax, orientation="vertical")

        # plot beam densities
        rho_beams = rho_beam1 + rho_beam2
        im = ax[0, 1].imshow(
            np.transpose(rho_beams),
            cmap="seismic",
            extent=[xmin, xmax, ymin, ymax],
            vmin=-np.max(np.abs(rho_beams)),
            vmax=np.max(np.abs(rho_beams)),
        )
        ax[0, 1].set_title(r"$\rho$ beams [C/m$^3$]")
        divider = make_axes_locatable(ax[0, 1])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax, orientation="vertical")

        # plot secondary densities
        rho2 = rho_ele1 + rho_pos1 + rho_ele2 + rho_pos2
        im = ax[1, 1].imshow(
            np.transpose(rho2),
            cmap="seismic",
            extent=[xmin, xmax, ymin, ymax],
            vmin=-np.max(np.abs(rho2)),
            vmax=np.max(np.abs(rho2)),
        )
        ax[1, 1].set_title(r"$\rho$ secondaries [C/m$^3$]")
        divider = make_axes_locatable(ax[1, 1])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax, orientation="vertical")

        for a in ax[-1, :].reshape(-1):
            a.set_xlabel(xlabel)
        for a in ax[:, 0].reshape(-1):
            a.set_ylabel(ylabel)

        fig.suptitle(f"Iteration = {n}, time [s] = {series.current_t}", fontsize=20)
        plt.tight_layout()

        image_file_name = "FIELDS_" + slice_axis + f"_{n:03d}.png"
        plt.savefig(image_file_name, dpi=100, bbox_inches="tight")
        plt.close()
```

![Slice across :math:`x` of different fields (:math:`|E|, |B|, \rho`) at timestep 45, in the middle of the collision.](https://gist.github.com/user-attachments/assets/04c9c0ec-b580-446f-a11a-437c1b244a41)

### Reduced Diagnostics

A similar script to the one below was used to produce the image showing the benchmark against Yakimenko *et al.* [[9](../../examples.md#id26)].

```python3
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.constants import c, nano, physical_constants

r_e = physical_constants["classical electron radius"][0]
my_dpi = 300
sigmaz = 10 * nano

fig, ax = plt.subplots(
    ncols=2, nrows=1, figsize=(2000.0 / my_dpi, 1000.0 / my_dpi), dpi=my_dpi
)

rdir = "./diags/reducedfiles/"

df_cr = pd.read_csv(f"{rdir}" + "ColliderRelevant_beam1_beam2.txt", sep=" ", header=0)
df_pn = pd.read_csv(f"{rdir}" + "ParticleNumber.txt", sep=" ", header=0)


times = df_cr[[col for col in df_cr.columns if "]time" in col]].to_numpy()
steps = df_cr[[col for col in df_cr.columns if "]step" in col]].to_numpy()

x = df_cr[[col for col in df_cr.columns if "]dL_dt" in col]].to_numpy()
coll_index = np.argmax(x)
coll_time = times[coll_index]

# number of photons per beam particle
np1 = df_pn[[col for col in df_pn.columns if "]pho1_weight" in col]].to_numpy()
np2 = df_pn[[col for col in df_pn.columns if "]pho2_weight" in col]].to_numpy()
Ne = df_pn[[col for col in df_pn.columns if "]beam1_weight" in col]].to_numpy()[0]
Np = df_pn[[col for col in df_pn.columns if "]beam2_weight" in col]].to_numpy()[0]

ax[0].plot((times - coll_time) / (sigmaz / c), (np1 + np2) / (Ne + Np), lw=2)
ax[0].set_title(r"photon number/beam particle")

# number of NLBW particles per beam particle
e1 = df_pn[[col for col in df_pn.columns if "]ele1_weight" in col]].to_numpy()
e2 = df_pn[[col for col in df_pn.columns if "]ele2_weight" in col]].to_numpy()

ax[1].plot((times - coll_time) / (sigmaz / c), (e1 + e2) / (Ne + Np), lw=2)
ax[1].set_title(r"NLBW particles/beam particle")

for a in ax.reshape(-1):
    a.set_xlabel(r"time [$\sigma_z/c$]")
image_file_name = "reduced.png"
plt.tight_layout()
plt.savefig(image_file_name, dpi=300, bbox_inches="tight")
plt.close("all")
```

![Photon and pair production rates in time throughout the collision.](https://gist.github.com/user-attachments/assets/c280490a-f1f2-4329-ad3c-46817d245dc1)
