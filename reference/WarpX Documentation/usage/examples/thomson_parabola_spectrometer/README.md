<a id="examples-thomson-parabola-spectrometer"></a>

# Thomson Parabola Spectrometer

This example simulates a Thomson parabola spectrometer (TPS) Rhee *et al.* [[13](../../examples.md#id43)].

A TPS is a type of detector that separates incoming ions according to their charge-to-mass ratio ($q/m$) and initial velocity (hence energy $E_0 = 1/2 m v_0^2$ if we assume non-relativistic dynamics).
TPSs are often used in laser-driven ion acceleration experiments, where different ion species are accelerated at once. To mimic this, we initialize a point-like source of 3 different ion species with different $q/m$ and $E_0$ (i.e. all ions have the same initial position, representative of a pinhole).

The ions propagate along $z$ through 4 subsequent regions:

> - a vacuum region, the distance between the pinhole and the TPS (0.1 m)
> - a region of constant electric field along $x$, (0.19 m, 1e5 V/m)
> - a region of constant magnetic field along $x$, (0.872 T, 0.12 m)
> - a vacuum region, the distance between the TPS and the screen of the detector (0.2 m)

The initial particle velocity $v_0$ is sampled from a uniform distribution in the range $[v_{min}, v_{max}]$ where $v_{min} = \sqrt{E_{max}/m}$, $v_{max} = \sqrt{2E_{max}/m}$, and $E_{max}$ is an input parameter for each species. We assume zero transverse momentum.

The ions are assumed to be test particles embedded in prescribed external fields, meaning that we neglect the self-field due to the ions’ motion and the ions do not interact with each other.

The detector is modeled using a `BoundaryScrapingDiagnostic` at the upper $z$ boundary of the domain, which stores the attributes of the particles when they exit the simulation box from the corresponding edge. Note that the transverse box size is large enough such that all particles exit the domain from the upper $z$ side.

## Run

The PICMI input file is not available for this example yet.

For [MPI-parallel](https://www.mpi-forum.org) runs, prefix these lines with `mpiexec -n 4 ...` or `srun -n 4 ...`, depending on the system.

```none
##############
#### CONSTANTS
##############
my_constants.MeV = 1e6*q_e

# distance between pinhole and electric field
my_constants.d1 = 0.1 # m
# length of the electric field region
my_constants.d2 = 0.19 # m
# length of the magnetic field region
my_constants.d3 = 0.12 # m
# distance between the magnetic field and the screen
my_constants.d4 = 0.2 # m

# constant fields in the TPS
my_constants.E0 = 1e5 # V/m
my_constants.B0 = 0.872 # T

# transverse domain
my_constants.xmin = -0.4 # m
my_constants.xmax =  0.4 # m
my_constants.ymin = -0.4 # m
my_constants.ymax =  0.4 # m

# longitudinal domain
my_constants.zmin= -1e-3 # m
my_constants.zmax = d1+d2+d3+d4

# each macroparticle corresponds to 1 real particle
my_constants.N_real_particles = 1e3
my_constants.N_macro_particles = 1e3

# maximum energy of the different species
# we assume that all the species have a
# uniform energy distribution in [0.5*Emax,Emax]
my_constants.Emax_hydrogen1_1 = 40*MeV
my_constants.Emax_carbon12_6 = 20*MeV
my_constants.Emax_carbon12_4 = 20*MeV

# velocity of a very slow particle
# used to estimate the simulation time
my_constants.vz = sqrt(2*1*MeV/(12*m_p))
my_constants.max_steps = 400
my_constants.max_time = (-zmin+d1+d2+d3+d4) / vz
my_constants.dt = max_time / max_steps

#############
#### NUMERICS
#############
algo.particle_shape = 1
algo.maxwell_solver = none
algo.particle_pusher = boris
amr.max_level = 0
warpx.verbose = 1

########
#### BOX
########
amr.n_cell = 8 8 8
geometry.dims = 3
geometry.prob_hi = xmax ymax zmax
geometry.prob_lo = xmin ymin zmin

#########
#### TIME
#########
stop_time = max_time
warpx.const_dt = dt

#############
#### BOUNDARY
#############
boundary.particle_hi = absorbing absorbing absorbing
boundary.particle_lo = absorbing absorbing absorbing

##############
#### PARTICLES
##############
particles.species_names = hydrogen1_1 carbon12_6 carbon12_4

hydrogen1_1.charge = q_e
hydrogen1_1.initialize_self_fields = 0
hydrogen1_1.injection_style = gaussian_beam
hydrogen1_1.mass = m_p
hydrogen1_1.momentum_distribution_type = uniform
hydrogen1_1.npart = N_macro_particles
hydrogen1_1.q_tot = N_real_particles*q_e
hydrogen1_1.ux_min = 0
hydrogen1_1.uy_min = 0
hydrogen1_1.uz_min = sqrt(Emax_hydrogen1_1/m_p)/clight
hydrogen1_1.ux_max = 0
hydrogen1_1.uy_max = 0
hydrogen1_1.uz_max = sqrt(2*Emax_hydrogen1_1/m_p)/clight
hydrogen1_1.x_m = 0
hydrogen1_1.x_rms = 0
hydrogen1_1.y_m = 0
hydrogen1_1.y_rms = 0
hydrogen1_1.z_m = 0
hydrogen1_1.z_rms = 0
hydrogen1_1.do_not_gather = 1
hydrogen1_1.do_not_deposit = 1

# carbon12_6 means carbon ions with 12 nucleons, of which 6 protons
carbon12_6.charge = 6*q_e
carbon12_6.initialize_self_fields = 0
carbon12_6.injection_style = gaussian_beam
carbon12_6.mass = 12*m_p
carbon12_6.momentum_distribution_type = uniform
carbon12_6.npart = N_macro_particles
carbon12_6.q_tot = N_real_particles*6*q_e
carbon12_6.ux_min = 0
carbon12_6.uy_min = 0
carbon12_6.uz_min = sqrt(Emax_carbon12_6/(12*m_p))/clight
carbon12_6.ux_max = 0
carbon12_6.uy_max = 0
carbon12_6.uz_max = sqrt(2*Emax_carbon12_6/(12*m_p))/clight
carbon12_6.x_m = 0
carbon12_6.x_rms = 0
carbon12_6.y_m = 0
carbon12_6.y_rms = 0
carbon12_6.z_m = 0
carbon12_6.z_rms = 0
carbon12_6.do_not_gather = 1
carbon12_6.do_not_deposit = 1

carbon12_4.charge = 4*q_e
carbon12_4.initialize_self_fields = 0
carbon12_4.injection_style = gaussian_beam
carbon12_4.mass = 12*m_p
carbon12_4.momentum_distribution_type = uniform
carbon12_4.npart = N_macro_particles
carbon12_4.q_tot = N_real_particles*4*q_e
carbon12_4.ux_min = 0
carbon12_4.uy_min = 0
carbon12_4.uz_min = sqrt(Emax_carbon12_4/(12*m_p))/clight
carbon12_4.ux_max = 0
carbon12_4.uy_max = 0
carbon12_4.uz_max = sqrt(2*Emax_carbon12_4/(12*m_p))/clight
carbon12_4.x_m = 0
carbon12_4.x_rms = 0
carbon12_4.y_m = 0
carbon12_4.y_rms = 0
carbon12_4.z_m = 0
carbon12_4.z_rms = 0
carbon12_4.do_not_gather = 1
carbon12_4.do_not_deposit = 1

###########
#### FIELDS
###########
particles.E_ext_particle_init_style = parse_E_ext_particle_function
particles.Ex_external_particle_function(x,y,z,t) = "E0*(z>d1)*(z<(d1+d2))"
particles.Ey_external_particle_function(x,y,z,t) = 0
particles.Ez_external_particle_function(x,y,z,t) = 0

particles.B_ext_particle_init_style = parse_B_ext_particle_function
particles.Bx_external_particle_function(x,y,z,t) = "B0*(z>d1+d2)*(z<(d1+d2+d3))"
particles.By_external_particle_function(x,y,z,t) = 0
particles.Bz_external_particle_function(x,y,z,t) = 0

################
#### DIAGNOSTICS
################
diagnostics.diags_names = diag0 screen diag1

diag0.diag_type = Full
diag0.fields_to_plot = none
diag0.format = openpmd
diag0.intervals = 0:0
diag0.write_species = 1
diag0.species = hydrogen1_1 carbon12_6 carbon12_4
diag0.dump_last_timestep = 0

# diagnostic that collects the particles at the detector's position,
# i.e. when a particle exits the domain from z_max = zhi
# we store it in the screen diagnostic
# we are assuming that most particles will exit the domain at z_max
# which requires a large enough transverse box
screen.diag_type = BoundaryScraping
screen.format = openpmd
screen.intervals = 1
hydrogen1_1.save_particles_at_zhi = 1
carbon12_6.save_particles_at_zhi = 1
carbon12_4.save_particles_at_zhi = 1

diag1.diag_type = Full
diag1.fields_to_plot = rho_hydrogen1_1 rho_carbon12_6 rho_carbon12_4
diag1.format = openpmd
diag1.intervals = 50:50
diag1.write_species = 1
diag1.species = hydrogen1_1 carbon12_6 carbon12_4
diag1.dump_last_timestep = 0
```

## Visualize

This figure below shows the ion trajectories starting from the pinhole (black star), entering the E and B field regions (purple box), up to the detector (gray plane).
The colors represent the different species: protons in blue, C <sub>+4</sub> in red, and C <sub>+6</sub> in green.
The particles are accelerated and deflected through the TPS.

![Ion trajectories through a synthetic TPS.](https://gist.github.com/assets/17280419/3e45e5aa-d1fc-46e3-aa24-d9e0d6a74d1a)

In our simulation, the virtual detector stores all the particle data once entering it (i.e. exiting the simulation box).
The figure below shows the ions colored according to their species (same as above) and shaded according to their initial energy.
The $x$ coordinate represents the electric deflection, while $y$ the magnetic deflection.

![Synthetic TPS screen.](https://gist.github.com/assets/17280419/4dd1adb7-b4ab-481d-bc24-8a7ca51471d9)
```none
#!/usr/bin/env python3

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from openpmd_viewer import OpenPMDTimeSeries
from scipy.constants import c, eV

mpl.use("Agg")
mpl.rcParams.update({"font.size": 18})

MeV = 1e6 * eV

# open the BoundaryScrapingDiagnostic that represents the detector
series = OpenPMDTimeSeries("./diags/screen/particles_at_zhi/")
# open the Full diagnostic at time zero
series0 = OpenPMDTimeSeries("./diags/diag0/")
# we use the data at time 0 to retrieve the initial energy
# of all the particles the boundary

# timesteps and real times
it = series.iterations
time = series.t  # s
N_iterations = len(it)

# list of species names
species = series.avail_species
N_species = len(species)

fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10, 8), dpi=300)

# some stuff for plotting
vmin = 0
vmax = 50
cmap = ["Reds", "Greens", "Blues"]

# loop through the species
for s in range(N_species):
    print(species[s])

    # arrays of positions and energies
    X, Y, E = [], [], []
    for i in range(N_iterations):
        # get particles at detector location
        x, y, z, ids = series.get_particle(
            ["x", "y", "z", "id"], iteration=it[i], species=species[s], plot=False
        )
        # get particles at initialization
        uz0, ids0, m = series0.get_particle(
            ["uz", "id", "mass"],
            iteration=series0.iterations[0],
            species=species[s],
            plot=False,
        )

        indeces = np.where(np.isin(ids0, ids))[0]

        E = np.append(E, 0.5 * m[indeces] * (uz0[indeces] * c) ** 2 / MeV)
        X = np.append(X, x)
        Y = np.append(Y, y)
    print(np.min(E), np.max(E))

    # sort particles according to energy for nicer plot
    sorted_indeces = np.argsort(E)
    ax.scatter(
        X[sorted_indeces],
        Y[sorted_indeces],
        c=E[sorted_indeces],
        vmin=vmin,
        vmax=vmax,
        cmap=cmap[s],
    )
    sorted_indeces = np.argsort(E)
    ax.scatter(
        X[sorted_indeces],
        Y[sorted_indeces],
        c=E[sorted_indeces],
        vmin=vmin,
        vmax=vmax,
        cmap=cmap[s],
    )

# dummy plot just to have a neutral colorbar
im = ax.scatter(np.nan, np.nan, c=np.nan, cmap="Greys_r", vmin=vmin, vmax=vmax)
plt.colorbar(im, label="E [MeV]")
ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")

plt.tight_layout()
fig.savefig("detect.png", dpi=300)
plt.close()
```
