<a id="examples-pierce-diode"></a>

# Pierce Diode at the Child–Langmuir Limit

This example shows how to simulate the physics of a 1D Pierce diode configuration operating
at the Child–Langmuir limit using WarpX. In this setup, an electron beam is injected
into a planar diode gap, consisting of by parallel conducting plates separated by the distance $d$ and powered by a voltage difference $V$ [Fig. 12](#fig-geom).
The injected current density is chosen to match the space-charge-limited current predicted by the Child–Langmuir law Zhang *et al.* [[16](../../examples.md#id253)].
The law predicts the maximum current density that can flow between two parallel plates due to space-charge effects.
This test demonstrates that WarpX correctly reproduces the Child–Langmuir law for a given voltage and gap length.

## Geometry

The figure below schematically illustrates the problem geometry described above.

<a id="fig-geom"></a>
![[fig:geom] Two parallel conducting plates separated by the distance :math:`d` and powered by a voltage difference :math:`V`. Given that the two plates are parallel, here we simulate the problem in 1D with WarpX.](https://gist.githubusercontent.com/oshapoval/aaafd8d131c3e1ed0fefe348bc8db28b/raw/92c4089e1b9eb23ae258f60c386e38e04f9499a2/geometry_pierce_diode.png)

## Сhild–Langmuir Limit

In steady state, the emitted current is limited by the Child–Langmuir law,
which defines the maximum current that can be transported across a planar diode for a given voltage and gap length Zhang *et al.* [[16](../../examples.md#id253)].
It can be shown that, at the Child-Langmuir limit (i.e. when this maximum current is reached), the potential and current density in the gap have the following expression:

<a id="equation-child-langmuir-phi"></a>
$$
\phi(z)=V\Big(\frac{z}{d}\Big)^{4/3},

$$

<a id="equation-child-langmuir-j"></a>
$$
J(z) = \frac{4}{9} \varepsilon_0 \sqrt{\frac{2 |q|}{m}} \frac{|V|^{3/2}}{d^2}.

$$

## Run

This example can be run with the WarpX executable using an input file: `warpx.1d inputs_test_1d_pierce_diode`.
For [MPI-parallel](https://www.mpi-forum.org) runs, prefix these lines with `mpiexec -n 4 ...` or `srun -n 4 ...`, depending on the system.

```none
# --- This input defines 1D Pierce Diode Example:
# --- This test simulates a classic Pierce diode configuration:
# --- two parallel conducting plates separated by a distance
# --- d_plate = 8 cm and powered by a voltage difference
# --- extractor_voltage = -93 kV. The injected current density is
# --- chosen to match the space-charge-limited current predicted by
# --- the Child–Langmuir law, that predicts the maximum current density
# --- that can flow between two parallel plates due to space-charge effects
# --- for a given voltage and gap length.

# constants
my_constants.ion_mass = 39*m_u # [kg] Assuming ion potassium
my_constants.kV = 1000
my_constants.cm = 0.01
my_constants.extractor_voltage  = -93.*kV
my_constants.d_plate = 8.*cm
my_constants.nz = 128
my_constants.dz = d_plate/nz
my_constants.vzfinal = sqrt(2.*abs(extractor_voltage)*q_e/ion_mass)
my_constants.dt = 0.4*(dz/vzfinal)
my_constants.ep0 = 8.8541878188e-12
# current density at Child-Langmuir limit
my_constants.J_CL = ( (4 / 9)* ep0 * sqrt(2 * abs(q_e) / ion_mass) * abs(extractor_voltage) ** (3 / 2) / d_plate**2 )
# algo
warpx.do_electrostatic = labframe
algo.particle_shape = 1

# amr
amr.n_cell = nz
amr.max_level = 0

# maxamx step
max_step = 5000

# timestep
warpx.const_dt = dt

# number of procs
warpx.numprocs = 2   # 2 MPI ranks

# geometry
geometry.dims = 1
geometry.prob_lo = 0.0
geometry.prob_hi = d_plate

# boundary
boundary.field_lo = pec
boundary.field_hi = pec
# set fixed potential at the plates (V difference)
boundary.potential_lo_z = 0.0      # cathode at x = 0 V
boundary.potential_hi_z = extractor_voltage  # anode at x = 1 kV
boundary.particle_lo = absorbing
boundary.particle_hi = absorbing

# ions
particles.species_names = ions
ions.charge = q_e
ions.mass = ion_mass
ions.do_continuous_injection = 1
ions.injection_style = "NFluxPerCell"
ions.num_particles_per_cell = 15
ions.surface_flux_pos = 0
ions.flux_normal_axis = z
ions.flux_direction = 1
ions.flux_profile = constant
ions.flux = J_CL/q_e # corresponds to Child-Langmuir limit
ions.momentum_distribution_type = constant

# diagnostics
diagnostics.diags_names = diag1
diag1.intervals = 5000
diag1.diag_type = Full
diag1.format=openpmd
diag1.fields_to_plot = Ez rho jz phi
```

## Visualize

The figure below shows the results of the simulation (orange curves), which agrees well with the analytical Child–Langmuir law (black curves) ([1](#equation-child-langmuir-phi), [2](#equation-child-langmuir-j)).

![Results of the WarpX Pierce Diode simulation.](https://gist.githubusercontent.com/oshapoval/aaafd8d131c3e1ed0fefe348bc8db28b/raw/fc76b371d323dbca4e1c43b45055405ff1fc6de4/Pierce_Diode.png)

This figure was obtained with the script below, which can be run with `python3 plot_sim.py`.

```none
import matplotlib.pyplot as plt
import numpy as np
from openpmd_viewer import OpenPMDTimeSeries
from scipy.constants import c, e, epsilon_0, m_u

ts = OpenPMDTimeSeries("./diags/diag1/")

kV = 1000
cm = 0.01
extractor_voltage = -93.0 * kV
d_plate = 8.0 * cm
ion_mass = 39 * m_u

it = ts.iterations[-1]
phi, meta = ts.get_field("phi", iteration=it, plot=False)
ez, _ = ts.get_field("E", "z", iteration=ts.iterations[-1], plot=False)
rho, _ = ts.get_field("rho", iteration=it, plot=False)
jz, _ = ts.get_field("j", "z", iteration=it, plot=False)
z, uz = ts.get_particle(["z", "uz"], iteration=it)
time_cur = ts.current_t

# Calculate theoretical Child-Langmuir limit for a given voltage
jz_CL_theory = (
    (4 / 9)
    * epsilon_0
    * np.sqrt(2 * abs(e) / ion_mass)
    * abs(extractor_voltage) ** (3 / 2)
    / d_plate**2
)
phi_CL_theory = extractor_voltage * (meta.z / d_plate) ** (4 / 3)
ez_CL_theory = -(4 / (3 * d_plate)) * extractor_voltage * (meta.z / d_plate) ** (1 / 3)
rho_CL_theory = (
    epsilon_0
    * (4 / (3 * d_plate) ** 2)
    * extractor_voltage
    * (meta.z / d_plate) ** (-2 / 3)
)
uz_CL_theory = -jz_CL_theory / rho_CL_theory / c

color = "orange"
title = r"$\Gamma_{ions}=38.79 \approx \Gamma_{CL}$"

fig, axs = plt.subplots(2, 2, figsize=(10, 8))
fig.suptitle(f"{title}\n time = {np.round(time_cur / 1e-6, 5)} $\\mu s$")
axs[0, 0].scatter(z, uz, color=color, label=title, s=0.2)
axs[0, 0].plot(meta.z, uz_CL_theory, ls=":", color="black")
axs[0, 0].set_title(r"$u_z$")
axs[0, 0].set_xlabel("z, mm")

axs[0, 1].plot(meta.z, ez, color=color, label="WarpX ")
axs[0, 1].set_title(r"$E_z$")
axs[0, 1].set_xlabel("z, mm")
axs[0, 1].plot(
    meta.z, ez_CL_theory, label="Child-Langmuir limit (theory)", ls=":", color="black"
)
axs[0, 1].legend()

axs[1, 0].plot(meta.z, jz, color=color)
axs[1, 0].set_title(r"$J_z$")
axs[1, 0].axhline(y=jz_CL_theory, color="black", linestyle="-")
axs[1, 0].set_xlabel("z, mm")

axs[1, 1].plot(meta.z, phi, color=color)
axs[1, 1].set_title(r"$\phi$")
axs[1, 1].set_xlabel("z, mm")
axs[1, 1].plot(meta.z, phi_CL_theory, label="theory", ls=":", color="black")

plt.tight_layout()

plt.savefig("Pierce_Diode.png")
```
