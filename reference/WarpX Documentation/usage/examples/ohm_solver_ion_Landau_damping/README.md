<a id="examples-ohm-solver-ion-landau-damping"></a>

# Ohm solver: Ion Landau Damping

Landau damping is a well known process in which electrostatic (acoustic) waves are
damped by transferring energy to particles satisfying a resonance condition.
The process can be simulated by seeding a plasma with a specific acoustic mode
(density perturbation) and tracking the strength of the mode as a function of time.

## Run

The same input script can be used for 1d, 2d or 3d simulations and to sweep different
temperature ratios.

### Script `inputs_test_2d_ohm_solver_landau_damping_picmi.py`

```python3
#!/usr/bin/env python3
#
# --- Test script for the kinetic-fluid hybrid model in WarpX wherein ions are
# --- treated as kinetic particles and electrons as an isothermal, inertialess
# --- background fluid. The script simulates ion Landau damping as described
# --- in section 4.5 of Munoz et al. (2018).

import argparse
import os
import sys
import time

import dill
import numpy as np
from mpi4py import MPI as mpi

from pywarpx import callbacks, libwarpx, picmi

constants = picmi.constants

comm = mpi.COMM_WORLD

simulation = picmi.Simulation(warpx_serialize_initial_conditions=True, verbose=0)


class IonLandauDamping(object):
    """This input is based on the ion Landau damping test as described by
    Munoz et al. (2018).
    """

    # Applied field parameters
    B0 = 0.1  # Initial magnetic field strength (T)
    beta = 2.0  # Plasma beta, used to calculate temperature

    # Plasma species parameters
    m_ion = 100.0  # Ion mass (electron masses)
    vA_over_c = 1e-3  # ratio of Alfven speed and the speed of light

    # Spatial domain
    Nz = 256  # number of cells in z direction
    Nx = 4  # number of cells in x (and y) direction for >1 dimensions

    # Temporal domain (if not run as a CI test)
    LT = 40.0  # Simulation temporal length (ion cyclotron periods)

    # Numerical parameters
    NPPC = [8192, 4096, 1024]  # Seed number of particles per cell
    DZ = 1.0 / 6.0  # Cell size (ion skin depths)
    DT = 1e-3  # Time step (ion cyclotron periods)

    # density perturbation strength
    epsilon = 0.03

    # Plasma resistivity - used to dampen the mode excitation
    eta = 1e-7
    # Number of substeps used to update B
    substeps = 20

    def __init__(self, test, dim, m, T_ratio, verbose):
        """Get input parameters for the specific case desired."""
        self.test = test
        self.dim = int(dim)
        self.m = m
        self.T_ratio = T_ratio
        self.verbose = verbose or self.test

        # sanity check
        assert dim > 0 and dim < 4, f"{dim}-dimensions not a valid input"

        # calculate various plasma parameters based on the simulation input
        self.get_plasma_quantities()

        self.dz = self.DZ * self.l_i
        self.Lz = self.Nz * self.dz
        self.Lx = self.Nx * self.dz

        diag_period = 1 / 16.0  # Output interval (ion cyclotron periods)
        self.diag_steps = int(diag_period / self.DT)

        self.total_steps = int(np.ceil(self.LT / self.DT))
        # if this is a test case run for only 100 steps
        if self.test:
            self.total_steps = 100

        self.dt = self.DT / self.w_ci  # self.DT * self.t_ci

        # dump all the current attributes to a dill pickle file
        if comm.rank == 0:
            with open("sim_parameters.dpkl", "wb") as f:
                dill.dump(self, f)

        # print out plasma parameters
        if comm.rank == 0:
            print(
                f"Initializing simulation with input parameters:\n"
                f"\tT = {self.T_plasma * 1e-3:.1f} keV\n"
                f"\tn = {self.n_plasma:.1e} m^-3\n"
                f"\tB0 = {self.B0:.2f} T\n"
                f"\tM/m = {self.m_ion:.0f}\n"
            )
            print(
                f"Plasma parameters:\n"
                f"\tl_i = {self.l_i:.1e} m\n"
                f"\tt_ci = {self.t_ci:.1e} s\n"
                f"\tv_ti = {self.v_ti:.1e} m/s\n"
                f"\tvA = {self.vA:.1e} m/s\n"
            )
            print(
                f"Numerical parameters:\n"
                f"\tdz = {self.dz:.1e} m\n"
                f"\tdt = {self.dt:.1e} s\n"
                f"\tdiag steps = {self.diag_steps:d}\n"
                f"\ttotal steps = {self.total_steps:d}\n"
            )

        self.setup_run()

    def get_plasma_quantities(self):
        """Calculate various plasma parameters based on the simulation input."""
        # Ion mass (kg)
        self.M = self.m_ion * constants.m_e

        # Cyclotron angular frequency (rad/s) and period (s)
        self.w_ci = constants.q_e * abs(self.B0) / self.M
        self.t_ci = 2.0 * np.pi / self.w_ci

        # Alfven speed (m/s): vA = B / sqrt(mu0 * n * (M + m)) = c * omega_ci / w_pi
        self.vA = self.vA_over_c * constants.c
        self.n_plasma = (self.B0 / self.vA) ** 2 / (
            constants.mu0 * (self.M + constants.m_e)
        )

        # Ion plasma frequency (Hz)
        self.w_pi = np.sqrt(constants.q_e**2 * self.n_plasma / (self.M * constants.ep0))

        # Skin depth (m)
        self.l_i = constants.c / self.w_pi

        # Ion thermal velocity (m/s) from beta = 2 * (v_ti / vA)**2
        self.v_ti = np.sqrt(self.beta / 2.0) * self.vA

        # Temperature (eV) from thermal speed: v_ti = sqrt(kT / M)
        self.T_plasma = self.v_ti**2 * self.M / constants.q_e  # eV

        # Larmor radius (m)
        self.rho_i = self.v_ti / self.w_ci

    def setup_run(self):
        """Setup simulation components."""

        #######################################################################
        # Set geometry and boundary conditions                                #
        #######################################################################

        if self.dim == 1:
            grid_object = picmi.Cartesian1DGrid
        elif self.dim == 2:
            grid_object = picmi.Cartesian2DGrid
        else:
            grid_object = picmi.Cartesian3DGrid

        self.grid = grid_object(
            number_of_cells=[self.Nx, self.Nx, self.Nz][-self.dim :],
            warpx_max_grid_size=self.Nz,
            lower_bound=[-self.Lx / 2.0, -self.Lx / 2.0, 0][-self.dim :],
            upper_bound=[self.Lx / 2.0, self.Lx / 2.0, self.Lz][-self.dim :],
            lower_boundary_conditions=["periodic"] * self.dim,
            upper_boundary_conditions=["periodic"] * self.dim,
            warpx_blocking_factor=4,
        )
        simulation.time_step_size = self.dt
        simulation.max_steps = self.total_steps
        simulation.current_deposition_algo = "direct"
        simulation.particle_shape = 1
        simulation.verbose = self.verbose

        #######################################################################
        # Field solver and external field                                     #
        #######################################################################

        self.solver = picmi.HybridPICSolver(
            grid=self.grid,
            gamma=1.0,
            Te=self.T_plasma / self.T_ratio,
            n0=self.n_plasma,
            plasma_resistivity=self.eta,
            substeps=self.substeps,
        )
        simulation.solver = self.solver

        #######################################################################
        # Particle types setup                                                #
        #######################################################################

        k_m = 2.0 * np.pi * self.m / self.Lz
        self.ions = picmi.Species(
            name="ions",
            charge="q_e",
            mass=self.M,
            initial_distribution=picmi.AnalyticDistribution(
                density_expression=f"{self.n_plasma}*(1+{self.epsilon}*cos({k_m}*z))",
                rms_velocity=[self.v_ti] * 3,
            ),
        )
        simulation.add_species(
            self.ions,
            layout=picmi.PseudoRandomLayout(
                grid=self.grid, n_macroparticles_per_cell=self.NPPC[self.dim - 1]
            ),
        )

        #######################################################################
        # Add diagnostics                                                     #
        #######################################################################

        callbacks.installafterstep(self.text_diag)

        if self.test:
            particle_diag = picmi.ParticleDiagnostic(
                name="diag1",
                period=100,
                species=[self.ions],
                data_list=["ux", "uy", "uz", "x", "z", "weighting"],
            )
            simulation.add_diagnostic(particle_diag)
            field_diag = picmi.FieldDiagnostic(
                name="diag1",
                grid=self.grid,
                period=100,
                data_list=["Bx", "By", "Bz", "Ex", "Ey", "Ez", "Jx", "Jy", "Jz"],
            )
            simulation.add_diagnostic(field_diag)

        self.output_file_name = "field_data.txt"
        # install a custom "reduced diagnostic" to save the average field
        callbacks.installafterEsolve(self._record_average_fields)
        try:
            os.mkdir("diags")
        except OSError:
            # diags directory already exists
            pass
        with open(f"diags/{self.output_file_name}", "w") as f:
            f.write("[0]step() [1]time(s) [2]z_coord(m) [3]Ez_lev0-(V/m)\n")

        self.prev_time = time.time()
        self.start_time = self.prev_time
        self.prev_step = 0

        #######################################################################
        # Initialize simulation                                               #
        #######################################################################

        simulation.initialize_inputs()
        simulation.initialize_warpx()

        # get ion particle container wrapper
        self.ions = simulation.particles.get("ions")

    def text_diag(self):
        """Diagnostic function to print out timing data and particle numbers."""
        step = simulation.extension.warpx.getistep(lev=0) - 1

        if step % (self.total_steps // 10) != 0:
            return

        wall_time = time.time() - self.prev_time
        steps = step - self.prev_step
        step_rate = steps / wall_time

        status_dict = {
            "step": step,
            "nplive ions": self.ions.size,
            "wall_time": wall_time,
            "step_rate": step_rate,
            "diag_steps": self.diag_steps,
            "iproc": None,
        }

        diag_string = (
            "Step #{step:6d}; "
            "{nplive ions} core ions; "
            "{wall_time:6.1f} s wall time; "
            "{step_rate:4.2f} steps/s"
        )

        if libwarpx.amr.ParallelDescriptor.MyProc() == 0:
            print(diag_string.format(**status_dict))

        self.prev_time = time.time()
        self.prev_step = step

    def _record_average_fields(self):
        """A custom reduced diagnostic to store the average E&M fields in a
        similar format as the reduced diagnostic so that the same analysis
        script can be used regardless of the simulation dimension.
        """
        step = simulation.extension.warpx.getistep(lev=0) - 1

        if step % self.diag_steps != 0:
            return

        Ez_warpx = simulation.fields.get("Efield_fp", dir="z", level=0)[...]

        if libwarpx.amr.ParallelDescriptor.MyProc() != 0:
            return

        t = step * self.dt
        z_vals = np.linspace(0, self.Lz, self.Nz, endpoint=False)

        if self.dim == 1:
            Ez = Ez_warpx
        elif self.dim == 2:
            Ez = np.mean(Ez_warpx, axis=0)
        else:
            Ez = np.mean(Ez_warpx, axis=(0, 1))

        with open(f"diags/{self.output_file_name}", "a") as f:
            for ii in range(self.Nz):
                f.write(f"{step:05d} {t:.10e} {z_vals[ii]:.10e} {Ez[ii]:+.10e}\n")


##########################
# parse input parameters
##########################

parser = argparse.ArgumentParser()
parser.add_argument(
    "-t",
    "--test",
    help="toggle whether this script is run as a short CI test",
    action="store_true",
)
parser.add_argument(
    "-d", "--dim", help="Simulation dimension", required=False, type=int, default=1
)
parser.add_argument(
    "-m", help="Mode number to excite", required=False, type=int, default=4
)
parser.add_argument(
    "--temp_ratio",
    help="Ratio of ion to electron temperature",
    required=False,
    type=float,
    default=1.0 / 3,
)
parser.add_argument(
    "-v",
    "--verbose",
    help="Verbose output",
    action="store_true",
)
args, left = parser.parse_known_args()
sys.argv = sys.argv[:1] + left

run = IonLandauDamping(
    test=args.test,
    dim=args.dim,
    m=args.m,
    T_ratio=args.temp_ratio,
    verbose=args.verbose,
)
simulation.step()
```

For [MPI-parallel](https://www.mpi-forum.org) runs, prefix these lines with `mpiexec -n 4 ...` or `srun -n 4 ...`, depending on the system.

> ```bash
> python3 inputs_test_2d_ohm_solver_landau_damping_picmi.py -dim {1/2/3} --temp_ratio {value}
> ```

## Analyze

The following script extracts the amplitude of the seeded mode as a function
of time and compares it to the theoretical damping rate.

### Script `analysis.py`

```python3
#!/usr/bin/env python3
#
# --- Analysis script for the hybrid-PIC example of ion Landau damping.

import dill
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from pywarpx import picmi

constants = picmi.constants

matplotlib.rcParams.update({"font.size": 20})

# load simulation parameters
with open("sim_parameters.dpkl", "rb") as f:
    sim = dill.load(f)

# theoretical damping rates were taken from Fig. 14b of Munoz et al.
theoretical_damping_rate = np.array(
    [
        [0.09456706, 0.05113443],
        [0.09864177, 0.05847507],
        [0.10339559, 0.0659153],
        [0.10747029, 0.07359366],
        [0.11290323, 0.08256106],
        [0.11833616, 0.09262114],
        [0.12580645, 0.10541121],
        [0.13327674, 0.11825558],
        [0.14006791, 0.13203098],
        [0.14889643, 0.14600538],
        [0.15772496, 0.16379615],
        [0.16791171, 0.18026693],
        [0.17606112, 0.19650209],
        [0.18828523, 0.21522808],
        [0.19983022, 0.23349062],
        [0.21273345, 0.25209216],
        [0.22835314, 0.27877403],
        [0.24465195, 0.30098317],
        [0.25959253, 0.32186286],
        [0.27657046, 0.34254601],
        [0.29626486, 0.36983567],
        [0.3139219, 0.38984826],
        [0.33157895, 0.40897973],
        [0.35195246, 0.43526107],
        [0.37368421, 0.45662113],
        [0.39745331, 0.47902942],
        [0.44974533, 0.52973074],
        [0.50747029, 0.57743925],
        [0.57334465, 0.63246726],
        [0.64193548, 0.67634255],
    ]
)

expected_gamma = np.interp(
    sim.T_ratio, theoretical_damping_rate[:, 0], theoretical_damping_rate[:, 1]
)

data = np.loadtxt("diags/field_data.txt", skiprows=1)
field_idx_dict = {"z": 2, "Ez": 3}

step = data[:, 0]

num_steps = len(np.unique(step))

# get the spatial resolution
resolution = len(np.where(step == 0)[0]) - 1

# reshape to separate spatial and time coordinates
sim_data = data.reshape((num_steps, resolution + 1, data.shape[1]))

z_grid = sim_data[1, :, field_idx_dict["z"]]
idx = np.argsort(z_grid)[1:]
dz = np.mean(np.diff(z_grid[idx]))
dt = np.mean(np.diff(sim_data[:, 0, 1]))

data = np.zeros((num_steps, resolution))
for i in range(num_steps):
    data[i, :] = sim_data[i, idx, field_idx_dict["Ez"]]

print(f"Data file contains {num_steps} time snapshots.")
print(f"Spatial resolution is {resolution}")

field_kt = np.fft.fft(data[:, :], axis=1)

t_norm = 2.0 * np.pi * sim.m / sim.Lz * sim.v_ti

# Plot the 4th Fourier mode
fig, ax1 = plt.subplots(1, 1, figsize=(10, 5))

t_points = np.arange(num_steps) * dt * t_norm
ax1.plot(
    t_points,
    np.abs(field_kt[:, sim.m] / field_kt[0, sim.m]),
    "r",
    label=f"$T_i/T_e$ = {sim.T_ratio:.2f}",
)

# Plot a line showing the expected damping rate
t_points = t_points[np.where(t_points < 8)]
ax1.plot(t_points, np.exp(-t_points * expected_gamma), "k--", lw=2)

ax1.grid()
ax1.legend()
ax1.set_yscale("log")
ax1.set_ylabel("$|E_z|/E_0$")
ax1.set_xlabel("t $(k_mv_{th,i})$")
ax1.set_xlim(0, 18)

ax1.set_title(f"Ion Landau damping - {sim.dim}d")
plt.tight_layout()
plt.savefig(f"diags/ion_Landau_damping_T_ratio_{sim.T_ratio}.png")
```

The figure below shows a set of such simulations with parameters matching those
described in section 4.5 of Muñoz *et al.* [[1](../../examples.md#id10)].
The straight lines show the theoretical damping rate for the given temperature ratios.

![Ion Landau damping](https://user-images.githubusercontent.com/40245517/230523935-3c8d63bd-ee69-4639-b111-f06dad5587f6.png)
