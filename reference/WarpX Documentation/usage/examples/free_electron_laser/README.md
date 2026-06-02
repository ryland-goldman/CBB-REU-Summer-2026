<a id="examples-free-electron-laser"></a>

# Free-electron laser

This example shows how to simulate the physics of a free-electron laser (FEL) using WarpX.
In this example, a relativistic electron beam is sent through an undulator (represented by an external,
oscillating magnetic field). The radiation emitted by the beam grows exponentially
as the beam travels through the undulator, due to the Free-Electron-Laser instability.

The parameters of the simulation are taken from section 5.1 of Fallahi [[10](../../examples.md#id39)].

The simulation is performed in 1D, and uses the boosted-frame technique as described in
Fawley and Vay [[11](../../examples.md#id40)] and Fawley and Vay [[12](../../examples.md#id41)] to reduce the computational cost (the Lorentz frame of the simulation is moving at the average speed of the beam in the undulator).
Even though the simulation is run in this boosted frame, the results are reconstructed in the
laboratory frame, using WarpX’s `BackTransformed` diagnostic.

The effect of space-charge is intentionally turned off in this example, as it may not be properly modeled in 1D.
This is achieved by initializing two species of opposite charge (electrons and positrons) to
represent the physical electron beam, as discussed in Fawley and Vay [[12](../../examples.md#id41)].

## Run

This example can be run with the WarpX executable using an input file: `warpx.1d inputs_test_1d_fel`. For [MPI-parallel](https://www.mpi-forum.org) runs, prefix these lines with `mpiexec -n 4 ...` or `srun -n 4 ...`, depending on the system.

```none
my_constants.gamma_bunch=100.6
my_constants.Bu = 0.5
my_constants.lambda_u = 3e-2
my_constants.k_u= 2*pi/lambda_u
my_constants.K = q_e*Bu/(m_e*clight*k_u) # Undulator parameter

warpx.gamma_boost = gamma_bunch/sqrt(1+K*K/2) # Lorentz factor of the ponderomotive frame
warpx.boost_direction = z
algo.maxwell_solver = yee
algo.particle_shape = 2
algo.particle_pusher = vay

# geometry
geometry.dims = 1
geometry.prob_hi = 0
geometry.prob_lo = -192e-6

amr.max_grid_size = 1024
amr.max_level = 0
amr.n_cell = 1024

# boundary
boundary.field_hi = absorbing_silver_mueller
boundary.field_lo = absorbing_silver_mueller
boundary.particle_hi = absorbing
boundary.particle_lo = absorbing

# diagnostics
diagnostics.diags_names = diag_labframe diag_boostedframe

# Diagnostic that show quantities in the frame
# of the simulation (boosted-frame)
diag_boostedframe.diag_type = Full
diag_boostedframe.format = openpmd
diag_boostedframe.intervals = 100

# Diagnostic that show quantities
# reconstructed in the lab frame
diag_labframe.diag_type = BackTransformed
diag_labframe.num_snapshots_lab = 25
diag_labframe.dz_snapshots_lab = 0.1
diag_labframe.format = openpmd
diag_labframe.buffer_size = 64

# Run the simulation long enough for
# all backtransformed diagnostic to be complete
warpx.compute_max_step_from_btd = 1

particles.species_names = electrons positrons
particles.rigid_injected_species= electrons positrons

electrons.charge = -q_e
electrons.injection_style = nuniformpercell
electrons.mass = m_e
electrons.momentum_distribution_type = constant
electrons.num_particles_per_cell_each_dim = 8
electrons.profile = constant
electrons.density = 2.7e19/2
electrons.ux = 0.0
electrons.uy = 0.0
electrons.uz = gamma_bunch
electrons.zmax = -25e-6
electrons.zmin = -125e-6
electrons.zinject_plane=0.0
electrons.rigid_advance=0

positrons.charge = q_e
positrons.injection_style = nuniformpercell
positrons.mass = m_e
positrons.momentum_distribution_type = constant
positrons.num_particles_per_cell_each_dim = 8
positrons.profile = constant
positrons.density = 2.7e19/2
positrons.ux = 0.0
positrons.uy = 0.0
positrons.uz = gamma_bunch
positrons.zmax = -25e-6
positrons.zmin = -125e-6
positrons.zinject_plane=0.0
positrons.rigid_advance=0

warpx.do_moving_window = 1
warpx.moving_window_dir = z
warpx.moving_window_v = sqrt(1-(1+K*K/2)/(gamma_bunch*gamma_bunch))

# Undulator field
particles.B_ext_particle_init_style = parse_B_ext_particle_function
particles.Bx_external_particle_function(x,y,z,t) = 0
particles.By_external_particle_function(x,y,z,t) = if( z>0, Bu*cos(k_u*z), 0 )
particles.Bz_external_particle_function(x,y,z,t) =0.0

warpx.cfl = 0.99
```

## Visualize

The figure below shows the results of the simulation. The left panel shows the exponential growth of the radiation along the undulator (note that the vertical axis is plotted in log scale). The right panel shows a snapshot of the simulation,
1.6 m into the undulator. Microbunching of the beam is visible in the electron density (blue). One can also see the
emitted FEL radiation (red) slipping ahead of the beam.

![Results of the WarpX FEL simulation.](https://gist.githubusercontent.com/RemiLehe/871a1e24c69e353c5dbb4625cd636cd1/raw/7f4e3da7e0001cff6c592190fee8622580bbe37a/FEL.png)

This figure was obtained with the script below, which can be run with `python3 plot_sim.py`.

```none
import matplotlib.pyplot as plt
from openpmd_viewer import OpenPMDTimeSeries

ts = OpenPMDTimeSeries("./diags/diag_labframe/")


def extract_peak_E(iteration):
    """
    Extract peak electric field and its position
    """
    Ex, info = ts.get_field("E", "x", iteration=iteration)
    Ex_max = abs(Ex).max()
    z_max = info.z[abs(Ex).argmax()]
    return z_max, Ex_max


# Loop through the lab-frame snapshots and extract the peak electric field
z_max, Ex_max = ts.iterate(extract_peak_E)

# Create a figure
plt.figure(figsize=(8, 4))

# Plot of the E field growth
plt.subplot(121)  # Span all rows in the first column
plt.semilogy(z_max, Ex_max)
plt.ylim(2e7, 2e9)
plt.xlabel("z (m)")
plt.ylabel("Peak $E_x$ (V/m)")
plt.title("Growth of the radiation field\n along the undulator")

# Plots of snapshot
iteration = 16
plt.subplot(122)  # Upper right panel


plt.ylabel("$E_x$ (V/m)")
plt.xlabel("")
ts.get_particle(["z"], iteration=iteration, nbins=300, species="electrons", plot=True)
plt.title("")
plt.ylim(0, 30e12)
plt.ylabel("Electron density (a. u.)", color="b")
plt.twinx()
Ex, info = ts.get_field("E", "x", iteration=iteration, plot=True)
plt.ylabel("$E_x$ (V/m)", color="r")
plt.plot(info.z, Ex, color="r")
plt.ylim(-0.6e9, 0.4e9)
plt.xlabel("z (m)")
plt.title("Snapshot 1.6 m into the undulator")

plt.tight_layout()

plt.savefig("FEL.png")
```
