<a id="examples-pwfa-boost"></a>

# In-Depth: PWFA

As described in the [Overview](../theory/intro.md), one of the key applications of the WarpX exascale computing platform is in modelling future, compact and economic plasma-based accelerators.
In this section we describe the simulation setup of a realistic *electron beam driven plasma wakefield accelerator* (PWFA) configuration.
For illustration purposes the setup can be explored with **WarpX** using the example input file [`PWFA`](../../../Examples/Physics_applications/plasma_acceleration/inputs_test_2d_plasma_acceleration_boosted).

The simulation setup consists of 4 particle species: drive beam (driver), witness beam (beam), plasma electrons (plasma_e), and plasma ions (plasma_p).
The species physical parameters are summarized in the following table.

| Species   | Parameters                                                             |
|-----------|------------------------------------------------------------------------|
| driver    | $\gamma$ = 48923; N = 2x10^8; $\sigma_z$ = 4.0 um; $\sigma_x$ = 2.0 um |
| beam      | $\gamma$ = 48923; N = 6x10^5; $\sigma_z$ = 1.0 mm; $\sigma_x$ = 0.5 um |
| plasma_e  | n = 1x10^23 m^-3; w = 70 um; lr = 8 mm; L = 200 mm                     |
| plasma_p  | n = 1x10^23 m^-3; w = 70 um; lr = 8 mm; L = 200 mm                     |

Where $\gamma$ is the beam relativistic Lorentz factor, N is the number of particles, and $\sigma_x$, $\sigma_y$, $\sigma_z$ are the beam widths (root-mean-squares of particle positions) in the transverse (x,y) and longitudinal directions.

The plasma, of total length L, has a density profile that consists of a lr long linear up-ramp, ranging from 0 to peak value n, is uniform within a transverse width of w and after the up-ramp.

With this configuration the driver excites a nonlinear plasma wake and drives the bubble depleted of plasma electrons where the beam accelerates, as can be seen in Fig. [[fig:PWFA]](#fig:PWFA).

![Driver, beam and plasma electron distribution in 2D PWFA simulation](usage/PWFA.png)

Listed below are the key arguments and best-practices relevant for choosing the pwfa simulation parameters used in the example.

## 2D Geometry

> 2D cartesian (with longitudinal direction z and transverse x) geometry simulations can give valuable physical and numerical insight into the simulation requirements and evolution.
> At the same time it is much less time consuming than the full 3D cartesian or cylindrical geometries.

## Finite Difference Time Domain

> For standard plasma wakefield configurations, it is possible to model the physics correctly using the [Particle-In-Cell (PIC)](../theory/models_algorithms/explicit_em_pic.md#theory-explicit-em-pic) Finite Difference Time Domain (FDTD) algorithms.
> If the simulation contains localised extremely high intensity fields, however, numerical instabilities might arise, such as the numerical Cherenkov instability ([Moving window and optimal Lorentz boosted frame](../theory/boosted_frame.md)).
> In that case, it is recommended to use the Pseudo Spectral Analytical Time Domain (PSATD) or the Pseudo-Spectral Time-Domain (PSTD) algorithms.
> In the example we are describing, it is sufficient to use FDTD.

## Cole-Karkkainen solver with Cowan coefficients

> There are two FDTD Maxwell field solvers that compute the field push implemented in WarpX: the Yee and Cole-Karkkainen solver with Cowan coefficients (CKC) solvers.
> The later includes a modification that allows the numerical dispersion of light in vacuum to be exact, and that is why we choose CKC for the example.

## Lorentz boosted frame

> WarpX simulations can be done in the laboratory or [Lorentz-boosted](https://warpx.readthedocs.io/en/latest/theory/boosted_frame.html) frames.
> In the laboratory frame, there is typically no need to model the plasma ions species, since they are mainly stationary during the short time scales associated with the motion of plasma electrons.
> In the boosted frame, that argument is no longer valid, as ions have relativistic velocities.
> The boosted frame still results in a substantial reduction to the simulation computational cost.

#### NOTE
Even if the simulations uses the boosted frame, most of its input file parameters are defined in respect to the laboratory frame.

We recommend that you design your numerical setup so that the width of the box is not significantly narrower than the distance from 0 to its right edge (done, for example, by setting the right edge equal to 0).

## Moving window

> To avoid having to simulate the whole 0.2 mm of plasma with the high resolution that is required to model the beam and plasma interaction correctly, we use the moving window.
> In this way we define a simulation box (grid) with a fixed size that travels at the speed-of-light ($c$), i.e. follows the beam.

> #### NOTE
> When using moving window the option of continuous injection needs to be active for all particles initialized outside of the simulation box.

## Resolution

> Longitudinal and transverse resolutions (i.e. number and dimensions of the PIC grid cells) should be chosen to accurately describe the physical processes taking place in the simulation.
> Convergence scans, where resolution in both directions is gradually increased, should be used to determine the optimal configuration.
> Multiple cells per beam length and width are recommended (our illustrative example resolution is coarse).

> #### NOTE
> To avoid spurious effects, in the boosted frame, we consider the constraint that the transverse cell size should be larger than the transverse one.
> Translating this condition to the cell transverse ($d_{x}$) and longitudinal dimensions ($d_{z}$) in the laboratory frame leads to: $d_{x} > (d_{z} (1+\beta_{b}) \gamma_{b})$, where $\beta_{b}$ is the boosted frame velocity in units of $c$.

## Time step

> The time step ($dt$) is used to iterated over the main PIC loop and is computed by WarpX differently depending on the Maxwell field FDTD solvers used:

> * **For Yee** is equal to the CFL parameter chosen in the input file ([Inputs: Parameter List](parameters.md)) times the Courant–Friedrichs–Lewy condition (CFL) that follows the analytical expression in [Explicit electromagnetic PIC](../theory/models_algorithms/explicit_em_pic.md#theory-explicit-em-pic)
> * **For CKC** is equal to CFL times the minimum between the boosted frame cell dimensions

> where CFL is chosen to be below unity and set an optimal trade-off between making the simulation faster and avoiding NCI and other spurious effects.

## Duration of the simulation

> To determine the total number of time steps of the simulation, we could either set the <zmax_plasma_to_compute_max_step> parameter to the end of the plasma ($z_{\textrm{end}}$), or compute it using:

> * boosted frame edge of the simulation box, $\textrm{corner} = l_{e}/ ((1-\beta_{b}) \gamma_{b})$
> * time of interaction in the boosted frame, $T = \frac{z_{\textrm{end}}/\gamma_{b}-\textrm{corner}}{c (1+\beta_{b})}$
> * total number of iterations, $i_{\textrm{max}} = T/dt$

> where $l_{e}$ is the position of the left edge of the simulation box (in respect to propagation direction).

## Plotfiles and snapshots

> WarpX allows the data to be stored in different formats, such as plotfiles (following the [yt guidelines](https://yt-project.org/doc/index.html)), hdf5 and openPMD (following its [standard](https://github.com/openPMD)).
> In the example, we are dumping plotfiles with boosted frame information on the simulation particles and fields.
> We are also requesting back transformed diagnostics that transform that information back to the laboratory frame.
> The diagnostics results are analysed and stored in snapshots at each time step and so it is best to make sure that the run does not end before filling the final snapshot.

## Maximum grid size and blocking factor

> These parameters are carfully chosen to improve the code parallelization, load-balancing and performance ([Inputs: Parameter List](parameters.md)) for each numerical configuration.
> They define the smallest and largest number of cells that can be contained in each simulation box and are carefully defined in the [AMReX](https://amrex-codes.github.io/amrex/docs_html/GridCreation.html?highlight=blocking_factor) documentation.
