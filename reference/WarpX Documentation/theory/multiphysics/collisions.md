<a id="multiphysics-collisions"></a>

# Collisions

WarpX includes several different models to capture collisional processes
including collisions between kinetic particles (Coulomb collisions, DSMC,
nuclear fusion) as well as collisions between kinetic particles and a fixed
(i.e. non-evolving) background species (MCC, background stopping).

<a id="multiphysics-collisions-mcc"></a>

## Background Monte Carlo Collisions (MCC)

Several types of collisions between simulation particles and a neutral
background gas are supported including elastic scattering, back scattering,
charge exchange, excitation collisions and impact ionization.

The so-called null collision strategy is used in order to minimize the
computational burden of the MCC module. This strategy is standard in PIC-MCC and
a detailed description can be found elsewhere, for example in Birdsall [[1](#id7)].
In short the maximum collision probability is found over a sensible range of
energies and is used to pre-select the appropriate number of macroparticles for
collision consideration. Only these pre-selected particles are then individually
considered for a collision based on their energy and the cross-sections of all
the different collisional processes included.

The MCC implementation assumes that the background neutral particles are **thermal**,
and are moving at non-relativistic velocities in the lab frame. For each
simulation particle considered for a collision, a velocity vector for a neutral
particle is randomly chosen given the user specified neutral temperature. The
particle velocity is then boosted to the stationary frame of the neutral through
a Galilean transformation. The energy of the collision is calculated using the
particle utility function, `ParticleUtils::getCollisionEnergy()`, as

> $$
> \begin{aligned}
>  E_{coll} &= \sqrt{(\gamma mc^2 + Mc^2)^2 - (mu)^2} - (mc^2 + Mc^2) \\
>           &= \frac{2Mmu^2}{M + m + \sqrt{M^2+m^2+2\gamma mM}}\frac{1}{\gamma + 1}
> \end{aligned}
> $$

where $u$ is the speed of the particle as tracked in WarpX (i.e.
$u = \gamma v$ with $v$ the particle speed), while $m$ and
$M$ are the rest masses of the simulation and background species,
respectively. The Lorentz factor is defined in the usual way,
$\gamma \equiv \sqrt{1 + u^2/c^2}$. Note that if $\gamma\to1$ the above
expression reduces to the classical equation
$E_{coll} = \frac{1}{2}\frac{Mm}{M+m} u^2$. The collision cross-sections
for all scattering processes are evaluated at the energy as calculated above.

Once a particle is selected for a specific collision process, that process determines how the particle is scattered as outlined below.

<a id="multiphysics-collisions-pulseddecay"></a>

## Pulsed Decay

This collision module can be used to have a parent species decay into two product
species with a user-defined decay rate. Mathematically, it solves the following
rate equations on a cell-by-cell basis:

> $$
> \begin{aligned}
>  \frac{dn_1}{dt} &= -\nu(t)n_1, \\
>  \frac{dn_A}{dt} &= +\nu(t)n_1 = \frac{dn_B}{dt}, \\
> \end{aligned}
> $$

where $n_1$ is the parent species density, $n_A$ and $n_B$ are the product species densities,
and $\nu(x,y,z,t)$ is the user-specified decay rate.

This can be used, for example, to represent ionization of a parent species by an externally applied laser pulse.

<a id="multiphysics-collisions-dsmc"></a>

## Direct Simulation Monte Carlo (DSMC)

The algorithm by which binary collisions are treated is outlined below. The
description assumes collisions between different species.

1. Particles from both species are sorted by grid-cells.
2. The order of the particles in each cell is shuffled.
3. Within each cell, particles are paired to form collision partners. Particles
   of the species with fewer members in a given cell is split in half so that
   each particle has exactly one partner of the other species.
4. Each collision pair is considered for a collision using the same logic as in
   the MCC description above.
5. Particles that are chosen for collision are scattered according to the
   selected collision process.

## Scattering processes

### Charge exchange

This process can occur when an ion and neutral (of the same species) collide
and results in the exchange of an electron. The ion velocity is simply replaced
with the neutral velocity and vice-versa.

### Elastic scattering

The `elastic` option uses isotropic scattering, i.e., with a differential
cross section that is independent of angle.
This scattering process as well as the ones below that relate to it, are all
performed in the center-of-momentum (COM) frame. Designating the COM velocity of
the particle as $\vec{u}_c$ and its labframe velocity as $\vec{u}_l$,
the transformation from lab frame to COM frame is done with a general Lorentz
boost (see function `ParticleUtils::doLorentzTransform()`):

> $$
> \begin{bmatrix}
>        \gamma_c c \\
>        u_{cx} \\
>        u_{cy} \\
>        u_{cz}
>    \end{bmatrix}
> = \begin{bmatrix}
>        \gamma & -\gamma\beta_x & -\gamma\beta_y & -\gamma\beta_z \\
>        -\gamma\beta_x & 1+(\gamma-1)\frac{\beta_x^2}{\beta^2} & (\gamma-1)\frac{\beta_x\beta_y}{\beta^2} & (\gamma-1)\frac{\beta_x\beta_z}{\beta^2} \\
>        -\gamma\beta_y & (\gamma-1)\frac{\beta_x\beta_y}{\beta^2} & 1 +(\gamma-1)\frac{\beta_y^2}{\beta^2} & (\gamma-1)\frac{\beta_y\beta_z}{\beta^2} \\
>        -\gamma\beta_z & (\gamma-1)\frac{\beta_x\beta_z}{\beta^2} & (\gamma-1)\frac{\beta_y\beta_z}{\beta^2} & 1+(\gamma-1)\frac{\beta_z^2}{\beta^2} \\
>    \end{bmatrix} \begin{bmatrix}
>        \gamma_l c \\
>        u_{lx} \\
>        u_{ly} \\
>        u_{lz}
>    \end{bmatrix}

> $$

where $\gamma$ is the Lorentz factor of the relative speed between the lab frame and the COM frame, $\beta_i = v^{COM}_i/c$ is the i’th component of the relative velocity between the lab frame and the COM frame with

> $$
> \vec{v}^{COM} = \frac{m \vec{u_c}}{\gamma_u m + M}
> $$

The particle velocity in the COM frame is then isotropically scattered using the function `ParticleUtils::RandomizeVelocity()`. After the direction of the velocity vector has been appropriately changed, it is transformed back to the lab frame with the reversed Lorentz transform as was done above followed by the reverse Galilean transformation using the starting neutral velocity.

### Back scattering

The process is the same as for elastic scattering above except the scattering angle is fixed at $\pi$, meaning the particle velocity in the COM frame is updated to $-\vec{u}_c$.

### Excitation

The process is also the same as for elastic scattering except the excitation energy cost is subtracted from the particle energy. This is done by reducing the velocity before a scattering angle is chosen.

### Forward scattering

This process operates in two ways:

1. If an excitation energy cost is provided, the energy cost is subtracted from the particle energy and no scattering is performed.
2. If an excitation energy cost is not provided, the particle is not scattered and the velocity is unchanged (corresponding to a scattering angle of $0$ in the elastic scattering process above).

See Janssen *et al.* [[2](#id8)] for a recommended use of this process.

## Benchmarks

See the [MCC example](../../usage/examples/capacitive_discharge/README.md#examples-capacitive-discharge) for a benchmark of the MCC
implementation against literature results.

### Particle cooling due to elastic collisions

It is straight forward to determine the energy a projectile loses during an elastic collision with another body, as a function of scattering angle, through energy and momentum conservation.
See for example Lim [[3](#id9)] for a derivation. The result is that given a projectile with mass $m$, a target with mass $M$, a scattering angle $\theta$, and collision energy $E$, the post collision energy of the projectile is given by

> $$
> \begin{aligned}
> E_{final} = E - &[(E + mc^2)\sin^2\theta + Mc^2 - \cos(\theta)\sqrt{M^2c^4 - m^2c^4\sin^2\theta}] \\
> &\times \frac{E(E+2mc^2)}{(E+mc^2+Mc^2)^2 - E(E+2mc^2)\cos^2\theta}
> \end{aligned}
> $$

The impact of incorporating relativistic effects in the MCC routine can be seen in the plots below where high energy collisions are considered with both a classical and relativistic implementation of MCC. It is observed that the classical version of MCC reproduces the classical limit of the above equation but especially for ions, this result differs substantially from the fully relativistic result.

![Comparison of classical and relativistic MCC collision results](https://user-images.githubusercontent.com/40245517/170900079-74e505a5-2790-44f5-ac84-5847eda954e6.png)
