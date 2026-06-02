<a id="theory-em-pic"></a>

# Electromagnetic PIC

In the *electromagnetic Particle-In-Cell method* [[2](explicit_em_pic.md#id257), [3](explicit_em_pic.md#id89)],
the fields are updated using Maxwell’s equations:

<a id="equation-faraday-1"></a>
$$
\frac{\partial \boldsymbol{B}}{\partial t} = -\nabla\times \boldsymbol{E}

$$

<a id="equation-ampere-1"></a>
$$
\frac{1}{c^2}\frac{\partial \boldsymbol{E}}{\partial t} = \nabla\times \boldsymbol{B}-\mu_0 \boldsymbol{j}

$$

where $\boldsymbol{E}$ and $\boldsymbol{B}$ are the electric and magnetic field
components, and $\boldsymbol{j}$ is the current density.

Because the electromagnetic PIC method retains the full Maxwell equations,
this method can capture the **physics of the electromagnetic waves**,
including their propagation and self-consistent interaction with particles.

The electromagnetic PIC method can be run either with an explicit or implicit time integration scheme:

> - In the **explicit integration scheme**, the particles and fields are updated sequentially at each time step
>   (see [Explicit electromagnetic PIC](explicit_em_pic.md#theory-explicit-em-pic)). This integration scheme is simple, but requires a small enough time step
>   size $\Delta t$ to ensure the stability of the simulation (e.g., CFL condition $c\Delta t \lessapprox \Delta x$,
>   need to resolve the plasma frequency $\omega_p \Delta t \leq 2$ [[2](explicit_em_pic.md#id257), [3](explicit_em_pic.md#id89)]).
> - In the **implicit integration scheme**, the particles and fields are updated simultaneously at each time step, using
>   an iterative solver (see [Implicit electromagnetic PIC](implicit_em_pic.md#theory-implicit-em-pic)). While this integration scheme is more complex, it can use
>   larger time step sizes $\Delta t$ and still retain the stability of the simulation. In addition, the implicit
>   integration scheme is exactly energy conserving.

For more details, see the sections below:

* [Explicit electromagnetic PIC](explicit_em_pic.md)
* [Implicit electromagnetic PIC](implicit_em_pic.md)
