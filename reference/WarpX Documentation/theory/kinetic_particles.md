<a id="theory-kinetic-particles"></a>

# Kinetic Particles

<a id="theory-kinetic-particles-push"></a>

## Particle push

A centered finite-difference discretization of the Newton-Lorentz
equations of motion is given by

<a id="equation-leapfrog-x"></a>
$$
\frac{\mathbf{x}^{i+1}-\mathbf{x}^{i}}{\Delta t} = \mathbf{v}^{i+1/2},

$$

<a id="equation-leapfrog-v"></a>
$$
\frac{\gamma^{i+1/2}\mathbf{v}^{i+1/2}-\gamma^{i-1/2}\mathbf{v}^{i-1/2}}{\Delta t} = \frac{q}{m}\left(\mathbf{E}^{i}+\mathbf{\bar{v}}^{i}\times\mathbf{B}^{i}\right).

$$

In order to close the system, $\bar{\mathbf{v}}^{i}$ must be
expressed as a function of the other quantities. The two implementations that have become the most popular are presented below.

<a id="theory-kinetic-particles-push-boris"></a>

### Boris relativistic velocity rotation

The solution proposed by Boris [[1](#id195)] is given by

<a id="equation-boris-v"></a>
$$
\mathbf{\bar{v}}^{i} = \frac{\gamma^{i+1/2}\mathbf{v}^{i+1/2}+\gamma^{i-1/2}\mathbf{v}^{i-1/2}}{2\bar{\gamma}^{i}}

$$

where $\bar{\gamma}^{i}$ is defined by $\bar{\gamma}^{i} \equiv (\gamma^{i+1/2}+\gamma^{i-1/2} )/2$.

The system ([19](#equation-leapfrog-v), [20](#equation-boris-v)) is solved very
efficiently following Boris’ method, where the electric field push
is decoupled from the magnetic push. Setting $\mathbf{u}=\gamma\mathbf{v}$, the
velocity is updated using the following sequence:

$$
\begin{aligned}
\mathbf{u^{-}}     & = \mathbf{u}^{i-1/2}+\left(q\Delta t/2m\right)\mathbf{E}^{i}
\\
\mathbf{u'}        & = \mathbf{u}^{-}+\mathbf{u}^{-}\times\mathbf{t}
\\
\mathbf{u}^{+}     & = \mathbf{u}^{-}+\mathbf{u'}\times2\mathbf{t}/(1+\mathbf{t}^{2})
\\
\mathbf{u}^{i+1/2} & = \mathbf{u}^{+}+\left(q\Delta t/2m\right)\mathbf{E}^{i}
\end{aligned}
$$

where $\mathbf{t}=\left(q\Delta t/2m\right)\mathbf{B}^{i}/\bar{\gamma}^{i}$ and where
$\bar{\gamma}^{i}$ can be calculated as $\bar{\gamma}^{i}=\sqrt{1+(\mathbf{u}^-/c)^2}$.

The Boris implementation is second-order accurate, time-reversible and fast. Its implementation is very widespread and used in the vast majority of PIC codes.

<a id="theory-kinetic-particles-push-vay"></a>

### Vay Lorentz-invariant formulation

It was shown in Vay [[2](#id121)] that the Boris formulation is
not Lorentz invariant and can lead to significant errors in the treatment
of relativistic dynamics. A Lorentz invariant formulation is obtained
by considering the following velocity average

<a id="equation-new-v"></a>
$$
\mathbf{\bar{v}}^{i} = \frac{\mathbf{v}^{i+1/2}+\mathbf{v}^{i-1/2}}{2}.

$$

This gives a system that is solvable analytically (see Vay [[2](#id121)]
for a detailed derivation), giving the following velocity update:

<a id="equation-pusher-gamma"></a>
$$
\mathbf{u^{*}} = \mathbf{u}^{i-1/2}+\frac{q\Delta t}{m}\left(\mathbf{E}^{i}+\frac{\mathbf{v}^{i-1/2}}{2}\times\mathbf{B}^{i}\right),

$$

<a id="equation-pusher-upr"></a>
$$
\mathbf{u}^{i+1/2} = \frac{\mathbf{u^{*}}+\left(\mathbf{u^{*}}\cdot\mathbf{t}\right)\mathbf{t}+\mathbf{u^{*}}\times\mathbf{t}}{1+\mathbf{t}^{2}},

$$

where

$$
\begin{align}
\mathbf{t} & = \boldsymbol{\tau}/\gamma^{i+1/2},
\\
\boldsymbol{\tau} & = \left(q\Delta t/2m\right)\mathbf{B}^{i},
\\
\gamma^{i+1/2} & = \sqrt{\sigma+\sqrt{\sigma^{2}+\left(\boldsymbol{\tau}^{2}+w^{2}\right)}},
\\
w & = \mathbf{u^{*}}\cdot\boldsymbol{\tau},
\\
\sigma & = \left(\gamma'^{2}-\boldsymbol{\tau}^{2}\right)/2,
\\
\gamma' & = \sqrt{1+(\mathbf{u}^{*}/c)^{2}}.
\end{align}
$$

This Lorentz invariant formulation
is particularly well suited for the modeling of ultra-relativistic
charged particle beams, where the accurate account of the cancellation
of the self-generated electric and magnetic fields is essential, as
shown in Vay [[2](#id121)].
