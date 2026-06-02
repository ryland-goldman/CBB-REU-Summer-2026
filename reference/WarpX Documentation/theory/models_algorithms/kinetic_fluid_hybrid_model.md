<a id="theory-kinetic-fluid-hybrid-model"></a>

# Ampere’s law coupled with Ohm’s law (a.k.a. “hybrid PIC”)

Many problems in plasma physics fall in a class where both electron kinetics and electromagnetic waves do not
play a critical role in the solution. Examples of such situations include the
study of collisionless magnetic reconnection and instabilities driven by ion
temperature anisotropy, to mention only two. For these kinds of problems the
computational cost of resolving the electron dynamics can be avoided by modeling
the electrons as a neutralizing fluid rather than kinetic particles. By further
using Ohm’s law to compute the electric field rather than evolving it with the
Maxwell-Faraday equation, light waves can be stepped over. The simulation resolution
can then be set by the ion time and length scales (commonly the ion cyclotron
period $1/\Omega_i$ and ion skin depth $l_i$, respectively), which
can reduce the total simulation time drastically compared to a simulation that
has to resolve the electron Debye length and CFL-condition based on the speed
of light.

Many authors have described variations of the kinetic ion & fluid electron model,
generally referred to as particle-fluid hybrid or just hybrid-PIC models. The
implementation in WarpX is described in detail in Groenewald *et al.* [[1](#id28)].
The “Model derivation” section below gives a detailed description of the model
that follows mostly from the above reference, but succinctly, the model
entails the following:

The magnetic field is advanced in time using Faraday’s law,

> $$
> \frac{\partial\vec{B}}{\partial t} = -\nabla\times\vec{E},
> $$

where the electric field is calculated from Ohm’s law which involves the currents,
the magnetic field, and the electron pressure (for which an additional closure is required,
see [here](#theory-hybrid-model-elec-temp)),

> $$
> \vec{E} = -\frac{1}{en_e}\left( \vec{J}_e\times\vec{B} + \nabla P_e \right)+\eta\vec{J}-\eta_h \nabla^2\vec{J}.
> $$

The electron current is in turn obtained by subtracting the ion current (obtained from
kinetic ion macro-particles) from the total current (obtained from Ampere’s law):

> $$
> \vec{J}_e = \vec{J} - \sum_{s\neq e}\vec{J}_s - \vec{J}_{ext}
> $$

where

> $$
> \mu_0\vec{J} = \vec{\nabla}\times\vec{B}.
> $$

## Algorithm details

#### NOTE
Various verification tests of the hybrid model implementation can be found in
the [examples section](../../usage/examples.md#examples-hybrid-model).

The kinetic-fluid hybrid extension mostly uses the same routines as the standard electromagnetic
PIC algorithm with the only exception that the E-field is calculated from Ohm’s law
rather than it being updated from the full Maxwell-Ampere equation. The E-field update occurs
after particle pushing and deposition (charge and current density) has been completed. Therefore, based
on the usual time-staggering in the PIC algorithm, when the E-field is updated
at timestep $t=t_n$, the quantities $\rho^n$, $\rho^{n+1}$, $\vec{J}_i^{n-1/2}$
and  $\vec{J}_i^{n+1/2}$ are all known.

### Field update

The field update is done in three steps as described below.

#### First half step

Firstly the E-field at $t=t_n$ is calculated for which the current density needs to
be interpolated to the correct time, using $\vec{J}_i^n = 1/2(\vec{J}_i^{n-1/2}+ \vec{J}_i^{n+1/2})$.
The electron pressure is simply calculated using $\rho^n$ and the B-field is also already
known at the correct time since it was calculated for $t=t_n$ at the end of the last step.
Once $\vec{E}^n$ is calculated, it is used to push $\vec{B}^n$ forward in time
(using the Maxwell-Faraday equation) to $\vec{B}^{n+1/2}$.

#### Second half step

Next, the E-field is recalculated to get $\vec{E}^{n+1/2}$. This is done
using the known fields $\vec{B}^{n+1/2}$, $\vec{J}_i^{n+1/2}$ and
interpolated charge density $\rho^{n+1/2}=1/2(\rho^n+\rho^{n+1})$ (which is
also used to calculate the electron pressure). Similarly as before, the B-field
is then pushed forward to get $\vec{B}^{n+1}$ using the newly calculated
$\vec{E}^{n+1/2}$ field.

#### Extrapolation step

Obtaining the E-field at timestep $t=t_{n+1}$ is a well documented issue for
the hybrid model. Currently the approach in WarpX is to simply extrapolate
$\vec{J}_i$ forward in time, using

> $$
> \vec{J}_i^{n+1} = \frac{3}{2}\vec{J}_i^{n+1/2} - \frac{1}{2}\vec{J}_i^{n-1/2}.
> $$

With this extrapolation all fields required to calculate $\vec{E}^{n+1}$
are known and the simulation can proceed.

### Sub-stepping

It is also well known that hybrid PIC routines require the B-field to be
updated with a smaller timestep than needed for the particles. A 4th order
Runge-Kutta scheme is used to update the B-field. The RK scheme is repeated a
number of times during each half-step outlined above. The number of sub-steps
used can be specified by the user through a runtime simulation parameter
(see [input parameters section](../../usage/parameters.md#running-cpp-parameters-hybrid-model)).

<a id="theory-hybrid-model-elec-temp"></a>

### Electron pressure

The electron pressure is assumed to be a scalar quantity and calculated using the given
input parameters, $T_{e0}$, $n_0$ and $\gamma$ using

> $$
> P_e = n_0T_{e0}\left( \frac{n_e}{n_0} \right)^\gamma.
> $$

The isothermal limit is given by $\gamma = 1$ while $\gamma = 5/3$
(default) produces the adiabatic limit.

### Electron current

WarpX’s displacement current diagnostic can be used to output the electron current in
the kinetic-fluid hybrid model since in the absence of kinetic electrons, and under
the assumption of zero displacement current, that diagnostic simply calculates the
hybrid model’s electron current.

## Model derivation

The basic justification for the hybrid model is that the system to which it is
applied is dominated by ion kinetics, with ions moving much slower than electrons
and photons. In this scenario two critical approximations can be made, namely,
neutrality ($n_e=n_i$) and the Maxwell-Ampere equation can be simplified by
neglecting the displacement current term [[2](#id10)], giving,

> $$
> \mu_0\vec{J} = \vec{\nabla}\times\vec{B},
> $$

where $\vec{J} = \sum_{s\neq e}\vec{J}_s + \vec{J}_e + \vec{J}_{ext}$ is the total electrical current,
i.e. the sum of electron and ion currents as well as any external current (not captured through plasma
particles). Since ions are treated in the regular
PIC manner, the ion current, $\sum_{s\neq e}\vec{J}_s$, is known during a simulation. Therefore,
given the magnetic field, the electron current can be calculated.

The electron momentum transport equation (obtained from multiplying the Vlasov equation by mass and
integrating over velocity), also called the generalized Ohm’s law, is given by:

> $$
> en_e\vec{E} = \frac{m}{e}\frac{\partial \vec{J}_e}{\partial t} + \frac{m}{e}\left( \vec{U}_e\cdot\nabla \right) \vec{J}_e - \nabla\cdot {\overleftrightarrow P}_e - \vec{J}_e\times\vec{B}+\vec{R}_e
> $$

where $\vec{U}_e = \vec{J}_e/(en_e)$ is the electron fluid velocity,
${\overleftrightarrow P}_e$ is the electron pressure tensor and
$\vec{R}_e$ is the drag force due to collisions between electrons and ions.
Applying the above momentum equation to the Maxwell-Faraday equation ($\frac{\partial\vec{B}}{\partial t} = -\nabla\times\vec{E}$)
and substituting in $\vec{J}$ calculated from the Maxwell-Ampere equation, gives,

> $$
> \frac{\partial\vec{J}_e}{\partial t} = -\frac{1}{\mu_0}\nabla\times\left(\nabla\times\vec{E}\right) - \frac{\partial\vec{J}_{ext}}{\partial t} - \sum_{s\neq e}\frac{\partial\vec{J}_s}{\partial t}.
> $$

Plugging this back into the generalized Ohm’s law gives:

> $$
> \left(en_e +\frac{m}{e\mu_0}\nabla\times\nabla\times\right)\vec{E} =&
> - \frac{m}{e}\left( \frac{\partial\vec{J}_{ext}}{\partial t} + \sum_{s\neq e}\frac{\partial\vec{J}_s}{\partial t} \right) \\
> &+ \frac{m}{e}\left( \vec{U}_e\cdot\nabla \right) \vec{J}_e - \nabla\cdot {\overleftrightarrow P}_e - \vec{J}_e\times\vec{B}+\vec{R}_e.
> $$

If we now further assume electrons are inertialess (i.e. $m=0$), the above equation simplifies to,

> $$
> en_e\vec{E} = -\vec{J}_e\times\vec{B}-\nabla\cdot{\overleftrightarrow P}_e+\vec{R}_e.
> $$

Making the further simplifying assumptions that the electron pressure is isotropic and that
the electron drag term can be written using a simple resistivity ($\eta$) and hyper-resistivity ($\eta_h$)
i.e. $\vec{R}_e = en_e(\eta-\eta_h \nabla^2)\vec{J}$, brings us to the implemented form of
Ohm’s law:

> $$
> \vec{E} = -\frac{1}{en_e}\left( \vec{J}_e\times\vec{B} + \nabla P_e \right)+\eta\vec{J}-\eta_h \nabla^2\vec{J}.
> $$

Lastly, if an electron temperature is given from which the electron pressure can
be calculated, the model is fully constrained and can be evolved given initial
conditions.
