<a id="theory-boostedframe"></a>

# Moving window and optimal Lorentz boosted frame

The simulations of plasma accelerators from first principles are extremely computationally intensive, due to the need to resolve the evolution of a driver (laser or particle beam) and an accelerated particle beam into a plasma structure that is orders of magnitude longer and wider than the accelerated beam. As is customary in the modeling of particle beam dynamics in standard particle accelerators, a moving window is commonly used to follow the driver, the wake and the accelerated beam. This results in huge savings, by avoiding the meshing of the entire plasma that is orders of magnitude longer than the other length scales of interest.

<a id="fig-boosted-frame"></a>
![Scale transformation in boosted frame simulation](theory/Boosted_frame.png)

Even using a moving window, however, a full PIC simulation of a plasma accelerator can be extraordinarily demanding computationally, as many time steps are needed to resolve the crossing of the short driver beam with the plasma column.
As it turns out, choosing an optimal frame of reference that travels close to the speed of light in the direction of the laser or particle beam (as opposed to the usual choice of the laboratory frame) enables speedups by orders of magnitude [[1](#id158), [2](#id132)].
This is a result of the properties of Lorentz contraction and dilation of space and time.
In the frame of the laboratory, a very short driver (laser or particle) beam propagates through a much longer plasma column, necessitating millions to tens of millions of time steps for parameters in the range of the BELLA or FACET-II experiments.
As sketched in [Fig. 43](#fig-boosted-frame), in a frame moving with the driver beam in the plasma at velocity $v=\beta c$ (where $c$ is the speed of light in vacuum), the beam length is now elongated by $\approx(1+\beta)\gamma$ while the plasma contracts by $\gamma$ (where $\gamma=1/\sqrt{1-\beta^2}$ is the relativistic factor associated with the frame velocity)
The number of time steps that is needed to simulate a “longer” beam through a “shorter” plasma is now reduced by up to $\approx(1+\beta) \gamma^2$ (a detailed derivation of the speedup is given below).

#### NOTE
For additional reading on inputs and outputs in boosted frame simulations, consider the following pages:

* [Inputs and Outputs](boosted_frame/input_output.md)

The modeling of a plasma acceleration stage in a boosted frame
involves the fully electromagnetic modeling of a plasma propagating at near the speed of light, for which Numerical Cerenkov
[[3](#id222), [4](#id183)] is a potential issue, as explained in more details below.
In addition, for a frame of reference moving in the direction of the accelerated beam (or equivalently the wake of the laser),
waves emitted by the plasma in the forward direction expand
while the ones emitted in the backward direction contract, following the properties of the Lorentz transformation.
If one had to resolve both forward and backward propagating
waves emitted from the plasma, there would be no gain in selecting a frame different from the laboratory frame. However,
the physics of interest for a laser wakefield is the laser driving the wake, the wake, and the accelerated beam.
Backscatter is weak in the short-pulse regime, and does not
interact as strongly with the beam as do the forward propagating waves
which stay in phase for a long period. It is thus often assumed that the backward propagating waves
can be neglected in the modeling of plasma accelerator stages. The accuracy of this assumption has been demonstrated by
comparison between explicit codes which include both forward and backward waves and envelope or quasistatic codes which neglect backward waves [[5](#id198), [6](#id106), [7](#id145)].

## Theoretical speedup dependency with the frame boost

The derivation that is given here reproduces the one given in Vay *et al.* [[2](#id132)], where the obtainable speedup is derived as an extension of the formula that was derived earlier [[1](#id158)], taking in addition into account the group velocity of the laser as it traverses the plasma.

Assuming that the simulation box is a fixed number of plasma periods long, which implies the use (which is standard) of a moving window following
the wake and accelerated beam, the speedup is given by the ratio of the time taken by the laser pulse and the plasma to cross each other, divided by the shortest time scale of interest, that is the laser period. To first order, the wake velocity $v_w$ is set by the 1D group velocity of the laser driver, which in the linear (low intensity) limit, is given by [[8](#id84)]:

$$
v_w/c=\beta_w=\left(1-\frac{\omega_p^2}{\omega^2}\right)^{1/2}

$$

where $\omega_p=\sqrt{(n_e e^2)/(\epsilon_0 m_e)}$ is the plasma frequency, $\omega=2\pi c/\lambda$ is the laser frequency, $n_e$ is the plasma density, $\lambda$ is the laser wavelength in vacuum, $\epsilon_0$ is the permittivity of vacuum, $c$ is the speed of light in vacuum, and $e$ and $m_e$ are respectively the charge and mass of the electron.

In practice, the runs are typically stopped when the last electron beam macro-particle exits the plasma, and a measure of the total time of the simulation is then given by

$$
T=\frac{L+\eta \lambda_p}{v_w-v_p}

$$

where $\lambda_p\approx 2\pi c/\omega_p$ is the wake wavelength, $L$ is the plasma length, $v_w$ and $v_p=\beta_p c$ are respectively the velocity of the wake and of the plasma relative to the frame of reference, and $\eta$ is an adjustable parameter for taking into account the fraction of the wake which exited the plasma at the end of the simulation.
For a beam injected into the $n^{th}$ bucket, $\eta$ would be set to $n-1/2$. If positrons were considered, they would be injected half a wake period ahead of the location of the electrons injection position for a given period, and one would have $\eta=n-1$. The numerical cost $R_t$ scales as the ratio of the total time to the shortest timescale of interest, which is the inverse of the laser frequency, and is thus given by

$$
R_t=\frac{T c}{\lambda}=\frac{\left(L+\eta \lambda_p\right)}{\left(\beta_w-\beta_p\right) \lambda}

$$

In the laboratory, $v_p=0$ and the expression simplifies to

$$
R_{lab}=\frac{T c}{\lambda}=\frac{\left(L+\eta \lambda_p\right)}{\beta_w \lambda}

$$

In a frame moving at $\beta c$, the quantities become

$$
\begin{aligned}
\lambda_p^* & = \lambda_p/\left[\gamma \left(1-\beta_w \beta\right)\right]
\\
L^* & = L/\gamma
\\
\lambda^* & = \gamma\left(1+\beta\right) \lambda
\\
\beta_w^* & = \left(\beta_w-\beta\right)/\left(1-\beta_w\beta\right)
\\
v_p^* & = -\beta c
\\
T^* & = \frac{L^*+\eta \lambda_p^*}{v_w^*-v_p^*}
\\
R_t^* & = \frac{T^* c}{\lambda^*} = \frac{\left(L^*+\eta \lambda_p^*\right)}{\left(\beta_w^*+\beta\right) \lambda^*}
\end{aligned}
$$

where $\gamma=1/\sqrt{1-\beta^2}$.

The expected speedup from performing the simulation in a boosted frame is given by the ratio of $R_{lab}$ and $R_t^*$

<a id="equation-eq-scaling1d0"></a>
$$
S=\frac{R_{lab}}{R_t^*}=\frac{\left(1+\beta\right)\left(L+\eta \lambda_p\right)}{\left(1-\beta\beta_w\right)L+\eta \lambda_p}

$$

We note that assuming that $\beta_w\approx1$ (which is a valid approximation for most practical cases of interest) and that $\gamma<<\gamma_w$, this expression is consistent with the expression derived earlier [[1](#id158)] for the laser-plasma acceleration case, which states that $R_t^*=\alpha R_t/\left(1+\beta\right)$ with $\alpha=\left(1-\beta+l/L\right)/\left(1+l/L\right)$, where $l$ is the laser length which is generally proportional to $\eta \lambda_p$, and $S=R_t/R_T^*$. However, higher values of $\gamma$ are of interest for maximum speedup, as shown below.

For intense lasers ($a\sim 1$) typically used for acceleration, the energy gain is limited by dephasing [[9](#id219)], which occurs over a scale length $L_d \sim \lambda_p^3/2\lambda^2$.
Acceleration is compromised beyond $L_d$ and in practice, the plasma length is proportional to the dephasing length, i.e. $L= \xi L_d$. In most cases, $\gamma_w^2>>1$, which allows the approximations $\beta_w\approx1-\lambda^2/2\lambda_p^2$, and $L=\xi \lambda_p^3/2\lambda^2\approx \xi \gamma_w^2 \lambda_p/2>>\eta \lambda_p$, so that Eq.([46](#equation-eq-scaling1d0)) becomes

<a id="equation-eq-scaling1d"></a>
$$
S=\left(1+\beta\right)^2\gamma^2\frac{\xi\gamma_w^2}{\xi\gamma_w^2+\left(1+\beta\right)\gamma^2\left(\xi\beta/2+2\eta\right)}

$$

For low values of $\gamma$, i.e. when $\gamma<<\gamma_w$, Eq.([47](#equation-eq-scaling1d)) reduces to

<a id="equation-eq-scaling1d-simpl2"></a>
$$
S_{\gamma<<\gamma_w}=\left(1+\beta\right)^2\gamma^2

$$

Conversely, if $\gamma\rightarrow\infty$, Eq.(Eq_scaling1d) becomes

<a id="equation-eq-scaling-gamma-inf"></a>
$$
S_{\gamma\rightarrow\infty}=\frac{4}{1+4\eta/\xi}\gamma_w^2

$$

Finally, in the frame of the wake, i.e. when $\gamma=\gamma_w$, assuming that $\beta_w\approx1$, Eq.([47](#equation-eq-scaling1d)) gives

<a id="equation-eq-scaling-gamma-wake"></a>
$$
S_{\gamma=\gamma_w}\approx\frac{2}{1+2\eta/\xi}\gamma_w^2

$$

Since $\eta$ and $\xi$ are of order unity, and the practical regimes of most interest satisfy $\gamma_w^2>>1$, the speedup that is obtained by using the frame of the wake will be near the maximum obtainable value given by Eq.([49](#equation-eq-scaling-gamma-inf)).

Note that without the use of a moving window, the relativistic effects that are at play in the time domain would also be at play in the spatial domain [[1](#id158)], and the $\gamma^2$ scaling would transform to $\gamma^4$. Hence, it is important to use a moving window even in simulations in a Lorentz boosted frame. For very high values of the boosted frame, the optimal velocity of the moving window may vanish (i.e. no moving window) or even reverse.

<a id="theory-boostedframe-galilean"></a>

## Numerical Stability and alternate formulation in a Galilean frame

The numerical Cherenkov instability (NCI) [[10](#id80)]
is the most serious numerical instability affecting multidimensional
PIC simulations of relativistic particle beams and streaming plasmas
[[11](#id184), [12](#id157), [13](#id88), [14](#id248), [15](#id122), [16](#id173)].
It arises from coupling between possibly numerically distorted electromagnetic modes and spurious
beam modes, the latter due to the mismatch between the Lagrangian
treatment of particles and the Eulerian treatment of fields [[17](#id108)].

In recent papers the electromagnetic dispersion
relations for the numerical Cherenkov instability were derived and solved for both FDTD [[15](#id122), [18](#id274)]
and PSATD [[19](#id186), [20](#id154)] algorithms.

Several solutions have been proposed to mitigate the NCI [[19](#id186), [20](#id154), [21](#id194), [22](#id179), [23](#id240), [24](#id188)]. Although
these solutions efficiently reduce the numerical instability,
they typically introduce either strong smoothing of the currents and
fields, or arbitrary numerical corrections, which are
tuned specifically against the NCI and go beyond the
natural discretization of the underlying physical equation. Therefore,
it is sometimes unclear to what extent these added corrections could impact the
physics at stake for a given resolution.

For instance, NCI-specific corrections include periodically smoothing
the electromagnetic field components [[11](#id184)],
using a special time step [[12](#id157), [13](#id88)] or
applying a wide-band smoothing of the current components [[12](#id157), [13](#id88), [25](#id276)]. Another set of mitigation methods
involve scaling the deposited
currents by a carefully-designed wavenumber-dependent factor
[[18](#id274), [20](#id154)] or slightly modifying the
ratio of electric and magnetic fields ($E/B$) before gathering their
value onto the macroparticles
[[19](#id186), [22](#id179)].
Yet another set of NCI-specific corrections
[[23](#id240), [24](#id188)] consists
in combining a small timestep $\Delta t$, a sharp low-pass spatial filter,
and a spectral or high-order scheme that is tuned so as to
create a small, artificial “bump” in the dispersion relation
[[23](#id240)]. While most mitigation methods have only been applied
to Cartesian geometry, this last
set of methods [[23](#id240), [24](#id188)]
has the remarkable property that it can be applied
[[24](#id188)] to both Cartesian geometry and
quasi-cylindrical geometry (i.e. cylindrical geometry with
azimuthal Fourier decomposition [[26](#id113), [27](#id193), [28](#id244)]). However,
the use of a small timestep proportionally slows down the progress of
the simulation, and the artificial “bump” is again an arbitrary correction
that departs from the underlying physics.

A new scheme was recently proposed, in Kirchen *et al.* [[29](#id277)], Lehe *et al.* [[30](#id278)], which
completely eliminates the NCI for a plasma drifting at a uniform relativistic velocity
– with no arbitrary correction – by simply integrating
the PIC equations in *Galilean coordinates* (also known as
*comoving coordinates*). More precisely, in the new
method, the Maxwell equations *in Galilean coordinates* are integrated
analytically, using only natural hypotheses, within the PSATD
framework (Pseudo-Spectral-Analytical-Time-Domain [[4](#id183), [31](#id159)]).

The idea of the proposed scheme is to perform a Galilean change of
coordinates, and to carry out the simulation in the new coordinates:

<a id="equation-change-var"></a>
$$
\boldsymbol{x}' = \boldsymbol{x} - \boldsymbol{v}_{gal}t

$$

where $\boldsymbol{x} = x\,\boldsymbol{u}_x + y\,\boldsymbol{u}_y + z\,\boldsymbol{u}_z$ and
$\boldsymbol{x}' = x'\,\boldsymbol{u}_x + y'\,\boldsymbol{u}_y + z'\,\boldsymbol{u}_z$ are the
position vectors in the standard and Galilean coordinates
respectively.

When choosing $\boldsymbol{v}_{gal}= \boldsymbol{v}_0$, where
$\boldsymbol{v}_0$ is the speed of the bulk of the relativistic
plasma, the plasma does not move with respect to the grid in the Galilean
coordinates $\boldsymbol{x}'$ – or, equivalently, in the standard
coordinates $\boldsymbol{x}$, the grid moves along with the plasma. The heuristic intuition behind this scheme
is that these coordinates should prevent the discrepancy between the Lagrangian and
Eulerian point of view, which gives rise to the NCI [[17](#id108)].

An important remark is that the Galilean change of
coordinates in Eq. ([51](#equation-change-var)) is a simple translation. Thus, when used in
the context of Lorentz-boosted simulations, it does
of course preserve the relativistic dilatation of space and time which gives rise to the
characteristic computational speedup of the boosted-frame technique.

Another important remark is that the Galilean scheme is *not*
equivalent to a moving window (and in fact the Galilean scheme can be
independently *combined* with a moving window). Whereas in a
moving window, gridpoints are added and removed so as to effectively
translate the boundaries, in the Galilean scheme the gridpoints
*themselves* are not only translated but in this case, the physical equations
are modified accordingly. Most importantly, the assumed time evolution of
the current $\boldsymbol{J}$ within one timestep is different in a standard PSATD scheme with moving
window and in a Galilean PSATD scheme [[30](#id278)].

In the Galilean coordinates $\boldsymbol{x}'$, the equations of particle
motion and the Maxwell equations take the form

<a id="equation-motion1"></a>
$$
\frac{d\boldsymbol{x}'}{dt} = \frac{\boldsymbol{p}}{\gamma m} - \boldsymbol{v}_{gal}

$$

<a id="equation-motion2"></a>
$$
\frac{d\boldsymbol{p}}{dt} = q \left( \boldsymbol{E} + \frac{\boldsymbol{p}}{\gamma m} \times \boldsymbol{B} \right)

$$

<a id="equation-maxwell1"></a>
$$
\left(  \frac{\partial \;}{\partial t} - \boldsymbol{v}_{gal}\cdot\boldsymbol{\nabla'}\right)\boldsymbol{B} = -\boldsymbol{\nabla'}\times\boldsymbol{E}

$$

<a id="equation-maxwell2"></a>
$$
\frac{1}{c^2}\left(  \frac{\partial \;}{\partial t} - \boldsymbol{v}_{gal}\cdot\boldsymbol{\nabla'}\right)\boldsymbol{E} = \boldsymbol{\nabla'}\times\boldsymbol{B} - \mu_0\boldsymbol{J}

$$

where $\boldsymbol{\nabla'}$ denotes a spatial derivative with respect to the
Galilean coordinates $\boldsymbol{x}'$.

Integrating these equations from $t=n\Delta
t$ to $t=(n+1)\Delta t$ results in the following update equations (see
Lehe *et al.* [[30](#id278)] for the details of the derivation):

<a id="equation-disc-maxwell1"></a>
$$
\begin{aligned}
\mathbf{\tilde{B}}^{n+1} & = \theta^2 C \mathbf{\tilde{B}}^n -\frac{\theta^2 S}{ck}i\boldsymbol{k}\times \mathbf{\tilde{E}}^n \nonumber
\\
                         & + \;\frac{\theta \chi_1}{\epsilon_0c^2k^2}\;i\boldsymbol{k} \times \mathbf{\tilde{J}}^{n+1/2}
\end{aligned}

$$

<a id="equation-disc-maxwell2"></a>
$$
\begin{aligned}
\mathbf{\tilde{E}}^{n+1} & = \theta^2 C  \mathbf{\tilde{E}}^n +\frac{\theta^2 S}{k} \,c i\boldsymbol{k}\times \mathbf{\tilde{B}}^n \nonumber
\\
                         & + \frac{i\nu \theta \chi_1 - \theta^2S}{\epsilon_0 ck} \; \mathbf{\tilde{J}}^{n+1/2}\nonumber
\\
                         & - \frac{1}{\epsilon_0k^2}\left(\; \chi_2\;\hat{\mathcal{\rho}}^{n+1} - \theta^2\chi_3\;\hat{\mathcal{\rho}}^{n} \;\right) i\boldsymbol{k}
\end{aligned}

$$

where we used the short-hand notations
$\mathbf{\tilde{E}}^n \equiv \mathbf{\tilde{E}}(\boldsymbol{k}, n\Delta t)$,
$\mathbf{\tilde{B}}^n \equiv \mathbf{\tilde{B}}(\boldsymbol{k}, n\Delta t)$ as well as:

<a id="equation-def-c-s"></a>
$$
C = \cos(ck\Delta t), \quad S = \sin(ck\Delta t), \quad k = |\boldsymbol{k}|,

$$

<a id="equation-def-nu-theta"></a>
$$
\nu = \frac{\boldsymbol{k}\cdot\boldsymbol{v}_{gal}}{ck}, \quad \theta = e^{i\boldsymbol{k}\cdot\boldsymbol{v}_{gal}\Delta t/2},

$$

<a id="equation-def-chi1"></a>
$$
\chi_1 = \frac{1}{1 -\nu^2} \left( \theta^* - C \theta + i \nu \theta S \right),

$$

<a id="equation-def-chi2"></a>
$$
\chi_2 = \frac{\chi_1 - \theta(1-C)}{\theta^*-\theta}

$$

<a id="equation-def-chi3"></a>
$$
\chi_3 = \frac{\chi_1-\theta^*(1-C)}{\theta^*-\theta}

$$

Note that, in the limit $\boldsymbol{v}_{gal}=\boldsymbol{0}$,
Eqs. ([56](#equation-disc-maxwell1)) and ([57](#equation-disc-maxwell2)) reduce to the standard PSATD
equations [[4](#id183)], as expected.
As shown in Kirchen *et al.* [[29](#id277)], Lehe *et al.* [[30](#id278)],
the elimination of the NCI with the new Galilean integration is verified empirically via PIC simulations of uniform drifting plasmas and laser-driven plasma acceleration stages, and confirmed by a theoretical analysis of the instability.
