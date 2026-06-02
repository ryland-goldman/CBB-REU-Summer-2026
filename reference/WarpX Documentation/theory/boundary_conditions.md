<a id="theory-bc"></a>

# Boundary conditions

<a id="theory-bc-pml"></a>

## Perfectly Matched Layer: open boundary condition for electromagnetic waves

For the transverse electric (TE) case, the original Berenger’s Perfectly Matched Layer (PML) paper [[1](#id158)] writes

<a id="equation-pml-def-1"></a>
$$
\varepsilon _{0}\frac{\partial E_{x}}{\partial t}+\sigma _{y}E_{x} = \frac{\partial H_{z}}{\partial y}

$$

<a id="equation-pml-def-2"></a>
$$
\varepsilon _{0}\frac{\partial E_{y}}{\partial t}+\sigma _{x}E_{y} = -\frac{\partial H_{z}}{\partial x}

$$

<a id="equation-pml-def-3"></a>
$$
\mu _{0}\frac{\partial H_{zx}}{\partial t}+\sigma ^{*}_{x}H_{zx} = -\frac{\partial E_{y}}{\partial x}

$$

<a id="equation-pml-def-4"></a>
$$
\mu _{0}\frac{\partial H_{zy}}{\partial t}+\sigma ^{*}_{y}H_{zy} = \frac{\partial E_{x}}{\partial y}

$$

<a id="equation-pml-def-5"></a>
$$
H_{z}  = H_{zx}+H_{zy}

$$

This can be generalized to

<a id="equation-apml-def-1"></a>
$$
\varepsilon _{0}\frac{\partial E_{x}}{\partial t}+\sigma _{y}E_{x} = \frac{c_{y}}{c}\frac{\partial H_{z}}{\partial y}+\overline{\sigma }_{y}H_{z}

$$

<a id="equation-apml-def-2"></a>
$$
\varepsilon _{0}\frac{\partial E_{y}}{\partial t}+\sigma _{x}E_{y} = -\frac{c_{x}}{c}\frac{\partial H_{z}}{\partial x}+\overline{\sigma }_{x}H_{z}

$$

<a id="equation-apml-def-3"></a>
$$
\mu _{0}\frac{\partial H_{zx}}{\partial t}+\sigma ^{*}_{x}H_{zx} = -\frac{c^{*}_{x}}{c}\frac{\partial E_{y}}{\partial x}+\overline{\sigma }_{x}^{*}E_{y}

$$

<a id="equation-apml-def-4"></a>
$$
\mu _{0}\frac{\partial H_{zy}}{\partial t}+\sigma ^{*}_{y}H_{zy} = \frac{c^{*}_{y}}{c}\frac{\partial E_{x}}{\partial y}+\overline{\sigma }_{y}^{*}E_{x}

$$

<a id="equation-apml-def-5"></a>
$$
H_{z} = H_{zx}+H_{zy}

$$

For $c_{x}=c_{y}=c^{*}_{x}=c^{*}_{y}=c$ and $\overline{\sigma }_{x}=\overline{\sigma }_{y}=\overline{\sigma }_{x}^{*}=\overline{\sigma }_{y}^{*}=0$,
this system reduces to the Berenger PML medium, while adding the additional
constraint $\sigma _{x}=\sigma _{y}=\sigma _{x}^{*}=\sigma _{y}^{*}=0$
leads to the system of Maxwell equations in vacuum.

<a id="theory-bc-propa-plane-wave"></a>

### Propagation of a Plane Wave in an APML Medium

We consider a plane wave of magnitude ($E_{0},H_{zx0},H_{zy0}$)
and pulsation $\omega$ propagating in the APML medium with an
angle $\varphi$ relative to the x axis

<a id="equation-plane-wave-apml-def-1"></a>
$$
E_{x} = -E_{0}\sin \varphi \: e^{i\omega \left( t-\alpha x-\beta y\right) }

$$

<a id="equation-plane-wave-apml-def-2"></a>
$$
E_{y} = E_{0}\cos \varphi \: e^{i\omega \left( t-\alpha x-\beta y\right) }

$$

<a id="equation-plane-wave-ampl-def-3"></a>
$$
H_{zx} = H_{zx0} \: e^{i\omega \left( t-\alpha x-\beta y\right) }

$$

<a id="equation-plane-wave-apml-def-4"></a>
$$
H_{zy} = H_{zy0} \: e^{i\omega \left( t-\alpha x-\beta y\right) }

$$

where $\alpha$ and $\beta$ are two complex constants to
be determined.

Introducing Eqs. ([34](#equation-plane-wave-apml-def-1)), ([35](#equation-plane-wave-apml-def-2)),
([36](#equation-plane-wave-ampl-def-3)) and ([37](#equation-plane-wave-apml-def-4))
into Eqs. ([29](#equation-apml-def-1)), ([30](#equation-apml-def-2)), ([31](#equation-apml-def-3))
and ([32](#equation-apml-def-4)) gives

<a id="equation-plane-wave-apml-1-1"></a>
$$
\varepsilon _{0}E_{0}\sin \varphi -i\frac{\sigma _{y}}{\omega }E_{0}\sin \varphi = \beta \frac{c_{y}}{c}\left( H_{zx0}+H_{zy0}\right) +i\frac{\overline{\sigma }_{y}}{\omega }\left( H_{zx0}+H_{zy0}\right)

$$

<a id="equation-plane-wave-apml-1-2"></a>
$$
\varepsilon _{0}E_{0}\cos \varphi -i\frac{\sigma _{x}}{\omega }E_{0}\cos \varphi = \alpha \frac{c_{x}}{c}\left( H_{zx0}+H_{zy0}\right) -i\frac{\overline{\sigma }_{x}}{\omega }\left( H_{zx0}+H_{zy0}\right)

$$

<a id="equation-plane-wave-apml-1-3"></a>
$$
\mu _{0}H_{zx0}-i\frac{\sigma ^{*}_{x}}{\omega }H_{zx0} = \alpha \frac{c^{*}_{x}}{c}E_{0}\cos \varphi -i\frac{\overline{\sigma }^{*}_{x}}{\omega }E_{0}\cos \varphi

$$

<a id="equation-plane-wave-apml-1-4"></a>
$$
\mu _{0}H_{zy0}-i\frac{\sigma ^{*}_{y}}{\omega }H_{zy0} = \beta \frac{c^{*}_{y}}{c}E_{0}\sin \varphi +i\frac{\overline{\sigma }^{*}_{y}}{\omega }E_{0}\sin \varphi

$$

Defining $Z=E_{0}/\left( H_{zx0}+H_{zy0}\right)$ and using Eqs. ([38](#equation-plane-wave-apml-1-1))
and ([39](#equation-plane-wave-apml-1-2)), we get

<a id="equation-plane-wave-apml-beta-of-g"></a>
$$
\beta = \left[ Z\left( \varepsilon _{0}-i\frac{\sigma _{y}}{\omega }\right) \sin \varphi -i\frac{\overline{\sigma }_{y}}{\omega }\right] \frac{c}{c_{y}}

$$

<a id="equation-plane-wave-apml-alpha-of-g"></a>
$$
\alpha = \left[ Z\left( \varepsilon _{0}-i\frac{\sigma _{x}}{\omega }\right) \cos \varphi +i\frac{\overline{\sigma }_{x}}{\omega }\right] \frac{c}{c_{x}}

$$

Adding $H_{zx0}$ and $H_{zy0}$ from Eqs. ([40](#equation-plane-wave-apml-1-3))
and ([41](#equation-plane-wave-apml-1-4)) and substituting the expressions
for $\alpha$ and $\beta$ from Eqs. ([42](#equation-plane-wave-apml-beta-of-g))
and ([43](#equation-plane-wave-apml-alpha-of-g)) yields

$$
\begin{aligned}
\frac{1}{Z} & = \frac{Z\left( \varepsilon _{0}-i\frac{\sigma _{x}}{\omega }\right) \cos \varphi \frac{c^{*}_{x}}{c_{x}}+i\frac{\overline{\sigma }_{x}}{\omega }\frac{c^{*}_{x}}{c_{x}}-i\frac{\overline{\sigma }^{*}_{x}}{\omega }}{\mu _{0}-i\frac{\sigma ^{*}_{x}}{\omega }}\cos \varphi \nonumber
\\
            & + \frac{Z\left( \varepsilon _{0}-i\frac{\sigma _{y}}{\omega }\right) \sin \varphi \frac{c^{*}_{y}}{c_{y}}-i\frac{\overline{\sigma }_{y}}{\omega }\frac{c^{*}_{y}}{c_{y}}+i\frac{\overline{\sigma }^{*}_{y}}{\omega }}{\mu _{0}-i\frac{\sigma ^{*}_{y}}{\omega }}\sin \varphi
\end{aligned}
$$

If $c_{x}=c^{*}_{x}$, $c_{y}=c^{*}_{y}$, $\overline{\sigma }_{x}=\overline{\sigma }^{*}_{x}$, $\overline{\sigma }_{y}=\overline{\sigma }^{*}_{y}$, $\frac{\sigma _{x}}{\varepsilon _{0}}=\frac{\sigma ^{*}_{x}}{\mu _{0}}$ and $\frac{\sigma _{y}}{\varepsilon _{0}}=\frac{\sigma ^{*}_{y}}{\mu _{0}}$ then

<a id="equation-apml-impedance"></a>
$$
Z = \pm \sqrt{\frac{\mu _{0}}{\varepsilon _{0}}}

$$

which is the impedance of vacuum. Hence, like the PML, given some
restrictions on the parameters, the APML does not generate any reflection
at any angle and any frequency. As for the PML, this property is not
retained after discretization, as shown subsequently.

Calling $\psi$ any component of the field and $\psi _{0}$
its magnitude, we get from Eqs. ([34](#equation-plane-wave-apml-def-1)), ([42](#equation-plane-wave-apml-beta-of-g)),
([43](#equation-plane-wave-apml-alpha-of-g)) and ([44](#equation-apml-impedance)) that

<a id="equation-plane-wave-absorption"></a>
$$
\psi =\psi _{0} \: e^{i\omega \left( t\mp x\cos \varphi /c_{x}\mp y\sin \varphi /c_{y}\right) }e^{-\left( \pm \frac{\sigma _{x}\cos \varphi }{\varepsilon _{0}c_{x}}+\overline{\sigma }_{x}\frac{c}{c_{x}}\right) x} e^{-\left( \pm \frac{\sigma _{y}\sin \varphi }{\varepsilon _{0}c_{y}}+\overline{\sigma }_{y}\frac{c}{c_{y}}\right) y}.

$$

We assume that we have an APML layer of thickness $\delta$ (measured
along $x$) and that $\sigma _{y}=\overline{\sigma }_{y}=0$
and $c_{y}=c.$ Using ([45](#equation-plane-wave-absorption)), we determine
that the coefficient of reflection given by this layer is

$$
\begin{aligned}
R_{\mathrm{APML}}\left( \theta \right) & = e^{-\left( \sigma _{x}\cos \varphi /\varepsilon _{0}c_{x}+\overline{\sigma }_{x}c/c_{x}\right) \delta }e^{-\left( \sigma _{x}\cos \varphi /\varepsilon _{0}c_{x}-\overline{\sigma }_{x}c/c_{x}\right) \delta },\nonumber
\\
                              & = e^{-2\left( \sigma _{x}\cos \varphi /\varepsilon _{0}c_{x}\right) \delta },
 \end{aligned}
$$

which happens to be the same as the PML theoretical coefficient of
reflection if we assume $c_{x}=c$. Hence, it follows that for
the purpose of wave absorption, the term $\overline{\sigma }_{x}$
seems to be of no interest. However, although this conclusion is true
at the infinitesimal limit, it does not hold for the discretized counterpart.

### Discretization

In the following we set $\varepsilon_0 = \mu_0 = 1$. We discretize Eqs. ([24](#equation-pml-def-1)), ([25](#equation-pml-def-2)), ([26](#equation-pml-def-3)), and ([27](#equation-pml-def-4)) to obtain

$$
\frac{E_x|^{n+1}_{j+1/2,k,l}-E_x|^{n}_{j+1/2,k,l}}{\Delta t} + \sigma_y \frac{E_x|^{n+1}_{j+1/2,k,l}+E_x|^{n}_{j+1/2,k,l}}{2} = \frac{H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j+1/2,k-1/2,l}}{\Delta y}

$$

$$
\frac{E_y|^{n+1}_{j,k+1/2,l}-E_y|^{n}_{j,k+1/2,l}}{\Delta t} + \sigma_x \frac{E_y|^{n+1}_{j,k+1/2,l}+E_y|^{n}_{j,k+1/2,l}}{2} = - \frac{H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j-1/2,k+1/2,l}}{\Delta x}

$$

$$
\frac{H_{zx}|^{n+3/2}_{j+1/2,k+1/2,l}-H_{zx}|^{n+1/2}_{j+1/2,k+1/2,l}}{\Delta t} + \sigma^*_x \frac{H_{zx}|^{n+3/2}_{j+1/2,k+1/2,l}+H_{zx}|^{n+1/2}_{j+1/2,k+1/2,l}}{2} = - \frac{E_y|^{n+1}_{j+1,k+1/2,l}-E_y|^{n+1}_{j,k+1/2,l}}{\Delta x}

$$

$$
\frac{H_{zy}|^{n+3/2}_{j+1/2,k+1/2,l}-H_{zy}|^{n+1/2}_{j+1/2,k+1/2,l}}{\Delta t} + \sigma^*_y \frac{H_{zy}|^{n+3/2}_{j+1/2,k+1/2,l}+H_{zy}|^{n+1/2}_{j+1/2,k+1/2,l}}{2} = \frac{E_x|^{n+1}_{j+1/2,k+1,l}-E_x|^{n+1}_{j+1/2,k,l}}{\Delta y}

$$

and this can be solved to obtain the following leapfrog integration equations

$$
\begin{aligned}
E_x|^{n+1}_{j+1/2,k,l} & = \left(\frac{1-\sigma_y \Delta t/2}{1+\sigma_y \Delta t/2}\right) E_x|^{n}_{j+1/2,k,l} + \frac{\Delta t/\Delta y}{1+\sigma_y \Delta t/2} \left(H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j+1/2,k-1/2,l}\right)
\\
E_y|^{n+1}_{j,k+1/2,l} & = \left(\frac{1-\sigma_x \Delta t/2}{1+\sigma_x \Delta t/2}\right) E_y|^{n}_{j,k+1/2,l} - \frac{\Delta t/\Delta x}{1+\sigma_x \Delta t/2} \left(H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j-1/2,k+1/2,l}\right)
\\
H_{zx}|^{n+3/2}_{j+1/2,k+1/2,l} & = \left(\frac{1-\sigma^*_x \Delta t/2}{1+\sigma^*_x \Delta t/2}\right) H_{zx}|^{n+1/2}_{j+1/2,k+1/2,l} - \frac{\Delta t/\Delta x}{1+\sigma^*_x \Delta t/2} \left(E_y|^{n+1}_{j+1,k+1/2,l}-E_y|^{n+1}_{j,k+1/2,l}\right)
\\
H_{zy}|^{n+3/2}_{j+1/2,k+1/2,l} & = \left(\frac{1-\sigma^*_y \Delta t/2}{1+\sigma^*_y \Delta t/2}\right) H_{zy}|^{n+1/2}_{j+1/2,k+1/2,l} + \frac{\Delta t/\Delta y}{1+\sigma^*_y \Delta t/2} \left(E_x|^{n+1}_{j+1/2,k+1,l}-E_x|^{n+1}_{j+1/2,k,l}\right)
\end{aligned}
$$

If we account for higher order $\Delta t$ terms, a better approximation is given by

$$
\begin{aligned}
E_x|^{n+1}_{j+1/2,k,l} & = e^{-\sigma_y\Delta t} E_x|^{n}_{j+1/2,k,l} + \frac{1-e^{-\sigma_y\Delta t}}{\sigma_y \Delta y} \left(H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j+1/2,k-1/2,l}\right)
\\
E_y|^{n+1}_{j,k+1/2,l} & = e^{-\sigma_x\Delta t} E_y|^{n}_{j,k+1/2,l} - \frac{1-e^{-\sigma_x\Delta t}}{\sigma_x \Delta x} \left(H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j-1/2,k+1/2,l}\right)
\\
H_{zx}|^{n+3/2}_{j+1/2,k+1/2,l} & = e^{-\sigma^*_x\Delta t} H_{zx}|^{n+1/2}_{j+1/2,k+1/2,l} - \frac{1-e^{-\sigma^*_x\Delta t}}{\sigma^*_x \Delta x} \left(E_y|^{n+1}_{j+1,k+1/2,l}-E_y|^{n+1}_{j,k+1/2,l}\right)
\\
H_{zy}|^{n+3/2}_{j+1/2,k+1/2,l} & = e^{-\sigma^*_y\Delta t} H_{zy}|^{n+1/2}_{j+1/2,k+1/2,l} + \frac{1-e^{-\sigma^*_y\Delta t}}{\sigma^*_y \Delta y} \left(E_x|^{n+1}_{j+1/2,k+1,l}-E_x|^{n+1}_{j+1/2,k,l}\right)
\end{aligned}
$$

More generally, this becomes

$$
\begin{aligned}
E_x|^{n+1}_{j+1/2,k,l} & = e^{-\sigma_y\Delta t} E_x|^{n}_{j+1/2,k,l} + \frac{1-e^{-\sigma_y\Delta t}}{\sigma_y \Delta y}\frac{c_y}{c} \left(H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j+1/2,k-1/2,l}\right)
\\
E_y|^{n+1}_{j,k+1/2,l} & = e^{-\sigma_x\Delta t} E_y|^{n}_{j,k+1/2,l} - \frac{1-e^{-\sigma_x\Delta t}}{\sigma_x \Delta x}\frac{c_x}{c} \left(H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j-1/2,k+1/2,l}\right)
\\
H_{zx}|^{n+3/2}_{j+1/2,k+1/2,l} & = e^{-\sigma^*_x\Delta t} H_{zx}|^{n+1/2}_{j+1/2,k+1/2,l} - \frac{1-e^{-\sigma^*_x\Delta t}}{\sigma^*_x \Delta x}\frac{c^*_x}{c} \left(E_y|^{n+1}_{j+1,k+1/2,l}-E_y|^{n+1}_{j,k+1/2,l}\right)
\\
H_{zy}|^{n+3/2}_{j+1/2,k+1/2,l} & = e^{-\sigma^*_y\Delta t} H_{zy}|^{n+1/2}_{j+1/2,k+1/2,l} + \frac{1-e^{-\sigma^*_y\Delta t}}{\sigma^*_y \Delta y}\frac{c^*_y}{c} \left(E_x|^{n+1}_{j+1/2,k+1,l}-E_x|^{n+1}_{j+1/2,k,l}\right)
\end{aligned}
$$

If we set

$$
\begin{aligned}
c_x & = c \: e^{-\sigma_x\Delta t} \frac{\sigma_x \Delta t}{1-e^{-\sigma_x\Delta t}} \\
c_y & = c \: e^{-\sigma_y\Delta t} \frac{\sigma_y \Delta t}{1-e^{-\sigma_y\Delta t}} \\
c^*_x & = c \: e^{-\sigma^*_x\Delta t} \frac{\sigma^*_x \Delta t}{1-e^{-\sigma^*_x\Delta t}} \\
c^*_y & = c \: e^{-\sigma^*_y\Delta t} \frac{\sigma^*_y \Delta t}{1-e^{-\sigma^*_y\Delta t}}\end{aligned}
$$

then this becomes

$$
\begin{aligned}
E_x|^{n+1}_{j+1/2,k,l} & = e^{-\sigma_y\Delta t} \left[ E_x|^{n}_{j+1/2,k,l} + \frac{\Delta t}{\Delta y} \left(H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j+1/2,k-1/2,l}\right) \right]
\\
E_y|^{n+1}_{j,k+1/2,l} & = e^{-\sigma_x\Delta t} \left[ E_y|^{n}_{j,k+1/2,l} - \frac{\Delta t}{\Delta x}  \left(H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j-1/2,k+1/2,l}\right) \right]
\\
H_{zx}|^{n+3/2}_{j+1/2,k+1/2,l} & = e^{-\sigma^*_x\Delta t} \left[ H_{zx}|^{n+1/2}_{j+1/2,k+1/2,l} - \frac{\Delta t}{\Delta x}  \left(E_y|^{n+1}_{j+1,k+1/2,l}-E_y|^{n+1}_{j,k+1/2,l}\right) \right]
\\
H_{zy}|^{n+3/2}_{j+1/2,k+1/2,l} & = e^{-\sigma^*_y\Delta t} \left[ H_{zy}|^{n+1/2}_{j+1/2,k+1/2,l} + \frac{\Delta t}{\Delta y}  \left(E_x|^{n+1}_{j+1/2,k+1,l}-E_x|^{n+1}_{j+1/2,k,l}\right) \right]
\end{aligned}
$$

When the generalized conductivities are zero, the update equations are

$$
\begin{aligned}
E_x|^{n+1}_{j+1/2,k,l} & = E_x|^{n}_{j+1/2,k,l} + \frac{\Delta t}{\Delta y} \left(H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j+1/2,k-1/2,l}\right)
\\
E_y|^{n+1}_{j,k+1/2,l} & = E_y|^{n}_{j,k+1/2,l} - \frac{\Delta t}{\Delta x} \left(H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j-1/2,k+1/2,l}\right)
\\
H_{zx}|^{n+3/2}_{j+1/2,k+1/2,l} & = H_{zx}|^{n+1/2}_{j+1/2,k+1/2,l} - \frac{\Delta t}{\Delta x} \left(E_y|^{n+1}_{j+1,k+1/2,l}-E_y|^{n+1}_{j,k+1/2,l}\right)
\\
H_{zy}|^{n+3/2}_{j+1/2,k+1/2,l} & = H_{zy}|^{n+1/2}_{j+1/2,k+1/2,l} + \frac{\Delta t}{\Delta y} \left(E_x|^{n+1}_{j+1/2,k+1,l}-E_x|^{n+1}_{j+1/2,k,l}\right)
\end{aligned}
$$

as expected.

<a id="theory-bc-pec"></a>

## Perfect Electrical Conductor

This boundary can be used to model a dielectric or metallic surface.
For the electromagnetic solve, at PEC, the tangential electric field and the normal magnetic
field are set to 0. In the guard-cell region, the tangential electric field is set equal and
opposite to the respective field component in the mirror location across the PEC
boundary, and the normal electric field is set equal to the field component in the
mirror location in the domain across the PEC boundary. Similarly, the tangential
(and normal) magnetic field components are set equal (and opposite) to the respective
magnetic field components in the mirror locations across the PEC boundary.

The PEC boundary condition also impacts the deposition of charge and current density.
On the boundary the charge density and parallel current density is set to zero. If
a reflecting boundary condition is used for the particles, density overlapping
with the PEC will be reflected back into the domain (for both charge and current
density). If absorbing boundaries are used, an image charge (equal weight but
opposite charge) is considered in the mirror location accross the boundary, and
the density from that charge is also deposited in the simulation domain. [Fig. 38](#fig-pec-boundary-deposition)
shows the effect of this. The left boundary is absorbing while
the right boundary is reflecting.

<a id="fig-pec-boundary-deposition"></a>
![Current deposition at absorbing and reflecting PEC boundaries](https://user-images.githubusercontent.com/40245517/221491318-b0a2bcbc-b04f-4b8c-8ec5-e9c92e55ee53.png)

<a id="theory-bc-pmc"></a>

## Perfect Magnetic Conductor

This boundary can be used to model a symmetric surface, where charges and current are
symmetric across the boundary.
This is equivalent to the Neumann (zero-derivative) boundary condition.
For the electromagnetic solve, at PMC, the tangential magnetic field and the normal electric
field are odd across the boundary and set to 0 on the boundary.
In the guard-cell region, those fields are set equal and
opposite to the respective field component in the mirror location across the PMC boundary.
The other components, the normal magnetic field and tangential electric field, are even
and set equal to the field component in the mirror location in the domain across the PMC boundary.

The PMC boundary condition also impacts the deposition of charge and current density.
The charge and current densities deposited into the guard cells are reflected back into
the domain, adding them to the mirror cells in the domain.
This represents the charge and current from the virtual symmetric particles in the guard cells.
