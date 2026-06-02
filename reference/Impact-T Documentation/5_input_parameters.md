# 5 Input Parameters Excluding Lattice

The following gives a line by line description of the input parameters appearing before the lattice layout used in the input file **ImpactT.in**.

> *Note: the comment line starting with `!` is not included in the line number.*

---

## Line 1: Npcol, Nprow

**Npcol** — Number of columns of processors, used to decompose domain along Y dimension.

**Nprow** — Number of rows of processors, used to decompose domain along Z dimension.

---

## Line 2: Dt, Ntstep, Nbunch

**Dt** — Time step size (secs).

**Ntstep** — Maximum number of time steps. IMPACT-T will stop after `Ntstep` time steps or when the center of the bunch goes beyond the end of the lattice, whichever is first.

**Nbunch** — During simulation of the rf photoinjector, when the laser pulse length is long and the rf acceleration gradient is high, the electron beam out of the cathode can have a significantly large energy spread. In this case, performing a single Lorentz transformation from the beam frame to the lab frame during the calculation of the space-charge forces is not sufficient. To model space-charge forces more accurately, the initial distribution of the bunch can be divided longitudinally into **Nbunch** slices. The Lorentz transformation can be done for each slice/bunch and the total space-charge forces are summed. More slices mean better accuracy, but at the cost of increased computation time. The computation time increases approximately linearly with `Nbunch`. When **Nbunch > 1**, more input files named **ImpactT2.in**, **ImpactT3.in**, ...**ImpactTNbunch.in** have to be provided. These input files can have independent definitions of the particle information. The total current is distributed among the different bunch files. However, the external lattice layout has to be the same. If each input file uses a read-in initial particle data file, **partcl2.in**, **partcl3.in**, ...**partclNbunch.in** have to be provided. This will be replaced as **partcl2.data**, **partcl3.data**, ...**partclNbunch.data** soon. Since each bin/bunch has independent parameters, the IMPACT-T code can be used to simulate a beam with multiple species.

---

## Line 3: Dim, Np, Flagmap, Flagerr, Flagdiag, Flagimg, Zimage

**Dim** — Random seed integer > 0.

**Np** — Number of macroparticles to track.

**Flagmap** — Type of integrator. Currently must be set to 1.

**Flagerr** — Error study flag. `0` - no misalignment and rotation errors; `1` - misalignment and rotation errors are allowed for Quadrupole, Multipole (Sextupole, Octupole, Decapole) and SolRF elements. This function can also be used to simulate the beam transport through rotated beam line elements such as skew quadrupole etc.

**Flagdiag** — Diagnostics flag: `1` - output the information at given time; `2` - output the information at the location of bunch centroid by drifting the particles to that location; `3` or more - no output.

**Flagimg** — Image charge flag. If set to 1 then the image charge forces due to the cathode are included. The cathode is always assumed to be at `z = 0`. To not include the image charge forces set `imchgF` to 0.

**Zimage** — z position beyond which image charge forces are neglected. Set z small to speed up the calculation but large enough so that the results are not affected.

---

## Line 4: Nx, Ny, Nz, Flagbc, Xrad, Yrad, Perdlen

**Nx, Ny, Nz** — Number of mesh points in x, y, and z. For the open boundary condition, currently they have to be a **power of 2** due to the FFT algorithm used. In the future, these will be generalized to arbitrary integer numbers.

**Flagbc** — Field boundary condition flag. Currently must be set to 1 which corresponds to an open boundary condition.

**Xrad, Yrad, Perdlen** — Size of computational domain. Xrad and Yrad define the transverse size, Perdlen defines the longitudinal size. Here, Perdlen has to be greater than the beam line lattice length.

---

## Line 5: Flagdist, Rstartflg, Flagsbstp, Nemission, Temission

**Flagdist** — Type of the initial distribution. This is a number between 1 and 27:

| Value | Name | Description |
|-------|------|-------------|
| 1 | Uniform | 6d uniform distribution |
| 2 | Gauss3 | 6d Gaussian distribution |
| 3 | Waterbag | 6d Waterbag distribution |
| 4 | Semigauss | 3d Waterbag distribution in spatial and 3d Gaussian distribution in momentum space |
| 5 | KV3d | Transverse KV distribution and longitudinal uniform distribution |
| 10 | ParobGauss | Transverse parabolic and longitudinal Gaussian distribution |
| 15 | SemicirGauss | Transverse semi-circle and longitudinal Gaussian distribution |
| 16 | Read | Read in an initial particle distribution from 6 column file **Partcl.data** (old version) |
| 166 | Read | Read in an initial particle distribution from 9 column file **Partcl.data** (new version including q/m, Q/Np, id) |
| 167 | Read | Read in an initial particle distribution from 9 column file **Partcl.data** of ImpactZ particle distribution output (first line of the file should contain: # of pts, gamma0, scaling length used in ImpactZ code) |
| 24 | readParmela | Read in Parmela particle format |
| 25 | readElegant | Read in Elegant particle format |
| 27 | CylcoldZSob | Uniform cylinder with longitudinal density modulation and Gaussian distribution in momentum space |
| ijk | Combine | Spatial transverse uniform or Gaussian (i=1,2), longitudinal flat-top with linear or Gaussian ramping or longitudinal Gaussian, or Tukey function (j=1,2,3,4), and 3D momentum distribution (k=1,2) |

The detailed input parameters for each type of initial distribution are defined in lines 6–8. See the description of Line 6–8 for more information.

**Rstartflg** — If `restartflag = 1`, restart the simulation from the previous check point. If `restartflag = 0`, start the simulation from the beginning.

**Flagsbstp** — Not used.

**Nemission** — There is a time period where the laser is shining on the cathode and electrons are being emitted. `Nemission` gives the number of numerical emission steps. More steps gives more accurate modeling but the computation time varies linearly with the number of steps. If `Nemission < 0`, there will be no cathode model. The particles are assumed to start in a vacuum.

**Temission** — Laser pulse emission time (sec). *Note: this time needs to be somewhat greater than the real emission time in the initial longitudinal distribution so that the time step size is changed after the whole beam is a few time steps out of the cathode.*

---

## Lines 6–8: Initial Distribution Parameters

These three lines give the initial distribution parameters in the x-px plane, y-py plane, and z-pz plane. These lines, along with `distType` as given in line 5, are used to form the initial distribution.

There are 21 parameters for the initial distribution except the `read` option, which reads in the particle data directly from the external file. The parameter names are:

```
sigx(m),  sigpx, muxpx, xscale,  pxscale, xmu1(m), xmu2,
sigy(m),  sigpy, muypy, yscale,  pyscale, xmu3(m), xmu4,
sigz(m),  sigpz, muzpz, zscale,  pzscale, xmu5(m), xmu6
```

**WARNING:**
1. `sigx`, `sigy`, `sigz` cannot be zeros even for read-in distribution.
2. `xmu6` is used to set beam energy.
3. `xscale` and `pxscale` are multipliers for `sigx` and `sigpx`. They are used when the `sigx` and `sigpx` numbers are derived from a match to the lattice and `xscale` and `pxscale` are then used to mismatch the initial distribution. If not used, `xscale` and `pxscale` should just be set to 1.
4. `sigz` is related to the laser pulse length through `v0 * t_laser` in some distributions, where `v0` is obtained from the kinetic energy `Bkenergy`.
5. For the `ijk` distribution, the `zscale` needs to be slightly greater than zero even for a hard edge cylinder beam.
6. `xmu1`–`xmu6` define the center offset of the initial distribution.

The phase space density is generally written as:

```
ρ(x̃, p̃x, ỹ, p̃y, z̃, p̃z) = ρx(x̃, p̃x) ρy(ỹ, p̃y) ρz(z̃, p̃z)
```

with the phase space coordinates:

```
x̃  = x  - xmu1
p̃x = px - xmu2
```

(with similar equations for the other four coordinates)

### Distribution Types

**1 Uniform**

```
ρx(x̃, p̃x) = const   if |x̃| < sqrt(3)*σx / sqrt(1 - muxpx²) and |p̃x| < sqrt(3)*σpx*(1 - muxpx/sqrt(1-muxpx²))
             0        otherwise
```

where `σ̃x = sigx` and `σ̃px = sigpx`. For `muxpx = 0`, in `(x̃, p̃x)` phase space, the distribution is a rectangle of constant density.

**2 Gauss3**

The form of the density function is:

```
ρx(x̃, p̃x) ∝ exp(-fx/2)
```

where:

```
fx(x̃, p̃x) = x̃²/σ̃x² + 2x̃p̃x*muxpx/(σ̃x*σpx) + p̃x²/σ̃px²
```

with similar forms for ρy and ρz.

**3 Waterbag**

The form of the density function is:

```
ρx(x̃, p̃x) = const   if fx(x̃, p̃x) < sqrt(8)
             0        otherwise
```

where `fx` is given in Gauss3. The density is a constant inside of an ellipse.

**4 Semigauss** — Essentially a Waterbag distribution in `(x̃, ỹ, z̃)` space and a Gaussian distribution in `(p̃x, p̃y, p̃z)` space.

**5 KV3d** — δ(fx, fy)-function transversely and uniform longitudinally. Here the transverse distribution is a self-consistent solution of the Poisson-Vlasov equation for a coasting beam. Not a realistic model but good for diagnostic checks.

**10 ParobGauss** — `ρ(r) = 1 - (r/A)²` transversely (`A = sqrt(6)*σx,y`) and uniform longitudinally. Gaussian distribution in momentum space with `pz > 0`.

**15 SemicirGauss** — `ρ(r) = sqrt(1 - (r/A)²)` transversely (`A = sqrt(5)*σx,y`) and uniform longitudinally. Gaussian distribution in momentum space with `pz > 0`.

**16 Read** — Read distribution from an external file called `partcl.data`. The file has the format:

```
nptot
x, px, y, py, z, pz
```

where `nptot` is number of particles, and `x(m)`, `px/mc`, `y(m)`, `py/mc`, `z(m)`, `pz/mc` are six phase space coordinates. Note: these coordinates will be shifted by the centroid parameters defined in lines 6–8.

**27 CylcoldZSob** — Uniform cylinder with longitudinal density modulation in spatial, Gaussian distribution in transverse momentum space and a semi-Gaussian distribution in longitudinal momentum space. Here, the relative modulation amplitude is "muzpz", the modulation wavelength is "zscale (m)".

**ijk Combine** — Generating initial particle distribution based on the combination of transverse spatial distribution (2 types), longitudinal spatial distribution (3 types) and 3D momentum distribution (4 types). For example, `flagdist = 111` denotes type 1 from the transverse spatial distribution, type 1 from the longitudinal spatial distribution, and type 1 from the 3D momentum distribution. In this case, it denotes a transverse uniform ellipse, longitudinal flat-top with linear ramping in spatial, and 3d full Gaussian distribution in momentum space.

For `flagdist = ijk`:
- `i = 1`: transverse uniform ellipse
- `i = 2`: Gaussian with cut-off set by xscale and yscale
- `j = 1`: flat-top set by sigz with linear ramping set by zscale
- `j = 2`: flat-top set by sigz with 2sigma Gaussian ramping set by zscale
- `j = 3`: Gaussian with cut-off set by zscale
- `k = 1`: 3D Gaussian momentum
- `k = 2`: transverse Gaussian momentum, longitudinal semi-Gaussian `pz*exp(-pz²/sigpz²)`
- `k = 3`: x-ray emission model `f(E,θ,φ) ∝ E/(E+Wk)²*sin(2θ)`, kinetic energy `0 ≤ E ≤ Emax (eV)`, Wk is the work function in eV. Emax is set by the input parameter `sigpz` and Wk is set by the input parameter `pzscale`.
- `k = 4`: 3 step model given by Dowell et al. [2]

---

## Line 9: Bcurr, Bkenergy, Bmass, Bcharge, Bfreq, Tini

**Bcurr** — Beam current in Amps with total charge Q = Bcurr/Bfreq.

**Bkenergy** — Initial beam pseudo-kinetic energy in eV. WARNING: this one is used to calculate the drift velocity that is needed to convert the laser pulse length in seconds into the longitudinal bunch length in meters in the initial distribution and to drift the particle out of the wall. The real initial beam energy needs to be input from "xmu6" in the initial distribution or the particle data file for the readin distribution.

**Bmass** — Mass of the particles in eV.

**Bcharge** — Particle charge in units of proton charge.

**Bfreq** — Reference frequency in Hz.

**Tini** — Initial reference time in seconds.
