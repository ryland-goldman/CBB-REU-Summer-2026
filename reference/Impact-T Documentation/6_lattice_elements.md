# 6 Lattice Beam Line Elements

Line 10 and beyond describe the lattice to track through. Each line of lattice input represents one element. Elements may overlap longitudinally. The general form of a lattice line is:

```
Blength, Bnseg, Bmpstp, Btype, V1 ... V23
```

> *NOTE: even though the IMPACT-T code allows the overlap of element field, the maximum length of the computational domain is set by the starting location of the last element plus the length of the last element. The beam line element has to be arranged so that `zmin(i+1) >= zmin(i)`, where `i = 1, ..., N`.*

**All elements should be laid out in sequence following the starting location of the element.**

---

**Blength** — The longitudinal length of the element.

**Bnseg, Bmpstp** — Not used except for Type less than 0, i.e. the BPM.

**Btype** — Type of element. An integer specifying the type of element. See below for more details.

**V1 ... V23** — Element parameters. See below for more details.

---

The RF E field is parameterized by:

```
field = f(x,y,z) * cos(2π*f*t + θ₀)
```

Some of the element parameters are:

| Parameter | Description |
|-----------|-------------|
| **zedge** | Longitudinal position of the entrance face of the element. |
| **scale** | This can be used to scale the field amplitude. Normally set to 1. |
| **RF frequency** | Frequency of the RF field in Hz. |
| **theta0** | Initial phase in degree. |
| **file ID** | ID number of the external file to open to read in field data. |
| **radius** | Aperture radius. |
| **x misalignment error** | Used in Quadrupole, Multipole (sextupole, octupole, decapole) and SolRF beam line elements. |
| **y misalignment error** | Used in Quadrupole, Multipole (sextupole, octupole, decapole) and SolRF beam line elements. |
| **rotation error x** | Used in Quadrupole, Multipole (sextupole, octupole, decapole) and SolRF beam line elements. |
| **rotation error y** | Used in Quadrupole, Multipole (sextupole, octupole, decapole) and SolRF beam line elements. Here, the notation of "+y" means clockwise y. |
| **rotation error z** | Used in Quadrupole, Multipole (sextupole, octupole, decapole) and SolRF beam line elements. |

---

## BPM / Diagnostic Elements (Btype < 0)

**-1, -2, -3, -4, -5, -6, -7, -8, -9, -11, -12, -13, -15, -16, -17, -99: BPM** — Beam position monitor and diagnostic elements. Parameters: blength, Bnseg, Bmpstp, btype, V1, V2, V3, V4, ... V8.

> *NOTE: the location of the first "-1,-2,-4,-7,-9,-11,-12" line should be put after the initial beam distribution and aligned in sequence!!!*

- **btype = -1**: Kick the transverse beam centroid at given location V2(m) by x offset V3(m), Px (γβx) offset V4, y offset V5(m), Py (γβy) offset V6, z offset V7(m), and Pz (γβz) offset V8.

- **btype = -2**: Output particle phase-space coordinate information at given location V3(m) into filename `fort.Bmpstp` with particle sample frequency Bnseg. Here, the maximum number of phase-space files which can be output is 100. Here, 40 and 50 should be avoided since these are used for initial and final phase space output.

- **btype = -3**: Output particle phase-space and prepare restart at given location V3(m) into filename `fort.(Bmpstp+myid)`. Here, myid is processor id. On single processor, it is 0. If there are multiple restart lines in the input file, only the last line matters.

- **btype = -4**: Change the time step size from the initial Dt (secs) into V4 (secs) after location V3(m). The maximum number of time step change is 100.

- **btype = -5**: Switch the simulation from azimuthal symmetry to fully 3d simulation after location V3(m). This location should be set as large negative number such as "-1000.0" in order to start the 3D simulation immediately after the electron emission. If there are multiple such lines in the input file, only the last line matters.

- **btype = -6**: Turn on the wake field effects between location V3(m) and V4(m). If `Bnseg` is greater than 0, the longitudinal and transverse wake function will be read in from file `"fort.Bmpstp"`. If `Bnseg ≤ 0`, the code will use analytical wake function described as follows. For analytical wake functions, the wake function parameters (iris radius) a = V5(m), (gap) g = V6(m), (period) L = V7(m). The definition of these parameters can be found from SLAC-PUB-9663, "Short-Range Dipole Wakefields in Accelerating Structures for the NLC," by Karl L.F. Bane. For the backward traveling wave structure, the iris radius "a" has to be set greater than 100, gap "g" set to the initialization location of BTW. For backward traveling wave structures, the wakes are hardwired inside the code following the report: P. Craievich, T. Weiland, I. Zagorodnov, "The short-range wakefields in the BTW accelerating structure of the ELETTRA linac," ST/M-04/02. For `-10 < a < 0`, it uses the analytical equation from the 1.3 GHz Tesla standing wave structure. For `a < -10`, it assumes the 3.9 GHz structure longitudinal wake function. For external supplied wake function, the maximum number of data points is 1000. The data points are assumed uniformly distributed between 0 and V7(m). The V6 has to be less than 0. Each line of the `fort.Bmpstp` contains longitudinal wake function (V/m) and two transverse wake functions x and y (V/m/m).

- **btype = -7**: Merge the multiple bins into only one bin at given location V3(m) in order to save computing time. If there are multiple such lines in the input file, only the last line matters.

- **btype = -8**: Switch on/off the space-charge calculation at given location V3(m) according to the sign of V2 (> 0 on, otherwise off).

- **btype = -9**: Output slice-based information at given location V3(m) into file `"fort.Bmpstp"` using "Bnseg" slices.

- **btype = -11**: Collimate particles at given location V2(m) with transverse aperture defined by V3(m) (xmin), V4(m) (xmax), V5(m) (ymin), V6(m) (ymax), V7 (switch for rectangular ≤ 10 or round > 10 aperture V3).

- **btype = -12**: Apply instant 6x6 linear matrix (read in from `linearmap.in`) kick at given location V1(m).

- **btype = -13**: Include dielectric wake field implemented by Dianel Mihalcea. V1 = Z start for the dielectric structure. V2 = Z end for the dielectric structure. V3 = 0.0 if geometry is rectangular and 1.0 if cylindrical. V4 = inner radius or height. V5 = outer radius or height. V6 = dielectric relative permittivity (eps). V7 = Lx (length in x-direction). This parameter is ignored if cylindrical geometry (V3=1.0). V8 = number of modes associated with x-axis. For cylindrical symmetry only monopole modes (V8=0) are supported. V9 = number of modes associated with y-axis (ref: PRST-AB, 15, 081304, (2012)). A typical line looks like: `0.0 0 0 -13 0.02 0.12 0.0 2.5e-3 5.0e-3 4.0 10.0e-3 6.0 10.0 /`

- **btype = -15**: Switch on the direct point-to-point N-body calculation of the space-charge forces. The cut-off radius of particle is given by V3(m).

- **btype = -16**: Heat the beam (γβz) at V1 location by rms size V2.

- **btype = -17**: Rotate the beam with respect to z-axis at V1 location by V2 radians.

- **btype = -99**: Stop the simulation at given location V3(m).

---

## 0: DriftTube — Drift space

```
V1: zedge
V2: radius
```

---

## 1: Quadrupole — Quadrupole

```
V1:  zedge
V2:  quad gradient (T/m)
V3:  file ID
     If > 0 & < 100, then include fringe field (using Enge function) and
     V3 = effective length of quadrupole.
     If >= 100, use read-in file for quadrupole. The file name is rfdataV3.
     This file contains discrete gradient, its first and 2nd derivative (w.r.p z).
     The format is:
       # of data points, starting location, ending location
       Bg, Bg', Bg'', Bg'''
       1 0.0 0.0
       0.0 0.0 0.0 0.0
V4:  radius (m)
V5:  x misalignment error (m)
V6:  y misalignment error (m)
V7:  rotation error x     (rad)
V8:  rotation error y     (rad)
V9:  rotation error z     (rad)
     If V9 != 0, skew quadrupole
V10: rf quadrupole frequency    (Hz)
V11: rf quadrupole phase        (degree)
```

---

## 2: ConstFoc — 3D constant focusing beam line element

```
V1: zedge
V2: x focusing gradient: kx0^2
V3: y focusing gradient: ky0^2
V4: z focusing gradient: kz0^2
V5: radius
```

---

## 3: Sol — Solenoid

```
V1: zedge
V2: scale of B field
V3: file ID
V4: radius
V5: x misalignment error  (Not used)
V6: y misalignment error  (Not used)
V7: rotation error x      (Not used)
V8: rotation error y      (Not used)
V9: rotation error z      (Not used)
```

The discrete magnetic field data is stored in `1Tv3.T7` file. The read-in format of `1Tv3.T7` is:

```fortran
! the input range units are cm
read(14,*,end=33) tmp1, tmp2, tmpint
this%RminRft = tmp1/100.0
this%RmaxRft = tmp2/100.0
this%NrIntvRft = tmpint
! the input range units are cm
read(14,*,end=33) tmp1, tmp2, tmpint
this%ZminRft = tmp1/100.0
this%ZmaxRft = tmp2/100.0
this%NzIntvRft = tmpint
n = 0
50  continue
    read(14,*,end=77) tmp1, tmp2
    nn = n+1
    j  = (nn-1)/(this%NrIntvRft+1) + 1
    i  = mod((nn-1),this%NrIntvRft+1) + 1
    this%brdatat(i,j) = tmp1
    this%bzdatat(i,j) = tmp2
    n = n + 1
goto 50
77  continue
```

Here, `ZminRft = 0.0`, `RminRft = 0.0`.

---

## 4: Dipole — Dipole bending magnet element

In the IMPACT-T code, the bending magnet is described by four linear functions that characterize the pole face and starting and ending face of the fringe field. The coordinate system consists of Cartesian coordinates X-Z with origin located before the bending magnet. X is the horizontal direction and Z is the longitudinal direction. The particles are tracked in this X-Z Cartesian coordinate system with external B field given in these coordinates. The space-charge forces are calculated by transforming the particle coordinates into rotated coordinates of the bunch centroid. The four functions which characterize the faces are: `z1 = k1*x1+b1`, `z2 = k2*x2+b2`, `z3 = k3*x3+b3`, `z4 = k4*x4+b4`. The region between z1 and z2 is the fringe field region for entrance, z2 and z3 is the region of constant By field, z3 and z4 is the region of constant By field for exit.

```
V1:  zedge
V2:  x field strength (T)
V3:  y field strength (T)
V4:  file ID          (file ID to contain the geometry information of bend)
V5:  half of gap width (m)
V6:  x misalignment error  (Not used)
V7:  y misalignment error  (Not used)
V8:  rotation error x      (Not used)
V9:  rotation error y      (Not used)
V10: rotation error z      (Not used)
```

The input file **rfdatav4** contains 22 lines:
- Line 1: Switch flag for 1D CSR wakefield. CSR wakefield will be included for a value greater than 0.
- Line 2: γ of the beam.
- Lines 3–10: k1, b1(m), k2, b2(m), k3, b3(m), k4, b4(m) — geometric description of pole faces.
- Line 11: Twice the shift z0 (m) at the entrance fringe field region (width of entrance region).
- Line 12: Twice the shift z0 (m) at the exit fringe field region (width of exit region).
- Lines 13–20: 8 coefficients of the Enge function (same for entrance and exit).
- Line 21: Effective starting location along the arc of the bend in meter.
- Line 22: Effective ending location along the arc of the bend in meter.

> *Note: the length (Blength) in the dipole element input line should contain the arc length of the reference particle inside the fringe fields and the dipole field region.*

Given the geometry of the bending magnetic field, the vertical field strength V3 should be adjusted so that the final bending angle attains the design value. The bending angle can be obtained from the output file `fort.38` that contains the reference particle trajectory. (That is `arctan(γβx/γβz)` in the last line of the `fort.38` in Output section.)

Also, see Section 8.9 for more description of the bending magnet model in the code.

---

## 5: Multipole — Multipole (see Quadrupole for Quadrupole)

```
V1:  zedge
V2:  pole type: sextupole (2), octupole (3), decapole (4)
V3:  field strength (T/m^n)
V4:  file ID        (If > 1e-5 then read in the fringe field)
V5:  radius
V6:  x misalignment error
V7:  y misalignment error
V8:  rotation error x
V9:  rotation error y
V10: rotation error z
```

---

## 101: DTL — Drift tube linac

```
V1:  zedge
V2:  scale
V3:  RF frequency
V4:  theta0
V5:  file ID
V6:  radius
V7:  quad 1 length
V8:  quad 1 gradient
V9:  quad 2 length
V10: quad 2 gradient
V11: x misalignment error for Quad 1
V12: y misalignment error for Quad 1
V13: rotation error x for Quad 1
V14: rotation error y for Quad 1
V15: rotation error z for Quad 1
V16: x misalignment error for Quad 2
V17: x misalignment error for Quad 2
V18: rotation error x for Quad 2
V19: rotation error y for Quad 2
V20: rotation error z for Quad 2
V21: x misalignment error for RF cavity
V22: y misalignment error for RF cavity
V23: rotation error x for RF cavity
V24: rotation error y for RF cavity
V25: rotation error z for RF cavity
```

---

## 102: CCDTL — Coupled-cavity-drift-tube-linac

See 104:SC for more details.

---

## 103: CCL — Coupled-cavity-linac

See 104:SC for more details.

---

## 104: SC — Superconducting cavity

```
V1:  zedge  (the real used field range in z is [zedge, zedge+Blength])
V2:  scale of RF field
V3:  RF frequency
V4:  theta0
V5:  file ID
V6:  radius
V7:  x misalignment error
V8:  y misalignment error
V9:  rotation error x
V10: rotation error y
V11: rotation error z
```

The `rfdataV5` file contains the Fourier coefficients for RF field. Example `rfdataV5` file:

```
1.2636782940838925  /Fourier coefficients of Ez on axis.
0.40391363899146518
-0.15415128679863189E-08
-0.14862013907801718
-0.46080547294426235E-0
0.86398401467898108E-01
0.36951164554240630E-09
```

---

## 105: SolRF — Solenoid with embedded RF field

```
V1:  zedge  (the real used field range in z is [zedge, zedge+Blength])
V2:  scale of RF field
V3:  RF frequency
V4:  theta0
V5:  file ID
V6:  radius
V7:  x misalignment error
V8:  y misalignment error
V9:  rotation error x
V10: rotation error y
V11: rotation error z
V12: scale of solenoid B field  [Only used with SolRF element]
```

The `rfdataV5` file contains Fourier coefficients for both E fields and B fields when `V5 ≤ 1000`. The first half contains E fields, the second half contains B fields. Example `rfdataV5` file:

```
7.0   /# of Fourier coef. of Ez on axis
-0.0616710037  /distance before the zedge.
0.0616710037   /distance after the zedge.
0.1233420074   /length of the Fourier expanded field.
 1.2636782940838925  /Fourier coefficients of Ez on axis.
 0.40391363889914653
 ...
7.0  /# of Fourier coef. of Bz on axis.
 0.9803477058479955E-04
 ...
```

Even though there is only one type of field (RF Ez or Bz), one still needs to keep five lines for the 0 field.

When `v5 > 1000`, the discrete field on-axis input format is:

```
# of data points for Ez, zstart, zend (in meters, w.r.p element edge)
Ez(z), Ez', Ez'', Ez'''
# of data points for Bz, zstart, zend
Bz(z), Bz', Bz'', Bz'''
```

---

## 110: EMfld — Read in discrete EM field data

Read in discrete EM field data as a function of (x,y,z) or (r,z) or analytical representation of EM field data. This element is not used in the IMPACT-T code but used in the z-based IMPACT code.

---

## 111: EMfldCart — Read in discrete EM field data on Cartesian grid

Read in discrete EM field data `Ex(V/m)`, `Ey(V/m)`, `Ez(V/m)`, `Bx(T)`, `By(T)` and `Bz(T)`, as a function of (x,y,z). Each field data point is a **complex number** to represent both traveling wave and standing wave fields. For standing wave field, the electric field is the real part of the complex number (E,0), the magnetic field is the imaginary part of the complex number (0,B). For a static electric field, it can be written as (E,0) with 0 frequency and phase. For a static magnetic field, it can be written as (0,B) with 0 frequency and -90 degree phase.

```
V1:  zedge (meter)
V2:  scale of RF field
V3:  RF frequency
V4:  theta0
V5:  file ID
V6:  radius
V7:  x misalignment error
V8:  y misalignment error
V9:  rotation error x
V10: rotation error y
V11: rotation error z
```

The discrete field data are stored in the file `1Tv3.T7`. Read-in FORTRAN format:

```fortran
! the input range units are m
read(14,*,end=33) tmp1, tmp2, tmpint
this%XminRfgt = tmp1
this%XmaxRfgt = tmp2
this%NxIntvRfgt = tmpint  !# of grid points -1
read(14,*,end=33) tmp1, tmp2, tmpint
this%YminRfgt = tmp1
this%YmaxRfgt = tmp2
this%NyIntvRfgt = tmpint  !# of grid points -1
read(14,*,end=33) tmp1, tmp2, tmpint
this%ZminRfgt = tmp1
this%ZmaxRfgt = tmp2
this%NzIntvRfgt = tmpint  !# of grid points -1
```

Data columns: `Ex`, `Ey`, `Ez` (complex V/m), `Bx`, `By`, `Bz` (complex T).

---

## 112: EMfldCyl — Read in discrete EM field data on cylindrical grid

Read in discrete EM field data `Ez(MV/m)`, `Er(MV/m)`, and `Hθ(A/m)` as a function of (r,z) from SUPERFISH output.

```
V1:  zedge
V2:  scale
V3:  RF frequency
V4:  theta0
V5:  file ID
V6:  radius
V7:  x misalignment error  (not used yet)
V8:  y misalignment error  (not used yet)
V9:  rotation error x      (not used yet)
V10: rotation error y      (not used yet)
V11: rotation error z      (not used yet)
```

The discrete field data is stored in `1Tv3.T7` file (units in cm). Data columns: `Ez`, `Er`, `|E| (MV/m)`, `H (A/m)`.

---

## 113: EMfldAna — EM field data as an analytical function

EM field data as an analytical function. The following three analytical functions are test functions.

```
V1:  zedge
V2:  escale
V3:  RF frequency
V4:  theta0
V5:  file ID
V6:
V7:
V8:
V9:
V10:
```

Analytical function options (selected by file ID):

**-1: Alpha magnet field:**
```
bx = 0.0d0
by = escale * pos(3)
bz = escale * pos(2)
```

**-2: Traveling wave field in meander plates:**
```
t0     = V3  !wave starting time
tt     = pos(4) - t0
tt1    = V4 - t0  !wave ending time
z0     = V6  !wave starting location
zslope = V7  !wave rising length
vzspeed = V8  !wave z speed
```

**-3: DC surface roughness field:**
```
ww  = V3 * 2 * pi
an  = V7
wkn = V8
ex  = escale*an*wkn*exp(-wkn*(pos(3)-zedge)) * sin(wkn*pos(1))
ey  = 0.0d0
ez  = escale + escale*an*wkn*exp(-wkn*(pos(3)-zedge)) * cos(wkn*pos(1))
```

**10: Dipole field for steering:**
- half gap = V6
- field profile in `rfdataV5` has the same format as the dipole element
- deflection switch = V7: vertical deflection for `int(V7) > 0`, otherwise horizontal deflection
