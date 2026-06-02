# BMAD Lattice Language Guide

BMAD (aka "Baby MAD" or "Better MAD" or just plain "BE MAD!") is a set of subroutines to read in lattice specification files that conform to the BMAD standard (similar to the MAD standard), compute Twiss parameters, track particles, etc. The subroutines are written in Fortran90 and have been developed to:

- Cut down on the time needed to develop programs,
- Cut down on programming errors, and
- Provide a standard input format for specifying lattices.

The following assumes familiarity with the MAD standard input format.

*Created by: DCS | Revision Date: March 3, 2003 | Topic revision: r2 - 23 Oct 2013*

---

## BMAD Elements

### Standard Elements from MAD recognized by BMAD

```
Element         Attributes (lowercase = attribute not in MAD)
-----------     --------------------------------------------
BEAMBEAM        TYPE, SIG_X, SIG_Y, CHARGE, sig_z, n_slice
DRIFT           TYPE, L
ECOLLIMATOR     TYPE, L, x_limit, y_limit
ELSEPARATOR     TYPE, L, gap
HKICKER         TYPE, L, KICK, TILT
INSTUMENT       TYPE, L
KICKER          TYPE, L, HKICK, VKICK, h_displace, v_displace
MARKER          TYPE
MONITOR         TYPE, L
MULTIPOLE       TYPE {KnL, Tn, n = 0 - 20}, tilt, radius
OCTUPOLE        TYPE, L, K3, TILT, b_gradient
QUADRUPOLE      TYPE, L, K1, TILT, b_gradient
RCOLLIMATOR     TYPE, L, x_limit, y_limit
RBEND           TYPE, L, ANGLE, E1, E2, K1, TILT, g, delta_g,
                         roll, rho, b_field
RFCAVITY        TYPE, L, VOLT, HARMON, phi0 (= LAG)
SBEND           TYPE, L, ANGLE, E1, E2, K1, TILT, g, delta_g,
                         roll, rho, b_field
SEXTUPOLE       TYPE, L, K2, TILT, b_gradient
SOLENOID        TYPE, L, KS, b_field
VKICKER         TYPE, L, KICK, TILT
```

### Nonstandard Elements in BMAD

```
AB_MULTIPOLE    TYPE, {An, Bn, n = 0 - 20}, radius, tilt
CUSTOM          TYPE, L, VAL1, VAL2, VAL3, ..., VAL12
GROUP           TYPE, COMMAND, OLD_COMMAND, COEF
LCAVITY         TYPE, L, GRADIENT, PHI0, RF_FREQUENCY,
                    delta_e, k_loss, e_loss, sr_wake_file, lr_wake_file
OVERLAY         TYPE, {See below}
SOL_QUAD        TYPE, L, K1, KS, TILT  ! Solenoid/Quad
TAYLOR          TYPE, {out_index: coef, e1 e2 e3 e4 e5 e6}
WIGGLER         TYPE, L, B_MAX, N_POLE, TILT,
                    {Bn, n = 2, 4, 6, 8}, RADIUS
WIGGLER         TYPE, L, polarity, z_patch,
                    term(i) = {coef, kx, ky, kz, phi_z}
PATCH           TYPE, X_OFFSET, X_PITCH, Y_OFFSET, Y_PITCH, Z_OFFSET,
                    DE_OFFSET, TILT
```

### Broken Elements (Can be resurrected if the need arises)

```
ACCEL_SOL       TYPE, L, VOLT, LAG, RF_WAVELENGTH, B_Z, B_X1, B_Y1,
                S_ST1, L_ST1, B_X2, B_Y2, S_ST2, L_ST2,
                X_BEG_LIMIT, Y_BEG_LIMIT  ! LINAC RF/Solenoid
```

### Nonstandard Elements in BMAD that cannot be defined in an input file

```
HYBRID          Created when elements are concatenated together.
                This is typically done for long-term tracking.
```

---

## Element Attribute Notes

### Multipole Attributes

The following elements can have `ab_multipole` attributes (A0, B0, etc.):

```
Element         Scale_Factor    Ref_Exponent
-----------     ------------    ------------
ELSEPARATOR     HKICK or VKICK  0
KICKER          HKICK or VKICK  0
OCTUPOLE        K3              L  3
QUADRUPOLE      K1              L  1
RBEND           G               L  0
SBEND           G               L  0
SEXTUPOLE       K2              L  2
SOLENOID        KS              L  1
SOL_QUAD        K1              L  1
```

Scaling formula: `An --> An * Scale_Factor * RADIUS^(ref_exp) / RADIUS^(n)`

If RADIUS is not present or is zero then RADIUS = 1 is used.

Example:
```
Q01W: QUAD, ... RADIUS = 0.06, B5 = 0.007
```

### Additional Attributes

All elements (except as noted) also have:

```
Attribute       Comment
---------       -------
ALIAS           16 character alternative name.
DESCRIP         200 character string for general use (use double quotes).
HKICK, VKICK    Transverse kicks. Except: Taylor, Patch, Marker.
X_OFFSET, Y_OFFSET  Transverse offsets. Except: Group, Marker, Overlay, Taylor.
S_OFFSET        Transverse offsets of the element. Except: Group, Marker, Overlay, Taylor.
X_PITCH, Y_PITCH    Yaw and pitch of the element. Except: Group, Marker, Overlay, Taylor.
X_LIMIT, Y_LIMIT    Aperture limits. Except: Group, Overlay, Patch.
```

Example:
```
my_quad: quadrupole, alias = "Type I", descrip = "any string", & hkick = 0.0023, ... etc.
```

- **X_OFFSET and Y_OFFSET**: Offset the element transversely with respect to the beam coordinate system.
- **S_OFFSET**: Moves the element longitudinally.
- **X_PITCH** (yaw): Rotates element about its center around the Y-axis (in the X-S plane).
- **Y_PITCH** (pitch): Rotates element about its center around the X-axis.
  - Sign convention: The point on the element (x, y, z) = (0, 0, 1) with respect to the center will, for a small Y_PITCH, end up at approximately (0, theta, 1) and for a small X_PITCH end up at approximately (theta, 0, 1).
  - Note: For MAD, the pivot point for misalignments is the start of the element, not the center. MAD defines positive rotation using a right hand rule which makes Y_PITCH the negative of the MAD rotation DPHI.
- **HKICK and VKICK**: Kicks along the X-axis and Y-axis respectively, independent of the value of TILT.
- **ALIAS**: An alternative name. Like TYPE, not used by any BMAD routine but can be used by a program to identify elements.
- **If X_LIMIT (or Y_LIMIT) is zero**: The convention is to assume that X_LIMIT (or Y_LIMIT) is infinite.
- **For ELSEPARATOR**: The sign of the kicks is for positrons. If the particle type is set to electrons the sign of the kicks will be reversed.

### Attribute Switches

All elements except groups and overlays have the following switches:

```
Attribute           Default         Comment
---------           -------         -------
MAT6_CALC_METHOD    BMAD_standard   For Jacobian calculations.
SYMPLECTIFY         False           Symplectify the Jacobian?
TRACKING_METHOD     BMAD_standard   For tracking particles.
NUM_STEPS           10              Number of steps for symplectic integration.
INTEGRATION_ORDER   2               Order of the PTC integrator (2, 4, or 6).
PTC_KIND            0               Do not set this unless you know what you are doing.
IS_ON               True            Set element On/Off (looks like a drift when off).
```

---

## Special Element Notes

### H_DISPLACE and V_DISPLACE

`H_DISPLACE` and `V_DISPLACE` displace a particle going through a kicker. Essentially this means that a kicker element can be used to change the coordinate system for the rest of the line.

### B_FIELD and B_GRADIENT Attribute

It is sometimes desirable to use the magnetic field strength itself as a variable parameter. This can happen where the beam is changing energy as it is in a LINAC. The `B_field` or `B_gradient` attribute is for this purpose.

### ELSEPARATOR Element

For an elseparator the kick is determined by HKICK and VKICK. GAP for an elseparator is used to compute the electric field for a given kick. The voltage is a dependent attribute determined by:
```
Voltage (V) = 10^-9 * kick * E [Gev] * gap [m] / L [m]
```

### BEAMBEAM Element

- In BMAD, X_OFFSET and Y_OFFSET are used to offset the BEAMBEAM element instead of the MAD standard attributes XMA and YMA.
- For a crossing angle use X_PITCH and Y_PITCH (full crossing angle, not the half-angle).
- N_SLICE is the number of equal charge chunks which the strong beam is sliced. Default is N_SLICE = 1. For slicing you need a non-zero SIG_Z.
- `CHARGE = -1`: Opposite beam has the opposite charge (default). `CHARGE = +1`: Opposite beam has the same charge.

### BEND Elements

- In MAD and BMAD, "L" for an RBEND is defined differently from "L" for an SBEND. See the MAD manual for more details.
- Internally all RBENDs are converted to SBENDs. This will be reflected in any listings of the elements.
- The TILT for a bend is as in MAD with a pi/2 TILT being equivalent to a downward curving bend. Warning: there may be some bugs associated with this so use with caution.
- The ROLL for a bend is considered a misalignment error. Thus ROLL does NOT change the reference orbit.
- G is the design 1/rho used to construct the layout of the ring. G + Delta_G is the actual 1/bending_radius felt by a particle.
- In an input file you can specify either G or Angle but not both. The one not specified will be calculated through the equation `Angle / G = Length`.
- Note: To specify a bend give L and either G, Rho, or Angle.

### MULTIPOLE Elements

- The TILT attribute rotates the multipole as a whole. The rotation angle for the nth multipole is Tn + TILT.
- The magnetic field for the nth multipole of a MULTIPOLE is:
  ```
  B_y + i B_x = KnL * e^(-i(n+1)Tn) * r^n * e^(i n theta) / n!
  ```
  where (r, theta) are the coordinates of the observation point. In terms of kick given to a particle:
  ```
  -F_x + i F_y = const * (B_y + i B_x)
  ```
- In terms of an AB_MULTIPOLE, the components of a MULTIPOLE are:
  ```
  KnL = n! * sqrt(An^2 + Bn^2)
  Tan[(n + 1) * Tn] = -a_n / b_n
  ```

### AB_MULTIPOLE Elements

- The AB_MULTIPOLE is similar to the MULTIPOLE element except that the fields are defined in terms of Bn (normal) and An (skew) components instead of KnL (magnitude) and Tn (tilt). There is also a factor of n! present. The relationship between the two:
  ```
  Bn + i*An = KnL * Exp(-i(n+1)Tn) / n!
  ```
- The magnetic field for the nth multipole of an AB_MULTIPOLE is:
  ```
  B_y + i B_x = (Bn + i An) * (x + i y)^n
  ```
- The TILT attribute is the same here as for a MULTIPOLE element.

### WIGGLER Elements

There are two types of wigglers:

**Map type** (described using a magnetic field map): The field is given by a sum of terms:
```
B(x,y,z) = Sum_i term(i)
```
each term specified in the input file as:
```
term(i) = {coef, kx, ky, kz, phi_z}
```
Note that z = 0 is the beginning of the wiggler. C = polarity * coef.

If `ky^2 = kx^2 + kz^2`:
```
B_x = -C * (kx/ky) * sin(kx*x) * sinh(ky*y) * cos(kz*z+phi_z)
B_y = C * cos(kx*x) * cosh(ky*y) * cos(kz*z+phi_z)
B_z = -C * (kz/ky) * cos(kx*x) * sinh(ky*y) * sin(kz*z+phi_z)
```

**Periodic type** (specifying the maximum field and number of poles):
- `B_MAX`: Maximum magnetic field on the wiggler centerline.
- `N_POLE`: Number of poles. The period is then `L / (2 * N_POLE)`.
- There is an attribute K1 (calculated from B_MAX) which is the vertical focusing strength of the wiggler (a wiggler does not have any focusing in the horizontal).
- `Rho_bend`: The bending radius at maximum field strength B_MAX. Calculated by BMAD.
- The multipole attributes B2, B4, B6, and B8 are for the individual poles. Since the strength of the poles alternate in sign, these attributes are not the same as the integrated multipole attributes of other elements.
- The tracking is designed so that, without multipole components, a particle entering at the origin (x=x'=y=y'=0) leaves at the origin.

### TAYLOR Elements

Taylor elements are specified by 6 Taylor series. The first gives x, the second Px, the third y, etc. Each term in a Taylor series is given in the input file by:
```
{out_index: coef, e1, e2, e3, e4, e5, e6}
```
`out_index` is a digit 1 through 6. The effect of the term is:
```
out = coef * x(in)^e1 * Px(in)^e2 * ... * Pz(in)^e6 + ... other terms
```
There are 6 default non-zero terms that make the Taylor series initially look like the unit matrix.

### CUSTOM Elements

`CUSTOM` is for programmer-defined elements. See the BMAD Programmers Guide for more details.

### PATCH

A PATCH element shifts the local reference coordinates while leaving the beam fixed.

### RFCAVITY Elements

- `PHI0` is the phase lag of the reference particle (in radians/2π). To preserve compatibility with MAD this can also be referred to as `LAG`.

### LCAVITY Elements

LCAVITY is a LINAC RF Accelerating Cavity.
- `GRADIENT`: Accelerating gradient in eV/m for a particle on crest (maximal acceleration).
- `PHI0`: Phase (in units of 2*pi) of the reference particle. Zero phase corresponds to maximal (on crest) acceleration.
- `RF_FREQUENCY`: RF Frequency in Hz.

### ACCEL_SOL Elements

An ACCEL_SOL element is a combination RF cavity and solenoid with two non-overlapping sets of steerings.
- Because a harmonic number is meaningless in a linac, `RF_WAVELENGTH` is specified rather than `HARMON`. VOLT and LAG are defined analogously to RFCAVITY.
- `B_Z` is the value of the solenoid magnetic field. If B_Z > 0, the field points downstream.
- `B_X1` and `B_Y1` are the transverse magnetic fields of the first steering set in the x- and y-directions. `B_X2` and `B_Y2` are the values for the second (further downstream) steering set.
- `S_ST1` (`S_ST2`): Position of the upstream end of the first (second) steering set relative to the beginning of the ACCEL_SOL.
- `L_ST1` (`L_ST2`): Length of the first (second) steering set.
- `X_BEG_LIMIT` and `Y_BEG_LIMIT` are the upstream aperture limits.

---

## BMAD Units

MKS units are used:

```
Angles              Radians     (Except phase angles)
Phase Angles        Radians/2π
Frequency           Hz          !!!  NOTE: Not standard MAD
Current             Amps
Kick                Radians
Length              meters
Magnetic Field      Tesla
Particle Energy     eV          !!!  NOTE: Not standard MAD
Voltage             Volts       !!!  NOTE: Not standard MAD
```

### TYPE and ALIAS Attributes

The TYPE and ALIAS attributes can be up to 16 characters, may contain spaces, and need to be put in double quotation marks. For Example:
```
q01: quadrupole, type = "mark1"
```

---

## Tracking & Jacobian (Mat6) Switches

There are several ways to track through an element or to get the linear 6×6 transfer matrix (Jacobian) about some orbit. This can be set in the input file:

```
TRACKING_METHOD = <tracking_switch>
MAT6_CALC_METHOD = <mat6_calc_switch>
```

Example:
```
Q02W: QUADRUPOLE, L= 0.6, K1 = 0.98, &
        TRACKING_METHOD = RUNGE_KUTTA, MAT6_CALC_METHOD = TAYLOR
```

### Possible TRACKING_METHOD Switches

**`Adaptive_Boris`**
Second order Boris integration with adaptive step size control. This is symplectic but slow.

**`BMAD_Standard`**
Quick and dirty (not necessarily symplectic). Appropriate when only interested in single-turn stuff (not long-term tracking). Does an exact calculation through sector bends using thin quads at either end for non-sector focusing. Sextupoles and octupoles are tracked using a single kick-drift-kick integration.

**`Custom`**
Calls a routine `track1_custom` which can be supplied by the user. The default `track1_custom` supplied with the BMAD release does Runge-Kutta tracking.

**`Linear`**
Just uses the 0th order vector with the 1st order 6×6 transfer matrix for an element. Very simple. May or may not be symplectic depending upon how the transfer matrix was generated.

**`Runge_Kutta`**
4th order Runge-Kutta integration algorithm with adaptive step size control (essentially ODEINT adopted from Numerical Recipes). May be slow but it should be accurate. This method is non-symplectic.

**`Symp_Lie_BMAD`**
Symplectic tracking using a Hamiltonian with Lie operation techniques. Similar to `Symp_Lie_PTC` except uses a BMAD routine. The difference: PTC tries to do things correctly while BMAD goes for speed by making approximations like the small angle approximation. Right now only implemented for Wigglers.

**`Symp_Lie_PTC`**
Symplectic tracking using a Hamiltonian with Lie operation techniques using Etienne's PTC software. This method is symplectic but can be slow. This method can only be used on elements that have a Hamiltonian (Quadrupoles, Solenoids, and most other element types do; Hybrid elements do not).

**`Symp_Map`**
Uses an implicit (partially inverted) Taylor map. Since the map is implicit a Newton search method must be used. This will slow things down from the Taylor method but this is guaranteed symplectic.

**`Taylor`**
Uses a Taylor map generated from Etienne's PTC package. Generating the map may take time but once you have it it should be very fast. One possible problem: you have to worry about the accuracy if you do tracking at points far from the point about which the series was made. This method is non-symplectic.

**`Wiedemann`**
Wiedemann's hard edge model of a wiggler.

### Possible MAT6_CALC_METHOD Switches

**`BMAD_Standard`**
Quick and dirty. Tries to be symplectic but this is not guaranteed. Sextupoles and octupoles are done using a simple kick-drift-kick model.

**`Custom`**
Calls a routine `make_mat6_custom` which may be supplied by the user. The default `make_mat6_custom` supplied by the BMAD release will use Runge-Kutta tracking.

**`Runge_Kutta`**
Uses a Runge-Kutta tracking algorithm adopted from Numerical Recipes. This method tracks 6 particles around the central orbit. Non-symplectic, susceptible to inaccuracies caused by nonlinearities. The advantage is that it relies only on the field map for the element.

**`Symp_Lie_BMAD`**
Symplectic integration using a Hamiltonian and Lie operators. Right now only implemented for Wigglers.

**`Symp_Lie_PTC`**
Symplectic integration using a Hamiltonian and Lie operators using Etienne's PTC software. Symplectic but can be slow. Only usable on elements that have a Hamiltonian.

**`Taylor`**
Uses a Taylor map generated from Etienne's PTC package. Non-symplectic. See notes on TAYLOR method above for tracking.

**`Tracking`**
Uses the tracking method set by TRACKING_METHOD to track 6 particles around the central orbit. Non-symplectic, susceptible to inaccuracies caused by nonlinearities. The advantage is that it is directly related to any tracking results.

### Symplectic Integration

Symplectic integration integrates the Hamiltonian H(y) where y here could be a 6-dimensional vector (for tracking) or be a Taylor series (for the mat6 calculation). The order at which a Taylor series is truncated is set by `TAYLOR_ORDER` (a global variable). In BMAD (or more precisely Etienne's PTC) you can use one of 3 methods, set by `INTEGRATION_ORDER`:

- `INTEGRATION_ORDER = n` (n = 2, 4, or 6) means that the error scales as `dz^n` where `dz` is the integration step size. The step size `dz` is set by the length of the element and the value of `NUM_STEPS`.
- Higher order does not necessarily imply higher accuracy.

`PTC_KIND` sets how the symplectic integrator divides up the Hamiltonian. Default `PTC_KIND = 0` means BMAD will choose what it thinks is best (see the routine `BMAD_ELE_TO_FIBRE` for more details).

---

## Superimpose/Overlay/Group

An element may be superimposed on a line using the `SUPERIMPOSE` switch. Additionally, `OVERLAY` and `GROUP` elements can be defined as explained below.

---

## SUPERIMPOSE

Possible switches to use with SUPERIMPOSE:

```
SUPERIMPOSE              Marks an element to be superimposed
REF = <element_name>     Reference for placement
OFFSET = <num>           Distance to offset superimposed element

REF_BEGINNING
REF_CENTER               Reference point to use on the reference element.
REF_END

ELE_BEGINNING
ELE_CENTER               Reference point to use on the superimposed element.
ELE_END

COMMON_LORD              Make 1 Lord element for multiple superpositions
```

The `SUPERIMPOSE` switch places an element in a line without having to worry that the element is "on top of" other elements. Example:

```
Q1: QUAD, L = 0.6, SUPERIMPOSE, REF = D03, REF_END, ELE_BEGINNING, OFFSET = 0.2
MY_LINE: LINE = ( ... , DO3, ... )
```

The element Q1 is superimposed on MY_LINE. The beginning of Q1 is positioned 0.2 meters further along than the end of D03.

An element may be superimposed simultaneously at different locations using wild card characters `"*"` and `"%"` ("multiple superposition"). Example:

```
MM: MULTIPOLE, ...., SUPERIMPOSE, COMMON_LORD, REF_END, REF = Q*
```

In the above example a multipole is put after each element that starts with the letter "Q". The `COMMON_LORD` switch is used to gang all the multipoles together using a single controller ("super_lord") element. `COMMON_LORD` can **ONLY** be used with MULTIPOLEs and AB_MULTIPOLEs.

### Notes

- A superimposed element cannot also be named in the line list.
- Defaults for a superimposed element:
  - `REF = {beginning of ring (s = 0)}`
  - `OFFSET = 0`
  - `REF_CENTER`
  - `ELE_CENTER`
- If a superimposed element extends past either end of the line then it is "wrapped around" (the line is always considered to be a closed ring).
- When an element is superimposed, BMAD forms composite ("slave") elements used for tracking. The values for the attributes of the slave elements are the sum of the attributes of the elements that are used in the forming process (the "lord" elements). Exceptions:
  - The LIMITs of a slave element are the minimum of the LIMITs of the lord elements whose LIMITs are non-zero.
  - The KICKs of a lord element are divided among the slave elements to keep the total kick the same. The division of the kicks is proportional to the slave elements' lengths.

---

## Overlay

An `OVERLAY` element is used to control a single attribute of other elements. Example:

```
OVER1: OVERLAY = {A_ELE, B_ELE/2.0}, HKICK = 0.003
OVER2: OVERLAY = {B_ELE}, HKICK
OVER2[HKICK] = 0.9
A_ELE: QUAD, HKICK = 0.05, ...
B_ELE: RBEND, ...
THIS_LINE: LINE = ( ... A_ELE, ... B_ELE, ... )
USE, THIS_LINE
```

In the example, OVER1 controls the HKICK attribute of A_ELE and B_ELE. OVER2 controls the HKICK attribute of B_ELE. OVER1 has a HKICK value of 0.003 and OVER2 has been assigned a value for HKICK of 0.9.

The resulting attribute values:
```
A_ELE[HKICK] = OVER1[HKICK] = 0.003
B_ELE[HKICK] = OVER2[HKICK] + 2 * OVER1[HKICK] = 0.906
```

There are coefficients associated with the control of a slave element. The default coefficient is 1.0. To specify a coefficient use a slash "/" after the element name followed by the coefficient. In the above example the coefficient for the control of B_ELE from OVER1 is 2.0.

An overlay will control all elements of a given name. Thus, in the above example, if there are multiple elements in THIS_LINE with the name B_ELE then the OVER1 and OVER2 overlays will control the HKICK attribute of all of them.

**Note**: Overlays completely determine the value of the attributes that are controlled by the overlay. In the above example, the HKICK of 0.05 assigned directly to A_ELE is overwritten by the overlay action of OVER1.

**Note**: The default value for an overlay is 0.

The `"\"` character can be used to distinguish between elements in a line list that have the same common specification:

```
H04: OVERLAY = {B05\H04}, HKICK
B05: RBEND, ...
THIS_LINE: LINE = ( ... B05, ... B05\H04, ... B05, ... )
```

The element definition for an element whose name has a `"\"` is obtained by neglecting the part of the name after the `"\"`. In this example the definition of B05\H04 is that of B05.

More than one type of attribute can be controlled by an overlay:
```
OV1: OVERLAY = {Q01W, ... , H_SEP_O8W[HKICK]/0.34, ...}, K1
```

In the above example the OV1 overlay generally controls the K1 attribute except for the element H_SEP_08W where the HKICK attribute is controlled.

---

## Group

`GROUP` is like `OVERLAY` in that a GROUP element controls the attribute values of other "slave" elements. Unlike an OVERLAY, however, a GROUP element is used to control *changes* (deltas), not the absolute value. In addition, a group element can control an element's position and length using the special attributes `ACCORDION_EDGE`, `SYMMETRIC_EDGE`, `START_EDGE`, and `END_EDGE`. Example:

```
Q10: QUAD, L = ...
D1: DRIFT, L = ...
D2: DRIFT, L = ...

GR1: GROUP = {Q10}, S_OFFSET
GR1[OLD_COMMAND] = 0.4
GR1[COMMAND] = 0.6
GR1[COEF] = 100

GR2: GROUP = {Q10}, START_EDGE = 0.1

GQ: GROUP = {Q10, Q11/-1}, K1

THIS_LINE: LINE = (... D1, Q10, D2, Q11, ...)
```

In this example the group element GR1 can be used by a program to control the longitudinal position of Q10. Similarly GR2 controls the placement of the starting edge of Q10 (the edge with the minimum s distance). In this case the lengths of D1 and Q10 are varied in such a way so that the total length of THIS_LINE is kept constant.

`ACCORDION_EDGE` varies the edges of an element so that the center of the element is fixed but the length varies. With `ACCORDION_EDGE` a change of, say, 0.1 in the group moves both edges of the element by 0.1 meters so that the length changes by 0.2 meters. `SYMMETRIC_EDGE`, `ACCORDION_EDGE`, `START_EDGE` and `END_EDGE` keep the total length of the lattice invariant.

### Group Notes

- Like overlays, coefficients can be specified for the individual elements under group control. In the above example, the GQ group controls the 2 quads in an anti-symmetric manner.
- Like overlays, groups can control more than one type of attribute.
- Unlike overlays, values are assigned to group elements using the `COMMAND` attribute. The `OLD_COMMAND` attribute sets the starting position for the group. In the above example, the effect of GR1 initially is to move the position of Q10 by 0.2 meters (= 0.6 - 0.4).
- The COEF attribute for a group has no meaning within BMAD and is used for communication with calling programs (for example, to define the conversion to Data Base computer units).
- When a lattice file is read in then COMMAND values for any groups are **always** applied last. This is independent of the order that they appear in the file. Example:

```
GR: GROUP = {Q1}, K1  GR[COMMAND] = 0.34  Q1[K1] = 0.57
```

In this example the value of Q1[K1] would be 0.91 = 0.34 + 0.57.

---

## Other MAD Stuff Recognized by BMAD

```
=                    ! variable assignment
!                    ! Comments!
LINE, -, *           ! Lines with reflection and repetition count,
                     !     also replacement lines are recognized.
LIST                 ! Replacement lists
USE, line_name       ! What line to use for generating the lattice
CALL, FILENAME = file_name  ! Include file
RETURN               ! return from a CALLed file
PI, TWOPI, E         ! Math constants
DEGRAD, RADDEG       ! Conversion factors
E_MASS, P_MASS       ! Mass in GeV
M_ELECTRON, M_PROTON ! Mass in eV (Not in MAD)
C_LIGHT              ! Speed of light.
R_E                  ! Classical electron radius.
E_CHARGE             ! Charge.
SIN, COS, TAN        ! Trig functions
ASIN, ACOS, ATAN     ! Inverse trig functions
ABS                  ! Absolute value (Not in MAD)
BEAM, ENERGY = GeV   ! Beam total energy (in GeV to be consistent with MAD)
BEAM, PARTICLE = name    ! name = electron or positron
BEAM, N_PART = real  ! number of particles in a bunch
```

---

## BMAD Stuff with No MAD Equivalent

The following parameters are recognized by BMAD:

```
lattice = name          ! the name of the lattice
parameter[lattice_type] = value  ! Type of lattice.
                         ! Possibilities are: CIRCULAR_LATTICE (default),
                         !   LINAC_LATTICE, LINEAR_LATTICE
parameter[symmetry] = value  ! symmetry of the ring, possibilities are:
                         !   NO_SYMMETRY (default),
                         !   MOBIUS_SYMMETRY, EW_ANTISYMMETRY, EW_SYMMETRY
parameter[taylor_order] = value  ! Order at which Taylor series are truncated.
parameter[beam_energy] = value   ! Beam total energy in eV.
```

Note: `"parameter[beam_energy]"` is the same as `"beam, energy"` except for the units.

### Starting Floor Position

```
beginning[x_position] = value     ! X position
beginning[y_position] = value     ! Y position
beginning[z_position] = value     ! Z position
beginning[theta_position] = value  ! Angle on floor
beginning[phi_position] = value    ! Angle of attack
```

### Starting Twiss Parameters

(Needed if the ring is not closed)

```
beginning[beta_x] = value     ! "a" mode beta
beginning[alpha_x] = value    ! "a" mode alpha
beginning[phi_x] = value      ! "a" mode phase
beginning[eta_x] = value      ! "a" mode dispersion
beginning[etap_x] = value     ! "a" mode dispersion derivative
beginning[beta_y] = value     ! "b" mode beta
beginning[alpha_y] = value    ! "b" mode alpha
beginning[phi_y] = value      ! "b" mode phase
beginning[eta_y] = value      ! "b" mode dispersion
beginning[etap_y] = value     ! "b" mode dispersion derivative
beginning[cij] = value        ! C coupling matrix. i,j = {"1" or "2"}
beginning[energy] = value     ! Total energy.
```

Note: `GAMMA_X`, `GAMMA_Y`, and `GAMMA_C` (the coupling gamma factor) will be kept consistent with the values set. If not set the default values are all zero.

### Other BMAD-Specific Syntax

- For DIMAD compatibility the `TITLE` command is recognized by BMAD. The next line after the `TITLE` command is used as a title.
- If a line begins with `END_FILE` then any lines thereafter in the file are ignored. This is exactly the same as using `RETURN`.
- If the input file contains the line `PARSER_DEBUG` then the parser will type out info on all variables, lines, and lists listed in the input file.
- Because of the time it takes to parse a file, `BMAD_PARSER` will save RING in a file with the name `'DIGESTED_' // IN_FILE`. For subsequent calls to the same `IN_FILE`, `BMAD_PARSER` will just read in the digested file. `BMAD_PARSER` will always check to see that the digested file is up-to-date.
- If the input file contains the line `NO_DIGESTED` then the parser will not make a digested file.
- The `"RETURN"` statement is optional in a called file. In fact `RETURN` and `END_FILE` statements are treated identically.

---

## Notes on BMAD

- **Variables and element attributes can be expressed using arithmetic expressions**:
  ```
  Q01W: Quadrupole, L = 2.3, K1 = 0.1
  X3 = 5.67 * Q01W[K1] / 4.5  ! Notice use of attribute value
  WIG1: WIGGLER, L = PI*X3 + 37  ! Define WIG1 and set its length
  WIG1[L] = 2.3^2  ! Redefine the length of WIG1
  Q01W[k1] = 0.3   ! This has no effect on the value of X3
  ```
  Equations are always "immediately evaluated" so the value of X3 is determined by the value of Q01W[k1] when the "X3 = ..." line is parsed.

- **Element classes** are recognized (Cf. MAD manual) so that the following example is acceptable in BMAD:
  ```
  GEN_BEND: RBEND, L=...   ! First define GEN_BEND
  B44E: GEN_BEND, E1=...   ! Then use GEN_BEND as a key name
  ```

- **No forward definitions of variables**. Do not use:
  ```
  b_var = 4.5 * a_var  a_var = 9
  ```
  Instead write:
  ```
  a_var = 9  b_var = 4.5 * a_var
  ```

- **No redefinition of lines, lists, or elements allowed**:
  ```
  abc: line = (...)  ! Define ABC
  abc: line = (...)  ! NO. Redef of ABC not allowed.
  ```

- **Element attributes, however, may be redefined**:
  ```
  q01w: quad, k1 = -.239, ...  ! define Q01W
  q01w[k1] = -0.493            ! OK. Redef of quad strength allowed.
  ```

- **Use of default TILT values** are accepted per the MAD standard. Example:
  ```
  sk05e: quad, tilt, k1 = ...  ! default tilt for quads = pi/4
  ```

- **Default BEAM parameters** are:
  ```
  particle = positron  energy = 0  ! GeV
  ```

- **BMAD is case insensitive**.

- **Use of "S" (longitudinal position) as an attribute** is allowed but only in "secondary" lattice files (read in after the "primary" lattice file is read in by BMAD_PARSER2):
  ```
  var = 2 * Q01W[S]  ! allowed in secondary lattice file only
  ```
