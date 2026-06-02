# BMAD Programming Guide

Helpful hints for using BMAD subroutines.

*Created by: David Sagan | Last updated: January 2003*

---

## Index

- [Basic Ring Structure](#basic-ring-structure)
- [Transfer Matrices Through Tracking](#transfer-matrices-through-tracking)
- [Custom & Runge Kutta Tracking](#custom--runge-kutta-tracking)
- [Etienne's FPP/PTC Software](#etiennes-fppppc-software)
- [Taylor Series Maps](#taylor-series-maps)
- [Custom Elements & Custom Calculations](#custom-elements--custom-calculations)
- [Finding Element Endpoints in BMAD Lattice](#finding-element-endpoints-in-bmad-lattice)
- [Runge Kutta Tracking](#runge-kutta-tracking)
- [Dependent Attributes](#dependent-attributes)
- [Allocation/Deallocation](#allocationdeallocation)
- [Module Dependencies](#module-dependencies)

---

## Basic Ring Structure *(06/2003 dcs)*

The ring structure is the basis upon which BMAD is built. It holds information about where elements of the ring (or transfer line or LINAC) are placed and their attribute values. Declare an instance via:

```fortran
type (ring_struct) ring
```

See `bmad_struct.f90` for definitions of all structures. The ring has an array of elements `ring%ele_(i)` where `i` runs from 0 through `ring%n_ele_max`.

- **Element 0**: Special element holding Twiss parameters at the ring start. Does not correspond to a physical entity.
- **Elements 1 to `ring%n_ele_ring`**: Regular "SLAVE and FREE" elements used in calculations. To track once around the ring, track from element 1 to element `ring%n_ele_ring`.
- **Elements `ring%n_ele_ring+1` to `ring%n_ele_max`**: LORD elements which control attributes of the regular elements.

The control type of an element is given by `ele%control_type`:

### Slave/Free Elements in the Regular Part of the Ring

| Type | May control | May be slave of |
|------|-------------|-----------------|
| `FREE` | None | `GROUP_LORD` |
| `OVERLAY_SLAVE` | None | `GROUP_LORD`, `OVERLAY_LORD` |
| `SUPER_SLAVE` | None | `SUPER_LORD` |

### Lord Elements in the Control Part of the Ring

| Type | May control | May be slave of |
|------|-------------|-----------------|
| `GROUP_LORD` | `GROUP_LORD`, `SUPER_LORD`, `FREE`, `OVERLAY_LORD`, `OVERLAY_SLAVE` | `GROUP_LORD` |
| `OVERLAY_LORD` | `OVERLAY_SLAVE`, `OVERLAY_LORD`, `SUPER_LORD` | `GROUP_LORD`, `OVERLAY_LORD` |
| `SUPER_LORD` | `SUPER_SLAVE` | `GROUP_LORD`, `OVERLAY_LORD` |

The control information is stored in the array `ring%control_(:)`. Each element of this array connects a lord with one of its slaves:

- `RING%CONTROL_(I)%IX_LORD` — Index of the Lord element.
- `RING%CONTROL_(I)%IX_SLAVE` — Index of the Slave element.
- `RING%CONTROL_(I)%IX_ATTRIB` — Index of the attribute controlled.
- `RING%CONTROL_(I)%COEF` — Coefficient.

`ring%control_(i)%ix_attrib` is used with overlay_lords and group_lords but not needed with super_lords (super_lords control all attributes of their slaves). For overlay_lords and group_lords, `ring%control_(i)%coef` gives how much the attribute of a slave changes per unit change in the lord attribute. For a super_lord, `ring%control_(i)%coef` is the ratio of the slave's length to the lord's length.

### Traversing Lord → Slave

```fortran
! Lord index is IX_LORD. Slave index is IX_SLAVE
do i = ring%ele_(ix_lord)%ix1_slave, ring%ele_(ix_lord)%ix2_slave
  ix_slave = ring%control_(i)%ix_slave  ! index of a slave
  ...
enddo
```

### Traversing Slave → Lord

```fortran
do i = ring%ele_(ix_slave)%ic1_lord, ring%ele_(ix_slave)%ic2_lord
  j = ring%ic_(i)
  ix_lord = ring%control_(j)%ix_lord
  ...
enddo
```

> **Note**: Elements controlled by a GROUP_LORD do not have the control information that they are so controlled. This is the only exception to the rule that control information between two elements is available starting from either element.

---

## Transfer Matrices Through Tracking *(03/2002 dcs)*

Subroutine chain for `mat6` when tracking is used for computation:

```
Note:  MCM == ele%mat6_calc_method
       TM  == ele%tracking_method
       Set: This is a temporary set. The original value
            is restored when the subroutine finishes.

                         make_mat6
        --------<------- [MCM == custom$]-------->--------
        |                [MCM == tracking$]               |
        |                [MCM == runge_kutta$]             |
        v                        v                        v
make_mat6_tracking        make_mat6_runge_kutta    make_mat6_custom
        |                (Set: TM = runge_kutta$)  (Set: TM = custom$)
        |                        v                        |
        ----------> transfer_matrix_through_tracking <----
                                 |
                                 v
                              track1
                                 |
                                 v
               See also: Section below on
               Custom & Runge Kutta Tracking
```

---

## Custom & Runge Kutta Tracking *(03/2002 dcs)*

Subroutine chain for `track1` for Custom and Runge Kutta tracking. See `track1.f90` for other tracking types. The default `track1_custom` routine supplied with the BMAD library uses Runge-Kutta tracking with the field supplied by `field_rk_custom`. `field_rk_custom` is not part of the BMAD library and must be supplied by you if needed.

```
Note: TM == ele%tracking_method

      track1
      [TM == custom$]-------->--------
      [TM == runge_kutta$]            |
              |                       v
              v                track1_custom
       track1_runge_kutta <---
              |
              v
          odeint_bmad
              |
              v
           rkqs_bmad
              |
              v
            rkck
              |
              v
         derivs_bmad
         [TM == custom$]------->-------
         [TM /= custom$]               |
              |                        v
              v               field_rk_custom
       field_rk_standard
```

---

## Etienne's FPP/PTC Software *(06/2003 dcs)*

The FPP/PTC (Full Polymorphic Package/Polymorphic Tracking Code) software package of Etienne Forest handles Taylor maps to any arbitrary order (also known as Truncated Power Series Algebra — TPSA) along with Lie Algebraic Operations. BMAD uses this software to compute Taylor maps and 6×6 transfer matrices.

FPP/PTC is a very general package and BMAD only makes use of a small part of its features. See the FPP/PTC documentation (starts on pg. 129 of the link) for more information.

---

## Taylor Series Maps *(03/2002 dcs)*

The order of the Taylor maps is set either from the BMAD input lattice file or by the program using the `set_taylor_order` routine.

Since it can take a while to compute the Taylor Maps, once a Taylor Map is made for an element it will never get discarded unless you kill it intentionally.

There are two Taylor series structures:
- **`taylor`**: Defined by Etienne's PTC/FPP
- **`taylor_struct`**: Defined by BMAD

The `taylor_struct` structure is:

```fortran
type taylor_term_struct
  real(rdef) :: coef
  integer :: exp(6)
end type

type taylor_struct
  real(rdef) ref
  type (taylor_term_struct), pointer :: term(:)
end type
```

In the `ele_struct` structure, a map is composed of 6 `taylor_struct`s, one for each coordinate (x, P_x, y, P_y, z, P_z):

```fortran
type ele_struct
  ...
  type (taylor_struct) :: taylor(6)
  ...
end type
```

> **Note**: In PTC/FPP the longitudinal components are reversed. Coordinates are (x, P_x, y, P_y, P_z, ct = -z). Conversion routines take this into account automatically.

The `term(:)` array defines the terms in the Taylor series. An individual term is:

```
out = coef * x^exp(1) * P_x^exp(2) * y^exp(3) ...
```

where "out" depends on which `taylor(i)` is chosen: i = 1 ==> "out" = x, etc.

The `ref` component indicates the reference point in phase space about which the Taylor series has been evaluated. For example, suppose the map through an element is quadratic:

```
out = a * x^2
```

The Taylor series evaluated to 1st order about `x_ref` is:

```
out = c1 + c2 * x
where:
  c1 = -a * x_ref^2
  c2 = 2 * a * x_ref
```

The Taylor series evaluated to 2nd or higher order about `x_ref` is:

```
out = a * x^2
```

Note: this is independent of `x_ref`. The coefficients of the Taylor series may depend on the reference point, but once evaluated the calculation of the output from the input phase space point is independent of the reference.

---

## Custom Elements & Custom Calculations *(03/2002 dcs)*

A Custom calculation involves linking in custom routines in a program. These routines are used when:

```
1) Dealing with a Custom element.
2) The TRACKING or MAT6_CALC switches for an element are set to CUSTOM.
```

There are 4 custom routines which are needed:

- `make_mat6_custom`
- `track1_custom`
- `custom_emitt_calc`
- `custom_radiation_integrals`

The `make_mat6_custom` and `track1_custom` routines supplied by BMAD use Runge-Kutta tracking with the field given by `field_rk_custom`. You need to supply this routine.

The `custom_emitt_calc` and `custom_radiation_integrals` routines supplied by BMAD will bomb the program if called.

See the dummy routines in the BMAD release for more details.

---

## Finding Element Endpoints in BMAD Lattice *(03/2000 dhr)*

The hierarchical structure of BMAD layouts (LORDS and SLAVES) can cause problems when finding the endpoint of a physical magnet. For example, permanent magnet quadrupoles (REQ's) end up split into two parts in the final ring layout.

In searching for a magnet or other physical device (by matching against `ring%ele_(i)%key` or `ring%ele_(i)%alias` or `ring%ele_(i)%name`) once a match has been made, check `ring%ele_(i)%control_type` for that member.

- If `ring%ele_(i)%control_type = free$` or `= overlay_slave$`: treat as an independent element and use its parameters directly.
- If `ring%ele_(i)%control_type = super_lord$`: use parameters in the element indexed by `iele = ring%control_(ix)%ix_slave` where `ix` runs from `ring%ele_(i)%ix1_slave` to `ring%ele_(i)%ix2_slave`.
- If `ring%ele_(i)%control_type = anything else`: ignore it and continue the search.

---

## Runge Kutta Tracking *(09/2001)*

Typical subroutine hierarchy for Runge-Kutta tracking:

```
Routine              Where From
---------            ----------
track_all            BMAD
  |
  v
track1               BMAD
  |
  v
track1_runge_kutta   Programmer supplied
  |
  v
track_runge_kutta    BMAD
  |
  v
odeint2              Numerical Recipes [Modified]
  |
  v
derivs               BMAD
  |
  v
field_rk             Programmer supplied
```

`Track1_runge_kutta` and `field_rk` must be supplied by the programmer. `Track1_runge_kutta` sets up the particular parameters for tracking through the current ring element and `field_rk` must be able to return the field as a function of position. See code files for the above routines for more details.

---

## Dependent Attributes *(DCS 01/2003)*

Some attributes of an element are designated as "dependent variables" which are dependent upon other independent variables:

```
               Dependent Variables        Independent Variables
               -------------------        --------------------
Rbend:         Rho, Angle, L_Cord         G, L
Sbend:         Rho, Angle, L_Cord         G, L
RFCavity:      RF_Wavelength              Harmon
BeamBeam:      BBI_Const                  Charge, Sig_x, Sig_y
Wiggler:       K1, Rho                    B_max
```

When `attribute_bookkeeper` is called (e.g., by `make_mat6`) the values of the dependent variables will be set based upon the values of the independent variables. Thus trying to vary the strength of a bend by varying, say, the Rho attribute is an exercise in futility.

---

## Allocation/Deallocation *(DCS 04/2002)*

Various structures have pointer elements which mean you, the programmer, have to be aware of allocation and deallocation. All BMAD pointers are solely used to point to arrays in the heap and must always be deallocated when going out of scope.

The `ele_struct` structure has pointers so the programmer is responsible for making sure there are no memory leaks. The following calls are used to allocate/deallocate:

```fortran
call init_ele (ele)
call deallocate_ele_pointers (ele)
```

Furthermore since the `ring_struct` holds an array of elements, the `ring_struct` also has pointers in it.

---

## Module Dependencies *(DCS 04/2002)*

*Created by: David Sagan | Maintained by: David Sagan | Last updated: January 2003*
*Topic revision: r3 - 23 Oct 2013, DraganaJusic*
