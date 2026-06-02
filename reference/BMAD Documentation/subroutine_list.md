# Library Subroutine List for BMAD

Subroutines sorted by functionality. Use `GETF` to get full header info from source:
```
$ GETF TWISS_PROP*
```
Skeleton templates: `[CESR.BMAD.CODE]SUBROUTINE_TEMPLATE.F90`

---

## Index

- [Reading/Writing a Lattice File](#readingwriting-a-lattice-file)
- [Choosing a Lattice](#choosing-a-lattice)
- [CESR Specific](#cesr-specific)
- [Twiss etc.](#twiss-etc)
- [Matrices](#matrices)
- [Routines called by MAKE_MAT6](#routines-called-by-make_mat6)
- [Low Level Matrix Routines](#low-level-matrix-routines)
- [Tracking, Closed Orbit](#tracking-closed-orbit)
- [Tracking Routines called by TRACK1](#tracking-routines-called-by-track1)
- [Low Level Tracking Routines](#low-level-tracking-routines)
- [Particle Coordinate Stuff](#particle-coordinate-stuff)
- [Ring Geometry](#ring-geometry)
- [Interface to PTC](#interface-to-ptc)
- [Taylor Map Routines](#taylor-map-routines)
- [Long Range Beam-Beam Interaction](#long-range-beam-beam-interaction)
- [Helper Subroutines: Informational](#helper-subroutines-informational)
- [Helper Subroutines: Elemental](#helper-subroutines-elemental)
- [Helper Subroutines: Transformational](#helper-subroutines-transformational)
- [Helper Subroutines: Multipolar](#helper-subroutines-multipolar)
- [Helper Subroutines: Miscellaneous](#helper-subroutines-miscellaneous)

---

## Reading/Writing a Lattice File

**`BMAD_PARSER (IN_FILE, RING)`**
Subroutine to parse a BMAD input file.

**`BMAD_PARSER2 (IN_FILE, RING)`**
Subroutine to parse (read in) a BMAD input file to modify an existing lattice.

**`READ_DIGESTED_BMAD_FILE (IN_FILE_NAME, RING, VERSION)`**
Subroutine to read in a digested file.

**`WRITE_DIGESTED_BMAD_FILE (DIGESTED_NAME, RING, N_FILES, FILE_NAMES)`**
Subroutine to write a digested file.

---

## Choosing a Lattice

**`CHOOSE_CESR_LATTICE (LATTICE, LAT_FILE, CURRENT_LAT, RING)`**
Subroutine to let the user choose a lattice. The subroutine will present a list to choose from.

**`GET_LATTICE_LIST (LAT_LIST, NUM_LATS, DIRECTORY)`**
Subroutine to get the names of the lattices of the form: `directory // BMAD_*.*`

---

## CESR Specific

**`BMAD_TO_CESR (RING, CESR)`**
Subroutine to transfer information from the RING structure returned from BMAD_PARSER to a structure for the CESR ring.

**`BMAD_TO_DB (RING, DB)`**
Subroutine to return information on the data base that pertains to CESR elements.

**`CHOOSE_CESR_LATTICE (LATTICE, LAT_FILE, CURRENT_LAT, RING)`**
Subroutine to let the user choose a lattice. The subroutine will present a list to choose from.

**`CREATE_VSP_VOLT_ELEMENTS (RING, ELE_TYPE)`**
Subroutine to create elements corresponding to the 6 data base elements in CSR VSP VOLT.

**`DB_GROUP_TO_BMAD (ING_NAME, ING_NUM, BIGGRP_SET, RING, CON_, N_CON, OK, TYPE_ERR)`**
Subroutine to take a data base group element and find the elements controlled along with the coefficients.

**`DB_GROUP_TO_BMAD_GROUP (GROUP_NAME, GROUP_NUM, I_BIGGRP, RING, IX_ELE, OK, TYPE_ERR)`**
Subroutine to set up a data base group knob in a bmad ring structure.

**`IDENTIFY_DB_NODE (DB_NAME, DB, DP_PTR, OK, TYPE_ERR)`**
Subroutine to find which array in DB is associated with DB_NAME.

**`K_TO_QUAD_CALIB (K_THEORY, ENERGY, CU_THEORY, K_BASE, DK_GEV_DCU, CU_PER_K_GEV)`**
Subroutine to return the calibration constants for the CESR quads.

**`LATTICE_TO_BMAD_FILE_NAME (LATTICE, BMAD_FILE_NAME)`**
Subroutine to convert a lattice name to the appropriate bmad file name.

**`QUAD_CALIB (LATTICE, K_THEORY, K_BASE, LEN_QUAD, CU_PER_K_GEV, QUAD_ROT, DK_GEV_DCU, CU_THEORY)`**
Subroutine to return the calibration constants for the CESR quads.

**`READ_BUTNS_FILE (BUTNS_NUM, BUTNS, DB, OK)`**
Subroutine to read in the information in a BUTNS.nnnnn file.

**`RING_TO_QUAD_CALIB (RING, CESR, K_THEORY, K_BASE, LEN_QUAD, CU_PER_K_GEV, QUAD_ROT, DK_GEV_DCU, CU_THEORY)`**
Subroutine to return the calibration constants for the CESR quads.

---

## Twiss etc.

**`CALC_Z_TUNE (RING)`**
Subroutine to calculate the synchrotron tune from the full 6×6 1-turn matrix.

**`CHROM_CALC (RING, DELTA_E, CHROM_X, CHROM_Y)`**
Subroutine to calculate the chromaticities by computing the tune change when the energy is changed.

**`CHROM_TUNE (RING, DELTA_E, TARGET_X, TARGET_Y, ERR_FLAG)`**
Subroutine to set the sextupole strengths so that the ring has the desired chromaticities.

**`EMITT_CALC (RING, WHAT, MODE)`**
Subroutine to calculate the emittance, energy spread, and synchrotron integrals.

**`MOBIUS_TWISS_CALC (ELE, V_MAT)`**
Subroutine to calculate the Mobius betas and etas which are effective projections of beta and eta in the X and Y planes.

**`QUAD_BETA_AVE (RING, IX_ELE, BETA_X_AVE, BETA_Y_AVE)`**
Subroutine to compute the average betas in a quad.

**`RADIATION_INTEGRALS (RING, ORB_, MODE)`**
Subroutine to calculate the synchrotron radiation integrals along with the emittance and energy spread.

**`RELATIVE_MODE_FLIP (ELE1, ELE2)`**
Function to see if the modes of ELE1 are flipped relative to ELE2.

**`SET_TUNE (PHI_X_SET, PHI_Y_SET, DK1, RING, ORB_, OK)`**
Subroutine to Q_tune a ring. Program will set the tunes to within 0.001 radian (0.06 deg).

**`SET_Z_TUNE (RING)`**
Subroutine to set the longitudinal tune by setting the RF voltages in the RF cavities.

**`TWISS_AND_TRACK (RING, ORB)`**
Subroutine to calculate the Twiss and orbit parameters. This is not necessarily the fastest routine.

**`TWISS_AND_TRACK_PARTIAL (ELE1, ELE2, PARAM, DEL_S, ELE3, START, END)`**
Subroutine to propagate partially through ELE2 the Twiss parameters and the orbit.

**`TWISS_AND_TRACK_BODY (ELE1, ELE2, PARAM, DEL_S, ELE3, START, END)`**
Subroutine to propagate partially through ELE2 the Twiss parameters and the orbit.

**`TWISS_AT_ELEMENT (RING, IX_ELE, START, END, AVERAGE)`**
Subroutine to return the Twiss parameters at the beginning, end or the average of an element.

**`TWISS_PROPAGATE_MANY (RING, IX_START, IX_END, DIRECTION)`**
Subroutine to propagate the Twiss parameters from one point in the ring to another.

**`TWISS_AT_S (RING, S, ELE)`**
Obsolete. Use `twiss_and_track_at_s` instead.

**`TWISS_AND_TRACK_AT_S (RING, S, ELE, ORB_, HERE)`**
Subroutine to calculate the Twiss parameters and orbit at a particular longitudinal position.

**`TWISS_AT_START (RING)`**
Subroutine to calculate the Twiss parameters at the start of the ring.

**`TWISS_FROM_TRACKING (RING, CLOSED_ORB_, D_ORB, ERROR)`**
Subroutine to compute from tracking, for every element in the ring, the Twiss parameters and the transfer matrices.

**`TWISS_PROPAGATE1 (ELE1, ELE2)`**
Subroutine to propagate the Twiss parameters from the end of ELE1 to the end of ELE2.

**`TWISS_PROPAGATE_ALL (RING)`**
Subroutine to propagate the Twiss parameters from the start to the end.

**`TWISS_TO_1_TURN_MAT (TWISS, PHI, MAT2)`**
Subroutine to form the 2×2 1-turn transfer matrix from the Twiss parameters.

---

## Matrices

**`C_TO_CBAR (ELE, CBAR_MAT)`**
Subroutine to compute Cbar from the C matrix and the Twiss parameters.

**`CLEAR_RING_1TURN_MATS (RING)`**
Clear the 1-turn matrices in the ring structure.

**`DO_MODE_FLIP (ELE, ELE_FLIP)`**
Subroutine to mode flip the twiss_parameters of an element.

**`MAKE_G2_MATS (TWISS, G_MAT, G_INV_MAT)`**
Subroutine to make the matrices needed to go from normal mode coords to coordinates with the beta function removed.

**`MAKE_G_MATS (ELE, G_MAT, G_INV_MAT)`**
Subroutine to make the matrices needed to go from normal mode coords to coordinates with the beta function removed.

**`MAKE_MAT6 (ELE, PARAM, C0, C1)`**
Subroutine to make the 6×6 transfer matrix for an element.

**`MAKE_MAT627 (ELE, PARAM, DIRECTION, MAT627)`**
Subroutine to make the 6×27 2nd order transfer matrix for an element.

**`MAKE_V_MATS (ELE, V_MAT, V_INV_MAT)`**
Subroutine to make the matrices needed to go from normal mode coords to X-Y coords and vice versa.

**`MAT6_TO_TAYLOR (MAT6, VEC0, BMAD_TAYLOR)`**
Subroutine to form a first order Taylor map from the 6×6 transfer matrix and the 0th order transfer vector.

**`MAT_INVERSE (MAT, MAT_INV)`**
Program to take the inverse of a square matrix.

**`MAT_SYMPLECTIFY (MAT_IN, MAT_SYMP)`**
Subroutine to form a symplectic matrix that is "close" to the input matrix.

**`MAT_SYMP_CHECK (MAT, ERROR)`**
Routine to check the symplecticity of a square matrix.

**`MAT_SYMP_DECOUPLE (T0, TOL, STAT, U, V, UBAR, VBAR, G, TWISS1, TWISS2, TYPE_OUT)`**
Subroutine to find the symplectic eigen modes of the one-turn 4×4 coupled transfer matrix T0.

**`MULTI_TURN_TRACKING_TO_MAT (TRACK, I_DIM, MAT1, TRACK0, CHI)`**
Subroutine to analyze 1-turn tracking data to find the 1-turn transfer matrix and the closed orbit offset.

**`ONE_TURN_MATRIX (RING, MAT6)`**
Subroutine to calculate the full 6×6 1-turn matrix.

**`ONE_TURN_MAT_AT_ELE (ELE, PHI_A, PHI_B, MAT4)`**
Subroutine to form the 4×4 1-turn coupled matrix with the reference point at the end of an element.

**`RING_MAKE_MAT6 (RING, IX_ELE, COORD_)`**
Subroutine to make the 6×6 linear transfer matrix for an element.

**`RING_MAKE_MAT627 (RING, IX_ELE, DIRECTION, MATS627)`**
Subroutine to make the 6×27 2nd order matrices for long term tracking. Used by, for example, TRACK_LONG.

**`TAYLOR_TO_MAT6 (A_TAYLOR, C0, MAT6, C1)`**
Subroutine to calculate the linear (Jacobian) matrix about some trajectory from a Taylor map.

**`TRANSFER_MAT_FROM_TRACKING (ELE, PARAM, ORB0, D_ORB, ERROR)`**
Subroutine to compute the transfer map for an element from tracking.

**`TRANSFER_MAT_FROM_TWISS (TWISS1, TWISS2, MAT)`**
Subroutine to make a 2×2 transfer matrix from the Twiss parameters at the end points.

**`TWISS_FROM_MAT2 (MAT, DET, TWISS, STAT, TOL, TYPE_OUT)`**
Subroutine to extract the Twiss parameters from one-turn 2×2 matrix.

**`TWISS_FROM_MAT6 (MAT6, ELE, STABLE, GROWTH_RATE)`**
Subroutine to extract the Twiss parameters from one-turn 6×6 matrix.

**`TWISS_TO_1_TURN_MAT (TWISS, PHI, MAT2)`**
Subroutine to form the 2×2 1-turn transfer matrix from the Twiss parameters.

---

## Routines called by MAKE_MAT6

**`MAKE_MAT6_BMAD (ELE, PARAM, C0, C1)`**
Subroutine to make the 6×6 transfer matrix for an element.

**`MAKE_MAT6_CUSTOM (ELE, PARAM, C0, C1)`**
Default routine for making the 6×6 transfer matrices for: 1) `ele%mat6_calc_method = custom$`  2) `ele%key = custom$`

**`MAKE_MAT6_RUNGE_KUTTA (ELE, PARAM, C0, C1)`**
Subroutine to make the 6×6 transfer matrix for an element using Runge-Kutta tracking.

**`MAKE_MAT6_SYMP_LIE_PTC (ELE, PARAM, C0, C1)`**
Subroutine to make the 6×6 transfer matrix for an element.

**`MAKE_MAT6_TAYLOR (ELE, PARAM, C0, C1)`**
Subroutine to make the 6×6 transfer matrix for an element.

**`MAKE_MAT6_TRACKING (ELE, PARAM, C0, C1)`**
Subroutine to make the 6×6 transfer matrix for an element using the Present tracking method.

---

## Low Level Matrix Routines

**`DRIFT_MAT6_CALC (MAT6, LENGTH, START, END)`**
Subroutine to calculate a drift transfer matrix with a possible kick.

**`MAT6_ADD_TILT_AT_END (MAT6, TILT)`**
Subroutine to add a tilt matrix to a transfer matrix: `mat6 -> tilt_mat * mat6`

**`MAT6_DISPERSION (MAT6, E_VEC)`**
Subroutine to put the dispersion into `ELE.MAT6` given the eta vector `E_VEC`.

**`QUAD_MAT_CALC (K1, LENGTH, MAT)`**
Subroutine to initialize the transfer matrix for a quad.

**`SOL_QUAD_MAT6_CALC (KS, K1, LENGTH, MAT6, ORB)`**
Subroutine to calculate the transfer matrix for a combination solenoid/quadrupole element.

**`TILT_MAT6 (MAT6, TILT)`**
Subroutine to transform a 6×6 transfer matrix to a new reference frame that is tilted in (x, Px, y, Py) with respect to the old reference frame.

---

## Tracking, Closed Orbit

**`CHECK_APERTURE_LIMIT (ORB, ELE, PARAM)`**
Subroutine to check if an orbit is outside the aperture.

**`CLOSED_ORBIT_AT_START (RING, CO, I_DIM, ITERATE)`**
Subroutine to calculate the closed orbit at the beginning of the ring.

**`CLOSED_ORBIT_FROM_TRACKING (RING, CLOSED_ORB_, I_DIM, EPS_REL, EPS_ABS, INIT_GUESS)`**
Subroutine to find the closed orbit via tracking.

**`DYNAMIC_APERTURE (RING, TRACK_INPUT, APERTURE)`**
Subroutine to determine the dynamic aperture of a lattice by tracking.

**`MULTI_TURN_TRACKING_ANALYSIS (TRACK, I_DIM, TRACK0, ELE, STABLE, GROWTH_RATE, CHI)`**
Subroutine to analyze multi-turn tracking data to get the Twiss parameters etc.

**`MULTI_TURN_TRACKING_TO_MAT (TRACK, I_DIM, MAT1, TRACK0, CHI)`**
Subroutine to analyze 1-turn tracking data to find the 1-turn transfer matrix and the closed orbit offset.

**`OFFSET_PARTICLE (ELE, PARAM, COORD, SET, SET_CANONICAL, SET_TILT, SET_MULTIPOLES, SET_HVKICKS, S_POS)`**
Subroutine to effectively offset an element by instead offsetting the particle position to correspond to the local element coordinates.

**`SETUP_RADIATION_TRACKING (RING, CLOSED_ORB, FLUCTUATIONS_ON, DAMPING_ON)`**
Subroutine to compute synchrotron radiation parameters prior to tracking.

**`TILT_COORDS (TILT_VAL, COORD, SET)`**
Subroutine to effectively tilt (rotate in the x-y plane) an element by instead rotating the particle position with negative the angle.

**`TRACK1 (START, ELE, PARAM, END)`**
Particle tracking through a single element.

**`TRACK_ALL (RING, ORBIT_)`**
Subroutine to track through the ring.

**`TRACK_BEND (START, ELE, END, IS_LOST)`**
Particle tracking through a bend element. Assumes no k1 quadrupole component. Used by track1.

**`TRACK_LONG (RING, ORBIT_, IX_START, DIRECTION, MATS627)`**
Subroutine to track for 1-turn using 2nd order transport matrices. Meant for long term tracking.

**`TRACK_MANY (RING, ORBIT_, IX_START, IX_END, DIRECTION)`**
Subroutine to track from one point in the ring to another.

**`TRACK_RUNGE_KUTTA (START, END, S_START, S_END, EPS, DEL_S_STEP, DEL_S_MIN, FUNC_TYPE, PARAM)`**
Subroutine to do tracking using Runge-Kutta integration.

**`TRACK_WIGGLER (START, ELE, PARAM, END)`**
Particle tracking through a wiggler. Used by track1.

**`TRANSFER_MAT_FROM_TRACKING (ELE, PARAM, ORB0, D_ORB, ERROR)`**
Subroutine to compute the transfer map for an element from tracking.

**`TWISS_AND_TRACK_AT_S (RING, S, ELE, ORB_, HERE)`**
Subroutine to calculate the Twiss parameters and orbit at a particular longitudinal position.

**`TWISS_FROM_TRACKING (RING, CLOSED_ORB_, D_ORB, ERROR)`**
Subroutine to compute from tracking, for every element in the ring, the Twiss parameters and the transfer matrices.

---

## Tracking Routines called by TRACK1

(Generally you don't call these routines directly.)

**`SYMP_LIE_BMAD (ELE, PARAM, START, END, CALC_MAT6)`**
Subroutine to track through an element (which gives the 0th order Taylor series) and optionally make the 6×6 transfer matrix.

**`TRACK1_ADAPTIVE_BORIS (START, ELE, PARAM, END, S_START, S_END)`**
Subroutine to do Boris tracking with adaptive step size control. Adapted from odeint in Numerical Recipes.

**`TRACK1_BORIS (START, ELE, PARAM, END, S_START, S_END)`**
Subroutine to do Boris tracking. See `boris_mod` documentation for more information.

**`TRACK1_BMAD (START, ELE, PARAM, END)`**
Particle tracking through a single element BMAD_standard style. NOT meant for long term tracking.

**`TRACK1_CUSTOM (START, ELE, PARAM, END)`**
Default routine for custom_tracking. This routine will do Runge-Kutta tracking.

**`TRACK1_LINEAR (START, ELE, PARAM, END)`**
Particle tracking through a single element assuming linearity (using `ele%mat6`).

**`TRACK1_RADIATION (START, ELE, PARAM, END, EDGE)`**
Subroutine to put in radiation damping and/or fluctuations.

**`TRACK1_RUNGE_KUTTA (START, ELE, PARAM, END)`**
Subroutine to do tracking using Runge-Kutta integration. Core routine is `odeint_bmad`.

**`TRACK1_SYMP_LIE_PTC (START, ELE, PARAM, END)`**
Particle tracking through a single element using a Hamiltonian and a symplectic integrator (Etienne's PTC).

**`TRACK1_SYMP_MAP (START, ELE, PARAM, END)`**
Particle tracking through a single element using a partially inverted Taylor map (genfield in PTC/FPP).

**`TRACK1_TAYLOR (START, ELE, PARAM, END)`**
Subroutine to track through an element using the element's Taylor series.

**`TRACK1_WIEDEMANN_WIGGLER (START, ELE, PARAM, END)`**
Subroutine to track through the body of a wiggler. Used by track1.

---

## Low Level Tracking Routines

**`ODEINT_BMAD (START, ELE, PARAM, END, S1, S2, REL_TOL, ABS_TOL, H1, HMIN)`**
Subroutine to do Runge-Kutta tracking.

**`TRACK_A_ACCEL_SOL (START, ELE, PARAM, END)`**
Subroutine to track through an accel_sol element.

**`TRACK1_BORIS_PARTIAL (START, ELE, PARAM, S, DS, END)`**
Subroutine to track 1 step using Boris tracking. Used by `track1_boris` and `track1_adaptive_boris`.

**`TRACK_A_DRIFT (ORB, LENGTH)`**
Subroutine to track through a drift.

**`TRACK_A_BEND (START, ELE, PARAM, END)`**
Particle tracking through a bend element.

---

## Particle Coordinate Stuff

**`CONVERT_COORDS (IN_TYPE_STR, COORD_IN, ELE, OUT_TYPE_STR, COORD_OUT)`**
Subroutine to convert between lab frame, normal mode, normalized normal mode, and action-angle coordinates.

**`TYPE_COORD (COORD)`**
Subroutine to type out a coordinate.

---

## Ring Geometry

**`RING_GEOMETRY (RING)`**
Subroutine to calculate the physical placement of all the elements in a ring (layout on the floor).

**`S_CALC (RING)`**
Subroutine to calculate the longitudinal distance S for the elements in a ring.

---

## Interface to PTC

**`CONCAT_REAL_8 (Y1, Y2, Y3)`**
Subroutine to concatenate two real_8 Taylor series.

**`ELE_TO_FIBRE (ELE, FIBER, PARAM, INTEG_ORDER, STEPS)`**
Subroutine to convert a BMAD element to a PTC fibre element. Allocates fresh storage for the fibre.

**`MAP_COEF (Y, I, J, K, L, STYLE)`**
Function to return the coefficient of the map `y(:)` up to 3rd order.

**`KILL_GEN_FIELD (GEN_FIELD)`**
Subroutine to kill a gen_field.

**`KIND_NAME (THIS_KIND)`**
Function to return the name of a PTC kind.

**`REAL_8_EQUAL_TAYLOR (Y8, BMAD_TAYLOR)`**
Subroutine to overload "=" in expressions `y8 = bmad_taylor`.

**`REAL_8_TO_TAYLOR (Y8, BMAD_TAYLOR, SWITCH_Z)`**
Subroutine to convert from a real_8 Taylor map in Etienne's PTC to a Taylor map in BMAD.

**`REAL_8_INIT (Y, SET_TAYLOR)`**
Subroutine to allocate a PTC real_8 variable.

**`REMOVE_CONSTANT_TAYLOR (TAYLOR_IN, TAYLOR_OUT, C0, remove_higher_order_terms)`**
Subroutine to remove the constant part of a Taylor series. Optionally removes terms higher order than `bmad_com%taylor_order` can handle.

**`RING_TO_LAYOUT (RING, PTC_LAYOUT)`**
Subroutine to create a PTC layout from a BMAD ring.

**`SET_PTC (PARAM, TAYLOR_ORDER, INTEG_ORDER, num_steps, no_cavity, exact_calc)`**
Subroutine to initialize PTC.

**`SET_TAYLOR_ORDER (ORDER, OVERRIDE_FLAG)`**
Subroutine to set the Taylor order.

**`SORT_UNIVERSAL_TERMS (UT_IN, UT_SORTED)`**
Subroutine to sort the Taylor terms from "lowest" to "highest". Needed since PTC output is not sorted.

**`TAYLOR_EQUAL_REAL_8 (BMAD_TAYLOR, Y8)`**
Subroutine to overload "=" in expressions `bmad_taylor = y8`.

**`TAYLOR_TO_REAL_8 (BMAD_TAYLOR, Y8, SWITCH_Z)`**
Subroutine to convert from a Taylor map in BMAD to a real_8 Taylor map in Etienne's PTC.

**`TYPE_LAYOUT (LAY)`**
Subroutine to print the global information in a PTC layout.

**`TYPE_MAP1 (Y, TYPE0, N_DIM, STYLE)`**
Subroutine to type the transfer map up to first order.

**`TYPE_FIBRE (FIB)`**
Subroutine to print the global information in a fibre.

**`TYPE_MAP (Y)`**
Subroutine to type the transfer maps of a real_8 array.

**`TYPE_REAL_8_TAYLORS (Y, SWITCH_Z)`**
Subroutine to type out the Taylor series from a real_8 array.

**`TAYLOR_TO_GENFIELD (BMAD_TAYLOR, GEN_FIELD, C0)`**
Subroutine to construct a genfield (partially inverted map) from a Taylor map.

**`UNIVERSAL_TO_BMAD_TAYLOR (U_TAYLOR, BMAD_TAYLOR, SWITCH_Z)`**
Subroutine to convert from a universal_taylor map in Etienne's PTC to a Taylor map in BMAD.

**`VEC_BMAD_TO_PTC (VEC_BMAD, VEC_PTC)`**
Subroutine to convert between BMAD and PTC coordinates.

**`VEC_PTC_TO_BMAD (VEC_PTC, VEC_BMAD)`**
Subroutine to convert between BMAD and PTC coordinates.

---

## Taylor Map Routines

**`CONCAT_TAYLOR (TAYLOR1, TAYLOR2, TAYLOR3)`**
Subroutine to concatenate two Taylor series: `taylor3(x) = taylor1(taylor2(x))`

**`ELE_TO_TAYLOR (ELE, ORB0, PARAM)`**
Subroutine to make a Taylor map for an element. Order set by `set_ptc`.

**`INIT_TAYLOR (BMAD_TAYLOR)`**
Subroutine to initialize (nullify) a BMAD Taylor map.

**`KILL_TAYLOR (BMAD_TAYLOR)`**
Subroutine to deallocate a BMAD Taylor map.

**`MAT6_TO_TAYLOR (MAT6, VEC0, BMAD_TAYLOR)`**
Subroutine to form a first order Taylor map from the 6×6 transfer matrix and the 0th order transfer vector.

**`SET_TAYLOR_ORDER (ORDER, OVERRIDE_FLAG)`**
Subroutine to set the Taylor order.

**`SORT_TAYLOR_TERMS (TAYLOR_IN, TAYLOR_SORTED)`**
Subroutine to sort the Taylor terms from "lowest" to "highest" of a Taylor series.

**`TAYLOR_EQUAL_TAYLOR (TAYLOR1, TAYLOR2)`**
Subroutine to transfer the values from one Taylor map to another: `Taylor1 <= Taylor2`

**`TAYLOR_TO_MAT6 (A_TAYLOR, C0, MAT6, C1)`**
Subroutine to calculate the linear (Jacobian) matrix about some trajectory from a Taylor map.

**`TAYLOR_INVERSE (TAYLOR_IN, TAYLOR_INV)`**
Subroutine to invert a Taylor map.

**`TAYLOR_PROPAGATE1 (TLR, ELE, PARAM)`**
Subroutine to track a real_8 Taylor map through an element. Alternative: `concat_taylor` if ele has a Taylor series.

**`TRACK_TAYLOR (START, BMAD_TAYLOR, END)`**
Subroutine to track using a Taylor map.

**`TYPE_TAYLORS (BMAD_TAYLOR)`**
Subroutine to print a BMAD Taylor map at the terminal in a nice format.

**`TYPE2_TAYLORS (BMAD_TAYLOR, LINES, N_LINES)`**
Subroutine to write a BMAD Taylor map in a nice format to a character array.

---

## Long Range Beam-Beam Interaction

**`INSERT_LRBBI (RING, RING_OPPOS, CROSS_POSITIONS, IX_LRBBI)`**
Uses a ring and a list of parasitic crossing points to create and insert beambeam elements at each crossing point.

**`LRBBI_CROSSINGS (N_BUCKET, OPPOS_BUCKETS, CROSS_POSITIONS)`**
Subroutine to calculate the location of the parasitic crossing points given a bunch and an array of positions of the bunches it will cross.

**`MAKE_LRBBI (MASTER_RING, MASTER_RING_OPPOS, RING, IX_LRBBI, MASTER_IX_LRBBI)`**
Subroutine to turn elements marking parasitic crossings into beam-beam elements.

**`MARK_LRBBI (MASTER_RING, MASTER_RING_OPPOS, RING, CROSSINGS)`**
Subroutine to insert named markers into the ring structure at the positions of parasitic crossings.

---

## Helper Subroutines: Informational

**`ATTRIBUTE_INDEX (KEY, NAME)`**
Function to return the index of an attribute for a given element type and the name of the attribute.

**`ATTRIBUTE_NAME (KEY, INDEX)`**
Function to return the name of an attribute for a particular type of element.

**`CHECK_RING_CONTROLS (RING, EXIT_ON_ERROR)`**
Subroutine to check if the control links in a ring structure are valid.

**`CHECK_ATTRIB_FREE (ELE, IX_ATTRIB, RING, ERR_FLAG, ERR_PRINT_FLAG)`**
Subroutine to check if an attribute is free to vary. Attributes that cannot be changed directly are super_slave attributes.

**`ELEMENTS_LOCATOR (KEY, RING, INDX)`**
Subroutine to locate all the elements of a certain kind in a ring. Note: super_slave elements are not included.

**`ELEMENT_LOCATOR (ELE_NAME, RING, IX_ELE)`**
Subroutine to locate an element in a ring.

**`FIND_ELEMENT_ENDS (RING, IX_ELE, IX_START, IX_END)`**
Subroutine to find the end points of an element.

**`TYPE_ELE (ELE, TYPE_ZERO_ATTRIB, TYPE_MAT6, TYPE_TWISS, TYPE_CONTROL)`**
Subroutine to type out the contents of an element.

**`TYPE2_ELE (ELE, TYPE_ZERO_ATTRIB, TYPE_MAT6, TYPE_TWISS, TYPE_CONTROL, LINES, N_LINES)`**
Like `TYPE_ELE` but the output is stored in the LINES array.

**`TYPE_TWISS (ELE, FREQUENCY_UNITS)`**
Subroutine to type out the Twiss parameters from an element.

**`TYPE2_TWISS (ELE, TYPE_ZERO_ATTRIB, TYPE_MAT6, TYPE_TWISS, TYPE_CONT, LINES, N_LINES)`**
Like `TYPE_TWISS` but the output is stored in the LINES array.

---

## Helper Subroutines: Elemental

(Adding elements, moving elements, etc.)

**`ADD_SUPERIMPOSE (RING, SUPER_ELE, IX_SUPER)`**
Subroutine to make a superimposed element.

**`ATTRIBUTE_BOOKKEEPER (ELE, PARAM)`**
Subroutine to make sure the attributes of an element are self-consistent.

**`CESR_CROSSINGS (I_TRAIN, J_CAR, SPECIES, N_TRAINS_TOT, N_CARS, CROSS_POSITIONS, N_CAR_SPACING, TRAIN_SPACING)`**
Subroutine to calculate all parasitic crossing points for a bunch at a given location.

**`CHANGE_BASIS (COORD, REF_ENERGY, REF_Z, TO_CART, TIME_DISP)`**
Subroutine to convert accelerator coordinates (x, x', y, y', z, z') to Cartesian coordinates and time derivatives, or vice versa.

**`CHECK_ELE_ATTRIBUTE_SET (RING, I_ELE, ATTRIB_NAME, IX_ATTRIB, ERR_FLAG, ERR_PRINT_FLAG)`**
Subroutine to check whether a particular attribute of an element can be changed directly.

**`CREATE_GROUP (RING, IX_ELE, N_CONTROL, CONTROL_)`**
Subroutine to create a group control element.

**`CREATE_OVERLAY (RING, IX_OVERLAY, IX_VALUE, N_SLAVE, CON_)`**
Subroutine to add the controller information to slave elements of an overlay_lord.

**`COMPRESS_RING (RING, OK)`**
Subroutine to compress the `ele_()` and `control_()` arrays to remove elements no longer used.

**`INSERT_ELEMENT (RING, INSERT_ELE, INSERT_INDEX)`**
Insert a new element into the regular part of the ring structure.

**`MAKE_HYBRID_RING (RING_IN, USE_ELE, REMOVE_MARKERS, RING_OUT, IX_OUT)`**
Subroutine to concatenate together elements to make a hybrid ring.

**`NEW_CONTROL (RING, IX_ELE)`**
Subroutine to create a new control element.

**`POINTER_TO_ATTRIBUTE (ELE, ATTRIB_NAME, DO_ALLOCATION, PTR_ATTRIB, IX_ATTRIB, ERR_FLAG, ERR_PRINT_FLAG)`**
Returns a pointer to an attribute of an element with name `attrib_name`.

**`RING_SET_ELE_ATTRIBUTE (RING, I_ELE, ATTRIB_NAME, ATTRIB_VALUE, ERR_FLAG, MAKE_MAT6_FLAG, ORBIT_)`**
Subroutine to set the attribute of an element, propagate the change to any slave elements, and optionally remake the 6×6 transfer matrix.

**`SET_ELE_ATTRIBUTE (RING, I_ELE, ATTRIB_NAME, ATTRIB_VALUE, ERR_FLAG, MAKE_MAT6_FLAG, ORBIT_)`**
Subroutine to set the attribute of an element and propagate the change to any slave elements.

**`SPLIT_RING (RING, S_SPLIT, IX_SPLIT, SPLIT_DONE)`**
Subroutine to split a ring at a point. Will not split if it would create a "runt" element shorter than 1 μm.

**`UPDATE_HYBRID_LIST (RING, N_IN, USE_ELE)`**
`USE_ELE` is a list of elements that should not be hybridized.

---

## Helper Subroutines: Transformational

**`ADJUST_CONTROL_STRUCT (RING, IX_ELE)`**
Subroutine to adjust the control structure of a ring so that extra control elements can be added.

**`CONTROL_BOOKKEEPER (RING, IX_ELE)`**
Subroutine to calculate the combined strength of the attributes for controlled elements.

**`INIT_ELE (ELE)`**
Subroutine to initialize an element. Element is initialized to be free.

**`REVERSE_ELE (ELE)`**
Subroutine to "reverse" an element for backward tracking.

**`RING_REVERSE (RING_IN, RING_REV)`**
Subroutine to construct a ring structure with the elements in reversed order (for backward tracking).

**`SET_DESIGN_LINEAR (RING)`**
Subroutine to set only those elements on that constitute the "design" lattice (quadrupoles, bends, and wigglers).

**`SET_ON (KEY, RING, ON_SWITCH, ORB_)`**
Subroutine to turn on or off a set of elements (quadrupoles, rfcavities, etc.) in a ring.

**`SET_SYMMETRY (SYMMETRY, RING)`**
Subroutine to set the symmetry of a ring.

**`TRANSFER_RING_PARAMETERS (RING_IN, RING_OUT)`**
Subroutine to transfer the ring parameters (such as `ring%name`, `ring%param`, etc.) from one ring to another.

---

## Helper Subroutines: Multipolar

**`MULTIPOLE_AB_TO_KT (AN, BN, KNL, TN)`**
Subroutine to convert ab type multipoles to kt (MAD standard) multipoles.

**`MULTIPOLE_ELE_TO_AB (ELE, PARTICLE, A, B, USE_ELE_TILT)`**
Subroutine to put the scaled element multipole components (normal and skew) into 2 vectors.

**`MULTIPOLE_ELE_TO_KT (ELE, PARTICLE, KNL, TILT, USE_ELE_TILT)`**
Subroutine to put the scaled element multipole components (strength and tilt) into 2 vectors.

**`MULTIPOLE_KICK (KNL, TILT, N, COORD)`**
Subroutine to put in the kick due to a multipole.

**`MULTIPOLE_KT_TO_AB (KNL, TN, AN, BN)`**
Subroutine to convert kt (MAD standard) multipoles to ab type multipoles.

---

## Helper Subroutines: Miscellaneous

**`C_MULTI (N, M)`**
Subroutine to compute multipole factors: `c_multi(n, m) = +/- ("n choose m") / n!`

**`COMPUTE_ELEMENT_ENERGY (RING)`**
Subroutine to compute the energy of the reference particle for each element in a ring structure.

---

*Created by: David Sagan | Last updated: January 2003 | Topic revision: r3 - 23 Oct 2013*
