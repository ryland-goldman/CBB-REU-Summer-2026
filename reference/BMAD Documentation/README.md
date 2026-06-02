# BMAD Documentation

BMAD (aka "Baby MAD" or "Better MAD" or "BE MAD!") is a set of Fortran90 subroutines to read in lattice specification files, compute Twiss parameters, track particles, etc. It conforms to the MAD input standard and was developed at Cornell (CESR/CLASSE) by David Sagan.

## Documents

| File | Description |
|------|-------------|
| [lattice_language_guide.md](lattice_language_guide.md) | BMAD Lattice Language Guide — elements, attributes, tracking switches, superimpose/overlay/group, input file syntax |
| [programming_guide.md](programming_guide.md) | BMAD Programming Guide — ring structure, transfer matrices, tracking subroutine chains, Taylor maps, custom elements, Runge-Kutta tracking |
| [subroutine_list.md](subroutine_list.md) | Library Subroutine List — all BMAD subroutines organized by category with signatures and descriptions |

## Quick Reference

- **Input format**: MAD-compatible lattice files; parsed by `BMAD_PARSER`
- **Core data structure**: `ring_struct` — array of `ele_(i)` from 0 to `ring%n_ele_max`
- **Particle coordinates**: (x, Px, y, Py, z, Pz)
- **Units**: MKS — meters, radians, Tesla, eV; phase angles in radians/2π; voltage in Volts
- **Tracking methods**: BMAD_Standard, Runge_Kutta, Symp_Lie_PTC, Taylor, Linear, Adaptive_Boris
- **Transfer matrix**: 6×6 Jacobian computed by `MAKE_MAT6`; Twiss by `TWISS_AND_TRACK`

## Key Concepts

- **LORD/SLAVE elements**: Controlling ("lord") elements in indices `n_ele_ring+1` to `n_ele_max`; regular ("slave/free") elements in indices 1 to `n_ele_ring`; element 0 holds ring-start Twiss
- **Control types**: `FREE`, `OVERLAY_SLAVE`, `SUPER_SLAVE` (regular); `GROUP_LORD`, `OVERLAY_LORD`, `SUPER_LORD` (lords)
- **SUPERIMPOSE**: Place an element on top of existing elements without explicit list position
- **OVERLAY**: Controls the absolute value of a single attribute across slave elements
- **GROUP**: Controls attribute changes (deltas) and can control position/length via ACCORDION_EDGE etc.
- **Taylor maps**: Computed via Etienne Forest's FPP/PTC; stored as 6 `taylor_struct`s (one per coordinate); order set globally by `SET_TAYLOR_ORDER`
