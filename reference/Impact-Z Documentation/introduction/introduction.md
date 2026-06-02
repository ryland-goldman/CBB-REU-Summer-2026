## 1. Introduction

IMPACT-Z is a 3D parallel/serial Particle-In-Cell (PIC) code based on multi-layer object-oriented design. The present version of IMPACT-Z can treat intense beams propagating through drifts, magnetic quadrupoles, magnetic solenoids, bending magnets, multipoles, and RF cavities, using map integrator or nonlinear Lorentz integrator.

> **Warning:** Some elements such as 3D EM field can be used **ONLY** for the Lorentz integrator.

It has a novel treatment of RF cavities, in which the gap transfer maps are computed during the simulations by reading in Superfish RF fields. The goal is to avoid time-consuming (and unnecessary) fine-scale integration of millions of particles through the highly z-dependent cavity fields. Instead, fine-scale integration is used to compute the maps (which involve a small number of terms), and the maps are applied to particles. If you are familiar with magnetic optics, then you will recognize that this is analogous to the technique used to simulate beam transport through magnets with fringe fields.

The version of IMPACT-Z (v2.2) currently has a 3D space-charge model that assumes 3D open boundary conditions. The other boundary conditions and more functions will be added later on in the new versions. The error studies can include the field errors, misalignment errors, and rotation errors.
