# 2 Recently Added Features

- Merge development version into master version. Each particle has 9 attributes now. In addition to the 6 phase space coordinates, charge over mass ratio, charge per macroparticle, and global id of each macroparticle are added.

- Nonlinear quadrupole element with input discrete gradient data.

- Element 105 to allow discrete data on axis when the input file name number > 1000.

- Beam loss check to enable element aperture radius: modify collimator element to allow for both round and rectangular aperture shape.

- Output 4x4 and x6 sigma matrices in fort.31 and fort.32.

- Three new distribution types.

- Global random seed number (Dim) for multiple simulations with different random numbers.

- A beam rotation with respect to longitudinal z axis (-17).

- A beam longitudinal heating (-16).

- A dielectric wakefield model from Dianel Mihalcea (-13).

- 3D field data on Cartesian grid (111 EMfldCart) has to be in the format of **complex number** to represent both traveling wave and standing wave fields.

- An alpha magnet field model, traveling wave in meander plates model, and DC surface roughness field model in (113) EMfldAna.

- New initial particle distribution `ijk`, which contains a total of 24 types of initial distribution behind the cathode.

- Python script code, `PhaseOpt.py`, automatically finds the RF initial driven phases of cavities used in the input file, `ImpactT.in`, based on the user specified design phases in the `ImpactT.in`.

- Switch (-12) for applying instant linear matrix kick to a beam at a given location.

- Modified switch (-1) for steering so that steering means an instant kick.

- Switch (-15) for mean-field grid based or point-to-point N-body space charge solver.

- Quadrupole element: If v9 (rotation angle with respect to z) is nonzero, this is a skew quadrupole. If v10 (rf frequency) and v11 (phase) are nonzero, this is an rf quadrupole.

- Collimation function at given z location by using "-11".

- Output slice-based information (current, uncorrelated energy spread, slice emittances, correlated energy spread) at given z location using "-9". The slice information for the initial distribution and the final output distribution is stored in file fort.60 and fort.70.

- Modified the structure wakefield calculation so that the code can use both analytical expressions and read-in transverse and longitudinal wake function from external files.

- Particle coordinates for initial read-in distribution are x(m), Px/mc, y(m), Py/mc, z(m), Pz/mc.

- Particle coordinates for output phase distribution at given location are x(m), Px/mc, y(m), Py/mc, z(m), Pz/mc.

- Switch flag for cathode (Nemission > 0, flagcathode = 1, cathode exists, otherwise, 0 no cathode model).

- 1D CSR wake module including the transient effects at the entrance and the exit of dipole bend magnet with integrated Green method.
