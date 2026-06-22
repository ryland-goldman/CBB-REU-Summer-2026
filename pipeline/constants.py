"""Physical constants (SI + accelerator-physics eV conventions), from scipy.

One definition shared by every stage so the rest-energy / charge / mass literals
cannot drift between files. See CLAUDE.md "Codebase Standards" for the unit
conventions (u = gamma*beta*c; momentum in eV/c). The four WarpX stages also have
picmi.constants available in-sim; these scipy values agree with them to machine
precision and are the single source for the pure-Python build/plot/IO code.
"""

import scipy.constants as _sc

C_LIGHT = _sc.c                          # m/s
E_CHARGE = _sc.e                         # C (elementary charge, positive)
M_E = _sc.m_e                            # kg (electron mass)
EPS0 = _sc.epsilon_0                     # F/m
K_B = _sc.k                              # J/K (Boltzmann)
K_B_EV = _sc.k / _sc.e                   # eV/K
MC2_EV = _sc.m_e * _sc.c**2 / _sc.e      # electron rest energy [eV] ~= 510998.95
