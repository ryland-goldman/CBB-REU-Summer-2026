# `ParticleGroup`

Container for particle data in the openPMD-beamphysics standard.

```python
ParticleGroup(h5=None, data=None)
```

Initialize from an HDF5 file/handle or a data dictionary.

**Parameters**
- `h5` (str, Path, h5py.File, or dict): source of particle data — filename, open HDF5
  handle, or dict-like openPMD particle data.
- `data` (dict): raw particle data with keys `x, px, y, py, z, pz, t, status, weight,
  species`.

**Required data keys**
- `x, y, z` — positions in meters
- `px, py, pz` — momenta in eV/c
- `t` — time in seconds
- `weight` — macro-charge weight in Coulombs
- `species` — species name (e.g. `'electron'`)

**Optional:** `id` — unique integer particle identifiers.

---

## Properties & derived attributes

**Position & momentum**
`x, y, z` (m); `px, py, pz` (eV/c); `t` (s); `xp, yp` (slopes px/pz, py/pz);
`r` (xy radius); `theta` (xy angle, rad); `pr, ptheta` (radial/angular momentum).

**Normalized coordinates**
`x_bar, px_bar, y_bar, py_bar` (√m); `Jx, Jy` (normalized amplitudes).

**Relativistic**
`gamma`; `beta, beta_x, beta_y, beta_z`; `p` (total momentum, eV/c); `energy` (total, eV);
`kinetic_energy` (eV); `mass` (rest mass, eV); `higher_order_energy` (energy with a
quadratic fit subtracted).

**Beam properties**
`charge` (C); `weight` (array, C); `species`; `species_charge` (C); `id`; `status`;
`n_particle`; `n_alive` (status == 1); `n_dead` (status ≠ 1); `Lz` (angular momentum
about z).

**Emittance & Twiss**
`norm_emit_x`, `norm_emit_y`, `norm_emit_4d`; `higher_order_energy_spread`.

**Coordinate-system flags**
`in_z_coordinates` (all particles at same z); `in_t_coordinates` (all at same t).

**Misc**
`average_current` (A); `data` (internal data dict).

---

## Methods

### Statistics
```python
min(key)   max(key)   ptp(key)   avg(key)   std(key)
cov(*keys) delta(key)
```
Weighted statistics on any particle property. `delta` returns the value minus its mean.

### Filtering & selection
```python
where(condition) -> ParticleGroup        # boolean-mask filter
__getitem__(key) -> ndarray | float | ParticleGroup   # pg['x'], pg['norm_emit_x'], ...
```

### Transformation
```python
copy() -> ParticleGroup
resample(n=0, equal_weights=False) -> ParticleGroup
split(n_chunks=100, key='z') -> list
fractional_split(fractions, key) -> list
assign_id() -> None
```

### Drift & rotation
```python
drift(delta_t)                  drift_to_z(z=None)        drift_to_t(t=None)
rotate(*, x_rot=0.0, y_rot=0.0, z_rot=0.0, order='zxy', xc=0.0, yc=0.0, zc=0.0)
rotate_x(theta, yc=0.0, zc=0.0) rotate_y(theta, xc=0.0, zc=0.0) rotate_z(theta, xc=0.0, yc=0.0)
linear_point_transform(mat3)
```

### Twiss & matching
```python
twiss(plane='x', fraction=1, p0c=None) -> dict
twiss_match(beta=None, alpha=None, plane='x', p0c=None, inplace=False) -> ParticleGroup
```

### Conversion (Bmad)
```python
to_bmad(p0c=None, tref=None) -> dict
from_bmad(bmad_dict) -> ParticleGroup
```

### I/O — export to accelerator-code formats
```python
write(h5)
write_bmad(filename)     write_elegant(filename)  write_astra(filename)
write_gpt(filename)      write_opal(filename)     write_impact(filename)
write_litrack(filename)  write_lucretia(filename) write_simion(filename)
write_genesis2_beam_file(filename)
write_genesis4_beam(filename)  write_genesis4_distribution(filename)
```

### Plotting
```python
plot(key1='x', key2=None, bins=None, xlim=None, ylim=None,
     return_figure=False, tex=True, nice=True, ellipse=False, **kwargs)
slice_plot(*keys, n_slice=100, slice_key=None, tex=True, nice=True,
           return_figure=False, xlim=None, ylim=None, **kwargs)
wakefield_plot(wake, key=None, nice=True, ax=None, xlim=None, ylim=None,
               tex=True, bins=None, **kwargs)
```

### Analysis
```python
slice_statistics(*keys, n_slice=100, slice_key=None) -> dict
bunching(wavelength) -> complex
higher_order_energy_calc(order=2) -> ndarray
histogramdd(*keys, bins=10, range=None) -> tuple
info(key) -> dict
units(key) -> pmd_unit
```

### Wakefield application
```python
apply_wakefield(wakefield, length, inplace=False, include_self_kick=True) -> ParticleGroup
```

### Magic methods
`__len__`, `__add__` (concatenate two groups), `__eq__`, `__contains__` (`'x' in pg`),
`__repr__`, `__str__`.
