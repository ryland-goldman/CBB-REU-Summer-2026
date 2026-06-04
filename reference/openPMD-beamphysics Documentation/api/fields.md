# `FieldMesh`

Container for external field-mesh data in the openPMD-beamphysics standard. This is the
class used in `impactt_cathode/` to build the DC gun field via `from_onaxis` and hand it
to IMPACT-T.

```python
FieldMesh(h5=None, data=None)
```

Initialize from an openPMD HDF5 file or a data dictionary.

**Parameters**
- `h5` (str or file handle): open HDF5 handle or file path.
- `data` (dict): raw data dictionary with `'attrs'` and `'components'` keys.

---

## Class-method constructors

### `from_onaxis()`
```python
@classmethod
from_onaxis(*, z=None, Bz=None, Ez=None, frequency=0, harmonic=None, eleAnchorPt='beginning')
```
Create a `FieldMesh` from on-axis field data.
- `z`: array of z-coordinates (**must be regularly spaced**).
- `Bz`: magnetic field at r=0 in Tesla (optional).
- `Ez`: electric field at r=0 in V/m (optional).
- `frequency`: fundamental frequency in Hz (default 0 → static).
- `harmonic`: harmonic number; defaults to 1 if `frequency ≠ 0`, else 0.
- `eleAnchorPt`: element anchor point — `'beginning'`, `'center'`, or `'end'`.
- **Returns:** `FieldMesh` instance.

### Other constructors
```python
from_superfish(filename, type=None, geometry='cylindrical')         # Superfish T7
from_ansys_ascii_3d(*, efile=None, hfile=None, frequency=None)      # ANSYS ASCII (E in V/m, H in A/m)
from_astra_3d(common_filename, frequency=0)                          # ASTRA 3D fieldmaps
from_impact_emfield_cartesian(filename, frequency=0, eleAnchorPt='beginning')  # Impact-T EMfldCart
```

---

## Properties

**Coordinates** — `coord_vecs`; `coord_vec(key)`; `axis_labels`; `axis_index(key)`.
**Field** — `factor` (complex scaling `scale·e^{iφ}`); `phase` (writable, `φ = -2π·RFphase`);
`meshgrid`.
**Classification** — `is_pure_electric`; `is_pure_magnetic`.

---

## Methods

### Interpolation
```python
interpolate(key, points)        # key e.g. 'Ez', 'magneticField/y'; points shape (3,) or (n,3)
interpolator(key)               # -> scipy RegularGridInterpolator
axis_values(axis_label, field_key, **kwargs)
```

### Components
```python
component_is_zero(key)          # True if the component is all zeros
scaled_component(key)           # component scaled by the complex factor
axis_points(axis_label)         # 3D points along an axis, for interpolation
```

### Visualization
```python
plot(component=None, *, cmap=None, nice=True, stream=False, mirror=None,
     density=2, linewidth=1, arrowsize=1, axes=None, return_figure=False, **kwargs)
```
- `component`: field to plot (defaults to `'B'` or `'E'`).
- `stream`: add streamlines; `mirror='r'` symmetrizes cylindrical data.
- `density, linewidth, arrowsize`: streamline styling. `axes`: matplotlib Axes.

### File I/O
```python
write(h5, name=None)                              # openPMD-beamphysics HDF5
write_gpt(filePath, asci2gdf_bin=None, verbose=True)
write_superfish(filePath, verbose=False)          # Poisson (static) / Fish (dynamic)
write_impact_emfield_cartesian(filename)          # Impact-T EMfldCart element
```

### Geometry / utility
```python
to_cylindrical()                # convert to cylindrical (uses y=0 slice if rectangular)
copy()                          # deep copy
units(key)                      # units for any field key
__eq__(other)                   # compare attrs + component data
__getitem__(key)                # component access with operators (re_, im_, abs_)
```
