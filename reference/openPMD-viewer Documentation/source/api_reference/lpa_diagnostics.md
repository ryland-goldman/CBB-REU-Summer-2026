# Domain-specific interface: `LpaDiagnostics`

Laser-plasma-acceleration diagnostics built on top of `OpenPMDTimeSeries`. Adds beam and
laser analysis methods (emittance, energy spread, current, a0, waist, …).

```python
class openpmd_viewer.addons.LpaDiagnostics(path_to_dir, check_all_files=True, backend=None)
```

**Parameters**

- `path_to_dir` (string): directory containing openPMD files. HDF5 only; files named with
  iteration numbers + `.h5` (e.g. `data0005000.h5`).
- `check_all_files` (bool, optional): verify consistency across all files. `False` for
  faster access. Default `True`.
- `backend` (string): `'openpmd-api'` or `'h5py'`. Defaults to openpmd-api if available,
  otherwise h5py.

Every method below accepts `t` (float, seconds) **or** `iteration` (int) — exactly one is
required — and most accept `species` (str) and `select` (dict of selection ranges, e.g.
`'x': [-4., 10.]`).

---

## Beam diagnostics

### `get_energy_spread()`
```python
get_energy_spread(t=None, iteration=None, species=None, select=None,
                  center='mean', width='std', property='energy')
```
Central energy and energy spread (weighted by particle weights).
- `center` (str): `'mean'` or `'median'`. Default `'mean'`.
- `width` (str): `'std'` or `'mad'`. Default `'std'`.
- `property` (str): `'energy'` (MeV) or `'gamma'` (Lorentz factor). Default `'energy'`.
- **Returns** `(central_energy, energy_spread)`; `NaN` if the selection is empty.

### `get_mean_gamma()`
```python
get_mean_gamma(t=None, iteration=None, species=None, select=None)
```
Weighted mean gamma and standard deviation.
- **Returns** `(mean_gamma, std)`; `NaN` if the selection is empty.

### `get_sigma_gamma_slice()`
```python
get_sigma_gamma_slice(dz, t=None, iteration=None, species=None, select=None, plot=False, **kw)
```
Gamma standard deviation within z-axis slices.
- `dz` (float): slice width in micrometers.
- `plot` (bool): default `False`. `**kw`: matplotlib options.
- **Returns** `(sigma_gamma_per_slice, central_z_positions)` (arrays).

### `get_charge()`
```python
get_charge(t=None, iteration=None, species=None, select=None)
```
Total electric charge of selected particles.
- **Returns** float, charge in Coulombs.

### `get_divergence()`
```python
get_divergence(t=None, iteration=None, species=None, select=None)
```
Divergence in the x and y planes.
- **Returns** `(divergence_x_rad, divergence_y_rad)`; `NaN` if empty.

### `get_emittance()`
```python
get_emittance(t=None, iteration=None, species=None, select=None, kind='normalized',
              description='projected', nslices=0, beam_length=None)
```
RMS emittance (Floettmann et al. methodology).
- `kind` (str): `'normalized'` or `'trace'`. Default `'normalized'`.
- `description` (str): `'projected'`, `'all-slices'`, or `'slice-averaged'`. Default
  `'projected'`.
- `nslices` (int): number of slices for slice-emittance calculations.
- `beam_length` (float): beam length in meters. Default 4× the z standard deviation.
- **Returns** — projected / slice-averaged: `(emittance_x, emittance_y)` (π·m·rad);
  all-slices: `(emittance_x_array, emittance_y_array, electrons_per_slice, slice_centers)`.

### `get_current()`
```python
get_current(t=None, iteration=None, species=None, select=None, bins=100, plot=False, **kw)
```
Electric current along the z-axis.
- `bins` (int): number of z-axis bins. Default 100. `plot` (bool): default `False`.
- **Returns** `(current_per_bin_ampere, FieldMetaInformation)`.

---

## Laser diagnostics

### `get_laser_envelope()`
```python
get_laser_envelope(t=None, iteration=None, pol=None, laser_propagation='z', m='all',
                   theta=0, slice_across=None, slice_relative_position=None, plot=False,
                   plot_range=[[None, None], [None, None]], **kw)
```
Laser envelope via high-frequency filtering.
- `pol` (str): field polarization `'x'`, `'y'`, or `'z'`.
- `laser_propagation` (str): propagation coordinate. Default `'z'`.
- `m`, `theta`, `slice_across`, `slice_relative_position`: as in `get_field`.
- **Returns** `(envelope_1d_or_2d_array, FieldMetaInformation)`.

### `get_main_frequency()`
```python
get_main_frequency(t=None, iteration=None, pol=None, m='all', method='max')
```
Laser angular frequency.
- `pol` (str): `'x'` or `'y'`. `method` (str): `'fit'` (Gaussian) or `'max'` (peak
  intensity). Default `'max'`.
- **Returns** float, mean angular frequency.

### `get_spectrum()`
```python
get_spectrum(t=None, iteration=None, pol=None, m='all', plot=False, **kw)
```
Laser spectrum via Fourier-transform magnitude.
- **Returns** `(spectrum_1d_array, FieldMetaInformation)`.

### `get_a0()`
```python
get_a0(t=None, iteration=None, pol=None)
```
Normalized vector potential `a0 = Emax·e / (me·c·ω)`.
- **Returns** float, `a0`.

### `get_ctau()`
```python
get_ctau(t=None, iteration=None, pol=None, method='fit')
```
Laser pulse length (longitudinal waist = √2·σz).
- `method` (str): `'fit'` (Gaussian) or `'rms'`. Default `'fit'`.
- **Returns** float, pulse length in meters.

### `get_laser_waist()`
```python
get_laser_waist(t=None, iteration=None, pol=None, theta=0, laser_propagation='z',
                method='fit', profile_method='peak')
```
Laser waist (√2·σr). In 3D, evaluates the x–z plane via a y-slice.
- `method` (str): `'fit'` or `'rms'`. Default `'fit'`.
- `profile_method` (str): `'peak'` or `'projection'`. Default `'peak'`.
- **Returns** float, waist in meters.

### `get_spectrogram()`
```python
get_spectrogram(t=None, iteration=None, pol=None, plot=False, **kw)
```
Laser spectrogram via FROG (Frequency-Resolved Optical Gating).
- **Returns** `(spectrogram_2d_array, FieldMetaInformation)`.
