# Generic interface: `OpenPMDTimeSeries`

The main entry point for loading an openPMD file series.

```python
class openpmd_viewer.OpenPMDTimeSeries(path_to_dir, check_all_files=True, backend=None)
```

Initialize an openPMD time series by scanning a directory to extract openPMD files and
relevant parameters.

**Parameters**

- `path_to_dir` (string): directory path containing openPMD files.
- `check_all_files` (bool, optional): validate consistency across all files (same
  fields, particles, metadata). Set `False` for faster access. Default `True`.
- `backend` (string): data-reading backend — `"openpmd-api"` or `"h5py"`. Defaults to
  openpmd-api if available, otherwise h5py.

---

## `get_field()`

```python
get_field(field=None, coord=None, t=None, iteration=None, m='all', theta=0.0,
          slice_across=None, slice_relative_position=None, plot=False,
          plot_range=[[None, None], [None, None]], **kw)
```

Extract a specified field from an openPMD file.

**Parameters**

- `field` (string, optional): field to extract.
- `coord` (string, optional): field component to extract.
- `m` (int or str, optional): for thetaMode geometry — `'all'` (sum all modes) or an
  integer (select a specific mode).
- `t` (float, seconds, optional): time for data retrieval; uses closest iteration if no
  exact match.
- `iteration` (int): iteration for data retrieval. Specify either `t` or `iteration`.
- `theta` (float or None, optional): for thetaMode, observation-plane angle relative to
  the x-axis. Returns a 2D array if specified; otherwise a 3D Cartesian array.
- `slice_across` (str or list of str, optional): slicing directions. Cartesian: `'z'`
  (1d); `'x'`/`'z'` (2d); `'x'`/`'y'`/`'z'` (3d). Cylindrical: `'r'`/`'z'`. Reduces
  dimension by one per slice.
- `slice_relative_position` (float or list of float, optional): values in [-1, 1] giving
  slice location (-1 lower edge, 0 middle, 1 upper edge). Default 0 for all directions.
- `plot` (bool, optional): display the requested quantity.
- `plot_range` (list of lists): clip plot values — two two-element lists for the 1st and
  2nd axes.
- `**kw` (dict, optional): additional matplotlib `imshow` options.

**Returns** — tuple `(F, info)`:
- `F` (2D array): field data.
- `info` (`FieldMetaInformation` object): associated metadata.

---

## `get_particle()`

```python
get_particle(var_list=None, species=None, t=None, iteration=None, select=None,
             plot=False, nbins=150, plot_range=[[None, None], [None, None]],
             use_field_mesh=True, histogram_deposition='cic', **kw)
```

Extract particle variables from an openPMD file, with optional histogram visualization.

**Parameters**

- `var_list` (list of string, optional): particle variables to extract. If omitted, the
  available quantities are printed.
- `species` (string): species name. Optional if only one species exists.
- `t` (float, seconds, optional): time for data retrieval; uses closest iteration if no
  exact match.
- `iteration` (int): iteration for data retrieval. Specify either `t` or `iteration`.
- `select` (dict or `ParticleTracker` object, optional): particle selection. Dict form
  specifies ranges, e.g. `'x': [-4., 10.]`, `'ux': [-0.1, 0.1]`, `'uz': [5., None]`.
- `plot` (bool, optional): display results. Available for one or two quantities only.
- `nbins` (int, optional): histogram bin count (when `plot=True`).
- `plot_range` (list of lists): histogram value bounds — two two-element lists.
- `use_field_mesh` (bool, optional): when `True`, auto-determine histogram extent and
  match bin spacing to the grid spacing, avoiding artifacts.
- `histogram_deposition` (string): `"ngp"` (Nearest Grid Point) or `"cic"` (Cloud-In-Cell)
  for particle deposition onto bins. Default `cic` (smoother).
- `**kw` (dict, optional): additional matplotlib `hist`/`hist2d` options.

**Returns** — list of 1D arrays corresponding to the requested variables, in `var_list`
order.

---

## `iterate()`

```python
iterate(called_method, *args, **kwargs)
```

Repeatedly invoke a method across every iteration in the time series.

**Parameters**

- `called_method`: method to call for each iteration.
- `*args`: positional arguments for the method (exclude `t`/`iteration`).
- `**kwargs`: keyword arguments for the method (exclude `t`/`iteration`).

**Returns** — results as a list or array (when possible), with iterations as the first
axis. If the method returns tuples/lists, `iterate` returns corresponding collections.

---

## `slider()`

```python
slider(figsize=(6, 5), fields_figure=0, particles_figure=1,
       exclude_particle_records=['charge', 'mass'], **kw)
```

Interactive navigation through simulation iterations via a slider control.

**Parameters**

- `figsize` (tuple): figure dimensions.
- `fields_figure` (int): matplotlib figure number for field display.
- `particles_figure` (int): matplotlib figure number for particle display.
- `exclude_particle_records` (list of strings): particle quantities to exclude from the
  slider display.
- `**kw` (dict): extra matplotlib `imshow` options (`cmap`, etc.) — sets initial plotting
  parameters, modifiable through the slider interface.
