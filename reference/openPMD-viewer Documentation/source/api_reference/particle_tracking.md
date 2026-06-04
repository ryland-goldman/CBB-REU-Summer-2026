# Particle tracking: `ParticleTracker`

```python
class openpmd_viewer.ParticleTracker(ts, species=None, t=None, iteration=None,
                                     select=None, preserve_particle_index=False)
```

Select particles at a given iteration so they can be retrieved (tracked by unique ID) at
a later iteration. A `ParticleTracker` instance can be passed as the `select` argument of
`OpenPMDTimeSeries.get_particle`.

**Parameters**

- `ts` (`OpenPMDTimeSeries` object): contains the particle data.
- `species` (string, optional): name of the particle species. Required only if multiple
  species exist.
- `t` (float, seconds, optional): time point for data retrieval. Uses closest iteration
  if no exact match. Either `t` or `iteration` is required.
- `iteration` (int, optional): iteration number for data retrieval. Either `t` or
  `iteration` is required.
- `select` (dict or 1D array of int, optional): particle selection. Dict form specifies
  ranges: `'x': [-4., 10.]` (position bounds), `'ux': [-0.1, 0.1]` (velocity),
  `'uz': [5., None]` (open range). Alternatively, pass a 1D integer array of particle IDs.
- `preserve_particle_index` (bool, optional): when `True`, particles keep consistent
  indices across iterations, with `NaN` for absent particles. When `False`, arrays shrink
  as particles vanish and indices are not preserved.
