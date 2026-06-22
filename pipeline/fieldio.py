"""Shared GDF→openPMD thetaMode (RZ, m=0) field-map IO for the build_*_field scripts.

One GDF flat-column parse (R-fastest reshape, optional r-pad, optional
descending-z row reversal) and one openPMD mesh writer, so the gun / injector /
linac_sec1 field builders share a single schema (axis order ["r","z"], m=0,
nodal position, V/m and T unit dimensions). Per-stage block names, field scales,
z-offsets, and physics asserts stay in each build_*_field.py.

See gun/README.md and CLAUDE.md for the axis-order / m-mode convention (a
deliberate, reader-validated deviation from WarpX's native RZ diagnostic schema).
"""

import os

import numpy as np
import easygdf
import openpmd_api as io

# openPMD unit_dimension dicts (SI base-unit exponents).
E_UNIT = {io.Unit_Dimension.M: 1.0, io.Unit_Dimension.L: 1.0,
          io.Unit_Dimension.T: -3.0, io.Unit_Dimension.I: -1.0}      # [V/m]
B_UNIT = {io.Unit_Dimension.M: 1.0,
          io.Unit_Dimension.T: -2.0, io.Unit_Dimension.I: -1.0}      # [T]


def load_cols(path, names):
    """Return the named flat columns from a GPT GDF field map."""
    d = easygdf.load(path)
    col = {b["name"]: np.asarray(b["value"]) for b in d["blocks"]}
    return [col[n] for n in names]


def to_grid(R, Z, *arrs, reverse_descending_z=False):
    """Reshape GDF flat columns (R fastest, Z slowest) to (nr, nz) grid arrays.

    With ``reverse_descending_z``, a GDF stored z-DESCENDING is row-reversed so the
    data ascends in z to match the ascending ``np.unique(Z)`` axis — otherwise the
    map is z-flipped vs its own axis, negating odd components (the prebuncher Er;
    see injector/README.md). Returns ``(r, z, *grids)``.
    """
    r = np.unique(R)
    z = np.unique(Z)
    nr, nz = r.size, z.size
    assert nr * nz == R.size, "field map is not a complete rectangular grid"
    grids = [a.reshape(nz, nr) for a in arrs]
    if reverse_descending_z and Z.reshape(nz, nr)[0, 0] > Z.reshape(nz, nr)[-1, 0]:
        grids = [a[::-1] for a in grids]
    grids = [a.T.copy() for a in grids]
    return (r, z, *grids)


def pad_r(r, rmax, *arrs):
    """Extend the (uniform-dr) r-grid with zero rows until it reaches ``rmax`` so
    r > the native map extent feels an exact zero field, not a WarpX extrapolation
    (no-op when the map already covers ``rmax``)."""
    dr = r[1] - r[0]
    if r[-1] >= rmax:
        return (r, *arrs)
    n_add = int(np.ceil((rmax - r[-1]) / dr))
    r_new = np.concatenate([r, r[-1] + dr * np.arange(1, n_add + 1)])
    out = [np.vstack([a, np.zeros((n_add, a.shape[1]))]) for a in arrs]
    return (r_new, *out)


def write_thetamode_series(out_file, r0, z0, dr, dz, meshes):
    """Write one openPMD thetaMode (m=0) RZ field file.

    ``meshes`` = ``[(name, [(component, (nr,nz) array), …], unit_dim_dict), …]``.
    ``r0``/``z0`` are the lab-frame coordinates of grid index 0 (grid_global_offset).
    """
    os.makedirs(os.path.dirname(out_file), exist_ok=True)
    series = io.Series(out_file, io.Access.create)
    it = series.iterations[0]
    for name, comps, unit_dim in meshes:
        m = it.meshes[name]
        m.geometry = io.Geometry.thetaMode
        m.geometry_parameters = "m=0;imag=+"
        m.axis_labels = ["r", "z"]
        m.grid_spacing = [dr, dz]
        m.grid_global_offset = [float(r0), float(z0)]
        m.grid_unit_SI = 1.0
        m.unit_dimension = unit_dim
        # thetaMode single (m=0) mode -> leading axis length 1 (np.newaxis). position
        # [0,0] is nodal centering for node-sampled GDF maps (WarpX's own RZ diags are
        # cell-centered [0.5,0.5]); deliberate — see gun/README.md.
        for cname, arr in comps:
            data = np.ascontiguousarray(arr[np.newaxis, :, :], dtype=np.float64)
            comp = m[cname]
            comp.position = [0.0, 0.0]
            comp.unit_SI = 1.0
            comp.reset_dataset(io.Dataset(data.dtype, data.shape))
            comp.store_chunk(data)
    series.flush()
    del series
