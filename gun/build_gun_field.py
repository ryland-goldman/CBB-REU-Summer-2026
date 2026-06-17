"""
Convert the CESR gun field map (`CESR_gun.gdf`) into an openPMD E-mesh that
WarpX loads as an externally applied electrode field (RZ, thetaMode m=0).

See gun/README.md for physics, voltage scaling, field-map layout, and gotchas.
"""

import os
import numpy as np
import easygdf
import openpmd_api as io

GDF_PATH = "fieldmaps/CESR_gun.gdf"
OUT_DIR = "gun/gun_field"
OUT_FILE = os.path.join(OUT_DIR, "gun_E.h5")

# Scale = -V_gun/V_map: NEGATIVE so electrons accelerate in +z (see README).
GUN_VOLTAGE = 150.0e3        # [V]
MAP_VOLTAGE = 1.0e3          # CESR_gun.gdf normalisation [V]


def load_gun_map(path):
    """Return regular-grid (r, z, Er, Ez) arrays from the GPT GDF field map."""
    d = easygdf.load(path)
    col = {b["name"]: np.asarray(b["value"]) for b in d["blocks"]}
    R, Z, Er, Ez = col["R"], col["Z"], col["Er"], col["Ez"]

    r = np.unique(R)
    z = np.unique(Z)
    nr, nz = r.size, z.size
    assert nr * nz == R.size, "field map is not a complete rectangular grid"
    # Origin must be (0,0): gun_sim.py/plot_gun.py assume axis + cathode plane.
    assert r[0] == 0.0 and z[0] == 0.0, (
        f"gun field map origin (r[0]={r[0]}, z[0]={z[0]}) must be (0, 0)")

    # GDF columns are R-fastest, Z-slowest -> reshape (nz, nr), transpose (nr, nz).
    Er = Er.reshape(nz, nr).T.copy()
    Ez = Ez.reshape(nz, nr).T.copy()
    return r, z, Er, Ez


def main():
    r, z, Er, Ez = load_gun_map(GDF_PATH)
    nr, nz = r.size, z.size
    dr = float(r[1] - r[0])
    dz = float(z[1] - z[0])

    # Compute here (not at import) so a config() GUN_VOLTAGE override applies.
    scale = -GUN_VOLTAGE / MAP_VOLTAGE
    Er = scale * Er
    Ez = scale * Ez
    Et = np.zeros_like(Er)

    print(f"Gun field map: nr={nr} (0–{r[-1]*1e3:.2f} mm), "
          f"nz={nz} (0–{z[-1]*1e3:.2f} mm)")
    ipk = np.argmax(np.abs(Ez[0]))
    print(f"Scaled by {scale:.0f}×  ->  -{GUN_VOLTAGE/1e3:.0f} kV cathode "
          f"(electrons accelerate in +z; Ez < 0 on axis)")
    print(f"On-axis Ez: cathode {Ez[0, 0]/1e6:.3f} MV/m, "
          f"peak {Ez[0, ipk]/1e6:.3f} MV/m at z={z[ipk]*1e3:.1f} mm")

    os.makedirs(OUT_DIR, exist_ok=True)
    series = io.Series(OUT_FILE, io.Access.create)
    it = series.iterations[0]

    E = it.meshes["E"]
    E.geometry = io.Geometry.thetaMode
    E.geometry_parameters = "m=0;imag=+"
    E.axis_labels = ["r", "z"]
    E.grid_spacing = [dr, dz]
    E.grid_global_offset = [float(r[0]), float(z[0])]
    E.grid_unit_SI = 1.0
    # E-field unit_dimension [V/m] = kg·m·s⁻³·A⁻¹
    E.unit_dimension = {
        io.Unit_Dimension.M: 1.0,
        io.Unit_Dimension.L: 1.0,
        io.Unit_Dimension.T: -3.0,
        io.Unit_Dimension.I: -1.0,
    }

    # thetaMode single m=0 mode -> leading axis of length 1 (the np.newaxis).
    for name, arr in (("r", Er), ("t", Et), ("z", Ez)):
        data = np.ascontiguousarray(arr[np.newaxis, :, :], dtype=np.float64)
        comp = E[name]
        comp.position = [0.0, 0.0]
        comp.unit_SI = 1.0
        comp.reset_dataset(io.Dataset(data.dtype, data.shape))
        comp.store_chunk(data)

    series.flush()
    del series
    print(f"\nWrote openPMD gun field -> {OUT_FILE}")


if __name__ == "__main__":
    main()
