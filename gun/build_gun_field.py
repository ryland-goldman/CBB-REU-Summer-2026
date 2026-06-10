"""
Convert the CESR gun Poisson–Superfish field map (`CESR_gun.gdf`) into an
openPMD file that WarpX can load as an externally applied electrode field.

`CESR_gun.gdf` (read with easygdf) is a 2D cylindrical (R, Z) map of the
electrostatic gun field, computed at a **1000 V** cathode→exit potential drop.
We scale it to the real CESR "Chili Gun Mk II" voltage (~150 kV) and write it
in the openPMD layout WarpX's `read_from_file` external-field reader expects for
RZ geometry:

  * geometry           = "thetaMode" with a single azimuthal mode (m = 0)
  * mesh record        = "E", components "r", "t", "z"
  * axisLabels         = ["r", "z"]   (theta is the leading, size-1 axis)
  * dataset shape      = (1, nr, nz)

The gun is purely electrostatic (the GDF carries no magnetic field), so only an
E mesh is written; `gun_sim.py` sets `particles.B_ext_particle_init_style = none`.

Run with:
    conda run -n CBB python gun/build_gun_field.py
"""

import os
import numpy as np
import easygdf
import openpmd_api as io

# ── Inputs / outputs ─────────────────────────────────────────────────────────
GDF_PATH = "fieldmaps/CESR_gun.gdf"
OUT_DIR = "gun/gun_field"
OUT_FILE = os.path.join(OUT_DIR, "gun_E.h5")

# The map is normalised to a +1 kV cathode (V = +1000 at the cathode, 0 at the
# exit), so its on-axis Ez = -dV/dz is POSITIVE — which would push electrons
# back into the cathode. The physical gun holds the cathode at NEGATIVE high
# voltage (anode grounded), so electrons accelerate in +z. We therefore scale by
# a negative factor: real field = map × (-V_gun / V_map).
GUN_VOLTAGE = 150.0e3        # CESR Chili Gun Mk II cathode potential magnitude [V]
MAP_VOLTAGE = 1.0e3          # normalisation of CESR_gun.gdf [V]


def load_gun_map(path):
    """Return regular-grid (r, z, Er, Ez) arrays from the GPT GDF field map.

    The GDF stores flat columns with R varying fastest, then Z. Er, Ez are in
    V/m for the 1 kV-normalised solution.
    """
    d = easygdf.load(path)
    col = {b["name"]: np.asarray(b["value"]) for b in d["blocks"]}
    R, Z, Er, Ez = col["R"], col["Z"], col["Er"], col["Ez"]

    r = np.unique(R)
    z = np.unique(Z)
    nr, nz = r.size, z.size
    assert nr * nz == R.size, "field map is not a complete rectangular grid"
    # gun_sim.py (RMAX/ZMAX, grid) and plot_gun.py assume the map starts at the
    # axis and the cathode plane; a swapped-in map with a nonzero native origin
    # would otherwise be silently mis-placed.
    assert r[0] == 0.0 and z[0] == 0.0, (
        f"gun field map origin (r[0]={r[0]}, z[0]={z[0]}) must be (0, 0)")

    # R fastest, Z slowest  ->  reshape to (nz, nr), then transpose to (nr, nz).
    Er = Er.reshape(nz, nr).T.copy()
    Ez = Ez.reshape(nz, nr).T.copy()
    return r, z, Er, Ez


def main():
    r, z, Er, Ez = load_gun_map(GDF_PATH)
    nr, nz = r.size, z.size
    dr = float(r[1] - r[0])
    dz = float(z[1] - z[0])

    # Recompute the scale here (not at import) so a config() override of
    # GUN_VOLTAGE takes effect on the written field map.  = -150 by default.
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
    E.grid_global_offset = [float(r[0]), float(z[0])]   # native origin (asserted (0, 0))
    E.grid_unit_SI = 1.0
    # Electric field  [V/m] = kg·m·s⁻³·A⁻¹
    E.unit_dimension = {
        io.Unit_Dimension.M: 1.0,
        io.Unit_Dimension.L: 1.0,
        io.Unit_Dimension.T: -3.0,
        io.Unit_Dimension.I: -1.0,
    }

    # thetaMode with a single (m = 0) mode -> leading axis of length 1.
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
