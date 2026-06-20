"""
Convert the CESR gun field map (`CESR_gun.gdf`) into an openPMD E-mesh that
WarpX loads as an externally applied electrode field (RZ, thetaMode m=0).

See gun/README.md for physics, voltage scaling, field-map layout, and gotchas.
"""

import os
import numpy as np

from pipeline.fieldio import E_UNIT, load_cols, to_grid, write_thetamode_series

GDF_PATH = "fieldmaps/CESR_gun.gdf"
OUT_DIR = "gun/gun_field"
OUT_FILE = os.path.join(OUT_DIR, "gun_E.h5")

# Scale = -V_gun/V_map: NEGATIVE so electrons accelerate in +z (see README).
GUN_VOLTAGE = 150.0e3        # [V]
MAP_VOLTAGE = 1.0e3          # CESR_gun.gdf normalisation [V]


def load_gun_map(path):
    """Return regular-grid (r, z, Er, Ez) arrays from the GPT GDF field map."""
    R, Z, Er, Ez = load_cols(path, ["R", "Z", "Er", "Ez"])
    r, z, Er, Ez = to_grid(R, Z, Er, Ez)
    # Origin must be (0,0): gun_sim.py/plot_gun.py assume axis + cathode plane.
    assert r[0] == 0.0 and z[0] == 0.0, (
        f"gun field map origin (r[0]={r[0]}, z[0]={z[0]}) must be (0, 0)")
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

    write_thetamode_series(OUT_FILE, float(r[0]), float(z[0]), dr, dz, [
        ("E", (("r", Er), ("t", Et), ("z", Ez)), E_UNIT),
    ])
    print(f"\nWrote openPMD gun field -> {OUT_FILE}")


if __name__ == "__main__":
    main()
