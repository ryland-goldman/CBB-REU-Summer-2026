"""
Convert the CESR Linac prebuncher field map (`prebuncher_25D.gdf`) into an
openPMD file that WarpX can load as an externally applied RF-cavity field.

`prebuncher_25D.gdf` (read with easygdf) is a 2.5-D, axisymmetric (R, Z) map of a
standing-wave TM cavity, normalised to **1 J** of stored energy. Its columns are
`R, Z, Er, Ez, H`, where `H` is the azimuthal magnetic field Bφ in **Tesla** (the
A/m interpretation gives an unphysically negligible B; with H-as-Tesla the peak
E/cB ≈ 15, sensible for a resonant cavity). This is the same map the reference
GPT model drives with

    Map25D_TM(..., "Er","Ez","H", scale, 0, phi, 2*pi*f_RF)

i.e. Er,Ez(t) = map·scale·cos(ωt+φ) and Bφ(t) = H·scale·sin(ωt+φ) — E and B 90°
out of phase. We store the **raw, 1-J-normalised** spatial map here; the runtime
`scale`, `cos(ωt+φ)` / `sin(ωt+φ)` modulation are applied in `prebuncher_sim.py`
via `picmi.LoadAppliedField`.

The map is written in the openPMD layout WarpX's `read_from_file` external-field
reader expects for RZ geometry, with TWO meshes (E and B):

  * geometry           = "thetaMode" with a single azimuthal mode (m = 0)
  * mesh records       = "E" (r,t,z) and "B" (r,t,z)
  * axisLabels         = ["r", "z"]   (theta is the leading, size-1 axis)
  * dataset shape      = (1, nr, nz)

The map's native z runs [-152.4, +152.4] mm about the cavity gap; we set
`grid_global_offset` so it lands at `Z_GAP_CENTER` in the simulation lab frame.

Run with:
    conda run -n CBB python prebuncher/build_prebuncher_field.py
"""

import os
import numpy as np
import easygdf
import openpmd_api as io

# ── Inputs / outputs ─────────────────────────────────────────────────────────
GDF_PATH = "fieldmaps/prebuncher_25D.gdf"
OUT_DIR = "prebuncher/prebuncher_field"
OUT_FILE = os.path.join(OUT_DIR, "prebuncher_EB.h5")

# Lab-frame z of the cavity gap centre (the map is gap-centred at its own z=0).
# Imported by prebuncher_sim.py so the field placement and the beam phasing agree.
Z_GAP_CENTER = 0.20          # [m]
MAP_HALF_Z = 0.1524          # [m] half-length of the map (±152.4 mm)

# On-axis 1-J effective gap voltage ∫|Ez(r=0,z)|dz of the committed map, in keV.
# Imported by prebuncher_sim.py / plot_prebuncher.py as the gap-voltage coefficient
# (V_gap = scale · V1J_KEV) so the run, the transit estimate, and the plots stay in
# sync with the map. Defined as a literal (not computed at import) to keep importing
# the module cheap; main() asserts it matches the integral of the loaded map.
V1J_KEV = 438.6              # [keV]

# RF-drive constants from reference/Linac Simulation Documentation/details.md. Defined
# here (this module is pywarpx-free) as the single source of truth, imported by BOTH
# prebuncher_sim.py (the run) and plot_prebuncher.py (which re-derives the RF scale and
# phase for its waveform figure) so the two cannot drift. Mirrors the V1J_KEV pattern.
F_RF = 499.7645e6 / 42 * 18  # 18 × master RF = 214.18 MHz
Q_L = 3000                   # loaded Q of prebuncher 1


def load_prebuncher_map(path):
    """Return regular-grid (r, z, Er, Ez, Bphi) arrays from the GPT GDF map.

    The GDF stores flat columns with R varying fastest, then Z. Er, Ez are in
    V/m and H (= Bφ) in Tesla, for the 1-J-normalised cavity solution.
    """
    d = easygdf.load(path)
    col = {b["name"]: np.asarray(b["value"]) for b in d["blocks"]}
    R, Z, Er, Ez, H = col["R"], col["Z"], col["Er"], col["Ez"], col["H"]

    r = np.unique(R)
    z = np.unique(Z)
    nr, nz = r.size, z.size
    assert nr * nz == R.size, "field map is not a complete rectangular grid"

    # R fastest, Z slowest  ->  reshape to (nz, nr), then transpose to (nr, nz).
    Er = Er.reshape(nz, nr).T.copy()
    Ez = Ez.reshape(nz, nr).T.copy()
    Bphi = H.reshape(nz, nr).T.copy()
    return r, z, Er, Ez, Bphi


def main():
    r, z, Er, Ez, Bphi = load_prebuncher_map(GDF_PATH)
    nr, nz = r.size, z.size
    dr = float(r[1] - r[0])
    dz = float(z[1] - z[0])

    # On-axis 1-J effective gap voltage: V1J = ∫ |Ez(r=0, z)| dz. The physical
    # peak gap voltage at a given case is V_gap = scale · V1J (scale = sqrt(E/1J)).
    ez_axis = Ez[0]
    v1j = float(np.trapezoid(np.abs(ez_axis), z))
    ipk = int(np.argmax(np.abs(ez_axis)))
    assert abs(v1j / 1e3 - V1J_KEV) < 0.5, (
        f"1-J gap voltage {v1j/1e3:.2f} kV drifted from V1J_KEV={V1J_KEV}; "
        "update the constant if the map changed")

    print(f"Prebuncher map: nr={nr} (0–{r[-1]*1e3:.2f} mm), "
          f"nz={nz} ({z[0]*1e3:.1f}–{z[-1]*1e3:.1f} mm)")
    print(f"Peak |Ez| {np.abs(Ez).max()/1e6:.3f} MV/m, "
          f"peak |Er| {np.abs(Er).max()/1e6:.3f} MV/m, "
          f"peak |Bφ| {np.abs(Bphi).max()*1e3:.3f} mT  (1 J normalisation)")
    print(f"On-axis peak |Ez| {np.abs(ez_axis[ipk])/1e6:.3f} MV/m at "
          f"z={z[ipk]*1e3:.1f} mm; 1-J gap voltage V1J = {v1j/1e3:.2f} kV")
    print(f"Cavity gap placed at lab z = {Z_GAP_CENTER*1e3:.1f} mm "
          f"(map spans {(Z_GAP_CENTER-MAP_HALF_Z)*1e3:.1f}–"
          f"{(Z_GAP_CENTER+MAP_HALF_Z)*1e3:.1f} mm)")

    os.makedirs(OUT_DIR, exist_ok=True)
    series = io.Series(OUT_FILE, io.Access.create)
    it = series.iterations[0]

    # Shift the gap-centred map (native z0 = -MAP_HALF_Z) to the lab frame.
    z_offset = Z_GAP_CENTER - MAP_HALF_Z

    def write_mesh(name, comps, unit_dim):
        m = it.meshes[name]
        m.geometry = io.Geometry.thetaMode
        m.geometry_parameters = "m=0;imag=+"
        m.axis_labels = ["r", "z"]
        m.grid_spacing = [dr, dz]
        m.grid_global_offset = [0.0, z_offset]
        m.grid_unit_SI = 1.0
        m.unit_dimension = unit_dim
        # thetaMode with a single (m = 0) mode -> leading axis of length 1.
        for cname, arr in comps:
            data = np.ascontiguousarray(arr[np.newaxis, :, :], dtype=np.float64)
            comp = m[cname]
            comp.position = [0.0, 0.0]
            comp.unit_SI = 1.0
            comp.reset_dataset(io.Dataset(data.dtype, data.shape))
            comp.store_chunk(data)

    zero = np.zeros_like(Er)
    # Electric field  [V/m] = kg·m·s⁻³·A⁻¹
    write_mesh("E", (("r", Er), ("t", zero), ("z", Ez)), {
        io.Unit_Dimension.M: 1.0, io.Unit_Dimension.L: 1.0,
        io.Unit_Dimension.T: -3.0, io.Unit_Dimension.I: -1.0,
    })
    # Magnetic field  [T] = kg·s⁻²·A⁻¹ ; only the azimuthal (t) component (Bφ).
    write_mesh("B", (("r", zero), ("t", Bphi), ("z", zero)), {
        io.Unit_Dimension.M: 1.0,
        io.Unit_Dimension.T: -2.0, io.Unit_Dimension.I: -1.0,
    })

    series.flush()
    del series
    print(f"\nWrote openPMD prebuncher field (E + B) -> {OUT_FILE}")


if __name__ == "__main__":
    main()
