"""
Convert the SLAC 3 m linac (Section 1) + solenoid GPT field maps into openPMD
files that WarpX loads as externally applied fields.

Three openPMD files are written (all thetaMode, single azimuthal mode m = 0):

  * ``linac_rf1.h5`` / ``linac_rf2.h5`` — the **two quadrature (Re/Im) components
    of one 86-cell, 2π/3 traveling-wave SLAC accelerating structure**, driven at
    ``Linac_RF`` = 2856 MHz. GPT builds the traveling wave as the sum of two
    standing waves 90° apart (``reference/Linac Simulation Documentation/details.md``):

        Map25D_TM(... "SLAC-3mLinac-field1.gdf", "ErRe","EzRe","HphiIm", scale, 0, φ,        2π·f);
        Map25D_TM(... "SLAC-3mLinac-field2.gdf", "ErIm","EzIm","HphiRe", scale, 0, φ+0.5π,  2π·f);

    i.e. for each map  E(t) = map·scale·cos(ωt+φ),  Bφ(t) = map·scale·sin(ωt+φ),
    with field2 offset by +π/2. Re[ (Ẽ_re + i Ẽ_im) e^{i(ωt+φ)} ] is a forward
    traveling wave. We store the **raw, 1-kW-normalised** spatial maps here; the
    runtime ``scale = sqrt(P_MW/0.001)`` and the cos/sin(ωt+φ) modulation are
    applied in ``linac_sec1_sim.py`` via two ``picmi.LoadAppliedField`` objects.

  * ``linac_sol.h5`` — a **static** solenoid/lens focusing map (``Br, Bz``),
    per-Ampere normalised. Applied with a constant time function = the chosen
    current ``I_SOL`` (and no E field). Selectable via ``SOL_MAP`` (SOL_0 or
    LENS_0A…0E). This is the transverse focusing that lets the diverging ~148 keV
    prebuncher beam be captured into the 9.5 mm structure bore.

Layout WarpX's ``read_from_file`` reader expects (per mesh):
  geometry "thetaMode" (m=0); records "E"(r,t,z) and/or "B"(r,t,z);
  axisLabels ["r","z"]; dataset shape (1, nr, nz).

The RF maps only reach r ≈ 9.55 mm (the structure bore); they are **zero-padded
in r out to the sim domain RMAX** so every applied field covers the domain (WarpX
then sees an explicit zero field in the bore shadow rather than relying on its
out-of-range behaviour). The solenoid map already reaches 40 mm. Each map is
placed in the lab frame via ``grid_global_offset`` (Z_STRUCT, SOL_Z).

Run with:
    conda run -n CBB python -c "import linac_sec1; linac_sec1.run(plots=False)"
"""

import os
import numpy as np
import easygdf
import openpmd_api as io

# ── Inputs / outputs ─────────────────────────────────────────────────────────
RF1_GDF = "fieldmaps/SLAC-3mLinac-field1.gdf"
RF2_GDF = "fieldmaps/SLAC-3mLinac-field2.gdf"
OUT_DIR = "linac_sec1/linac_sec1_field"
RF1_FILE = os.path.join(OUT_DIR, "linac_rf1.h5")
RF2_FILE = os.path.join(OUT_DIR, "linac_rf2.h5")
SOL_FILE = os.path.join(OUT_DIR, "linac_sol.h5")

# Which solenoid/lens map to convert (SOL_0 or LENS_0A…0E). The .gdf path is
# rebuilt in main() so a config(SOL_MAP=...) override lands.
SOL_MAP = "SOL_0"

# RF operating point used only for the build-time gradient/gain report (the maps
# themselves are power-independent, 1-kW-normalised). Mirrors linac_sec1_sim.py so
# a config(POWER_MW=...) override makes the report track the actual run.
RF_NORM_MW = 0.001           # field-map power normalisation (1 kW)
POWER_MW = 15.0              # RF input power [MW]  (~37 MeV on crest)

# Traveling-wave 1-kW synchronous voltage ∫|Ez|dz of the committed SLAC maps
# (RF1_GDF + RF2_GDF), reported by the build below. Imported by linac_sec1_sim.py
# as the on-crest gain coefficient (gain = sqrt(P_MW/1e-3)·V1KW_KEV) so the sim's
# transit estimate stays in sync with the maps. Defined as a literal (not computed
# at import) to keep importing the module cheap; main() asserts it matches.
V1KW_KEV = 331.2             # [keV] = on-axis ∫|Ez|dz of the 1-kW maps

# ── Shared geometry (imported by linac_sec1_sim.py so field/phasing/domain agree) ─
# Lab-frame z of grid index 0 of each map (openPMD grid_global_offset). The SLAC
# map's own z runs −3.3…3012 mm; placing index 0 at Z_STRUCT puts the structure
# entrance there. Z_STRUCT also anchors the RF phase reference in the sim.
Z_STRUCT = 0.10              # [m] structure entrance (after a short injection drift)
# Solenoid map index-0 z. SOL_0 peaks 813 mm into its own grid; a negative offset
# slides that peak to ≈ lab z 0.21 m so the strongest focusing sits in the
# low-energy capture region just inside the structure entrance (the beam is most
# rigid-limited at ~148 keV there). The strongly-focusing peak sits inside the
# domain (lab z ≈ 0.21 m); the rising edge below lab z = 0 is clipped — note the
# clip plane is already at ~98% of peak |Bz|, so it is the map's far grid edge
# (native z = 0), not this clip plane, where Bz → 0.
SOL_Z = -0.60                # [m]
RMAX = 0.012                 # [m] sim radial domain; RF maps are zero-padded in r to here

# Electron-volt-free SI unit dimensions for the openPMD meshes.
E_UNIT = {io.Unit_Dimension.M: 1.0, io.Unit_Dimension.L: 1.0,
          io.Unit_Dimension.T: -3.0, io.Unit_Dimension.I: -1.0}      # [V/m]
B_UNIT = {io.Unit_Dimension.M: 1.0,
          io.Unit_Dimension.T: -2.0, io.Unit_Dimension.I: -1.0}      # [T]


def load_cols(path, names):
    """Return the named flat columns from a GPT GDF field map."""
    d = easygdf.load(path)
    col = {b["name"]: np.asarray(b["value"]) for b in d["blocks"]}
    return [col[n] for n in names]


def to_grid(R, Z, *arrs):
    """Reshape GDF flat columns (R fastest, Z slowest) to (nr, nz) grid arrays."""
    r = np.unique(R)
    z = np.unique(Z)
    nr, nz = r.size, z.size
    assert nr * nz == R.size, "field map is not a complete rectangular grid"
    out = [a.reshape(nz, nr).T.copy() for a in arrs]
    return (r, z, *out)


def pad_r(r, rmax, *arrs):
    """Extend the (uniform-dr) r-grid with zero rows until it reaches ``rmax``.

    Guarantees every RF applied field explicitly covers the sim domain, so a
    particle in the bore shadow (r > structure bore) feels an exact zero RF field
    rather than whatever WarpX would extrapolate past the map edge.
    """
    dr = r[1] - r[0]
    if r[-1] >= rmax:
        return (r, *arrs)
    n_add = int(np.ceil((rmax - r[-1]) / dr))
    r_new = np.concatenate([r, r[-1] + dr * np.arange(1, n_add + 1)])
    out = [np.vstack([a, np.zeros((n_add, a.shape[1]))]) for a in arrs]
    return (r_new, *out)


def write_series(out_file, z_offset, dr, dz, meshes):
    """Write one openPMD field file. ``meshes`` = [(name, [(comp, arr)], unit_dim)]."""
    os.makedirs(OUT_DIR, exist_ok=True)
    series = io.Series(out_file, io.Access.create)
    it = series.iterations[0]
    for name, comps, unit_dim in meshes:
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
    series.flush()
    del series


def _build_rf(gdf, ez_name, er_name, h_name, out_file):
    """Build one quadrature RF file; return (r, z, Ez_on_axis) for reporting."""
    R, Z, Er, Ez, Hphi = load_cols(gdf, ["R", "Z", er_name, ez_name, h_name])
    r, z, Er, Ez, Hphi = to_grid(R, Z, Er, Ez, Hphi)
    r, Er, Ez, Hphi = pad_r(r, RMAX, Er, Ez, Hphi)
    dr, dz = float(r[1] - r[0]), float(z[1] - z[0])
    zero = np.zeros_like(Er)
    # E uses cos(ωt+φ); Bφ (the H column) uses sin(ωt+φ) — supplied at runtime.
    write_series(out_file, Z_STRUCT, dr, dz, [
        ("E", (("r", Er), ("t", zero), ("z", Ez)), E_UNIT),
        ("B", (("r", zero), ("t", Hphi), ("z", zero)), B_UNIT),
    ])
    return r, z, Ez[0]


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # ── RF: two quadrature traveling-wave maps ────────────────────────────────
    r, z, ez1 = _build_rf(RF1_GDF, "EzRe", "ErRe", "HphiIm", RF1_FILE)
    _, _, ez2 = _build_rf(RF2_GDF, "EzIm", "ErIm", "HphiRe", RF2_FILE)
    nr, nz = r.size, z.size
    L = float(z[-1] - z[0])
    # Traveling-wave on-axis amplitude |Ez| = |EzRe + i EzIm|; its z-integral is the
    # 1-kW synchronous voltage. Physical gain (on crest) = sqrt(P_MW/1e-3)·this.
    env = np.sqrt(ez1**2 + ez2**2)
    v1kW = float(np.trapezoid(env, z))
    assert abs(v1kW / 1e3 - V1KW_KEV) < 0.5, (
        f"1-kW voltage {v1kW/1e3:.2f} kV drifted from V1KW_KEV={V1KW_KEV}; "
        "update the constant if the SLAC maps changed")
    print(f"SLAC Section 1 RF: nr={nr} (0–{r[-1]*1e3:.2f} mm, padded to RMAX="
          f"{RMAX*1e3:.0f} mm), nz={nz}, L={L:.3f} m, entrance at lab z={Z_STRUCT*1e3:.0f} mm")
    print(f"  peak on-axis |Ez| {env.max()/1e3:.2f} kV/m (1 kW); traveling-wave "
          f"1-kW voltage ∫|Ez|dz = {v1kW/1e3:.1f} kV")
    sc = np.sqrt(POWER_MW / RF_NORM_MW)
    print(f"  → at P={POWER_MW:g} MW (scale={sc:.1f}): peak gradient "
          f"{env.max()*sc/1e6:.2f} MV/m, on-crest gain ≈ {sc*v1kW/1e6:.1f} MeV")

    # ── Solenoid: static focusing map (Br, Bz), per-Ampere ────────────────────
    sol_gdf = f"fieldmaps/{SOL_MAP}.gdf"
    Rs, Zs, Br, Bz = load_cols(sol_gdf, ["R", "Z", "Br", "Bz"])
    rs, zs, Br, Bz = to_grid(Rs, Zs, Br, Bz)
    dr_s, dz_s = float(rs[1] - rs[0]), float(zs[1] - zs[0])
    zero_s = np.zeros_like(Br)
    write_series(SOL_FILE, SOL_Z, dr_s, dz_s, [
        ("B", (("r", Br), ("t", zero_s), ("z", Bz)), B_UNIT),
    ])
    bz_axis = Bz[0]
    ipk = int(np.argmax(np.abs(bz_axis)))
    print(f"Solenoid map '{SOL_MAP}': nr={rs.size} (0–{rs[-1]*1e3:.0f} mm), "
          f"nz={zs.size} (0–{zs[-1]*1e3:.0f} mm native), index-0 at lab z={SOL_Z*1e3:.0f} mm")
    print(f"  peak on-axis |Bz| {np.abs(bz_axis[ipk])*1e3:.4f} mT/A at native "
          f"z={zs[ipk]*1e3:.0f} mm; edge |Bz| {abs(bz_axis[0])*1e3:.4f}/"
          f"{abs(bz_axis[-1])*1e3:.4f} mT/A")
    print(f"\nWrote openPMD linac fields → {RF1_FILE}, {RF2_FILE}, {SOL_FILE}")


if __name__ == "__main__":
    main()
