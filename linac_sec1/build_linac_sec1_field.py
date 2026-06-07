"""
Convert the SLAC 3 m linac (Section 1) GPT field maps into openPMD files that WarpX
loads as externally applied fields.

Two openPMD files are written (both thetaMode, single azimuthal mode m = 0):

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

**The in-linac solenoid map (``linac_sol.h5``) was removed in the injector upgrade.**
Transverse focusing is now applied in the ``injector`` stage by the three real lenses
(Lens 0A / Sol 0 / Lens 0E) at their true lab z, and the injector hands the linac a
beam already focused and collimated to the 9.547 mm iris. The linac owns only the two
SLAC RF maps; it no longer carries ``SOL_FILE``/``SOL_MAP``/``SOL_Z``/``I_SOL``.

Layout WarpX's ``read_from_file`` reader expects (per mesh):
  geometry "thetaMode" (m=0); records "E"(r,t,z) and/or "B"(r,t,z);
  axisLabels ["r","z"]; dataset shape (1, nr, nz).

The RF maps only reach r ≈ 9.55 mm (the structure bore); they are **zero-padded in r
out to the sim domain RMAX (now 9.547 mm = the bore/iris)** so every applied field
covers the domain. Each map is placed in the lab frame via ``grid_global_offset``
(Z_STRUCT).

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
# NOTE: the in-linac solenoid map (linac_sol.h5) was removed in the injector upgrade.
# Transverse focusing now lives in the injector stage (Lens 0A / Sol 0 / Lens 0E at
# their true lab z ≈ 0.23/1.90/1.91 m); the linac no longer owns a solenoid field.

# RF operating point used only for the build-time gradient/gain report (the maps
# themselves are power-independent, 1-kW-normalised). Mirrors linac_sec1_sim.py so
# a config(POWER_MW=...) override makes the report track the actual run.
RF_NORM_MW = 0.001           # field-map power normalisation (1 kW)
POWER_MW = 11.0              # RF input power [MW] (sec1_input_power in the original LinacSim gpt_master.in)

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
# RMAX is now the SLAC bore / injector→linac collimator radius (9.547 mm). The injector
# delivers a beam already collimated to this iris (gpt scatteriris 9.547 mm at z=1.922 m
# + pipe), and the structure bore is ~9.55 mm, so the radial domain IS the aperture — a
# particle outside it is scraped at injection, exactly as the real machine does. (Was
# 12 mm, sized for the old blown-up, unfocused beam; the focused+collimated injector beam
# makes the faithful 9.547 mm correct. Do NOT widen it to contain a re-expanded envelope —
# that would accept charge the real iris scrapes and inflate the capture number.)
RMAX = 0.009547              # [m] sim radial domain = SLAC bore / collimator iris
BORE_R = 0.00955             # [m] structure bore radius (native r-extent of the SLAC maps);
                             # particles beyond this feel zero RF field. ≈ RMAX (the iris),
                             # so the bore and the aperture coincide. Single source of truth.

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

    # (The in-linac solenoid map was removed in the injector upgrade — transverse
    # focusing now lives in the injector stage at the lenses' true lab z. The linac
    # owns only the two SLAC quadrature RF maps.)
    print(f"\nWrote openPMD linac fields → {RF1_FILE}, {RF2_FILE}")


if __name__ == "__main__":
    main()
