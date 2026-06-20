"""
Build the injector's openPMD external-field maps from GPT .gdf sources: two
prebuncher RF cavities (preb1/preb2_EB.h5, same forward field, different lab-z) and
six static B-only solenoid lenses.

See injector/README.md for physics, the reversed-Preb-2 / PREB2_REV_PHASE reasoning,
the solenoid placement, and gotchas.

Run with:
    conda run -n CBB python injector/build_injector_field.py
"""

import os
import numpy as np
import easygdf
import openpmd_api as io

GDF_PATH = "fieldmaps/prebuncher_25D.gdf"
OUT_DIR = "injector/injector_field"
# Both cavities use the FORWARD field; Preb 2 differs only in lab-z (grid_global_offset)
# and the run-time +π. read_from_file reads spatial position from the file, so each gap z
# needs its own file even with identical field VALUES.
OUT_FILE_1 = os.path.join(OUT_DIR, "preb1_EB.h5")
OUT_FILE_2 = os.path.join(OUT_DIR, "preb2_EB.h5")

# Solenoid lenses (static B-only, per-Ampere maps), in z-order. The grids differ
# (LENS_0A nr=189/nz=16, others nr=16/nz~601) so each is its own file.
SOL_NAMES = ("LENS_0A", "LENS_0B", "LENS_0C", "LENS_0D", "SOL_0", "LENS_0E")
SOL_GDF = {n: f"fieldmaps/{n}.gdf" for n in SOL_NAMES}
SOL_FILES = {n: os.path.join(OUT_DIR, n.lower().replace("_", "") + ".h5")
             for n in SOL_NAMES}
# GUI lab-z annotation of each lens (gpt_master.in "@GUI element" lines). For the thin
# lenses this is the field peak (cross-checked below); for the flat-top SOL_0 it is a
# center/edge label, NOT the argmax — placement is native absolute z, not GUI-aligned.
SOL_GUI_Z = {"LENS_0A": 0.225, "LENS_0B": 1.603, "LENS_0C": 1.692,
             "LENS_0D": 1.838, "SOL_0": 1.897, "LENS_0E": 1.914}
SOL_TOL = 0.001              # [m] floor for the |lab_peak − GUI_z| narrow-lens cross-check

Z_HANDOFF = 2.03             # [m] linac handoff plane; every focusing peak must be upstream
ZMAX = 2.10                  # [m] injector domain end
RMAX = 0.036                 # [m] sim radial domain

# Lab-frame z of each cavity gap centre (the map is gap-centred at its own z=0).
Z_GAP_CENTER_1 = 0.534       # [m] Prebuncher 1
Z_GAP_CENTER_2 = 1.318       # [m] Prebuncher 2
MAP_HALF_Z = 0.1524          # [m] half-length of the map (±152.4 mm)

Z_GAP_CENTER = Z_GAP_CENTER_1   # back-compat alias

# On-axis 1-J gap voltage ∫|Ez(r=0,z)|dz, keV (V_gap = scale·V1J_KEV); literal kept
# cheap, main() asserts it matches the loaded map's integral.
V1J_KEV = 438.6              # [keV]

# RF-drive constants (this module is pywarpx-free; single source of truth imported by
# injector_sim.py and plot_injector.py).
F_RF = 499.7645e6 / 42 * 18  # 18 × master RF = 214.18 MHz
Q_L_1 = 3000                 # loaded Q of prebuncher 1
Q_L_2 = 4300                 # loaded Q of prebuncher 2
Q_L = Q_L_1                  # back-compat alias


def load_prebuncher_map(path):
    """Return regular-grid (r, z, Er, Ez, Bphi) arrays from the GPT GDF map.

    GDF flat columns: R fastest, then Z. Er, Ez in V/m; the H column is Bφ in
    Tesla (NOT A/m — the A/m reading gives a negligible B; see README).
    """
    d = easygdf.load(path)
    col = {b["name"]: np.asarray(b["value"]) for b in d["blocks"]}
    R, Z, Er, Ez, H = col["R"], col["Z"], col["Er"], col["Ez"], col["H"]

    r = np.unique(R)
    z = np.unique(Z)
    nr, nz = r.size, z.size
    assert nr * nz == R.size, "field map is not a complete rectangular grid"

    # R fastest, Z slowest  ->  reshape to (nz, nr).  prebuncher_25D.gdf stores z
    # DESCENDING (+152.4 → −152.4 mm), so the flat rows must be reversed to match
    # the ascending `np.unique(Z)` axis BEFORE transposing — otherwise the map is
    # z-flipped relative to its axis, which negates the odd Er and corrupts the
    # off-axis transverse RF force (Ez/Bphi are even, so they would be unaffected).
    Er = Er.reshape(nz, nr)
    Ez = Ez.reshape(nz, nr)
    Bphi = H.reshape(nz, nr)
    if Z.reshape(nz, nr)[0, 0] > Z.reshape(nz, nr)[-1, 0]:
        Er, Ez, Bphi = Er[::-1], Ez[::-1], Bphi[::-1]
    Er = Er.T.copy()
    Ez = Ez.T.copy()
    Bphi = Bphi.T.copy()
    return r, z, Er, Ez, Bphi


def write_field(out_file, r, z, Er, Ez, Bphi, z_gap):
    """Write one prebuncher openPMD field file (E + B meshes), placed at ``z_gap``.

    Forward 1-J field written verbatim; Preb 2's reversal is a run-time +π, not a
    negation here (see README). z_offset = z_gap − MAP_HALF_Z lands the gap (native
    z=0) at lab z_gap; the assertion enforces the required ±MAP_HALF_Z symmetry.
    """
    nr, nz = r.size, z.size
    dr = float(r[1] - r[0])
    dz = float(z[1] - z[0])

    assert abs(z[0] + z[-1]) < 1e-6, (
        f"map z-extent must be symmetric about the gap (z spans "
        f"[{z[0]*1e3:.3f}, {z[-1]*1e3:.3f}] mm); z_offset=z_gap−MAP_HALF_Z "
        f"assumes a ±MAP_HALF_Z map")

    z_offset = z_gap - MAP_HALF_Z

    os.makedirs(OUT_DIR, exist_ok=True)
    series = io.Series(out_file, io.Access.create)
    it = series.iterations[0]

    def write_mesh(name, comps, unit_dim):
        m = it.meshes[name]
        m.geometry = io.Geometry.thetaMode
        m.geometry_parameters = "m=0;imag=+"
        m.axis_labels = ["r", "z"]
        m.grid_spacing = [dr, dz]
        m.grid_global_offset = [0.0, z_offset]
        m.grid_unit_SI = 1.0
        m.unit_dimension = unit_dim
        # thetaMode single (m=0) mode -> leading axis length 1; axisLabels ["r","z"].
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


B_UNIT = {io.Unit_Dimension.M: 1.0,
          io.Unit_Dimension.T: -2.0, io.Unit_Dimension.I: -1.0}      # [T]


def _sol_to_grid(R, Z, *arrs):
    """Reshape GDF flat columns (R fastest, Z slowest) to (nr, nz) grid arrays."""
    r = np.unique(R)
    z = np.unique(Z)
    nr, nz = r.size, z.size
    assert nr * nz == R.size, "solenoid map is not a complete rectangular grid"
    out = [a.reshape(nz, nr).T.copy() for a in arrs]
    return (r, z, *out)


def _sol_pad_r(r, rmax, *arrs):
    """Extend the (uniform-dr) r-grid with zero rows until it reaches ``rmax``
    (no-op when the map already covers it; robust against a future smaller map)."""
    dr = r[1] - r[0]
    if r[-1] >= rmax:
        return (r, *arrs)
    n_add = int(np.ceil((rmax - r[-1]) / dr))
    r_new = np.concatenate([r, r[-1] + dr * np.arange(1, n_add + 1)])
    out = [np.vstack([a, np.zeros((n_add, a.shape[1]))]) for a in arrs]
    return (r_new, *out)


def _write_b_series(out_file, z_offset, dr, dz, Br, Bz):
    """Write one single-mesh (B only) openPMD field file in the WarpX RZ layout."""
    os.makedirs(OUT_DIR, exist_ok=True)
    series = io.Series(out_file, io.Access.create)
    it = series.iterations[0]
    m = it.meshes["B"]
    m.geometry = io.Geometry.thetaMode
    m.geometry_parameters = "m=0;imag=+"
    m.axis_labels = ["r", "z"]
    m.grid_spacing = [dr, dz]
    m.grid_global_offset = [0.0, z_offset]
    m.grid_unit_SI = 1.0
    m.unit_dimension = B_UNIT
    zero = np.zeros_like(Br)
    for cname, arr in (("r", Br), ("t", zero), ("z", Bz)):
        data = np.ascontiguousarray(arr[np.newaxis, :, :], dtype=np.float64)
        comp = m[cname]
        comp.position = [0.0, 0.0]
        comp.unit_SI = 1.0
        comp.reset_dataset(io.Dataset(data.dtype, data.shape))
        comp.store_chunk(data)
    series.flush()
    del series


def build_solenoids():
    """Build the per-Ampere B-only solenoid maps, each placed in the lab frame.

    Placement is NATIVE absolute machine-z, matching gpt_master.in, which installs every
    solenoid with `Map2D_B("wcs", "z", 0.0, …)` — i.e. the GDF's stored Z column IS the
    absolute machine z, with NO peak-alignment shift. grid_global_offset is the lab-z of
    grid INDEX 0, so it is simply the native origin z[0] (0 for LENS_0A/SOL_0, 0.8 m for
    the pre-shifted LENS_0B…0E grids).

    Do NOT align argmax→GUI z: that is wrong for the flat-top SOL_0, whose argmax (0.813 m)
    is an arbitrary point on its [0.35, 1.87] m plateau — forcing it to GUI 1.897 m shifts
    the whole channel +1.08 m and leaves PB1/PB2 unfocused. The thin lenses' native peaks
    already match their GUI z (LENS_0A 0.233≈0.225, LENS_0E 1.915≈1.914), confirming native
    z is absolute. See README → Solenoid lenses.

    The READ-BACK assertions below re-read the WRITTEN placement (not an input recompute,
    which can't catch a bad offset).
    """
    for name in SOL_NAMES:
        d = easygdf.load(SOL_GDF[name])
        col = {b["name"]: np.asarray(b["value"]) for b in d["blocks"]}
        R, Z, Br, Bz = col["R"], col["Z"], col["Br"], col["Bz"]
        r, z, Br, Bz = _sol_to_grid(R, Z, Br, Bz)
        r, Br, Bz = _sol_pad_r(r, RMAX, Br, Bz)
        dr, dz = float(r[1] - r[0]), float(z[1] - z[0])
        ipk = int(np.argmax(np.abs(Bz[0])))             # peak on the axis (r=0 row)
        z_peak_native = float(z[ipk])
        gui_z = SOL_GUI_Z[name]
        # native absolute placement: grid index 0 sits at the GDF's own z[0]
        offset = float(z[0])

        _write_b_series(SOL_FILES[name], offset, dr, dz, Br, Bz)

        # READ-BACK guard: STORED peak + FWHM from grid_global_offset + argmax·dz.
        chk = io.Series(SOL_FILES[name], io.Access.read_only)
        mB = chk.iterations[0].meshes["B"]
        off_z, ddz = float(mB.grid_global_offset[1]), float(mB.grid_spacing[1])
        bz_stored = np.asarray(mB["z"].load_chunk())
        chk.flush()
        bz_axis = bz_stored[0][0]                        # (1,nr,nz) -> r=0 row
        lab_peak = off_z + int(np.argmax(np.abs(bz_axis))) * ddz
        half = np.flatnonzero(np.abs(bz_axis) >= 0.5 * np.abs(bz_axis).max())
        fwhm = float((half[-1] - half[0]) * ddz)
        del chk

        assert 0.0 <= lab_peak <= ZMAX, (
            f"{name} STORED lab-z peak {lab_peak*1e3:.1f} mm outside the [0, {ZMAX*1e3:.0f}] mm domain")
        assert lab_peak < Z_HANDOFF, (
            f"{name} STORED lab-z peak {lab_peak*1e3:.1f} mm is NOT upstream of the "
            f"{Z_HANDOFF*1e3:.0f} mm handoff plane — the linac would inherit a beam still "
            f"inside this lens")
        # Peak≈GUI z is a loose sanity ONLY for narrow (peaked) lenses; the flat-top SOL_0's
        # argmax is an arbitrary plateau point and its GUI z is a center/edge annotation, so
        # the cross-check is skipped there. Tol = half the FWHM: the GUI annotation is only
        # accurate to the field's own width (native vs GUI differ up to ~14 mm here), while
        # an 800 mm-class placement bug is still well outside.
        if fwhm < 0.30:
            tol = max(SOL_TOL, 0.5 * fwhm)
            assert abs(lab_peak - gui_z) < tol, (
                f"{name} STORED lab-z peak {lab_peak*1e3:.2f} mm differs from GUI z {gui_z*1e3:.1f} mm "
                f"by more than {tol*1e3:.1f} mm (grid_global_offset bug?)")

        # Report: per-Ampere peak |Bz| (mT/A), native placement, FWHM (flags flat-top).
        print(f"Solenoid {name}: nr={r.size} nz={z.size}, native z=[{z[0]*1e3:.0f},{z[-1]*1e3:.0f}] mm, "
              f"offset={offset*1e3:+.1f} mm -> STORED lab-z peak {lab_peak*1e3:.1f} mm "
              f"(GUI {gui_z*1e3:.1f}, FWHM {fwhm*1e3:.0f} mm), "
              f"peak |Bz| {abs(Bz[0][ipk])*1e3:.4f} mT/A -> {SOL_FILES[name]}")


def main():
    r, z, Er, Ez, Bphi = load_prebuncher_map(GDF_PATH)
    nr, nz = r.size, z.size
    dr = float(r[1] - r[0])
    dz = float(z[1] - z[0])

    # On-axis 1-J gap voltage V1J = ∫|Ez(r=0,z)|dz; V_gap = scale·V1J.
    ez_axis = Ez[0]
    v1j = float(np.trapezoid(np.abs(ez_axis), z))
    ipk = int(np.argmax(np.abs(ez_axis)))
    assert abs(v1j / 1e3 - V1J_KEV) < 0.5, (
        f"1-J gap voltage {v1j/1e3:.2f} kV drifted from V1J_KEV={V1J_KEV}; "
        "update the constant if the map changed")

    # Gap parity sanity check on the source .gdf: the reversed-Preb-2 reasoning assumes
    # Ez EVEN, Er ODD, Bφ EVEN about the gap (180° rotation flips all 3 = global E,B sign
    # flip ≡ +π absolute phase). Whether that +π is applied as PREB2_REV_PHASE is
    # convention-dependent: at the zc default it is applied explicitly (PREB2_REV_PHASE=π);
    # under the legacy crest+GUI convention it was absorbed by crest-referencing (rev=0). See
    # README -> Reversed install. A future map breaking parity invalidates this.
    def _parity(arr):
        f = arr[::-1]
        denom = float(np.sqrt((arr * arr).sum() * (f * f).sum()))
        return float((arr * f).sum() / denom) if denom > 0 else 0.0
    jr = int(np.argmax(np.abs(Er).max(axis=1)))      # row with largest |Er|
    p_ez, p_er, p_bphi = _parity(Ez[jr]), _parity(Er[jr]), _parity(Bphi[jr])
    assert p_ez > 0.99, f"Ez not EVEN about the gap (corr {p_ez:+.4f}); the reversed-install reasoning assumes Ez EVEN"
    assert p_er < -0.99, f"Er not ODD about the gap (corr {p_er:+.4f}); the reversed-install reasoning assumes Er ODD"
    assert p_bphi > 0.99, f"Bφ not EVEN about the gap (corr {p_bphi:+.4f}); the reversed-install reasoning assumes Bφ EVEN (TM0)"

    # Raw-GDF orientation check: parity (above) is invariant under a z-flip, so it
    # CANNOT catch a row-order/axis mismatch that negates the odd Er. Compare the
    # stored map against the raw GDF column at a fixed off-axis +z point: the sign
    # and magnitude must agree, guarding the descending-z row reversal in
    # load_prebuncher_map.
    _raw = easygdf.load(GDF_PATH)
    _rawcol = {b["name"]: np.asarray(b["value"]) for b in _raw["blocks"]}
    _Rr, _Zr, _Err = _rawcol["R"], _rawcol["Z"], _rawcol["Er"]
    _iz = int(np.argmax(z))                          # +z end of the ascending axis
    _ir = nr // 2                                    # an off-axis row (Er ≠ 0)
    _m = (np.isclose(_Rr, r[_ir]) & np.isclose(_Zr, z[_iz]))
    assert _m.sum() == 1, "could not locate the (r,+z) sample in the raw GDF"
    _er_raw = float(_Err[_m][0])
    assert np.sign(Er[_ir, _iz]) == np.sign(_er_raw) and \
        abs(Er[_ir, _iz] - _er_raw) < 1e-3 * max(abs(_er_raw), 1.0), (
        f"stored Er at (r={r[_ir]*1e3:.1f}mm, z=+{z[_iz]*1e3:.1f}mm)={Er[_ir, _iz]:.3f} "
        f"disagrees with raw GDF {_er_raw:.3f} — z row order / axis mismatch (Er sign flip)")

    print(f"Prebuncher map: nr={nr} (0–{r[-1]*1e3:.2f} mm), "
          f"nz={nz} ({z[0]*1e3:.1f}–{z[-1]*1e3:.1f} mm)")
    print(f"Peak |Ez| {np.abs(Ez).max()/1e6:.3f} MV/m, "
          f"peak |Er| {np.abs(Er).max()/1e6:.3f} MV/m, "
          f"peak |Bφ| {np.abs(Bphi).max()*1e3:.3f} mT  (1 J normalisation)")
    print(f"On-axis peak |Ez| {np.abs(ez_axis[ipk])/1e6:.3f} MV/m at "
          f"z={z[ipk]*1e3:.1f} mm; 1-J gap voltage V1J = {v1j/1e3:.2f} kV")
    print(f"Gap parity (z-flip corr, peak-|Er| row r={r[jr]*1e3:.1f}mm): "
          f"Ez {p_ez:+.4f} (EVEN), Er {p_er:+.4f} (ODD), Bφ {p_bphi:+.4f} (EVEN) "
          f"→ 180° rotation flips all 3 = +π in ABSOLUTE phase (applied explicitly as "
          f"PREB2_REV_PHASE=π at the zc default; absorbed by crest-ref under legacy crest)")

    write_field(OUT_FILE_1, r, z, Er, Ez, Bphi, Z_GAP_CENTER_1)
    print(f"Prebuncher 1 gap at lab z = {Z_GAP_CENTER_1*1e3:.1f} mm "
          f"(field spans {(Z_GAP_CENTER_1-MAP_HALF_Z)*1e3:.1f}–"
          f"{(Z_GAP_CENTER_1+MAP_HALF_Z)*1e3:.1f} mm) -> {OUT_FILE_1}")

    # Preb 2 reuses the SAME forward field; reversal applied at run time.
    write_field(OUT_FILE_2, r, z, Er, Ez, Bphi, Z_GAP_CENTER_2)
    print(f"Prebuncher 2 gap at lab z = {Z_GAP_CENTER_2*1e3:.1f} mm "
          f"(field spans {(Z_GAP_CENTER_2-MAP_HALF_Z)*1e3:.1f}–"
          f"{(Z_GAP_CENTER_2+MAP_HALF_Z)*1e3:.1f} mm; reversal = run-time phase) -> {OUT_FILE_2}")

    print(f"\nWrote openPMD injector prebuncher fields (E + B) -> "
          f"{OUT_FILE_1}, {OUT_FILE_2}")

    # ── Solenoid lenses (static B-only focusing maps) ─────────────────────────
    build_solenoids()
    print(f"Wrote openPMD injector solenoid fields (B) -> "
          f"{', '.join(SOL_FILES.values())}")


if __name__ == "__main__":
    main()
