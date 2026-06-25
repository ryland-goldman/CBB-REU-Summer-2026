"""GDF -> openPMD field-map builders for every stage that loads an applied field.

Consolidates the old pipeline/fieldio.py primitives and the four per-stage
build_*_field.py scripts. Reads GPT .gdf maps from fieldmaps/gdf/ and writes openPMD
thetaMode (RZ, m=0) .h5 maps to fieldmaps/h5/. The SLAC traveling-wave maps are built
once and shared by linac sections 1-3.

Map-geometry constants (gap centres, bore, 1-J/1-kW voltages) live here as facts of the
fixed GDF inputs; tunable operating-point values (gun voltage, RF power, frequency, Q) live
in the config YAMLs and are passed in or read by the stage drivers.

Axis-order / m-mode convention (axis ["r","z"], m=0, nodal position, V/m & T unit
dimensions) is a deliberate, reader-validated deviation from WarpX's native RZ diag schema.
"""

import os

import numpy as np
import easygdf
import openpmd_api as io

GDF_DIR = "fieldmaps/gdf"
H5_DIR = "fieldmaps/h5"

# openPMD unit_dimension dicts (SI base-unit exponents).
E_UNIT = {io.Unit_Dimension.M: 1.0, io.Unit_Dimension.L: 1.0,
          io.Unit_Dimension.T: -3.0, io.Unit_Dimension.I: -1.0}      # [V/m]
B_UNIT = {io.Unit_Dimension.M: 1.0,
          io.Unit_Dimension.T: -2.0, io.Unit_Dimension.I: -1.0}      # [T]


# ── GDF parse / grid / write primitives ──────────────────────────────────────────
def load_cols(path, names):
    """Return the named flat columns from a GPT GDF field map."""
    d = easygdf.load(path)
    col = {b["name"]: np.asarray(b["value"]) for b in d["blocks"]}
    return [col[n] for n in names]


def to_grid(R, Z, *arrs, reverse_descending_z=False):
    """Reshape GDF flat columns (R fastest, Z slowest) to (nr, nz) grid arrays.

    With reverse_descending_z, a GDF stored z-DESCENDING is row-reversed so the data
    ascends in z to match the ascending np.unique(Z) axis (else odd components like the
    prebuncher Er are negated). Returns (r, z, *grids).
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
    """Extend the (uniform-dr) r-grid with zero rows to rmax so r > the native map extent
    feels an exact zero field, not a WarpX extrapolation (no-op if already covers rmax)."""
    dr = r[1] - r[0]
    if r[-1] >= rmax:
        return (r, *arrs)
    n_add = int(np.ceil((rmax - r[-1]) / dr))
    r_new = np.concatenate([r, r[-1] + dr * np.arange(1, n_add + 1)])
    out = [np.vstack([a, np.zeros((n_add, a.shape[1]))]) for a in arrs]
    return (r_new, *out)


def write_thetamode_series(out_file, r0, z0, dr, dz, meshes):
    """Write one openPMD thetaMode (m=0) RZ field file.

    meshes = [(name, [(component, (nr,nz) array), ...], unit_dim_dict), ...].
    r0/z0 are the lab-frame coordinates of grid index 0 (grid_global_offset).
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
        for cname, arr in comps:
            data = np.ascontiguousarray(arr[np.newaxis, :, :], dtype=np.float64)
            comp = m[cname]
            comp.position = [0.0, 0.0]              # nodal centering for node-sampled GDF maps
            comp.unit_SI = 1.0
            comp.reset_dataset(io.Dataset(data.dtype, data.shape))
            comp.store_chunk(data)
    series.flush()
    del series


# ── Gun electrode field ──────────────────────────────────────────────────────────
GUN_GDF = os.path.join(GDF_DIR, "CESR_gun.gdf")
GUN_OUT = os.path.join(H5_DIR, "gun_E.h5")
GUN_MAP_VOLTAGE = 1.0e3                              # CESR_gun.gdf normalisation [V]


def build_gun_field(gun_voltage):
    """Build the gun applied-E map. scale = -gun_voltage/MAP_VOLTAGE (NEGATIVE so electrons
    accelerate in +z). `gun_voltage` [V] comes from config/gun.yaml."""
    R, Z, Er, Ez = load_cols(GUN_GDF, ["R", "Z", "Er", "Ez"])
    r, z, Er, Ez = to_grid(R, Z, Er, Ez)
    assert r[0] == 0.0 and z[0] == 0.0, (
        f"gun field map origin (r[0]={r[0]}, z[0]={z[0]}) must be (0, 0)")
    dr, dz = float(r[1] - r[0]), float(z[1] - z[0])
    scale = -gun_voltage / GUN_MAP_VOLTAGE
    Er, Ez = scale * Er, scale * Ez
    Et = np.zeros_like(Er)
    write_thetamode_series(GUN_OUT, float(r[0]), float(z[0]), dr, dz, [
        ("E", (("r", Er), ("t", Et), ("z", Ez)), E_UNIT)])
    ipk = int(np.argmax(np.abs(Ez[0])))
    print(f"Gun field: nr={r.size} nz={z.size}, scaled {scale:.0f}x -> "
          f"-{gun_voltage/1e3:.0f} kV cathode; on-axis peak Ez "
          f"{Ez[0, ipk]/1e6:.3f} MV/m at z={z[ipk]*1e3:.1f} mm -> {GUN_OUT}", flush=True)


# ── Injector prebunchers + solenoids ─────────────────────────────────────────────
PREB_GDF = os.path.join(GDF_DIR, "prebuncher_25D.gdf")
PREB1_OUT = os.path.join(H5_DIR, "preb1_EB.h5")
PREB2_OUT = os.path.join(H5_DIR, "preb2_EB.h5")

Z_GAP_CENTER_1 = 0.534          # [m] Prebuncher 1 gap centre (lab z)
Z_GAP_CENTER_2 = 1.318          # [m] Prebuncher 2 gap centre (lab z)
MAP_HALF_Z = 0.1524             # [m] half-length of the prebuncher map (+-152.4 mm)
V1J_KEV = 438.6                 # [keV] 1-J on-axis gap voltage int|Ez|dz (asserted in build)
INJ_Z_HANDOFF = 2.03            # [m] linac handoff plane; every focusing peak must be upstream
INJ_ZMAX = 2.10                 # [m] injector domain end
INJ_RMAX = 0.036                # [m] injector radial domain

SOL_NAMES = ("LENS_0A", "LENS_0B", "LENS_0C", "LENS_0D", "SOL_0", "LENS_0E")
SOL_GDF = {n: os.path.join(GDF_DIR, f"{n}.gdf") for n in SOL_NAMES}
SOL_FILES = {n: os.path.join(H5_DIR, n.lower().replace("_", "") + ".h5") for n in SOL_NAMES}
# GUI lab-z annotation per lens. For thin lenses this is the field peak (cross-checked);
# for the flat-top SOL_0 it is a centre/edge label -- placement is native absolute z.
SOL_GUI_Z = {"LENS_0A": 0.225, "LENS_0B": 1.603, "LENS_0C": 1.692,
             "LENS_0D": 1.838, "SOL_0": 1.897, "LENS_0E": 1.914}
SOL_TOL = 0.001                 # [m] floor for the |lab_peak - GUI_z| narrow-lens cross-check


def _load_prebuncher_map():
    """Regular-grid (r, z, Er, Ez, Bphi) from the prebuncher GDF. The map stores z DESCENDING
    (+152.4 -> -152.4 mm); reverse_descending_z row-reverses so the odd Er is not negated."""
    R, Z, Er, Ez, H = load_cols(PREB_GDF, ["R", "Z", "Er", "Ez", "H"])  # H is Bphi [T]
    return to_grid(R, Z, Er, Ez, H, reverse_descending_z=True)


def _write_prebuncher(out_file, r, z, Er, Ez, Bphi, z_gap):
    """Write one prebuncher openPMD file (E + B), gap (native z=0) placed at lab z_gap."""
    dr, dz = float(r[1] - r[0]), float(z[1] - z[0])
    assert abs(z[0] + z[-1]) < 1e-6, (
        f"prebuncher map z-extent must be symmetric about the gap (z spans "
        f"[{z[0]*1e3:.3f}, {z[-1]*1e3:.3f}] mm)")
    z_offset = z_gap - MAP_HALF_Z
    zero = np.zeros_like(Er)
    write_thetamode_series(out_file, 0.0, z_offset, dr, dz, [
        # E uses cos(wt+phi); Bphi (the H column) uses sin(wt+phi) -- supplied at runtime.
        ("E", (("r", Er), ("t", zero), ("z", Ez)), E_UNIT),
        ("B", (("r", zero), ("t", Bphi), ("z", zero)), B_UNIT)])


def _build_solenoids():
    """Build the per-Ampere B-only solenoid maps, each at its NATIVE absolute machine-z
    (matching gpt_master.in: no peak-alignment shift). grid_global_offset = native z[0]."""
    for name in SOL_NAMES:
        R, Z, Br, Bz = load_cols(SOL_GDF[name], ["R", "Z", "Br", "Bz"])
        r, z, Br, Bz = to_grid(R, Z, Br, Bz)
        r, Br, Bz = pad_r(r, INJ_RMAX, Br, Bz)
        dr, dz = float(r[1] - r[0]), float(z[1] - z[0])
        ipk = int(np.argmax(np.abs(Bz[0])))
        offset = float(z[0])                        # native absolute placement
        zero = np.zeros_like(Br)
        write_thetamode_series(SOL_FILES[name], 0.0, offset, dr, dz, [
            ("B", (("r", Br), ("t", zero), ("z", Bz)), B_UNIT)])

        # READ-BACK guard: STORED peak + FWHM from grid_global_offset + argmax*dz.
        chk = io.Series(SOL_FILES[name], io.Access.read_only)
        mB = chk.iterations[0].meshes["B"]
        off_z, ddz = float(mB.grid_global_offset[1]), float(mB.grid_spacing[1])
        bz_stored = np.asarray(mB["z"].load_chunk())
        chk.flush()
        bz_axis = bz_stored[0][0]
        lab_peak = off_z + int(np.argmax(np.abs(bz_axis))) * ddz
        half = np.flatnonzero(np.abs(bz_axis) >= 0.5 * np.abs(bz_axis).max())
        fwhm = float((half[-1] - half[0]) * ddz)
        del chk

        assert 0.0 <= lab_peak <= INJ_ZMAX, (
            f"{name} STORED lab-z peak {lab_peak*1e3:.1f} mm outside [0, {INJ_ZMAX*1e3:.0f}] mm")
        assert lab_peak < INJ_Z_HANDOFF, (
            f"{name} STORED lab-z peak {lab_peak*1e3:.1f} mm not upstream of the "
            f"{INJ_Z_HANDOFF*1e3:.0f} mm handoff plane")
        # Peak~=GUI z is a loose sanity ONLY for narrow (FWHM < 0.30 m) lenses; broad/flat-top
        # maps (SOL_0, LENS_0B) are exempt -- their argmax is an arbitrary plateau point.
        if fwhm < 0.30:
            tol = max(SOL_TOL, 0.5 * fwhm)
            assert abs(lab_peak - SOL_GUI_Z[name]) < tol, (
                f"{name} STORED lab-z peak {lab_peak*1e3:.2f} mm differs from GUI z "
                f"{SOL_GUI_Z[name]*1e3:.1f} mm by > {tol*1e3:.1f} mm (offset bug?)")
        print(f"Solenoid {name}: lab-z peak {lab_peak*1e3:.1f} mm (FWHM {fwhm*1e3:.0f} mm), "
              f"peak |Bz| {abs(Bz[0][ipk])*1e3:.4f} mT/A -> {SOL_FILES[name]}", flush=True)


def build_injector_fields():
    """Build the two prebuncher RF maps (same forward field, different lab-z) and the six
    per-Ampere solenoid maps. Parity / orientation asserts guard the reversed-Preb-2 reasoning."""
    r, z, Er, Ez, Bphi = _load_prebuncher_map()
    nr, nz = r.size, z.size

    ez_axis = Ez[0]
    v1j = float(np.trapezoid(np.abs(ez_axis), z))
    assert abs(v1j / 1e3 - V1J_KEV) < 0.5, (
        f"1-J gap voltage {v1j/1e3:.2f} kV drifted from V1J_KEV={V1J_KEV}")

    # Gap parity: Ez EVEN, Er ODD, Bphi EVEN about the gap -- so reversing this symmetric gap is a
    # +pi phase flip (the injector's PREB2_REV_PHASE=pi), NOT a z-flip (which inverts only odd Er).
    def _parity(arr):
        f = arr[::-1]
        denom = float(np.sqrt((arr * arr).sum() * (f * f).sum()))
        return float((arr * f).sum() / denom) if denom > 0 else 0.0
    jr = int(np.argmax(np.abs(Er).max(axis=1)))
    p_ez, p_er, p_bphi = _parity(Ez[jr]), _parity(Er[jr]), _parity(Bphi[jr])
    assert p_ez > 0.99, f"Ez not EVEN about the gap (corr {p_ez:+.4f})"
    assert p_er < -0.99, f"Er not ODD about the gap (corr {p_er:+.4f})"
    assert p_bphi > 0.99, f"Bphi not EVEN about the gap (corr {p_bphi:+.4f})"

    # Raw-GDF orientation check (parity is z-flip invariant): stored Er sign/magnitude at a
    # fixed off-axis +z point must match the raw GDF column (guards the descending-z reversal).
    _Rr, _Zr, _Err = load_cols(PREB_GDF, ["R", "Z", "Er"])
    _iz, _ir = int(np.argmax(z)), nr // 2
    _m = (np.isclose(_Rr, r[_ir]) & np.isclose(_Zr, z[_iz]))
    assert _m.sum() == 1, "could not locate the (r,+z) sample in the raw GDF"
    _er_raw = float(_Err[_m][0])
    assert np.sign(Er[_ir, _iz]) == np.sign(_er_raw) and \
        abs(Er[_ir, _iz] - _er_raw) < 1e-3 * max(abs(_er_raw), 1.0), (
        "stored Er disagrees with raw GDF -- z row order / axis mismatch (Er sign flip)")

    print(f"Prebuncher map: nr={nr} nz={nz}, 1-J gap voltage {v1j/1e3:.2f} kV; "
          f"parity Ez {p_ez:+.3f} Er {p_er:+.3f} Bphi {p_bphi:+.3f}", flush=True)
    _write_prebuncher(PREB1_OUT, r, z, Er, Ez, Bphi, Z_GAP_CENTER_1)
    _write_prebuncher(PREB2_OUT, r, z, Er, Ez, Bphi, Z_GAP_CENTER_2)   # reversal applied at runtime
    print(f"Prebunchers at lab z = {Z_GAP_CENTER_1*1e3:.0f}, {Z_GAP_CENTER_2*1e3:.0f} mm "
          f"-> {PREB1_OUT}, {PREB2_OUT}", flush=True)
    _build_solenoids()


# ── SLAC traveling-wave maps (linac sections 1-3, shared) ────────────────────────
SLAC_RF1_GDF = os.path.join(GDF_DIR, "SLAC-3mLinac-field1.gdf")
SLAC_RF2_GDF = os.path.join(GDF_DIR, "SLAC-3mLinac-field2.gdf")
SLAC_RF1_OUT = os.path.join(H5_DIR, "linac_rf1.h5")
SLAC_RF2_OUT = os.path.join(H5_DIR, "linac_rf2.h5")

Z_STRUCT = 0.10                 # [m] lab-z of grid index 0 = structure entrance; anchors RF phase
RMAX = 0.009547                 # [m] sim radial domain = SLAC bore / collimator iris (IS the aperture)
BORE_R = 0.00955                # [m] native map r-extent; r > this feels zero RF
V1KW_KEV = 331.2                # [keV] on-axis int|Ez|dz of the 1-kW maps (asserted in build)


def onaxis_quadrature_ez():
    """(z_maplocal [m], Ez1_onaxis, Ez2_onaxis) [V/m] of the two SLAC quadrature halves (r=0),
    z from the structure entrance. Used by the frozen-setpoint derivation (scratchpad)."""
    R, Z, Ez1 = load_cols(SLAC_RF1_GDF, ["R", "Z", "EzRe"])
    r, z, Ez1 = to_grid(R, Z, Ez1)
    R, Z, Ez2 = load_cols(SLAC_RF2_GDF, ["R", "Z", "EzIm"])
    _, _, Ez2 = to_grid(R, Z, Ez2)
    return z - z[0], Ez1[0], Ez2[0]


def _build_slac_rf(gdf, ez_name, er_name, h_name, out_file):
    """Build one SLAC quadrature RF file; return (r, z, Ez_on_axis) for reporting."""
    R, Z, Er, Ez, Hphi = load_cols(gdf, ["R", "Z", er_name, ez_name, h_name])
    r, z, Er, Ez, Hphi = to_grid(R, Z, Er, Ez, Hphi)
    r, Er, Ez, Hphi = pad_r(r, RMAX, Er, Ez, Hphi)
    dr, dz = float(r[1] - r[0]), float(z[1] - z[0])
    zero = np.zeros_like(Er)
    write_thetamode_series(out_file, 0.0, Z_STRUCT, dr, dz, [
        ("E", (("r", Er), ("t", zero), ("z", Ez)), E_UNIT),
        ("B", (("r", zero), ("t", Hphi), ("z", zero)), B_UNIT)])
    return r, z, Ez[0]


def build_linac_slac_fields():
    """Build the two SLAC quadrature RF files (1-kW normalised) shared by linac sections 1-3."""
    r, z, ez1 = _build_slac_rf(SLAC_RF1_GDF, "EzRe", "ErRe", "HphiIm", SLAC_RF1_OUT)
    _, _, ez2 = _build_slac_rf(SLAC_RF2_GDF, "EzIm", "ErIm", "HphiRe", SLAC_RF2_OUT)
    env = np.sqrt(ez1 ** 2 + ez2 ** 2)
    v1kW = float(np.trapezoid(env, z))
    assert abs(v1kW / 1e3 - V1KW_KEV) < 0.5, (
        f"1-kW voltage {v1kW/1e3:.2f} kV drifted from V1KW_KEV={V1KW_KEV}")
    print(f"SLAC RF: nr={r.size} nz={z.size}, entrance at lab z={Z_STRUCT*1e3:.0f} mm, "
          f"1-kW voltage int|Ez|dz = {v1kW/1e3:.1f} kV -> {SLAC_RF1_OUT}, {SLAC_RF2_OUT}",
          flush=True)
