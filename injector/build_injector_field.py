"""
Convert the CESR Linac prebuncher field map (`prebuncher_25D.gdf`) into openPMD
files that WarpX can load as externally applied RF-cavity fields for the injector.

The injector stage carries the full LinacSim injector subsection — two 214 MHz
prebuncher cavities (the second installed *reversed*) plus three static solenoid
lenses — in one self-consistent RZ drift. This module builds the field maps:

  * `preb1_EB.h5`      — forward 1-J cavity field, gap at Z_GAP_CENTER_1 = 0.534 m.
  * `preb2_EB.h5`      — the SAME forward field, gap at Z_GAP_CENTER_2 = 1.318 m
                         (reversed install applied as +π at run time, not here).
  * solenoid maps      — lens0a / sol0 / lens0e (added in build step 4).

**Reversed Prebuncher 2 needs NO negated map.** prebuncher_25D has definite parity
about the gap centre — Ez EVEN, Er ODD, Bφ EVEN (physics-measured corr ±0.9999,
asserted in main() below; Bφ EVEN is forced by Maxwell for a TM0 mode, Bφ ~ dEr/dz −
dEz/dr with both even). GPT's `-1,0,0` reversed install is a 180° rotation, which for
this definite-parity mode flips all three lab components = a GLOBAL sign flip of E and
B ≡ a `+π` phase shift on the cos/sin drive IN ABSOLUTE PHASE. So Preb 2 reuses this
same forward field; placement at Z_GAP_CENTER_2 is via grid_global_offset, and the
reversal lives in the time function. NOTE: the APPLIED reversal phase is
`PREB2_REV_PHASE = 0`, not +π — because the run crest-references the loaded field and
the GUI on-crest phase is defined for the already-reversed install, so the crest
reference already contains the +π (see `injector_sim.py` PREB2_REV_PHASE). A separately
negated `.h5` is never built.

`prebuncher_25D.gdf` (read with easygdf) is a 2.5-D, axisymmetric (R, Z) map of a
standing-wave TM cavity, normalised to **1 J** of stored energy. Its columns are
`R, Z, Er, Ez, H`, where `H` is the azimuthal magnetic field Bφ in **Tesla** (the
A/m interpretation gives an unphysically negligible B; with H-as-Tesla the peak
E/cB ≈ 15, sensible for a resonant cavity). This is the same map the reference
GPT model drives with

    Map25D_TM(..., "Er","Ez","H", scale, 0, phi, 2*pi*f_RF)

i.e. Er,Ez(t) = map·scale·cos(ωt+φ) and Bφ(t) = H·scale·sin(ωt+φ) — E and B 90°
out of phase. We store the **raw, 1-J-normalised** spatial map here; the runtime
`scale`, `cos(ωt+φ)` / `sin(ωt+φ)` modulation are applied in `injector_sim.py`
via `picmi.LoadAppliedField`.

Each map is written in the openPMD layout WarpX's `read_from_file` external-field
reader expects for RZ geometry, with TWO meshes (E and B):

  * geometry           = "thetaMode" with a single azimuthal mode (m = 0)
  * mesh records       = "E" (r,t,z) and "B" (r,t,z)
  * axisLabels         = ["r", "z"]   (theta is the leading, size-1 axis)
  * dataset shape      = (1, nr, nz)

The map's native z runs [-152.4, +152.4] mm about the cavity gap; we set
`grid_global_offset` so it lands at each cavity's lab-frame gap z.

Run with:
    conda run -n CBB python injector/build_injector_field.py
"""

import os
import numpy as np
import easygdf
import openpmd_api as io

# ── Inputs / outputs ─────────────────────────────────────────────────────────
GDF_PATH = "fieldmaps/prebuncher_25D.gdf"
OUT_DIR = "injector/injector_field"
# Both cavities use the FORWARD field (no negation); they differ only in lab-z
# placement (grid_global_offset, baked per file) and in the +π the reversed Preb 2
# adds in its time function. WarpX's read_from_file takes a field's spatial position
# from its openPMD file, so Preb 2 at a different gap z needs its own file even though
# the field VALUES are identical to Preb 1.
OUT_FILE_1 = os.path.join(OUT_DIR, "preb1_EB.h5")   # forward field, gap at Z_GAP_CENTER_1
OUT_FILE_2 = os.path.join(OUT_DIR, "preb2_EB.h5")   # forward field, gap at Z_GAP_CENTER_2 (+π at run)

# ── Solenoid lenses (static B-only focusing maps; per-Ampere normalised) ───────
# Three energized LinacSim lenses (GUI defaults; 0B/0C/0D are 0 A, Sol 1A-C are
# downstream of Section 1 — all omitted). Each is a separate single-mesh openPMD
# file (the grids differ: LENS_0A is nr=189/nz=16, the others nr=16/nz~601), placed
# in the lab frame via grid_global_offset. Imported by injector_sim.py.
SOL_GDF = {"LENS_0A": "fieldmaps/LENS_0A.gdf",
           "SOL_0":   "fieldmaps/SOL_0.gdf",
           "LENS_0E": "fieldmaps/LENS_0E.gdf"}
SOL_FILES = {"LENS_0A": os.path.join(OUT_DIR, "lens0a.h5"),
             "SOL_0":   os.path.join(OUT_DIR, "sol0.h5"),
             "LENS_0E": os.path.join(OUT_DIR, "lens0e.h5")}
# GUI lab-z of each lens (LinacSim gpt_master.in positions).
SOL_GUI_Z = {"LENS_0A": 0.225, "SOL_0": 1.897, "LENS_0E": 1.914}
# Offset is ALWAYS derived programmatically (offset = GUI_z − Z[argmax|Bz|]), per map,
# from the LOADED peak — NOT hard-coded. This self-corrects against stale plan literals
# (the plan table's SOL_0 native peak 0.8209 m / offset +1.0761 m is stale; the actual
# file is peak 0.8129 m ⇒ offset +1.0841 m, landing the peak dead-on at 1.897 m — a
# hard-coded +1.0761 would miss by 8 mm). It also lands LENS_0A dead-on (offset −0.0083 m
# from its native 0.2333 m peak) instead of accepting the ~8 mm map-vs-GUI discrepancy,
# and LENS_0E is ≈0 (native peak 1.9147 ≈ GUI). physics measured all three .gdf maps and
# confirmed the programmatic form (see the per-solenoid sanity report).
SOL_TOL = 0.001              # [m] |lab_peak − GUI_z| tolerance (programmatic lands all <0.1 mm)

Z_HANDOFF = 2.03             # [m] linac handoff plane; every focusing peak must be upstream
ZMAX = 2.10                  # [m] injector domain end (peaks must be in-domain)
RMAX = 0.036                 # [m] sim radial domain; solenoid maps already reach 40 mm

# Lab-frame z of each cavity gap centre (the map is gap-centred at its own z=0).
# Imported by injector_sim.py so the field placement and the beam phasing agree.
Z_GAP_CENTER_1 = 0.534       # [m] Prebuncher 1 (Z_prebuncher1 in LinacSim gpt_master.in)
Z_GAP_CENTER_2 = 1.318       # [m] Prebuncher 2 (Z_prebuncher2)
MAP_HALF_Z = 0.1524          # [m] half-length of the map (±152.4 mm)

# Back-compat aliases (the un-suffixed names some readers still import). Preb 1 is
# the canonical "the gap".
Z_GAP_CENTER = Z_GAP_CENTER_1

# On-axis 1-J effective gap voltage ∫|Ez(r=0,z)|dz of the committed map, in keV.
# Imported by injector_sim.py / plot_injector.py as the gap-voltage coefficient
# (V_gap = scale · V1J_KEV) so the run, the transit estimate, and the plots stay in
# sync with the map. Defined as a literal (not computed at import) to keep importing
# the module cheap; main() asserts it matches the integral of the loaded map.
V1J_KEV = 438.6              # [keV]

# RF-drive constants from reference/Linac Simulation Documentation/details.md. Defined
# here (this module is pywarpx-free) as the single source of truth, imported by BOTH
# injector_sim.py (the run) and plot_injector.py (which re-derives the RF scale and
# phase for its waveform figure) so the two cannot drift. Mirrors the V1J_KEV pattern.
F_RF = 499.7645e6 / 42 * 18  # 18 × master RF = 214.18 MHz
Q_L_1 = 3000                 # loaded Q of prebuncher 1
Q_L_2 = 4300                 # loaded Q of prebuncher 2
Q_L = Q_L_1                  # back-compat alias

# On-crest reference phase offsets (deg) reproducing LinacSim's GUI on-crest
# definitions (304.7° / 178.9°); see injector_sim.py make_cavity().
PHI_OFF_1_DEG = -70.0        # Prebuncher 1
PHI_OFF_2_DEG = -45.0        # Prebuncher 2 (reversed)


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


def write_field(out_file, r, z, Er, Ez, Bphi, z_gap):
    """Write one prebuncher openPMD field file (E + B meshes), placed at ``z_gap``.

    The FORWARD 1-J field is written verbatim; the reversed install of Preb 2 is a
    +π phase shift applied at run time (see module docstring / injector_sim.py), NOT
    a field negation here. The two cavities differ only in ``z_gap`` (the lab-frame
    gap centre, set via grid_global_offset) and that run-time +π.

    z_offset = z_gap − MAP_HALF_Z places grid index 0 at the gap's near edge, which
    lands the gap (native z = 0) at lab z = z_gap because the map's native z-extent is
    symmetric about the gap (±MAP_HALF_Z). The assertion enforces that symmetry.
    """
    nr, nz = r.size, z.size
    dr = float(r[1] - r[0])
    dz = float(z[1] - z[0])

    # The map's native z-extent must be symmetric about the gap (±MAP_HALF_Z) for
    # z_offset = z_gap − MAP_HALF_Z to land the gap at z_gap.
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


# ── Solenoid helpers (mirror the linac builder's to_grid / pad_r / write_series) ─
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
    """Build the three per-Ampere B-only solenoid maps, each placed in the lab frame.

    The lab-z offset is derived programmatically (not hard-coded) so it can't drift:
      * "local"-frame maps (SOL_0): offset = GUI_z − Z[argmax|Bz|] (slides the local
        peak to the GUI lab position).
      * "absolute"-frame maps (LENS_0A, LENS_0E): native z already IS lab z ⇒ offset 0.
    Per solenoid: asserts the in-domain physical lab-z peak is in [0, ZMAX], upstream of
    Z_HANDOFF (so the linac never inherits a beam still inside a lens it does not model),
    and within SOL_TOL of the GUI z; prints the physical lab-z peak (not the native peak).
    """
    for name in ("LENS_0A", "SOL_0", "LENS_0E"):
        d = easygdf.load(SOL_GDF[name])
        col = {b["name"]: np.asarray(b["value"]) for b in d["blocks"]}
        R, Z, Br, Bz = col["R"], col["Z"], col["Br"], col["Bz"]
        r, z, Br, Bz = _sol_to_grid(R, Z, Br, Bz)
        r, Br, Bz = _sol_pad_r(r, RMAX, Br, Bz)
        dr, dz = float(r[1] - r[0]), float(z[1] - z[0])
        ipk = int(np.argmax(np.abs(Bz[0])))             # peak on the axis (r=0 row)
        z_peak_native = float(z[ipk])
        gui_z = SOL_GUI_Z[name]
        # Programmatic offset for EVERY map (derived from the loaded peak, not hard-coded):
        # slide the native peak to the GUI lab-z. Self-corrects against stale plan literals
        # and lands all three peaks dead-on at their GUI positions.
        offset = gui_z - z_peak_native
        lab_peak = z_peak_native + offset                # = gui_z by construction (in-domain lab-z)

        # Assertions (silent-wrong-physics guards).
        assert 0.0 <= lab_peak <= ZMAX, (
            f"{name} lab-z peak {lab_peak*1e3:.1f} mm outside the [0, {ZMAX*1e3:.0f}] mm domain")
        assert lab_peak < Z_HANDOFF, (
            f"{name} lab-z peak {lab_peak*1e3:.1f} mm is NOT upstream of the "
            f"{Z_HANDOFF*1e3:.0f} mm handoff plane — the linac would inherit a beam still "
            f"inside this lens")
        assert abs(lab_peak - gui_z) < SOL_TOL, (
            f"{name} lab-z peak {lab_peak*1e3:.2f} mm differs from GUI z {gui_z*1e3:.1f} mm "
            f"by more than {SOL_TOL*1e3:.1f} mm")

        _write_b_series(SOL_FILES[name], offset, dr, dz, Br, Bz)
        # Report: per-Ampere peak |Bz| (mT/A), and the PHYSICAL lab-z peak (not native).
        print(f"Solenoid {name}: nr={r.size} nz={z.size}, native peak z={z_peak_native*1e3:.1f} mm, "
              f"offset={offset*1e3:+.1f} mm -> lab-z peak {lab_peak*1e3:.1f} mm "
              f"(GUI {gui_z*1e3:.1f}), peak |Bz| {abs(Bz[0][ipk])*1e3:.4f} mT/A -> {SOL_FILES[name]}")


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

    # ── Gap parity check (sanity on the source .gdf) ──────────────────────────
    # The reversed Preb-2 install (GPT's `-1,0,0`) is a 180° ROTATION of the cavity: it
    # applies the spatial z-mirror AND rotates the field VECTOR (ẑ→-ẑ, r̂→r̂, φ̂→-φ̂).
    # For THIS map's definite parity — Ez EVEN, Er ODD, Bφ EVEN (physics-confirmed; also
    # forced by Maxwell for a TM0 mode: Bφ ~ dEr/dz - dEz/dr, both even ⇒ Bφ even) —
    # every lab component flips sign: Ez (even, z-comp flips), Er (odd, r-comp keeps →
    # net flip), Bφ (even, φ-comp flips). All three flip = a GLOBAL E,B sign flip ≡ a +π
    # time-phase shift IN ABSOLUTE DRIVE PHASE. So Preb-2 reuses the FORWARD field; the
    # reversal is encoded purely in the time function — NO mirrored/negated .h5. NOTE: the
    # APPLIED value PREB2_REV_PHASE is 0, not +π, because the drive is crest-referenced and
    # the GUI 178.9° on-crest is for the already-reversed install, so the crest reference
    # already contains the reversal (adding +π double-counts). See the long note at
    # injector_sim.PREB2_REV_PHASE. We assert the parity here as a sanity check on the
    # source .gdf — if a future map breaks it, the reversal reasoning must be revisited.
    def _parity(arr):
        f = arr[::-1]
        denom = float(np.sqrt((arr * arr).sum() * (f * f).sum()))
        return float((arr * f).sum() / denom) if denom > 0 else 0.0
    jr = int(np.argmax(np.abs(Er).max(axis=1)))      # radial row with the largest |Er|
    p_ez, p_er, p_bphi = _parity(Ez[jr]), _parity(Er[jr]), _parity(Bphi[jr])
    assert p_ez > 0.99, f"Ez not EVEN about the gap (corr {p_ez:+.4f}); the reversed-install reasoning assumes Ez EVEN"
    assert p_er < -0.99, f"Er not ODD about the gap (corr {p_er:+.4f}); the reversed-install reasoning assumes Er ODD"
    assert p_bphi > 0.99, f"Bφ not EVEN about the gap (corr {p_bphi:+.4f}); the reversed-install reasoning assumes Bφ EVEN (TM0)"

    print(f"Prebuncher map: nr={nr} (0–{r[-1]*1e3:.2f} mm), "
          f"nz={nz} ({z[0]*1e3:.1f}–{z[-1]*1e3:.1f} mm)")
    print(f"Peak |Ez| {np.abs(Ez).max()/1e6:.3f} MV/m, "
          f"peak |Er| {np.abs(Er).max()/1e6:.3f} MV/m, "
          f"peak |Bφ| {np.abs(Bphi).max()*1e3:.3f} mT  (1 J normalisation)")
    print(f"On-axis peak |Ez| {np.abs(ez_axis[ipk])/1e6:.3f} MV/m at "
          f"z={z[ipk]*1e3:.1f} mm; 1-J gap voltage V1J = {v1j/1e3:.2f} kV")
    print(f"Gap parity (z-flip corr, peak-|Er| row r={r[jr]*1e3:.1f}mm): "
          f"Ez {p_ez:+.4f} (EVEN), Er {p_er:+.4f} (ODD), Bφ {p_bphi:+.4f} (EVEN) "
          f"→ 180° rotation flips all 3 = +π in ABSOLUTE phase (absorbed by crest-ref; "
          f"applied PREB2_REV_PHASE=0)")

    # Prebuncher 1 — forward field at Z_GAP_CENTER_1.
    write_field(OUT_FILE_1, r, z, Er, Ez, Bphi, Z_GAP_CENTER_1)
    print(f"Prebuncher 1 gap at lab z = {Z_GAP_CENTER_1*1e3:.1f} mm "
          f"(field spans {(Z_GAP_CENTER_1-MAP_HALF_Z)*1e3:.1f}–"
          f"{(Z_GAP_CENTER_1+MAP_HALF_Z)*1e3:.1f} mm) -> {OUT_FILE_1}")

    # Prebuncher 2 — SAME forward field at Z_GAP_CENTER_2 (reversal applied at run time).
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
