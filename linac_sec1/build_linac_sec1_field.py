"""
Convert the SLAC 3 m linac (Section 1) GPT field maps into the two openPMD RF files
(both thetaMode, m=0) that WarpX loads as externally applied fields: linac_rf1.h5 /
linac_rf2.h5, the Re/Im quadrature halves of one traveling-wave structure. The raw
1-kW-normalised spatial maps are stored here; the runtime power scale and cos/sin(ωt+φ)
modulation are applied in linac_sec1_sim.py.

See linac_sec1/README.md for physics, parameters, field-map provenance, and gotchas.
"""

import os
import numpy as np

from pipeline.fieldio import (B_UNIT, E_UNIT, load_cols, pad_r, to_grid,
                              write_thetamode_series)

# ── Inputs / outputs ─────────────────────────────────────────────────────────
RF1_GDF = "fieldmaps/SLAC-3mLinac-field1.gdf"
RF2_GDF = "fieldmaps/SLAC-3mLinac-field2.gdf"
OUT_DIR = "linac_sec1/linac_sec1_field"
RF1_FILE = os.path.join(OUT_DIR, "linac_rf1.h5")
RF2_FILE = os.path.join(OUT_DIR, "linac_rf2.h5")

# RF operating point: build-time gradient/gain report only (maps are 1-kW-normalised).
RF_NORM_MW = 0.001           # field-map power normalisation (1 kW)
POWER_MW = 11.0              # RF input power [MW]

V1KW_KEV = 331.2             # [keV] on-axis ∫|Ez|dz of the 1-kW maps; literal so import
                             # stays cheap (main() asserts it matches the built maps).

# Shared geometry (imported by linac_sec1_sim.py so field/phasing/domain agree):
Z_STRUCT = 0.10              # [m] lab-frame z of grid index 0 = structure entrance; anchors RF phase
# RMAX IS the aperture (SLAC bore / iris); do NOT widen to contain a re-expanded
# envelope — that accepts charge the real iris scrapes and inflates capture.
RMAX = 0.009547              # [m] sim radial domain = SLAC bore / collimator iris
BORE_R = 0.00955             # [m] structure bore (native map r-extent); r>this feels zero RF

def _build_rf(gdf, ez_name, er_name, h_name, out_file):
    """Build one quadrature RF file; return (r, z, Ez_on_axis) for reporting."""
    R, Z, Er, Ez, Hphi = load_cols(gdf, ["R", "Z", er_name, ez_name, h_name])
    r, z, Er, Ez, Hphi = to_grid(R, Z, Er, Ez, Hphi)
    r, Er, Ez, Hphi = pad_r(r, RMAX, Er, Ez, Hphi)
    dr, dz = float(r[1] - r[0]), float(z[1] - z[0])
    zero = np.zeros_like(Er)
    # E uses cos(ωt+φ); Bφ (the H column) uses sin(ωt+φ) — supplied at runtime.
    write_thetamode_series(out_file, 0.0, Z_STRUCT, dr, dz, [
        ("E", (("r", Er), ("t", zero), ("z", Ez)), E_UNIT),
        ("B", (("r", zero), ("t", Hphi), ("z", zero)), B_UNIT),
    ])
    return r, z, Ez[0]


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    r, z, ez1 = _build_rf(RF1_GDF, "EzRe", "ErRe", "HphiIm", RF1_FILE)
    _, _, ez2 = _build_rf(RF2_GDF, "EzIm", "ErIm", "HphiRe", RF2_FILE)
    nr, nz = r.size, z.size
    L = float(z[-1] - z[0])
    # |Ez| = |EzRe + i EzIm|; its z-integral is the 1-kW synchronous voltage.
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

    print(f"\nWrote openPMD linac fields → {RF1_FILE}, {RF2_FILE}")


if __name__ == "__main__":
    main()
