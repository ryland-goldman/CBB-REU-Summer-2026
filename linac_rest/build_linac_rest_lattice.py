"""
Per-section physics table + √P gradient scaling, and the chained 7-section Impact-T deck
builder (``build_impact``) for Cornell Linac Sections 2–8 (the ``linac_rest`` stage).
Sections 2–8 have no field maps; the vendored S-band TW shape (rfdata4–7) is reused and all
per-section physics lives in the calibrated field scale.

See linac_rest/README.md for the field model, energy budget, calibration, and gotchas.
"""

import copy
import math
import os

from impact import Impact
from impact.fieldmaps import read_fieldmap_rfdata

# ── RF / operating point (shared across all 7 sections) ──────────────────────
RF_FREQ_HZ = 2856e6          # S-band drive frequency [Hz]
POWER_MW = 11.0             # RF input power [MW] per section; √P-scaled below
TABLE_POWER_MW = 15.0       # power the details.md ΔE column is quoted at (the √P reference)
PHASE_DEG = 0.0             # on-crest synchronous phase

CELL_LENGTH_M = 0.035       # S-band cell length d = 3.5 cm
BETA0_D = 2.0 * math.pi * RF_FREQ_HZ * CELL_LENGTH_M / 299792458.0  # β₀ d = ω d / c
SIN_BETA0_D = math.sin(BETA0_D)   # ≈ 0.8657; body scale = entrance / this

IN_TO_M = 0.0254            # inch → metre


# ── Per-section table (details.md) ───────────────────────────────────────────
# Fields: name (CEA/CU label); length_m (TW structure length L); de15_mev (ΔE @ 15 MW);
# bore_in ((entrance, exit) bore *diameter* taper [inches] → radii via section_bore_radii);
# quad_in (drift-quad length after the section [inches]); quad_label (source quad Q2…Q8).
SECTIONS = (
    {"name": "CEA 2", "length_m": 2.94, "de15_mev": 33.0,
     "bore_in": (0.99, 0.78), "quad_in": 11.0, "quad_label": "Q2"},
    {"name": "CEA 3", "length_m": 2.94, "de15_mev": 33.0,
     "bore_in": (0.99, 0.78), "quad_in": 18.0, "quad_label": "Q3"},
    {"name": "CU 5", "length_m": 4.97, "de15_mev": 51.0,
     "bore_in": (1.16, 0.92), "quad_in": 25.0, "quad_label": "Q4"},
    {"name": "CEA 4", "length_m": 5.15, "de15_mev": 55.0,
     "bore_in": (1.16, 0.92), "quad_in": 16.2, "quad_label": "Q5"},
    {"name": "CEA 5", "length_m": 5.15, "de15_mev": 55.0,
     "bore_in": (1.16, 0.92), "quad_in": 22.0, "quad_label": "Q6"},
    {"name": "CU 3", "length_m": 4.97, "de15_mev": 51.0,
     "bore_in": (1.16, 0.92), "quad_in": 22.0, "quad_label": "Q7"},
    {"name": "CU 4", "length_m": 4.97, "de15_mev": 51.0,
     "bore_in": (1.16, 0.92), "quad_in": 20.9, "quad_label": "Q8"},
)
N_SECTIONS = len(SECTIONS)   # 7 (sections 2–8)

DRIFT_M = 0.4                # placeholder inter-section drift [m] (girder gaps not in details.md)


# ── √P scaling helpers ───────────────────────────────────────────────────────
def power_scale(power_mw=None):
    """√P field-amplitude scale factor relative to the 15 MW table: sqrt(P_op / 15)."""
    p = POWER_MW if power_mw is None else power_mw
    return math.sqrt(p / TABLE_POWER_MW)


def section_de_target(index, power_mw=None):
    """Per-section energy-gain target ΔE_target [MeV] at the operating power.

    ΔE_target(P_op) = ΔE_table × sqrt(P_op / 15); the value the per-section field scale is
    calibrated to (the gain is NOT computed analytically).
    """
    return SECTIONS[index]["de15_mev"] * power_scale(power_mw)


def section_gradient(index=None, power_mw=None):
    """Average accelerating gradient [MV/m] at the operating power.

    With ``index`` given, returns G_i(P_op) = (ΔE_table,i / L_i)·√(P/15); with ``index=None``,
    the full 7-section tuple (lab order).
    """
    if index is None:
        return tuple(section_gradient(i, power_mw) for i in range(N_SECTIONS))
    sec = SECTIONS[index]
    g15 = sec["de15_mev"] / sec["length_m"]      # @15 MW reference gradient [MV/m]
    return g15 * power_scale(power_mw)


def section_bore_radii(index):
    """(entrance, exit) bore *radius* [m] for a section (details.md diameters → radii)."""
    d_in, d_out = SECTIONS[index]["bore_in"]
    return (d_in * IN_TO_M / 2.0, d_out * IN_TO_M / 2.0)


def section_quad_length_m(index):
    """Real tabulated quad length [m] after a section (details.md, inches → m)."""
    return SECTIONS[index]["quad_in"] * IN_TO_M


def _quad_transfer_2x2(k1, length):
    """2×2 single-plane transfer matrix (x, x') of a thick quadrupole [SI], as a flat tuple
    ``(m11, m12, m21, m22)``. ``k1`` is the geometric focusing strength K1 [1/m²] in THIS plane
    (k1>0 focusing → trig matrix; k1<0 defocusing → hyperbolic; k1=0 → drift)."""
    if k1 > 0.0:
        s = math.sqrt(k1)
        c_, sn = math.cos(s * length), math.sin(s * length)
        return (c_, sn / s, -s * sn, c_)
    if k1 < 0.0:
        s = math.sqrt(-k1)
        ch, sh = math.cosh(s * length), math.sinh(s * length)
        return (ch, sh / s, s * sh, ch)
    return (1.0, length, 0.0, 1.0)


def _mat2_mul(a, b):
    """Multiply two flat 2×2 matrices (a·b), each ``(m11, m12, m21, m22)``."""
    return (a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3],
            a[2] * b[0] + a[3] * b[2], a[2] * b[1] + a[3] * b[3])


def _doublet_cell_half_trace(k1, l_q, l_drift):
    """½·Tr of one FODO-doublet cell's transfer matrix (cos μ = ½·Tr; |½·Tr| < 1 ⇒ stable).

    Cell: half-gap drift → +K1 half-quad (L_q/2) → −K1 half-quad (L_q/2) → half-gap drift → the
    following RF section as a field-free drift (``l_drift`` = L_section(i+1)).
    """
    half_q = l_q / 2.0
    half_gap = DRIFT_M / 2.0
    m = _quad_transfer_2x2(0.0, half_gap)                       # half inter-section gap
    m = _mat2_mul(_quad_transfer_2x2(k1, half_q), m)            # +K1 lead half-quad
    m = _mat2_mul(_quad_transfer_2x2(-k1, half_q), m)           # −K1 trailing half-quad
    m = _mat2_mul(_quad_transfer_2x2(0.0, half_gap), m)         # half inter-section gap
    m = _mat2_mul(_quad_transfer_2x2(0.0, l_drift), m)          # RF section as a drift
    return 0.5 * (m[0] + m[3])


def _solve_doublet_k1(mu_deg, l_q, l_drift, k1_max):
    """Geometric K1 [1/m²] of a ± doublet giving per-cell phase advance ``mu_deg`` (bisection on
    ``cos μ = ½·Tr``; ``½·Tr`` is monotone-decreasing from +1 at K1=0). If μ is unreachable within
    ``k1_max``, returns ``k1_max`` (the strongest stable focusing the cell supports)."""
    target = math.cos(math.radians(mu_deg))
    lo, hi = 1e-4, k1_max
    if _doublet_cell_half_trace(hi, l_q, l_drift) > target:
        return hi                      # cell can't reach this μ within k1_max → use the ceiling
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _doublet_cell_half_trace(mid, l_q, l_drift) > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def fodo_quad_gradients(*, phase_adv_deg=50.0, k1_max=14.0, mc2_mev=0.510998950,
                        ke_in_mev=25.0):
    """Energy-scaled FODO doublet base gradients [T/m] (PLACEHOLDER optics — guessed K1, A→T
    undocumented; exploratory, see README -> "Space charge & quads").

    Returns a length-``N_SECTIONS`` (7) list of signed lead-pole base gradients ``g_i`` [T/m],
    one per placed gap (``quad2``…``quad7``); the 7th entry (Q8, after section 8) is never placed
    and is fixed at ``0.0`` so the shape matches the ``[0.0]*N_SECTIONS`` default. The lead-pole
    sign alternation gap-to-gap (``(-1)**i``) is baked into ``g_i`` here.

    Recipe per gap (exact thick-lens cell matrix — the thin-lens formula is NOT used; the
    half-quad phase √K1·(L_q/2) ≈ 0.46 rad is not small):
        K1_i : solve cos μ = ½·Tr(cell_i)   [1/m²]  (both planes; symmetric ± doublet)
        Bρ_i = √(KE·(KE + 2·mc²)) / c       [T·m]   (per-section EXIT energy)
        g_i  = (-1)**i · K1_i · Bρ_i        [T/m]

    ``ke_in_mev`` defaults to ≈25 MeV for standalone calls; the sim passes the measured sec-1
    handoff ⟨KE⟩ so the Bρ energy-scaling tracks the actual beam.
    """
    c = 299792458.0                       # speed of light [m/s]
    mc2_ev = mc2_mev * 1e6                 # electron rest energy [eV]

    grads = []
    ke_ev = ke_in_mev * 1e6               # running cumulative KE at section EXITs [eV]
    for i in range(N_SECTIONS - 1):       # 6 placed gaps (Q2…Q7); Q8 appended as 0.0 below
        ke_ev += section_de_target(i) * 1e6             # per-section EXIT energy (table sum)
        # Bρ = p/e: with energies in eV, p[eV/c] = √(KE·(KE+2mc²)); /c (m/s) → SI T·m.
        b_rho = math.sqrt(ke_ev * (ke_ev + 2.0 * mc2_ev)) / c
        l_q = section_quad_length_m(i)                  # this gap's real (full) quad length [m]
        l_drift = SECTIONS[i + 1]["length_m"]           # following RF section, treated as a drift
        k1 = _solve_doublet_k1(phase_adv_deg, l_q, l_drift, k1_max)   # [1/m²], both planes
        grads.append(((-1) ** i) * k1 * b_rho)          # signed lead-pole base gradient g_i [T/m]
    grads.append(0.0)                     # Q8 (after the last section) is never placed
    return grads


def total_rf_length_m():
    """Σ active TW structure length over sections 2–8 [m]."""
    return sum(s["length_m"] for s in SECTIONS)


def total_lattice_length_m(n_drifts=None):
    """Σ (section length) + inter-section spacing [m].

    Each inter-section spacing is a ``DRIFT_M`` field-free margin PLUS a real-length quadrupole
    (gap/2 drift, quad, gap/2 drift — see ``build_impact``); one after every section except the
    last. With ``n_drifts`` given, only that many spacings (and quad lengths) are counted.
    """
    n = (N_SECTIONS - 1) if n_drifts is None else n_drifts
    quads = sum(section_quad_length_m(i) for i in range(n))
    return total_rf_length_m() + n * DRIFT_M + quads


def expected_exit_ke_mev(ke_in_mev, power_mw=None):
    """Validation gate: exit ⟨KE⟩ = measured ⟨KE⟩_in + Σ ΔE_target,i(P_op) [MeV].

    ``ke_in_mev`` is the MEASURED mean KE of the read-in sec-1 exit dump (NOT hardcoded 25).
    """
    return ke_in_mev + sum(section_de_target(i, power_mw) for i in range(N_SECTIONS))


# ═════════════════════════════════════════════════════════════════════════════
# LATTICE ASSEMBLY — reuses the vendored rfdata4–7 TW field shape (see README ->
# "Field model"). Per-section 4-line solrf superposition (SLAC-PUB-2295):
#   entrance rfdata4 θ₀=base+0°  scale=S
#   body_1   rfdata5 θ₀=base+30° scale=S/sin(β₀d)
#   body_2   rfdata6 θ₀=base+90° scale=S/sin(β₀d)
#   exit     rfdata7 θ₀=base+0°  scale=S
# The rfdata Fourier reconstruction uses the period stored INSIDE the file (~0.105 m 3-cell
# block), NOT the element `L` — `L` only sets the integrated z-range, so a longer section is
# more cells of the same field. PLACEHOLDER_SCALE S is calibrated per section by the sim.
# ═════════════════════════════════════════════════════════════════════════════

RFDATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "rfdata")
RFDATA_FILES = ("rfdata4", "rfdata5", "rfdata6", "rfdata7")   # entrance, body_1, body_2, exit
FILE_ID = {"entrance": 4, "body_1": 5, "body_2": 6, "exit": 7}

# Template coupler-cell lengths (rfdata4/rfdata7 z-extent); body length per section is (L − these).
L_ENTRANCE = 0.052464        # template solrf_entrance L [m] (rfdata4)
L_EXIT = 0.05246             # template solrf_exit L [m] (rfdata7)

# 4-line TW phase offsets [deg] relative to the section base phase (SLAC two-SW decomposition).
LINE_PHASE_OFFSET = {"entrance": 0.0, "body_1": 30.0, "body_2": 90.0, "exit": 0.0}

# Placeholder per-section field scale [V/m] before calibration; sim overwrites it.
PLACEHOLDER_SCALE = 11.5e6

# Impact-T deck header knobs (overridable by the sim via config()). SC OFF by default (Bcurr=0).
DECK_NP = 4000               # default macroparticle count (sim overrides via Np)
DECK_NXYZ = 16               # SC mesh per axis (unused at Bcurr=0; active under SPACE_CHARGE; power of 2)
DECK_DT = 2.0e-12            # time step [s]
DECK_NTSTEP = 80000          # step cap (sized for ~33 m at Dt=2e-12; sim asserts mean_z reached)
# Transverse computational-domain half-width [m] — the box wall; kept at a physical beam-pipe
# scale (NOT widened to fake transmission). See README -> "Space charge & quads".
XYRAD_M = 0.02               # 20 mm domain half-width

# solrf `radius` bore aperture: the REAL tapered section bore. Default True — a physically-
# anchored no-focusing lower bound. See README -> "Space charge & quads".
BORE_APERTURE_ON = True


def _section_subelements(index, zedge, scale, base_phase_deg, name_prefix, bore_aperture_on):
    """Return the 4 `solrf` sub-element dicts for one TW section, placed at `zedge`.

    `scale` is the entrance/exit field scale S; the body lines get S / sin(β₀d). The
    inter-line phase pattern (+0/+30/+90/+0) is added to `base_phase_deg`. The entrance/exit
    coupler cells keep the template short length; the body carries (L − L_entrance − L_exit).

    `bore_aperture_on` gates the solrf `radius`: the real tapered bore when True, else 0 ⇒ no
    scrape.
    """
    L = SECTIONS[index]["length_m"]
    r_in = section_bore_radii(index)[0] if bore_aperture_on else 0.0
    L_body = L - L_ENTRANCE - L_EXIT
    if L_body <= 0:
        raise ValueError(f"section {index} length {L} m too short for the coupler cells")
    geom = (
        ("entrance", zedge,                       L_ENTRANCE, scale),
        ("body_1",   zedge + L_ENTRANCE,          L_body,     scale / SIN_BETA0_D),
        ("body_2",   zedge + L_ENTRANCE,          L_body,     scale / SIN_BETA0_D),
        ("exit",     zedge + L_ENTRANCE + L_body, L_EXIT,     scale),
    )
    eles = []
    for line, ze, length, sc in geom:
        eles.append({
            "type": "solrf",
            "name": f"{name_prefix}_{line}",
            "L": length,
            "zedge": ze,
            "rf_field_scale": sc,
            "rf_frequency": RF_FREQ_HZ,
            "theta0_deg": base_phase_deg + LINE_PHASE_OFFSET[line],
            "filename": f"rfdata{FILE_ID[line]}",
            "radius": r_in,
            "solenoid_field_scale": 0.0,
        })
    return eles


def section_group_names(index, name_prefix=None):
    """The 4 sub-element names of a section (for the Task-5 scale ControlGroup)."""
    prefix = name_prefix or f"sec{index + 2}"      # sections are labelled 2..8
    return [f"{prefix}_{line}" for line in ("entrance", "body_1", "body_2", "exit")]


def _load_vendored_fieldmaps():
    """Read the vendored rfdata4–7 into the lume-impact fieldmap dict, keyed by `rfdataN`."""
    fieldmaps = {}
    for fname in RFDATA_FILES:
        path = os.path.join(RFDATA_DIR, fname)
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"vendored field shape missing: {path} — rfdata4–7 must be committed in "
                f"{RFDATA_DIR} (see linac_rest/README.md).")
        fieldmaps[fname] = read_fieldmap_rfdata(path)
    return fieldmaps


def build_impact(power_mw=None, phase_deg=None, drift_m=None, np_particles=None,
                 dt=None, ntstep=None, nxyz=None, scales=None, quads_on=False,
                 quad_k=None, bcurr=None, verbose=False):
    """Assemble the chained 7-section Impact-T deck and return a configured `Impact`.

    The sections are placed at increasing `zedge` with a `drift` of `drift_m` after every
    section except the last, and a `quadrupole` (real tabulated length, K1 default 0) inside
    each inter-section drift. The field scale per section is `scales[i]` if given, else
    PLACEHOLDER_SCALE. All sections reuse the vendored rfdata4–7 shapes.

    No `write_beam` slice dumps — they break the repeated `track_to_s` calls the calibration
    makes (fort.10N stat-column conflict in `load_many_fort`); per-section vs-z evolution comes
    from `I.stat(...)`. The final beam is `I.particles["final_particles"]` for the handoff OUT.

    Returns an `Impact` object with `.input['lattice']`, `.input['fieldmaps']`, `.header`
    populated and `.configure()` called (no `initial_particles` yet — the sim sets those).
    """
    # power_mw is accepted for API symmetry; the deck encodes power only via `scales`.
    base_phase = PHASE_DEG if phase_deg is None else phase_deg
    gap = DRIFT_M if drift_m is None else drift_m
    # `is None` (NOT `quad_k or …`): an explicit all-zero override must not be discarded by `or`.
    if quad_k is None:
        quad_k = fodo_quad_gradients() if quads_on else [0.0] * N_SECTIONS
    assert len(quad_k) >= N_SECTIONS - 1, (
        f"quad_k needs ≥{N_SECTIONS - 1} entries (one per placed quad, Q2..Q7); "
        f"got {len(quad_k)}")
    bore_aperture_on = bool(BORE_APERTURE_ON or quads_on)

    I = Impact(verbose=verbose)
    I.input["fieldmaps"] = _load_vendored_fieldmaps()

    lattice = []
    z = 0.0
    for i in range(N_SECTIONS):
        prefix = f"sec{i + 2}"                       # sec2 .. sec8
        scale = (scales[i] if scales is not None else PLACEHOLDER_SCALE)
        lattice += _section_subelements(i, z, scale, base_phase, prefix, bore_aperture_on)
        z += SECTIONS[i]["length_m"]
        if i < N_SECTIONS - 1:
            # Inter-section spacing: gap/2 drift, real-length quadrupole, gap/2 drift. The quad
            # length is NOT subtracted from `gap` (several real quads exceed the 0.4 m margin).
            qL = section_quad_length_m(i)
            half = gap / 2.0
            # Quad/drift bore radius gated on `quads_on` (NOT the already-True `bore_aperture_on`)
            # so the quads-OFF headline stays byte-identical (radius 0.0, no extra scrape plane).
            # Uses the section EXIT radius [1] (quad sits downstream of the exit taper; the solrf
            # body uses the ENTRANCE radius [0] — the real tapered bore). See README.
            r_exit = section_bore_radii(i)[1] if quads_on else 0.0
            lattice.append({"type": "drift", "name": f"drift{i + 2}a",
                            "L": half, "zedge": z, "radius": r_exit})
            z += half
            if quads_on:
                # H/V doublet: split the tabulated quad into two opposite-sign qL/2 halves (lead
                # sign quad_k[i], trailing its negation) so the ± pair net-focuses BOTH planes;
                # halves sum to qL so downstream zedges are unchanged. See README.
                g_lead = quad_k[i]
                qhalf = qL / 2.0
                lattice.append({
                    "type": "quadrupole", "name": f"quad{i + 2}a", "L": qhalf, "zedge": z,
                    "b1_gradient": g_lead,
                    "file_id": 0,                   # 0 ⇒ hard-edge (no Enge fringe field)
                    "radius": r_exit})
                z += qhalf
                lattice.append({
                    "type": "quadrupole", "name": f"quad{i + 2}b", "L": qhalf, "zedge": z,
                    "b1_gradient": -g_lead,
                    "file_id": 0,
                    "radius": r_exit})
                z += qhalf
            else:
                # Quads-OFF: single zero-K1 quad (optically a drift of its length); do NOT split,
                # so the element list / transmission is byte-identical.
                lattice.append({
                    "type": "quadrupole", "name": f"quad{i + 2}", "L": qL, "zedge": z,
                    "b1_gradient": 0.0,
                    "file_id": 0,                   # 0 ⇒ hard-edge (no Enge fringe field)
                    "radius": r_exit})
                z += qL
            lattice.append({"type": "drift", "name": f"drift{i + 2}b",
                            "L": half, "zedge": z, "radius": r_exit})
            z += half
    total_len = z

    I.input["lattice"] = lattice
    I.ele = {e["name"]: e for e in lattice}

    h = I.header
    h["Npcol"], h["Nprow"] = 1, 1
    h["Bcurr"] = 0.0 if bcurr is None else bcurr     # 0 ⇒ space charge OFF; >0 ⇒ SC current [A]
    h["Flagimg"] = 0                                 # no image charge (no cathode)
    h["Dt"] = DECK_DT if dt is None else dt
    h["Ntstep"] = DECK_NTSTEP if ntstep is None else ntstep
    h["Np"] = DECK_NP if np_particles is None else np_particles
    n = DECK_NXYZ if nxyz is None else nxyz
    h["Nx"], h["Ny"], h["Nz"] = n, n, n
    # Transverse domain bound (NOT a physical pipe; the bore aperture is the solrf `radius`).
    h["Xrad"], h["Yrad"] = XYRAD_M, XYRAD_M
    h["Perdlen"] = total_len + 1.0                   # > total lattice length
    h["Bkenergy"] = 25.0e6                           # ref energy [eV] (~sec-1 exit; sim resets)
    h["Bfreq"] = RF_FREQ_HZ
    h["Bmass"] = 0.51099895e6
    h["Bcharge"] = -1.0

    I.configure()
    return I, total_len


def main():
    """Stage `build.main()` contract — no artifact to write (deck is in-memory); prints a
    section-table + total-length sanity banner. The sim calls `build_impact()` directly."""
    _, total_len = build_impact()
    print(f"linac_rest lattice: {N_SECTIONS} TW sections (2–8), "
          f"Σ RF {total_rf_length_m():.2f} m, total {total_len:.2f} m "
          f"(rfdata4–7 reused; SC off; quads K1=0).", flush=True)


# ── Self-check: helpers reproduce the details.md @11 MW table within rounding ──
if __name__ == "__main__":
    print(f"sin(beta0 d) = {SIN_BETA0_D:.4f}  (expect ~0.8657)")
    print(f"power_scale(11) = {power_scale(11.0):.4f}  (sqrt(11/15) ~ 0.8563)")
    print(f"\n{'Sec':<7}{'L[m]':>7}{'dE@15':>8}{'G@15':>8}{'dE@11':>8}{'G@11':>8}"
          f"{'bore_r[mm]':>14}")
    # Expected @11 MW column from plan §4 (for the rounding check).
    de11_expected = (28.3, 28.3, 43.7, 47.1, 47.1, 43.7, 43.7)
    g11_expected = (9.6, 9.6, 8.8, 9.1, 9.1, 8.8, 8.8)
    for i, sec in enumerate(SECTIONS):
        g15 = sec["de15_mev"] / sec["length_m"]
        de11 = section_de_target(i, 11.0)
        g11 = section_gradient(i, 11.0)
        r_in, r_out = section_bore_radii(i)
        print(f"{sec['name']:<7}{sec['length_m']:>7.2f}{sec['de15_mev']:>8.0f}{g15:>8.2f}"
              f"{de11:>8.1f}{g11:>8.2f}{r_in*1e3:>7.2f}->{r_out*1e3:.2f}")
        # G@15 ≈ ΔE_table / L (definitional) and √P scaling reproduces the @11 MW column.
        assert abs(g15 - sec["de15_mev"] / sec["length_m"]) < 1e-9
        assert abs(de11 - de11_expected[i]) < 0.05, (sec["name"], de11, de11_expected[i])
        assert abs(g11 - g11_expected[i]) < 0.05, (sec["name"], g11, g11_expected[i])
    print(f"\nN_SECTIONS = {N_SECTIONS}")
    print(f"Sigma RF length      = {total_rf_length_m():.2f} m")
    print(f"Sigma lattice length = {total_lattice_length_m():.2f} m "
          f"(+{N_SECTIONS - 1} x [DRIFT_M={DRIFT_M} m + real quad length])")
    print(f"Sigma dE_target @11  = {sum(section_de_target(i, 11.0) for i in range(N_SECTIONS)):.1f}"
          f" MeV (details.md: 329 @15 -> 282 @11)")
    print(f"exit KE from 25 MeV  = {expected_exit_ke_mev(25.0, 11.0):.1f} MeV (~307 expected)")
    print("\nT1 self-check passed.")
