"""
Build the chained Impact-T lattice for Cornell Linac **Sections 2–8** (the ``linac_rest``
stage): seven electron travelling-wave (TW) accelerating sections, run as ONE Impact-T
deck downstream of ``linac_sec1``.

Sections 2–8 have **no field maps** (none exist; LinacSim/BMAD model them with the generic
constant-gradient linac function — ``reference/Linac Simulation Documentation/details.md``
lines 151–153). We reuse the shipped S-band TW template field shape (``rfdata4–7`` from
lume-impact's ``traveling_wave_cavity`` template) verbatim as the *shape only* and put **all**
per-section physics in the calibrated field scale. The rfdata carries no R/τ/shunt-impedance
information — those are already embedded in the per-section energy-gain table below, so encoding
them in the field profile would double-count.

This module owns the **per-section physics constants + the √P gradient-scaling helper** (Task 1).
The lattice-assembly code (rfdata reuse, the 4-line ``solrf`` TW superposition, drifts, quads,
deck header) is appended **below** the ``SECTION-TABLE / SCALING (Task 1)`` block by the
lattice-build tasks (T2/T3).

Per-section energy budget (ONE power convention, plan §4):

    G_section(15 MW)   = ΔE_table / L                       # @15 MW reference gradient
    ΔE_target(P_op)    = ΔE_table × sqrt(P_op / 15)         # √P scaling to the operating power
    section_gradient(P_op) = G_section(15) × sqrt(P_op / 15)

The default operating power is ``POWER_MW = 11`` (the section-1 faithful klystron point), one
convention for the whole linac. The @15 MW column in ``details.md`` corresponds to a *different*
(15 MW) section-1 run with a different, higher captured input — it is NOT a co-equal target for
this beam, so the faithfulness gate is the measured ⟨KE⟩_in (from the read-in sec-1 exit dump)
plus Σ ΔE_target,i(P_op), computed from actuals (see linac_rest/README.md and plan §4).

Frequency 2856 MHz (S-band), on-crest (θ₀ = 0), space charge off by default (the sim's
`SPACE_CHARGE=True` opts into an exploratory SC run via `bcurr`), quads OFF (K1 = 0) for the
headline beam (the A→T quad calibrations are undocumented — details.md line 178).
"""

import copy
import math
import os

from impact import Impact
from impact.fieldmaps import read_fieldmap_rfdata

# ── RF / operating point (shared across all 7 sections) ──────────────────────
RF_FREQ_HZ = 2856e6          # S-band drive frequency [Hz] (details.md; matches linac_sec1)
POWER_MW = 11.0             # RF input power [MW] per section — ONE convention for the whole
                           # linac (the section-1 faithful klystron point). √P-scaled below.
TABLE_POWER_MW = 15.0       # power the details.md ΔE column is quoted at (the √P reference)
PHASE_DEG = 0.0             # on-crest synchronous phase; β > 0.999 ⇒ no per-section re-phasing

# Travelling-wave cell geometry, reused verbatim from the SLAC/lume-impact template.
# The body-cell field scale is entrance / sin(β₀ d) and the inter-line phases are fixed by the
# two-standing-wave superposition (entrance, body+30°, body+90°, exit) — see T2/T3 below.
CELL_LENGTH_M = 0.035       # S-band cell length d = 3.5 cm
BETA0_D = 2.0 * math.pi * RF_FREQ_HZ * CELL_LENGTH_M / 299792458.0  # β₀ d = ω d / c
SIN_BETA0_D = math.sin(BETA0_D)   # ≈ 0.8657 (S-band, d = 3.5 cm); body scale = entrance/this

IN_TO_M = 0.0254            # inch → metre


# ── Per-section table (plan §4 / details.md lines 159–167) ───────────────────
# Sections 2–8 of the Cornell linac (the straight electron line to CHESS). The e+ compressor
# (CU 2) is OUT OF SCOPE — its lattice role is not established in details.md, so no topology is
# asserted here (deferred as future work).
#
# Fields per section:
#   name      : human label (CEA/CU design designation from details.md)
#   length_m  : active TW structure length L [m]
#   de15_mev  : energy gain @ 15 MW [MeV] (details.md "Energy gain @ 15 MW" column)
#   bore_in   : (entrance, exit) bore *diameter* taper [inches] (details.md "Bore (in)" column)
#   quad_in   : drift-quad length after the section [inches] (details.md "Used in simulation")
#   quad_label: which quad (Q2…Q8) the length is sourced from (details.md quad table)
#
# Bore diameters are converted to *radii* in metres by ``section_bore_radii`` below; the radius
# is what Impact-T's solrf ``radius`` aperture and any tapered-bore scrape would use.
#
# Caveats propagated from details.md: sec 5/6 cell counts are guesses (marked * in the source)
# and their bore taper is "taken from Dan Fromowitz's code, not independently verified" (marked
# †). The quad A→T (current→field) calibrations are unknown for EVERY quad (details.md line 178),
# so quads are OFF (K1 = 0) for the headline beam — the lengths here only place real-length
# drift-equivalent gaps; the K1 values live with the exploratory FODO knob in the sim module.
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

# Placeholder inter-section drift [m] (girder gaps are not in details.md — flagged placeholder).
DRIFT_M = 0.4


# ── √P scaling helpers (Task 1) ──────────────────────────────────────────────
def power_scale(power_mw=None):
    """√P field-amplitude scale factor relative to the 15 MW table: sqrt(P_op / 15).

    Energy gain in a constant-gradient TW section scales as the gradient, which scales as
    sqrt(input power). At the default 11 MW this is ≈ 0.856.
    """
    p = POWER_MW if power_mw is None else power_mw
    return math.sqrt(p / TABLE_POWER_MW)


def section_de_target(index, power_mw=None):
    """Per-section energy-gain target ΔE_target [MeV] at the operating power.

    ΔE_target(P_op) = ΔE_table × sqrt(P_op / 15). This is the value the per-section field scale
    is calibrated to (Task 5 via ``autophase_and_scale``); the gain is NOT computed analytically.
    """
    return SECTIONS[index]["de15_mev"] * power_scale(power_mw)


def section_gradient(index=None, power_mw=None):
    """Average accelerating gradient [MV/m] at the operating power.

    With ``index`` given, returns that section's gradient G_i(P_op) = (ΔE_table,i / L_i)·√(P/15).
    With ``index=None``, returns the full tuple for all 7 sections (lab order, sections 2–8).
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
    """½·Tr of one FODO-doublet cell's transfer matrix in the lead-focusing plane.

    The cell is: half-gap drift → +K1 half-quad (L_q/2) → −K1 half-quad (L_q/2) → half-gap drift
    → the following RF SECTION treated as a field-free drift (``l_drift`` = L_section(i+1)). The
    per-cell phase advance is ``cos μ = ½·Tr(M_cell)``; ``|½·Tr| < 1`` ⇒ stable (μ real). The
    OTHER (defocusing-lead) plane has the same |½·Tr| for a symmetric ± doublet — that plane
    equality is exactly why the doublet net-focuses BOTH planes (unlike a single thick quad)."""
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
    ``cos μ = ½·Tr``). ``½·Tr`` is monotone-decreasing from +1 at K1=0 down through the target as
    K1 rises, so a simple bracket [≈0, k1_max] converges. If the target μ is unreachable within
    ``k1_max`` (a too-weak short cell — e.g. gap 2 at high μ), returns ``k1_max`` (the strongest
    stable focusing the cell supports), which the caller energy-scales as usual."""
    target = math.cos(math.radians(mu_deg))
    lo, hi = 1e-4, k1_max
    if _doublet_cell_half_trace(hi, l_q, l_drift) > target:
        return hi                      # cell can't reach this μ within k1_max → use the ceiling
    for _ in range(200):               # ~2e-60 bracket; pure-math, no scipy dependency
        mid = 0.5 * (lo + hi)
        if _doublet_cell_half_trace(mid, l_q, l_drift) > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def fodo_quad_gradients(*, phase_adv_deg=50.0, k1_max=14.0, mc2_mev=0.510998950,
                        ke_in_mev=25.0):
    """Energy-scaled FODO **doublet** base gradients [T/m] (PLACEHOLDER optics).

    Returns a length-``N_SECTIONS`` (7) list of signed LEAD-POLE base gradients ``g_i`` [T/m],
    one per **placed gap** (``quad2``…``quad7``, i.e. one after every section except the last);
    the 7th entry (Q8, after the final section 8) is **never placed** and is fixed at ``0.0`` so
    the shape matches the ``[0.0]*N_SECTIONS`` default and no downstream length asymmetry /
    IndexError can occur (plan R7). ``build_impact`` places each gap as an **H/V doublet**: a
    ``+g_i`` half-quad (L=L_q/2) immediately followed by a ``−g_i`` half-quad (L=L_q/2). The
    sign alternation gap-to-gap (``(-1)**i``) is baked into ``g_i`` here (the lead-pole sign).

    These are **NOT** measured quad strengths — the A→T (current→field) calibration is
    undocumented for every quad (``details.md`` line 178). They are derived from accelerator
    optics: a constant per-cell phase-advance FODO of ± doublets, energy-scaled by the local
    beam momentum (Bρ). Caveats stated on every QUADS_ON output: the inter-doublet "drift" is
    really a multi-metre accelerating RF section, so the realized μ is **nominal** (the lattice
    is non-periodic and the beam accelerates through it) — boundedness of σ_x AND σ_y is the
    acceptance, not a measured μ/cell.

    **Why a doublet, not a single thick quad.** A single-sign thick quad per gap focuses one
    plane and DEFOCUSES the other over the 1.6–2.8 m half-cell; the defocused plane diverges and
    scrapes the real bore, so transmission came out WORSE than no-focusing (49–58 % vs the 78.5 %
    baseline) at every μ (verified). An H/V doublet (± pair) net-focuses BOTH planes — for a
    symmetric ± doublet the cell's ``|½·Tr|`` is identical in the two planes, so both get the
    same phase advance. That plane symmetry is the whole point of going to the doublet.

    **Recipe (exact thick-lens cell matrix — the thin-lens (K1·l)²·d doublet formula is NOT used,
    it is unreliable here because the half-quad phase √K1·(L_q/2) ≈ 0.46 rad (~27°) is not small):**

        cell_i   = drift(gap/2) · (+K1 half-quad, L_q/2) · (−K1 half-quad, L_q/2) · drift(gap/2)
                   · drift(L_section(i+1))                  (RF section as a field-free drift)
        K1_i     : solve  cos μ = ½·Tr(cell_i)             [1/m²]  (per-gap, bisection; both planes)
        Bρ_i     = √(KE·(KE + 2·mc²)) / c                  [T·m]   (per-section EXIT energy)
        g_i      = (-1)**i · K1_i · Bρ_i                   [T/m]   (energy-scaling enters via Bρ)

    **μ defaults to 50°.** The uniform design μ is bounded by the WEAKEST cell (gap 2 — the short
    2.94 m CEA-2 section — tops out near 67°), so 50° is reachable by every gap and sits mid-band
    (0–180° stable). A μ sweep (45/50/55°, doublet) showed the result is robust: exit ⟨KE⟩ stays
    ≈308–310 MeV (gate 2 PASS) and transmission clusters at 77.6–78.2 % with a bounded oscillating
    σ_x/σ_y (max RMS ≈4.4 mm, well inside the 9.9 mm exit bore) — 50° gives the best transmission.
    The ~78 % ceiling is set by the injected sec-1 beam expanding before the first quad (after
    section 2) plus the added quad-location apertures, NOT by the focusing strength. ``k1_max``
    caps the bisection bracket; a cell that cannot reach the target μ within it falls back to
    ``k1_max`` (strongest stable focusing it supports). μ and ``k1_max`` are documented design
    knobs, **not** reverse-fit to a transmission number.

    ``ke_in_mev`` defaults to the ≈25 MeV nominal for standalone calls; the sim passes the
    **measured** sec-1 handoff ⟨KE⟩ so the Bρ energy-scaling tracks the actual beam.
    """
    c = 299792458.0                       # speed of light [m/s]
    mc2_ev = mc2_mev * 1e6                 # electron rest energy [eV]

    grads = []
    ke_ev = ke_in_mev * 1e6               # running cumulative KE at section EXITs [eV]
    for i in range(N_SECTIONS - 1):       # 6 placed gaps (Q2…Q7); Q8 appended as 0.0 below
        # Per-section EXIT energy: KE_in + Σ_{j≤i} ΔE_target,j (deterministic from the table).
        ke_ev += section_de_target(i) * 1e6
        # Magnetic rigidity Bρ = p/e. With energies in eV, p[eV/c] = √(KE·(KE+2mc²)); dividing by
        # c (m/s) converts eV/c → SI momentum/charge in T·m (the eV cancels e, leaving J·s/(C·m)).
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

    Each inter-section spacing is a ``DRIFT_M`` field-free margin PLUS a real-length
    quadrupole (gap/2 drift, quad, gap/2 drift — see ``build_impact``); there is one such
    spacing after every section except the last (6 for 7 sections, matching the plan §4
    "Drift after" column where section 8 has none). With ``n_drifts`` given, only that many
    spacings are counted (and the first ``n_drifts`` quad lengths).
    """
    n = (N_SECTIONS - 1) if n_drifts is None else n_drifts
    quads = sum(section_quad_length_m(i) for i in range(n))
    return total_rf_length_m() + n * DRIFT_M + quads


def expected_exit_ke_mev(ke_in_mev, power_mw=None):
    """Validation gate: exit ⟨KE⟩ = measured ⟨KE⟩_in + Σ ΔE_target,i(P_op) [MeV].

    ``ke_in_mev`` is the MEASURED mean KE of the read-in sec-1 exit dump (NOT hardcoded 25).
    At 11 MW from a ~25 MeV input this PREDICTED (table-sum) value is ≈ 307 MeV; the achieved
    calibrated run reaches ≈ 308 MeV (307.97 survivors through the real bore / 309.2 full-beam).
    """
    return ke_in_mev + sum(section_de_target(i, power_mw) for i in range(N_SECTIONS))


# ═════════════════════════════════════════════════════════════════════════════
# LATTICE ASSEMBLY (Task 2 — rfdata reuse; Task 3 — chained deck)
# ═════════════════════════════════════════════════════════════════════════════
#
# Sections 2–8 have no field maps, so we reuse the shipped S-band TW template field
# *shape* (rfdata4–7, vendored into linac_rest/rfdata/) verbatim. Each section is the
# template's 4-line `solrf` superposition of two standing-wave maps (entrance coupler
# cell + two body lines 90° apart + exit coupler cell — G. A. Loew et al., SLAC-PUB-2295):
#
#   entrance : rfdata4, short coupler cell, θ₀ = base+0°,  scale = S
#   body_1   : rfdata5, the bulk length,    θ₀ = base+30°, scale = S / sin(β₀d)
#   body_2   : rfdata6, the bulk length,    θ₀ = base+90°, scale = S / sin(β₀d)
#   exit     : rfdata7, short coupler cell, θ₀ = base+0°,  scale = S
#
# The rfdata Fourier reconstruction uses the PERIOD stored INSIDE the file (a single
# ~0.105 m 3-cell block) as its wavelength — NOT the lattice element's `L`. The element
# `L` only sets the active z-range [zedge, zedge+L] that Impact-T integrates the periodic
# field over, so a longer section is simply more cells of the same per-cell field. That is
# why "rescale length" is just setting the body element `L` per section — the field shape
# is reused unchanged (verified: a 2.94 m section reaches the section length and accelerates
# on-crest at θ₀ = 0). The entrance/exit coupler cells keep their own short template length;
# the body lines carry the remaining (L − L_entrance − L_exit).
#
# The per-section field scale S is a placeholder here (PLACEHOLDER_SCALE); the sim module
# calibrates it per section via autophase_and_scale to hit ΔE_target (Task 5). The relative
# body/coupler ratio (1/sin β₀d) and the +0/+30/+90/+0 inter-line phase pattern are preserved
# across the calibration (one ControlGroup scales all four lines together).

RFDATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "rfdata")
RFDATA_FILES = ("rfdata4", "rfdata5", "rfdata6", "rfdata7")   # entrance, body_1, body_2, exit
FILE_ID = {"entrance": 4, "body_1": 5, "body_2": 6, "exit": 7}

# Template coupler-cell lengths (rfdata4/rfdata7 z-extent) — the short entrance/exit cells.
# Reused verbatim; the body length per section is (L − these two).
L_ENTRANCE = 0.052464        # template solrf_entrance L [m] (rfdata4)
L_EXIT = 0.05246             # template solrf_exit L [m] (rfdata7)

# 4-line TW phase offsets [deg] relative to the section base phase (the SLAC two-SW
# decomposition: entrance 0, body +30, body +90, exit 0). Reused verbatim.
LINE_PHASE_OFFSET = {"entrance": 0.0, "body_1": 30.0, "body_2": 90.0, "exit": 0.0}

# Placeholder per-section field scale [V/m] before Task-5 calibration. Picked near the
# template's 25.5e6 entrance scale so an uncalibrated deck is sane (and so a single-section
# smoke run accelerates), but the headline run overwrites it via autophase_and_scale.
PLACEHOLDER_SCALE = 11.5e6

# Impact-T deck header knobs (overridable by the sim via config()). SC is OFF by default (Bcurr=0,
# Npcol=Nprow=1, Flagimg=0); the SC mesh (Nx/Ny/Nz) is unused at Bcurr=0 but must be a valid power
# of 2 (and is the active SC grid when the sim passes bcurr>0 via SPACE_CHARGE=True).
DECK_NP = 4000               # default macroparticle count (sim overrides via Np)
DECK_NXYZ = 16               # SC mesh per axis (unused at Bcurr=0; active under SPACE_CHARGE; power of 2)
DECK_DT = 2.0e-12            # time step [s] (template value; sim may override)
DECK_NTSTEP = 80000          # step cap (sized for ~33 m at Dt=2e-12; sim asserts mean_z reached)
# Transverse computational-domain half-width Xrad=Yrad [m] — the Impact-T grid bound; particles
# outside it are lost at the box wall. Kept at a PHYSICALLY MEANINGFUL scale (~beam-pipe radius),
# NOT widened to fake high transmission: counting particles at tens of cm as "transmitted" would
# be meaningless (no real pipe is that wide). The quads-OFF headline reports the HONEST count-
# based transmission through this box as a no-focusing artifact (the real machine's FODO quads
# keep the beam in, so the true transmission is higher) — the headline DELIVERABLE is the energy
# (≈308 MeV) and per-section calibration, not the transmission. Physics-neutral at the headline:
# SC off (Bcurr=0) ⇒ the box has no field-solve effect, pure containment. lume-impact's 0.015 m default is the
# bore scale; 0.02 m here is a hair above the widest section bore radius (14.7 mm).
XYRAD_M = 0.02               # 20 mm domain half-width (beam-pipe scale; honest no-focusing loss)

# Bore aperture: ON for the headline — the REAL tapered section bore (section_bore_radii,
# 12.6→9.9 / 14.7→11.7 mm) is the solrf `radius` aperture. Transmission against the real bore is
# a physically-anchored, un-tunable number: "with the real aperture and NO transverse focusing,
# X% survives." That's an honest PESSIMISTIC LOWER BOUND — the only missing ingredient (the real
# machine's FODO quads) can only INCREASE it. This is cleaner than radius=0 (which falsely reports
# ~100% — the beam isn't really contained) and cleaner than an arbitrary wide box (whose width
# isn't a physical aperture). The headline DELIVERABLE is the energy (≈308 MeV) + per-section
# calibration; transmission is reported as this no-focusing lower bound. Set False only for an
# energy-only study where the transverse divergence shouldn't scrape (then the box XYRAD_M binds).
BORE_APERTURE_ON = True


def _section_subelements(index, zedge, scale, base_phase_deg, name_prefix, bore_aperture_on):
    """Return the 4 `solrf` sub-element dicts for one TW section, placed at `zedge`.

    `scale` is the entrance/exit field scale S; the body lines get S / sin(β₀d). The
    inter-line phase pattern (+0/+30/+90/+0) is added to `base_phase_deg`. The entrance/exit
    coupler cells keep the template short length; the body carries (L − L_entrance − L_exit).

    `bore_aperture_on` gates the solrf `radius`: the real tapered bore when True (the
    exploratory QUADS_ON FODO case, which has the focusing to keep the beam off the bore), else
    0 ⇒ no scrape (the quads-OFF headline — the upstream iris already collimates well inside the
    bore, so a bore scrape with zero focusing over 36 m is a no-optics artifact).
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
    """Read the vendored rfdata4–7 into the lume-impact fieldmap dict layout.

    Keyed by `rfdataN` (the solrf `filename`/`file_id` Impact-T expects). lume-impact
    writes these into the run workdir on configure(); we never depend on ~/Downloads.
    """
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
    each inter-section drift. The field scale per section is `scales[i]` if given (the
    calibrated values from Task 5), else PLACEHOLDER_SCALE — the sim calibrates before the
    headline run. All sections reuse the vendored rfdata4–7 shapes.

    No `write_beam` slice dumps: the per-section vs-z evolution comes from Impact-T's continuous
    `I.stat(...)` arrays (energy/σ/emittance vs z). `write_beam` elements were tried and dropped
    — they break the repeated `track_to_s` calls the calibration makes (the fort.10N stat
    columns conflict in `load_many_fort`). The final beam is `I.particles["final_particles"]`
    for the openPMD handoff OUT; a single exit dump is sufficient for plot_chain (sorts by ⟨z⟩).

    Returns an `Impact` object with `.input['lattice']`, `.input['fieldmaps']`, `.header`
    populated and `.configure()` called (no `initial_particles` yet — the sim sets those).
    """
    # power_mw is accepted for API symmetry but the deck encodes power only via the per-section
    # `scales` (calibrated by the sim); the raw value isn't used in lattice assembly.
    base_phase = PHASE_DEG if phase_deg is None else phase_deg
    gap = DRIFT_M if drift_m is None else drift_m
    # Quad gradients [T/m]. Use `is None` (NOT `quad_k or …`): an explicit all-zero override
    # must NOT be silently discarded by `or`. If unset and quads are on, fall back to the derived
    # energy-scaled FODO gradients; otherwise zeros (quads-OFF headline ⇒ optically a drift).
    if quad_k is None:
        quad_k = fodo_quad_gradients() if quads_on else [0.0] * N_SECTIONS
    assert len(quad_k) >= N_SECTIONS - 1, (
        f"quad_k needs ≥{N_SECTIONS - 1} entries (one per placed quad, Q2..Q7); "
        f"got {len(quad_k)}")
    # Bore aperture follows the quads: ON for the exploratory QUADS_ON FODO study (real focusing
    # keeps the beam off the bore, so the scrape is meaningful), OFF for the quads-OFF headline
    # (a bore scrape with zero focusing over 36 m is a no-optics artifact). BORE_APERTURE_ON
    # forces it on for a standalone bore-loss study even with quads off.
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
            # Inter-section spacing = a field-free drift margin (`gap`) split around a
            # `quadrupole` at its REAL tabulated length: gap/2 drift, quad, gap/2 drift.
            # The quad length is NOT subtracted from `gap` (several real quads, e.g. Q4 =
            # 0.635 m, exceed the 0.4 m placeholder margin — subtracting would give a
            # negative, non-monotonic drift). K1 defaults to 0 (quads OFF for the headline —
            # A→T calibration is unknown), so a zero-strength quad is optically a drift of
            # its own length; the total inter-section optical length is gap + quad_L.
            qL = section_quad_length_m(i)
            half = gap / 2.0
            # NEW quad/drift bore aperture, gated on `quads_on` (NOT `bore_aperture_on`): the
            # quad sits downstream of the section EXIT taper, so it uses the section EXIT radius
            # `[1]` (the solrf body already takes the ENTRANCE radius `[0]` — entrance-on-body /
            # exit-on-quad is the real tapered bore, not an inconsistency). Gating on `quads_on`
            # (not the already-True `bore_aperture_on`) keeps the quads-OFF headline byte-identical:
            # `radius` stays 0.0 there, so no new scrape plane is added and the published 78.5%
            # is unchanged. (Impact-T ignores drift `radius`, so the drift change is forward-looking;
            # the quad `radius>0` is the load-bearing real loss aperture on the focused path.)
            r_exit = section_bore_radii(i)[1] if quads_on else 0.0
            lattice.append({"type": "drift", "name": f"drift{i + 2}a",
                            "L": half, "zedge": z, "radius": r_exit})
            z += half
            if quads_on:
                # H/V DOUBLET: a single thick quad of one sign focuses one plane but DEFOCUSES the
                # other over the multi-metre half-cell (→ that plane scrapes, transmission worse than
                # no-focus). Split the tabulated quad into two opposite-sign halves (qL/2 each, back-
                # to-back): the +g/−g pair nets focusing in BOTH planes (1/f_net ∝ (K1·l)²·d). The
                # lead-pole sign is `quad_k[i]` (already alternating gap-to-gap via the helper's
                # (-1)**i), and the trailing half is its negation. Halves sum to qL, so the lattice
                # length and downstream zedges are unchanged vs. the single-quad placement.
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
                # Quads-OFF headline: keep the SINGLE zero-strength quad (byte-identical lattice —
                # do NOT split, so the exact element list / transmission is preserved). A zero-K1
                # quad is optically a drift of its own length.
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
    # Transverse domain bound (NOT a physical pipe). Wide enough to contain the unfocused
    # quads-OFF beam so particles aren't lost at the box wall (the no-optics artifact); the
    # physical bore aperture is the solrf `radius` (BORE_APERTURE_ON), separate from this.
    h["Xrad"], h["Yrad"] = XYRAD_M, XYRAD_M
    h["Perdlen"] = total_len + 1.0                   # > total lattice length
    h["Bkenergy"] = 25.0e6                           # ref energy [eV] (~sec-1 exit; sim resets)
    h["Bfreq"] = RF_FREQ_HZ
    h["Bmass"] = 0.51099895e6
    h["Bcharge"] = -1.0
    # (Flagimg already 0 above — no image charge; the sim sets a non-cathode coasting beam
    #  via initial_particles, so no Flagdist/Nemission cathode-emission settings are needed.)

    I.configure()
    return I, total_len


def main():
    """No build artifact to write — the deck is assembled in-memory by `build_impact`.

    Kept for the Stage `build.main()` contract (mirrors `build_*_field.py`); the sim calls
    `build_impact()` directly. Prints the section table + total length as a sanity banner.
    """
    _, total_len = build_impact()
    print(f"linac_rest lattice: {N_SECTIONS} TW sections (2–8), "
          f"Σ RF {total_rf_length_m():.2f} m, total {total_len:.2f} m "
          f"(rfdata4–7 reused; SC off; quads K1=0).", flush=True)


# ── Self-check: helper reproduces the plan §4 table within rounding (Task 1 acceptance) ──
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
