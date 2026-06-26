"""
Cornell Linac sections 5-8 (Impact-T). main(): handoff IN (the positron core from the e+/e-
converter, which sits after the WarpX section 4) -> build the chained 4-section traveling-wave
Impact-T deck from the vendored rfdata4-7 field shapes -> apply the FROZEN per-section field scale
+ crest phase -> I.run() (space charge OFF; quad doublets + sec5/6 cavity-solenoids ON) -> openPMD
handoff OUT + injection_summary.json.

Four S-band (2856 MHz) TW sections (CEA 4/5 + CU 3/4) chained into ONE Impact-T deck and
integrated as one time-ordered beam. Sections 5-8 have no field maps; the vendored S-band TW shape
(rfdata4-7) is reused verbatim and all per-section physics lives in the per-section field scale.

CALIBRATION IS FROZEN. The old stage ran a per-section brentq scale-fit + crest-phase scan + §5
validation gates each run; here the field_scale + crest_phase_deg are hardcoded in
config/linac5-8.yaml (derived once) and applied directly. See docs/linac5-8.md.

Run as `python sim/linac5-8.py` (hyphenated name is not importable). sim/plot/linac5-8.py makes
the figures. Do NOT call build/run from import -- main() does everything.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# Must precede `import numpy`: OpenMP latches OMP_NUM_THREADS at load, so prepare_env()'s later
# set is ignored and the grid oversubscribes.
os.environ.setdefault("OMP_NUM_THREADS", os.environ.get("OMP_THREADS", "1"))

import json
import math
import shutil

import numpy as np
import yaml

from sim.helpers.tools import MC2_EV, C_LIGHT, M_E, E_CHARGE as Q_E, prepare_env
from sim.helpers.loadparticles import (
    read_warpx_dump, write_openpmd_particles, upstream_exit_lab_z, upsample_smeared)
from sim.helpers.metrics import beam_quality
from sim.helpers.tqdmwrapper import impact_progress

CONFIG = "config/linac5-8.yaml"
MC2_MEV = MC2_EV / 1e6                  # electron rest energy [MeV]

# 4-line TW phase offsets [deg] relative to the section base phase (SLAC-PUB-2295 two-SW
# decomposition): entrance +0, body_1 +30, body_2 +90, exit +0.
LINE_PHASE_OFFSET = {"entrance": 0.0, "body_1": 30.0, "body_2": 90.0, "exit": 0.0}
FILE_ID = {"entrance": 4, "body_1": 5, "body_2": 6, "exit": 7}   # rfdata file per line
RFDATA_FILES = ("rfdata4", "rfdata5", "rfdata6", "rfdata7")
IN_TO_M = 0.0254                        # inch -> metre

# Capture-optics calibration (Fromowitz, capture_optics_specs.md):
GRAD_PER_K1 = 0.534                     # Table 6.3: k1 = 1 m^-2 <-> 0.534 T/m gradient (160 MeV ref)
SOL_COIL_RADIUS_M = 0.1304             # cavity-solenoid mean coil radius (inner 8.93 + outer 17.15 cm)/2
SOL_END_MARGIN_M = 0.1                  # coil inset from each section end. Small on purpose: the coil
                                        # must span nearly the full cavity so the capture solenoid grabs
                                        # the divergent beam at the sec5 ENTRANCE and across the gaps -- a
                                        # large inset leaves the low-energy core unfocused there and it
                                        # scrapes the bore (a 0.6 m inset cut transmission ~18x). At 0.1 m
                                        # the on-axis Bz is ~20% of peak and EQUAL at both ends, so the
                                        # periodic Fourier extension stays continuous (slope kink only).
SOL_N_FOURIER = 40                      # Fourier terms for the on-axis Bz expansion
SOL_FILE_ID = 50                        # cavity-solenoid fieldmap file IDs: SOL_FILE_ID + section index


# ── Config helpers ───────────────────────────────────────────────────────────────
def load_config(path=CONFIG):
    """Load the YAML config, coercing the scientific-notation numeric fields to float.

    PyYAML's safe loader only recognises sci-notation as a float when the exponent carries a sign
    (`1.0e-12` ok; `2856.0e6`, `2.019127e7` parse as STRINGS). Rather than litter the config with
    `e+` signs, coerce the known-numeric fields here so the config stays readable.
    """
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    cfg["rf"]["rf_freq_hz"] = float(cfg["rf"]["rf_freq_hz"])
    for sec in cfg["sections"]:
        sec["field_scale"] = float(sec["field_scale"])
        sec["crest_phase_deg"] = float(sec["crest_phase_deg"])
        if "crest_offset_deg" in sec:                # optimizer knob (added on top of the crest)
            sec["crest_offset_deg"] = float(sec["crest_offset_deg"])
    return cfg


def _beta0_d(cfg):
    """beta0*d = omega*d/c for the S-band cell (sets the body-line scale 1/sin(beta0 d))."""
    return 2.0 * math.pi * cfg["rf"]["rf_freq_hz"] * cfg["rf"]["cell_length_m"] / C_LIGHT


def section_bore_radii(sec):
    """(entrance, exit) bore RADIUS [m] for a section (config bore_in diameters [in] -> radii)."""
    d_in, d_out = sec["bore_in"]
    return (d_in * IN_TO_M / 2.0, d_out * IN_TO_M / 2.0)


def section_quad_length_m(sec):
    """Real tabulated quad length [m] after a section (config quad_in inches -> m)."""
    return sec["quad_in"] * IN_TO_M


def section_quad_gradient(sec):
    """Pole-tip gradient [T/m] for a section's trailing quad: machine k1 [m^-2] x 0.534 T/m
    (Table 6.3, 160 MeV reference). Placed as a +-grad doublet in build_impact (net focusing)."""
    return float(sec.get("quad_k1", 0.0)) * GRAD_PER_K1


def _solenoid_fieldmap(L_sec, file_id):
    """Build a solrf rfdata fieldmap (Ez block = 0, Bz block = normalised finite-solenoid on-axis
    profile, peak 1) for a cavity-solenoid spanning one section. Returns (filename, fieldmap_dict).

    The on-axis Bz of a thin-shell solenoid (mean radius SOL_COIL_RADIUS_M, inset SOL_END_MARGIN_M
    from each section end -- small, so the coil spans nearly the full cavity and focuses the divergent
    beam at the section entrance and across the inter-section gaps; the residual end field is equal at
    both boundaries so the periodic Fourier extension stays continuous) is decomposed into Impact-T
    Fourier coefficients. The element's `solenoid_field_scale` then sets the
    PEAK Bz [T]; Impact-T type-105 expands it paraxially as a STATIC field (Bz = B0 - B''0 r^2/4)."""
    from beamphysics.interfaces.impact import create_fourier_coefficients
    a = SOL_COIL_RADIUS_M
    margin = min(SOL_END_MARGIN_M, 0.4 * L_sec)
    b = (L_sec - 2.0 * margin) / 2.0                 # coil half-length
    zc = L_sec / 2.0
    zg = np.linspace(0.0, L_sec, 1201)
    zr = zg - zc
    bz = 0.5 * ((zr + b) / np.sqrt((zr + b) ** 2 + a ** 2)
                - (zr - b) / np.sqrt((zr - b) ** 2 + a ** 2))
    bz = bz / bz.max()                               # peak 1 -> solenoid_field_scale carries [T]
    bz_coefs = create_fourier_coefficients(zg, bz, n=SOL_N_FOURIER)
    ez_block = np.array([1.0, 0.0, L_sec, L_sec, 0.0])           # n=1, z0, z1, L, single 0 coef -> Ez=0
    bz_block = np.concatenate(([float(len(bz_coefs)), 0.0, L_sec, L_sec], bz_coefs))
    data = np.concatenate((ez_block, bz_block))
    return f"rfdata{file_id}", {"info": {"format": "rfdata"}, "data": data}


def section_de_target(sec, cfg):
    """Per-section energy-gain target ΔE_target [MeV] = ΔE_table * sqrt(P_op / table_power).

    Recorded in the frozen-calibration table for the section_gains plot; not used to drive the
    deck (the field scale is frozen, not fit to this).
    """
    return sec["de15_mev"] * math.sqrt(cfg["rf"]["power_mw"] / cfg["rfdata"]["table_power_mw"])


def section_ele_names(cfg, index):
    """(entrance, body_1, body_2, exit) solrf element names for section `index` (0-based)."""
    prefix = f"sec{index + cfg['lattice']['first_section']}"   # sections labelled 5..8
    return tuple(f"{prefix}_{line}" for line in ("entrance", "body_1", "body_2", "exit"))


def section_group_name(cfg, index):
    return f"sec{index + cfg['lattice']['first_section']}_scale"


# ── Frozen-calibration APPLY helpers (ported from the old calibration.py apply path; the
# search/validation is dropped). The rf_field_scale ControlGroup is absolute=True with factors
# [1, 1/sin(beta0 d), 1/sin(beta0 d), 1], so its value S sets entrance/exit=S, body=S/sin(beta0 d),
# preserving the template body ratio. theta0_deg is ABSOLUTE per solrf sub-element. ──────────────
def _ensure_section_group(I, cfg, index):
    """Create (idempotently) the rf_field_scale ControlGroup over a section's 4 solrf cells.

    Returns the group name. Must be called BEFORE _set_group_scale (an absolute group defaulting 0
    would otherwise zero the field).
    """
    gname = section_group_name(cfg, index)
    if gname in getattr(I, "group", {}):
        return gname
    inv_sin = 1.0 / math.sin(_beta0_d(cfg))
    I.add_group(
        gname,
        ele_names=list(section_ele_names(cfg, index)),
        var_name="rf_field_scale",
        factors=[1.0, inv_sin, inv_sin, 1.0],
        absolute=True,
    )
    return gname


def _set_group_scale(I, gname, value):
    """Set a section's ControlGroup field-scale value. Impact.__setitem__ needs the
    "name:attribute" form (a bare I[gname]=value splits on ':' and raises)."""
    I[f"{gname}:rf_field_scale"] = float(value)


def _set_section_phase(I, cfg, index, phase_deg):
    """Pin every solrf sub-element of a section to its on-crest driven base phase, with the
    template's fixed inter-line offsets (entrance +0, body_1 +30, body_2 +90, exit +0) added on."""
    entrance, body1, body2, exit_ = section_ele_names(cfg, index)
    I[entrance]["theta0_deg"] = phase_deg + LINE_PHASE_OFFSET["entrance"]
    I[body1]["theta0_deg"] = phase_deg + LINE_PHASE_OFFSET["body_1"]
    I[body2]["theta0_deg"] = phase_deg + LINE_PHASE_OFFSET["body_2"]
    I[exit_]["theta0_deg"] = phase_deg + LINE_PHASE_OFFSET["exit"]


# ── Deck assembly (ported from build_linac_rest_lattice; SC-free, quad doublets + sec5/6 sol.) ────
def _section_subelements(cfg, index, zedge, scale, base_phase_deg, name_prefix, bore_on):
    """The 4 `solrf` sub-element dicts for one TW section, placed at `zedge`.

    `scale` is the entrance/exit field scale S; the body lines get S/sin(beta0 d). The inter-line
    phase pattern (+0/+30/+90/+0) is added to `base_phase_deg`. The entrance/exit coupler cells keep
    the template short length; the body carries (L - l_entrance - l_exit). `bore_on` gates the solrf
    `radius` to the real ENTRANCE bore (else 0 => no scrape).
    """
    sec = cfg["sections"][index]
    L = sec["length_m"]
    l_entrance = cfg["rfdata"]["l_entrance"]
    l_exit = cfg["rfdata"]["l_exit"]
    inv_sin = 1.0 / math.sin(_beta0_d(cfg))
    r_in = section_bore_radii(sec)[0] if bore_on else 0.0
    L_body = L - l_entrance - l_exit
    if L_body <= 0:
        raise ValueError(f"section {index} length {L} m too short for the coupler cells")
    geom = (
        ("entrance", zedge,                       l_entrance, scale),
        ("body_1",   zedge + l_entrance,          L_body,     scale * inv_sin),
        ("body_2",   zedge + l_entrance,          L_body,     scale * inv_sin),
        ("exit",     zedge + l_entrance + L_body, l_exit,     scale),
    )
    eles = []
    for line, ze, length, sc in geom:
        eles.append({
            "type": "solrf",
            "name": f"{name_prefix}_{line}",
            "L": length,
            "zedge": ze,
            "rf_field_scale": sc,
            "rf_frequency": cfg["rf"]["rf_freq_hz"],
            "theta0_deg": base_phase_deg + LINE_PHASE_OFFSET[line],
            "filename": f"rfdata{FILE_ID[line]}",
            "radius": r_in,
            "solenoid_field_scale": 0.0,
        })
    return eles


def _load_vendored_fieldmaps(cfg):
    """Read the vendored rfdata4-7 into the lume-impact fieldmap dict, keyed by `rfdataN`."""
    from impact.fieldmaps import read_fieldmap_rfdata
    rfdata_dir = cfg["rfdata"]["dir"]
    fieldmaps = {}
    for fname in RFDATA_FILES:
        path = os.path.join(rfdata_dir, fname)
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"vendored field shape missing: {path} -- rfdata4-7 must be present in "
                f"{rfdata_dir} (see docs/linac5-8.md).")
        fieldmaps[fname] = read_fieldmap_rfdata(path)
    return fieldmaps


def build_impact(cfg, workdir=None):
    """Assemble the chained 4-section Impact-T deck (quad doublets + sec5/6 solenoids, SC OFF) from the vendored
    rfdata4-7 shapes and return (configured `Impact`, total_lattice_length_m, section_bounds),
    where section_bounds is the (z_entry, z_exit) [m] of each TW section in deck z (real geometry,
    used by the section_gains figure instead of an even split).

    Each section is placed at increasing `zedge`; a `drift`/zero-K1 `quadrupole`/`drift` spacing
    follows every section except the last (the quad at its real tabulated length is optically a
    drift). The per-section field scale starts at the frozen `field_scale`; the run also applies it
    via the rf_field_scale ControlGroup so the body ratio is exact (see main()).

    `workdir` (with use_temp_dir=False) runs Impact-T IN PLACE there, so fort.18 lands at a known
    <workdir>/fort.18 the progress bar can poll (the default temp dir is non-deterministic). No
    write_beam slice dumps -- per-section vs-z evolution comes from I.stat(...). The final beam is
    I.particles["final_particles"] for the handoff OUT.
    """
    from impact import Impact

    base_phase = cfg["rf"]["phase_deg"]
    gap = cfg["lattice"]["drift_m"]
    bore_on = bool(cfg["lattice"]["bore_aperture_on"])
    first = cfg["lattice"]["first_section"]
    sections = cfg["sections"]
    n_sec = len(sections)

    if workdir is not None:
        I = Impact(verbose=False, use_temp_dir=False, workdir=workdir)
    else:
        I = Impact(verbose=False)
    I.input["fieldmaps"] = _load_vendored_fieldmaps(cfg)

    lattice = []
    section_bounds = []                                  # (z_entry, z_exit) [m] per TW section
    z = 0.0
    for i, sec in enumerate(sections):
        prefix = f"sec{i + first}"                       # sec5 .. sec8
        z_entry = z
        lattice += _section_subelements(cfg, i, z, sec["field_scale"], base_phase, prefix, bore_on)
        # Cavity-solenoid (sec 5/6 only): a solenoid-only solrf (rf_field_scale 0, static Bz)
        # OVERLAPPING the cavity over [z_entry, z_entry+L] -- the continued transverse focusing of
        # the divergent positrons. solenoid_field_scale = peak Bz [T]; the Bz Fourier map is
        # normalised to peak 1. It does not advance z (overlaps the RF), so the absolute crest
        # geometry is unchanged.
        b_sol = float(sec.get("solenoid_b_tesla", 0.0))
        if b_sol > 0.0:
            fname, fmap = _solenoid_fieldmap(sec["length_m"], SOL_FILE_ID + i)
            I.input["fieldmaps"][fname] = fmap
            lattice.append({
                "type": "solrf", "name": f"{prefix}_solenoid",
                "L": sec["length_m"], "zedge": z_entry,
                "rf_field_scale": 0.0, "rf_frequency": cfg["rf"]["rf_freq_hz"],
                "theta0_deg": 0.0, "filename": fname,
                "radius": 0.0, "solenoid_field_scale": b_sol,
            })
        z += sec["length_m"]
        section_bounds.append((z_entry, z))
        if i < n_sec - 1:
            # Inter-section spacing: gap/2 drift, REAL-LENGTH quad DOUBLET, gap/2 drift -- matching
            # the deck geometry the crests were derived on. The doublet (+grad over qL/2 then -grad
            # over qL/2) is the thesis QH/QV per section (Fig 6.32 machine k1 -> gradient k1*0.534
            # T/m); sec5's k1~0 leaves it optically a drift (the solenoid focuses there). The total
            # gap+qL length MUST match: the crest phases are ABSOLUTE (Impact-T theta0, t=0 ref), so
            # the cumulative path length -- the bunch arrival time at sections 6-8 -- must hold or
            # those sections fall off-crest. radius 0.0 => no extra scrape plane.
            qL = section_quad_length_m(sec)
            grad = section_quad_gradient(sec)
            half = gap / 2.0
            qh = qL / 2.0
            lattice.append({"type": "drift", "name": f"drift{i + first}a",
                            "L": half, "zedge": z, "radius": 0.0})
            z += half
            lattice.append({"type": "quadrupole", "name": f"quad{i + first}a", "L": qh,
                            "zedge": z, "b1_gradient": +grad, "file_id": 0, "radius": 0.0})
            z += qh
            lattice.append({"type": "quadrupole", "name": f"quad{i + first}b", "L": qh,
                            "zedge": z, "b1_gradient": -grad, "file_id": 0, "radius": 0.0})
            z += qh
            lattice.append({"type": "drift", "name": f"drift{i + first}b",
                            "L": half, "zedge": z, "radius": 0.0})
            z += half
    total_len = z

    I.input["lattice"] = lattice
    I.ele = {e["name"]: e for e in lattice}

    h = I.header
    h["Npcol"], h["Nprow"] = 1, 1
    h["Bcurr"] = 0.0                                     # space charge OFF (overridden if SC on)
    h["Flagimg"] = 0                                     # no image charge (no cathode)
    h["Dt"] = cfg["deck"]["dt"]
    h["Ntstep"] = cfg["deck"]["ntstep"]
    h["Np"] = cfg["beam"]["np"]
    n = cfg["deck"]["nxyz"]
    h["Nx"], h["Ny"], h["Nz"] = n, n, n
    h["Xrad"], h["Yrad"] = cfg["deck"]["xyrad_m"], cfg["deck"]["xyrad_m"]
    h["Perdlen"] = total_len + 1.0                       # > total lattice length
    h["Bkenergy"] = 78.0e6                               # placeholder [eV]; lume-impact resets it from
                                                         # initial_particles -- NOT the theta0 phase ref
    h["Bfreq"] = cfg["rf"]["rf_freq_hz"]
    h["Bmass"] = MC2_EV                                  # same mass for e- and e+
    # Bcharge sign follows the beam species: +1 positrons, -1 electrons.
    h["Bcharge"] = 1.0 if str(cfg["beam"].get("species", "electrons")).startswith("positron") else -1.0

    I.configure()
    return I, total_len, section_bounds


# ── Handoff IN: the positron core from the converter (which sits after the WarpX section 4) ──────
def load_converter_core(cfg):
    """Read the converter positron beam, keep the captured core (KE >= MIN_KE_MEV), match the core
    to Np (downsample if richer; smeared upsample to Np if `beam.upsample` and the converter yield
    left fewer -- so the ~1% survivors are enough for statistics), drift to mean t + zero z for
    Impact-T injection. Returns (ParticleGroup, info dict). The ParticleGroup carries the
    captured-core charge (no renormalisation).
    """
    diag = cfg["io"]["conv_particles"]
    summary = cfg["io"]["conv_summary"]
    min_ke_mev = cfg["beam"]["min_ke_mev"]
    np_keep = cfg["beam"]["np"]
    rng_seed = cfg["beam"]["rng_seed"]

    species = cfg["beam"].get("species", "electrons")
    P = read_warpx_dump(diag, species=species)           # configurable species, t-coords, last dump
    n_all = P.n_particle
    z_local = float(P["mean_z"])                         # converter LOCAL frame
    z_inject_lab = upstream_exit_lab_z(summary, z_local)  # chain local->lab so the segment places right
    q_exit = float(P["charge"])                          # all converter-beam charge (denominator)

    ke_mev = (P.energy - MC2_MEV * 1e6) / 1e6
    core = ke_mev >= min_ke_mev
    if core.sum() < 50:
        raise RuntimeError(
            f"only {int(core.sum())} converter positrons above MIN_KE_MEV={min_ke_mev} MeV -- "
            f"capture cut too aggressive or converter yield too low")
    Pc = P[core]
    q_core = float(Pc.charge)
    n_core_raw = int(Pc.n_particle)                          # genuine converter-core macro count

    if Pc.n_particle > np_keep:
        rng = np.random.default_rng(rng_seed)
        sel = rng.choice(Pc.n_particle, np_keep, replace=False)
        Pc = Pc[sel]
        Pc.weight = Pc.weight * (q_core / float(Pc.charge))   # restore total core charge
    elif Pc.n_particle < np_keep and cfg["beam"].get("upsample", False):
        # The converter yield caps the core macro count, so the surviving handful is too few for
        # downstream statistics. SC is OFF, so each macro is dynamically independent: split the core
        # to Np with local phase-space smear (NOT plain bootstrap -- coincident duplicates track
        # identically and add nothing). Survival fraction + moments are preserved; only the SAMPLING
        # density rises. The genuine resolution is still set by n_core_raw (recorded for honesty).
        Pc = upsample_smeared(Pc, np_keep, rng_seed=rng_seed,
                              smear=float(cfg["beam"].get("upsample_smear", 0.2)))
        # Re-impose the KE floor: smearing the transverse momentum of a low-pz / large-angle
        # parent (kept only by its transverse p in the TOTAL-energy cut) can pull a clone below
        # MIN_KE_MEV. Drop those so the relativistic-core guarantee (beta_min, no-slip crest) holds.
        keep = (Pc.energy - MC2_MEV * 1e6) / 1e6 >= min_ke_mev
        if not keep.all():
            Pc = Pc[keep]
            Pc.weight = Pc.weight * (q_core / float(Pc.charge))   # restore total core charge

    # Impact-T injects at a common time with z == 0: drift to mean t, then translate z to 0.
    Pc.drift_to_t(Pc["mean_t"])
    Pc.z = Pc.z - Pc["mean_z"]

    ke_in = float(Pc["mean_energy"] / 1e6 - MC2_MEV)
    ke_min = float(Pc.energy.min() / 1e6 - MC2_MEV)
    g_min = 1.0 + ke_min / MC2_MEV
    beta_min = math.sqrt(max(0.0, 1.0 - 1.0 / (g_min * g_min)))
    info = dict(
        n_conv_in=int(n_all), n_core=int(Pc.n_particle), n_core_raw=n_core_raw,
        q_conv_in_C=q_exit, q_core_C=q_core,
        core_charge_frac=(q_core / q_exit if q_exit else 0.0),
        min_ke_mev_cut=float(min_ke_mev), ke_in_mev=ke_in,
        ke_min_core_mev=ke_min, beta_min_core=beta_min,
        z_inject_lab_m=z_inject_lab,
    )
    print(f"converter positron beam: {n_all} parts, {q_exit*1e12:.1f} pC; captured core "
          f"(KE>={min_ke_mev} MeV): {Pc.n_particle} parts, {q_core*1e12:.1f} pC "
          f"({info['core_charge_frac']*100:.1f}% of exit charge). <KE>_in {ke_in:.2f} MeV, "
          f"min-core KE {ke_min:.2f} MeV (beta_min={beta_min:.5f}), inject lab-z "
          f"{z_inject_lab:.3f} m", flush=True)
    return Pc, info


def _stat_vs_z(I, n=200, q_core_C=None, n_in=None):
    """Thin Impact-T's I.stat(...) z-arrays to ~`n` samples for the vs-z plots (write_beam dumps
    are off). sigma_KE uses sigma_gamma*mc2 (sigma_energy is not a stat key). With `q_core_C`/`n_in`
    the surviving macro count (fort.28 `n_particle`, aligned on mean_z) becomes a charge[pC] vs z
    column -- the aperture-loss curve."""
    zc = I.stat("mean_z")
    if len(zc) == 0:
        return {}
    idx = np.unique(np.linspace(0, len(zc) - 1, min(n, len(zc))).astype(int))
    out = {
        "z_m": zc[idx].tolist(),
        "ke_mev": (I.stat("mean_kinetic_energy")[idx] / 1e6).tolist(),
        "sigma_ke_mev": (I.stat("sigma_gamma")[idx] * MC2_MEV).tolist(),
        "sigma_x_m": I.stat("sigma_x")[idx].tolist(),
        "sigma_y_m": I.stat("sigma_y")[idx].tolist(),
        "norm_emit_x": I.stat("norm_emit_x")[idx].tolist(),
        "norm_emit_y": I.stat("norm_emit_y")[idx].tolist(),
    }
    npart = I.stat("n_particle")
    if len(npart) == len(zc):                            # fort.28 cadence matches fort.18
        out["n_particle"] = npart[idx].tolist()
        if q_core_C is not None and n_in:
            out["charge_pc"] = (npart[idx] / n_in * q_core_C * 1e12).tolist()
    return out


def _write_outputs(I, outdir, inj, species="electrons", charge=-Q_E):
    """Write the surviving ParticleGroups as WarpX-layout openPMD slices (sorted by <z>) plus
    injection_summary.json. Group charges were already re-imposed in main(). `species`/`charge`
    must match the beam (positron handoff => species="positrons", charge=+Q_E)."""
    part_dir = os.path.join(outdir, "particles")
    os.makedirs(part_dir, exist_ok=True)

    slices = []
    for _name, pg in I.particles.items():
        if _name == "initial_particles":                 # the injected beam (z~0), not an exit slice
            continue
        if pg is None or pg.n_particle < 50:
            continue
        slices.append((float(pg["mean_z"]), pg))
    slices.sort(key=lambda t: t[0])
    if not slices:                                       # e.g. uncaptured positrons fully scraped
        print("  no surviving particle group >=50 macroparticles -- writing summary only", flush=True)

    for it, (_zc, pg) in enumerate(slices):
        write_openpmd_particles(pg, part_dir, iteration=it, time=float(pg["mean_t"]),
                                species=species, charge=charge)

    with open(os.path.join(outdir, "injection_summary.json"), "w") as fh:
        json.dump(inj, fh, indent=2)


def main():
    prepare_env()
    cfg = load_config()
    outdir = cfg["io"]["outdir"]
    workdir = cfg["io"]["workdir"]

    # Fresh diags (regenerated, git-ignored): clear so a rerun doesn't mix old iterations.
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(workdir, exist_ok=True)

    # ── Handoff IN: positron core from the converter (after the WarpX section 4) ──────────
    P_in, core_info = load_converter_core(cfg)

    # ── Build the deck (quad doublets + sec5/6 solenoids, SC off) and apply per-section scale + crest ──
    # Run in-place under workdir so fort.18 lands at <workdir>/fort.18 for the progress poll.
    I, total_len, section_bounds = build_impact(cfg, workdir=workdir)
    n_sec = len(cfg["sections"])
    first = cfg["lattice"]["first_section"]
    power_mw = cfg["rf"]["power_mw"]

    # Apply the frozen calibration via the rf_field_scale ControlGroup (NOT just the build-time
    # element scales): the group is absolute=True defaulting 0, so adding it + configure() would
    # overwrite the baked-in scales with 0 (silent no-acceleration). Set the group scale AND the
    # absolute crest phase per section (theta0 is ABSOLUTE), then configure once.
    calib = []
    for i, sec in enumerate(cfg["sections"]):
        gname = _ensure_section_group(I, cfg, i)
        # crest_offset_deg is the optimizer's phase knob, added on top of the autophased crest
        # (autophase owns crest_phase_deg, so it cannot write that key). Absent => 0 (standalone runs).
        offset = float(sec.get("crest_offset_deg", 0.0))
        _set_section_phase(I, cfg, i, sec["crest_phase_deg"] + offset)
        _set_group_scale(I, gname, sec["field_scale"])
        calib.append({
            "index": i, "name": sec["name"],
            "scale": float(sec["field_scale"]),
            "crest_phase_deg": float(sec["crest_phase_deg"]),
            "crest_offset_deg": offset,
            "target_de_mev": float(section_de_target(sec, cfg)),
            "z_entry_m": float(section_bounds[i][0]),
            "z_exit_m": float(section_bounds[i][1]),
        })
    I.initial_particles = P_in
    I.configure()
    print(f"Deck: {n_sec} TW sections ({first}-{first + n_sec - 1}), Sigma {total_len:.2f} m, "
          f"P={power_mw:g} MW, frozen per-section scale+crest applied, SC off, quad doublets + "
          f"sec5/6 cavity-solenoids ON -> {outdir}/", flush=True)

    # ── Full run, bar driven from fort.18 (col 1 = reference z [m]) ────────────────────────
    # build_impact(workdir=) ran configure() with use_temp_dir=False, so I.path == workdir and
    # fort.18 lands at <workdir>/fort.18 (I.path is authoritative; poll it not the config value).
    print(f"Running Impact-T ({n_sec} sections, Ntstep={cfg['deck']['ntstep']})...", flush=True)
    fort18 = os.path.join(I.path, "fort.18")
    with impact_progress(fort18, total_len, desc="linac5-8"):
        I.run()
    if not I.finished or I.error:
        raise RuntimeError(f"Impact-T did not finish cleanly (finished={I.finished}, "
                           f"error={I.error})")

    # ── Transmission from MACRO COUNT, measured BEFORE re-imposing charge ──────────────────
    # n_out/n_in on the macro count (uniform per-macro weight) is the only honest transmission;
    # from charge AFTER the re-impose below it would force 1.0 and mask aperture loss.
    n_in = int(P_in.n_particle)
    q_core = float(P_in["charge"])
    P_out = I.particles["final_particles"] if "final_particles" in I.particles else None
    if P_out is None or P_out.n_particle == 0:
        # Whole bunch lost: the positron beam scraped on the bore before the deck end (the sec5/6
        # cavity-solenoids + quad doublets focus it, but a too-divergent core can still be lost).
        # Report 0% transmission.
        n_out, transmission, q_out, ke_out = 0, 0.0, 0.0, None
    else:
        n_out = int(P_out.n_particle)
        transmission = (n_out / n_in) if n_in else 0.0
        q_out = q_core * transmission                   # physically transmitted core charge
        # ── Re-impose physical charge for the openPMD `weighting` (SC-OFF loses it) ────────
        # Impact-T returns a default 1 C normalisation; rescale each group to q_core*(group n/n_in).
        for _name, _pg in I.particles.items():
            if _pg is not None and _pg.n_particle > 0:
                _pg.charge = q_core * (_pg.n_particle / n_in)
        ke_out = float(P_out["mean_energy"] / 1e6 - MC2_MEV)
    _mz = I.stat("mean_z")
    mean_z_reached = float(_mz[-1]) if len(_mz) else 0.0

    # ── Handoff OUT: openPMD + summary ────────────────────────────────────────────────────
    inj = dict(
        # HONEST capture denominator: the FULL converter positron-beam charge, NOT the
        # post-cut core -- so q_out/q_injected counts the dropped tail + in-run loss.
        q_injected_C=core_info["q_conv_in_C"],
        q_core_injected_C=core_info["q_core_C"],
        z_inject_lab_m=core_info["z_inject_lab_m"],
        z_inject_local_m=0.0,
        total_lattice_length_m=float(total_len),
        power_mw=float(power_mw), phase_deg=float(cfg["rf"]["phase_deg"]),
        quads_on=True,
        sec56_solenoid_b_tesla=[float(s.get("solenoid_b_tesla", 0.0)) for s in cfg["sections"]],
        quad_k1=[float(s.get("quad_k1", 0.0)) for s in cfg["sections"]],
        bore_aperture_on=bool(cfg["lattice"]["bore_aperture_on"]),
        xyrad_m=float(cfg["deck"]["xyrad_m"]),
        ke_in_mev=core_info["ke_in_mev"],
        ke_out_mev=ke_out,
        mean_z_reached_m=mean_z_reached,
        beta_min_core=core_info["beta_min_core"],
        # Transmission from MACRO COUNT (n_out/n_in), measured before re-imposing charge.
        n_core_in=n_in, n_out=n_out,
        transmission_core=transmission,
        q_out_C=q_out,
        core_charge_frac=core_info["core_charge_frac"],
        n_conv_in=core_info["n_conv_in"], n_core=core_info["n_core"],
        n_core_raw=core_info["n_core_raw"],
        min_ke_mev_cut=core_info["min_ke_mev_cut"],
        # The frozen per-section calibration table (scale, crest phase, ΔE target).
        calibration=calib,
        stat_vs_z=_stat_vs_z(I, q_core_C=q_core, n_in=n_in),
        # Exit beam-quality metrics (emittance, energy spread, spot); all-NaN on the no-survivor
        # case (P_out None) so the optimizer's read never KeyErrors. See sim/helpers/metrics.py.
        **beam_quality(P_out),
    )
    sp = str(cfg["beam"].get("species", "electrons"))
    _write_outputs(I, outdir, inj, species=sp,
                   charge=(Q_E if sp.startswith("positron") else -Q_E))
    ke_str = f"{ke_out:.1f} MeV" if ke_out is not None else "n/a (no survivors)"
    print(f"\nDone. Exit <KE> {ke_str} (in {core_info['ke_in_mev']:.1f} MeV); beam reached "
          f"{mean_z_reached:.2f}/{total_len:.2f} m; transmission {transmission*100:.1f}% "
          f"(quad doublets + sec5/6 solenoids). -> {outdir}/", flush=True)


if __name__ == "__main__":
    main()
