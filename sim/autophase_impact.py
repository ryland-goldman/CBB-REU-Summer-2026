"""
Auto-phase the linac 5-8 (Impact-T) RF sections for the POSITRON beam and rewrite the YAML.

`sim/main.py` runs this **before the linac5-8 stage**; it also works standalone. The per-section
crest_phase_deg are absolute Impact-T theta0 values that depend on both the species (positrons:
q=+e crest ~180 deg from electrons, lower injection energy shifts it) and the deck geometry, so they
must be re-found on the real deck with the positron core injected.

NUMERICAL model -- no Impact-T is launched (the analog of sim/autophase.py for the WarpX sections).
The Impact-T crest scan is already a pure 1D longitudinal problem (pencil-ised core, SC off, quads
K1=0), so a longitudinal RK4 integration through the exact on-axis Ez loses no physics and replaces
the old minutes-to-tens-of-minutes Impact-T scan with ~30 s of pure numpy.

Method: reconstruct each TW section's on-axis Ez from the vendored rfdata4-7 Fourier shapes as
Impact-T does (matches lume-impact's `impact.fieldmaps.ele_field` to float noise: the 4-line
+0/+30/+90/+0 superposition, the 1/sin(beta0 d) body scale, the absolute phase cos(2 pi f t +
theta0)). Place the sections at the real deck z (sections + real-length zero-K1 quad + drift spacing)
so the cumulative arrival time -- which the ABSOLUTE theta0 crest depends on -- matches the deck.
Pencil-ise the converter positron core (momentum onto +z, on-axis) and RK4-integrate the WHOLE bunch
through a base-phase scan; the crest maximises bunch-averaged exit KE (coarse -> fine -> parabolic, as
sim/autophase.find_crest). Earlier sections are pinned to their already-found crest so each section's
crest is found in their actual field, with the correct entering energy and time.

CREST vs ENERGY. The crests this finds match Impact-T's (re-derived ones agree with the previously
Impact-T-derived YAML values to ~1 deg). The reported per-section <KE> is a MODEL energy used only to
LOCATE the crest -- its magnitude is ~1.7x Impact-T's actual gain, a constant offset between
lume-impact's `ele_field` field reconstruction and the Fortran Impact-T field. A constant field-scale
offset does not move the argmax over phase, so the crest is unaffected; do NOT read the <KE> as the
stage energy (the real beam energy is what sim/linac5-8.py reports). --verify confirms the crest
directly with a tight Impact-T phase scan (robust to the amplitude offset), not the energy magnitude.

  python sim/autophase_impact.py            # phase sections 5 6 7 8, rewrite the YAML
  python sim/autophase_impact.py 5 6        # only sections 5-6
  python sim/autophase_impact.py --dry-run  # scan + report, write nothing
  python sim/autophase_impact.py --verify   # after scanning, run Impact-T at crest+-delta on the last
                                            # phased section and confirm its exit KE peaks at the crest
"""

import importlib.util
import math
import os
import re
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np

from sim.helpers.tools import C_LIGHT, E_CHARGE, M_E, MC2_EV, prepare_env

CONFIG = "config/linac5-8.yaml"
MC2_MEV = MC2_EV / 1e6

SUBSAMPLE = 1024             # macroparticles for the fine scan (bunch-averaged crest converges)
COARSE_SUB = 384             # lighter bunch for the wide scan (the peak location converges faster)
STEPS_PER_PERIOD = 40        # RF-cycle time resolution of the RK4 push (crest insensitive beyond this)
COARSE_STEP_DEG = 15.0       # wide 0..360 scan resolution (the crest is a broad cosine peak)
FINE_HALF_DEG = 10.0         # +-window about the coarse max
FINE_STEP_DEG = 1.0          # fine scan resolution


def _load_driver():
    """Import sim/linac5-8.py (hyphenated, not a normal module name) for its config loader, deck
    geometry helpers, and the positron-core handoff so the field shapes, section placement, and
    injected beam are byte-identical to what the driver runs."""
    path = os.path.join(os.path.dirname(__file__), "linac5-8.py")
    spec = importlib.util.spec_from_file_location("linac58_driver", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ── On-axis Ez reconstruction (matches impact.fieldmaps.ele_field to float noise) ────────────────
def _fourier_ez(z_local, coef, z0, zlen):
    """Impact-T Fourier on-axis field at element-local z (vectorised over z_local). Reproduces
    impact.fieldmaps.fourier_field_reconsruction: fz = c0/2 + Re[sum_n e^{i phi_n} (cos_n + i sin_n)],
    phi_n = -2 pi n ((z_local - z0)/zlen - 1/2)."""
    fcomplex = coef[1::2] + 1j * coef[2::2]
    n = np.arange(len(fcomplex)) + 1
    phi = -2.0 * np.pi * (np.outer(z_local, n) / zlen - (z0 / zlen) - 0.5)
    return coef[0] / 2.0 + np.real(np.exp(1j * phi) @ fcomplex)


def _build_section_fields(drv, cfg):
    """One entry per TW section: its (z_entry, z_exit) in deck z and the 4 precomputed field lines
    (entrance/body_1/body_2/exit). Each line carries a dense on-axis shape array (rf_field_scale *
    Fourier(z)) over [zedge, zedge+L] and its fixed inter-line phase offset; the time factor
    cos(omega t + base + offset) is applied per RK4 step. Section placement (sections + gap/2 +
    real-length quad + gap/2) mirrors linac5-8.build_impact so the absolute arrival time matches."""
    from impact.fieldmaps import read_fieldmap_rfdata, solrf_field_from_data

    fields = {f: solrf_field_from_data(
        read_fieldmap_rfdata(os.path.join(cfg["rfdata"]["dir"], f))["data"])
        for f in drv.RFDATA_FILES}

    l_entrance = cfg["rfdata"]["l_entrance"]
    l_exit = cfg["rfdata"]["l_exit"]
    inv_sin = 1.0 / math.sin(drv._beta0_d(cfg))
    gap = cfg["lattice"]["drift_m"]
    sections = cfg["sections"]
    n_sec = len(sections)

    def _line(zedge, length, scale, line):
        ez = fields[f"rfdata{drv.FILE_ID[line]}"]["Ez"]
        ng = max(400, int(length / ez["L"] * 60) + 1)        # ~60 samples per field period
        zg = np.linspace(zedge, zedge + length, ng)
        shape = scale * _fourier_ez(zg - zedge, ez["fourier_coefficients"], ez["z0"], ez["L"])
        return dict(zgrid=zg, shape=shape, offset_rad=math.radians(drv.LINE_PHASE_OFFSET[line]))

    out = []
    z = 0.0
    for i, sec in enumerate(sections):
        z_entry = z
        L, S = sec["length_m"], sec["field_scale"]
        L_body = L - l_entrance - l_exit
        lines = [
            _line(z_entry, l_entrance, S, "entrance"),
            _line(z_entry + l_entrance, L_body, S * inv_sin, "body_1"),
            _line(z_entry + l_entrance, L_body, S * inv_sin, "body_2"),
            _line(z_entry + l_entrance + L_body, l_exit, S, "exit"),
        ]
        z += L
        out.append(dict(z_entry=z_entry, z_exit=z, lines=lines))
        if i < n_sec - 1:                                    # gap/2 + real-length zero-K1 quad + gap/2
            z += gap + drv.section_quad_length_m(sec)
    return out


def _ez_section(zarr, t, base_rad, lines, omega):
    """Sum the 4 lines' on-axis Ez [V/m] at the particle z-positions for one section at time t."""
    total = np.zeros_like(zarr)
    for ln in lines:
        s = np.interp(zarr, ln["zgrid"], ln["shape"], left=0.0, right=0.0)   # 0 outside the element
        total += s * math.cos(omega * t + base_rad + ln["offset_rad"])
    return total


def _push(z, u, t, lines, base_rad, z_stop, omega, kq, dt):
    """RK4-integrate the pencil bunch (dz/dt = c u/gamma, du/dt = kq Ez) through one section's field
    until every particle clears z_stop. kq = +e/(m c) for positrons. Returns (z, u, t)."""
    z, u = z.copy(), u.copy()

    def dstate(zz, uu, tt):
        g = np.sqrt(1.0 + uu * uu)
        return C_LIGHT * uu / g, kq * _ez_section(zz, tt, base_rad, lines, omega)

    n_max = int(1.5 * (z_stop - z.min()) / (C_LIGHT * dt)) + 100
    for _ in range(n_max):
        dz1, du1 = dstate(z, u, t)
        dz2, du2 = dstate(z + 0.5 * dt * dz1, u + 0.5 * dt * du1, t + 0.5 * dt)
        dz3, du3 = dstate(z + 0.5 * dt * dz2, u + 0.5 * dt * du2, t + 0.5 * dt)
        dz4, du4 = dstate(z + dt * dz3, u + dt * du3, t + dt)
        z = z + (dt / 6.0) * (dz1 + 2 * dz2 + 2 * dz3 + dz4)
        u = u + (dt / 6.0) * (du1 + 2 * du2 + 2 * du3 + du4)
        t += dt
        if z.min() >= z_stop:
            break
    return z, u, t


def _mean_ke_mev(u, w):
    return float(np.average((np.sqrt(1.0 + u * u) - 1.0) * MC2_MEV, weights=w))


def find_crest(scan_coarse, scan_fine):
    """Crest base phase [deg, in [0,360)] from two phase->KE callables (a light bunch for the wide
    0..360 scan, the full bunch for the +-FINE_HALF fine scan about the max), then a parabolic refine
    of the top 3 (mirrors sim/autophase.find_crest). Returns (crest_deg, exit_ke_mev_at_crest)."""
    coarse = np.arange(0.0, 360.0, COARSE_STEP_DEG)
    c0 = coarse[int(np.argmax([scan_coarse(b) for b in coarse]))]
    fine = c0 + np.arange(-FINE_HALF_DEG, FINE_HALF_DEG + 1e-9, FINE_STEP_DEG)
    g = np.array([scan_fine(b) for b in fine])
    i = int(np.argmax(g))
    if 0 < i < len(fine) - 1:                                # parabolic vertex of the top 3
        y0, y1, y2 = g[i - 1], g[i], g[i + 1]
        denom = y0 - 2 * y1 + y2
        d = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
        crest = fine[i] + d * (fine[1] - fine[0])
    else:
        crest = fine[i]
    crest = float(crest % 360.0)
    return crest, scan_fine(crest)


def _phase_yaml_lines(path, found, indices):
    """Section-aware writeback: the sections are inline-flow mappings
    `  - {name: "CEA 4", ..., crest_phase_deg: 179.5768}    # sec 5`, so a top-level `key:` regex does
    NOT match. Walk the file; for each `^\\s*- \\{` line carrying a `crest_phase_deg`, if its running
    section index is in `indices`, substitute the number (keeping every other field + the trailing
    `# sec N` comment). Returns the rewritten text and a list of (index, old_str, new_str)."""
    with open(path) as fh:
        lines = fh.readlines()
    pat = re.compile(r"(crest_phase_deg:\s*)(\S+?)(\s*[},])")
    new_to = dict(zip(indices, found))
    changes = []
    sec_i = 0
    for li, line in enumerate(lines):
        if not re.match(r"^\s*- \{", line) or "crest_phase_deg" not in line:
            continue
        if sec_i in new_to:
            m = pat.search(line)
            old = m.group(2) if m else "?"
            new = f"{new_to[sec_i]:.4f}"
            lines[li] = pat.sub(lambda mm: mm.group(1) + new + mm.group(3), line, count=1)
            changes.append((sec_i, old, new))
        sec_i += 1
    return "".join(lines), changes


VERIFY_DELTA_DEG = 12.0      # +-offset for the --verify crest-confirmation phase scan


def _verify_impact(drv, cfg, found, last_idx):
    """Confirm the found crest with Impact-T directly: pin every section to its found crest, then run
    the deck at the LAST phased section's crest and crest +- VERIFY_DELTA (3 Impact-T runs, pencil
    core + aperture off). The found crest is confirmed if the centre run's exit KE is the maximum.
    This tests the crest -- the tool's deliverable -- and is immune to the constant field-amplitude
    offset between `ele_field` and Fortran Impact-T (a uniform scale does not move the argmax)."""
    import copy
    import shutil
    workdir = "logs/diags/autophase_impact_verify"

    P_in, _ = drv.load_converter_core(cfg)
    pmag = np.sqrt(np.asarray(P_in.px) ** 2 + np.asarray(P_in.py) ** 2 + np.asarray(P_in.pz) ** 2)
    n = P_in.n_particle
    P_in.x, P_in.y, P_in.px, P_in.py, P_in.pz = (np.zeros(n), np.zeros(n), np.zeros(n),
                                                 np.zeros(n), pmag)
    cfg_noap = copy.deepcopy(cfg)
    cfg_noap["lattice"]["bore_aperture_on"] = False          # no scrape: isolate the longitudinal physics

    print(f"\n── --verify: Impact-T crest check, section {last_idx + cfg['lattice']['first_section']} "
          f"at crest +- {VERIFY_DELTA_DEG:g}° (pencil core, aperture off) ──────")
    results = []
    for off in (-VERIFY_DELTA_DEG, 0.0, VERIFY_DELTA_DEG):
        if os.path.isdir(workdir):
            shutil.rmtree(workdir)
        os.makedirs(workdir, exist_ok=True)
        I, _, _ = drv.build_impact(cfg_noap, workdir=workdir)
        for i, sec in enumerate(cfg_noap["sections"]):
            gname = drv._ensure_section_group(I, cfg_noap, i)
            drv._set_section_phase(I, cfg_noap, i, found[i] + (off if i == last_idx else 0.0))
            drv._set_group_scale(I, gname, sec["field_scale"])
        I.initial_particles = P_in
        I.configure()
        I.run()
        P = I.particles.get("final_particles") if I.finished and not I.error else None
        ke = (float(P["mean_energy"]) / 1e6 - MC2_MEV) if P is not None and P.n_particle else float("nan")
        results.append((off, ke))
        print(f"  crest {off:+5.1f}° → Impact-T exit ⟨KE⟩ {ke:.3f} MeV", flush=True)

    kes = [r[1] for r in results]
    if np.argmax(kes) == 1:
        print("  ✓ crest confirmed: Impact-T exit KE peaks at the found crest.")
    else:
        print("  ✗ crest NOT at the Impact-T peak — investigate the field/phase reproduction.")


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    dry = "--dry-run" in sys.argv or "-n" in sys.argv
    verify = "--verify" in sys.argv

    prepare_env()
    drv = _load_driver()
    cfg = drv.load_config(CONFIG)
    first = cfg["lattice"]["first_section"]
    n_sec = len(cfg["sections"])
    omega = 2.0 * np.pi * cfg["rf"]["rf_freq_hz"]
    dt = 1.0 / (STEPS_PER_PERIOD * cfg["rf"]["rf_freq_hz"])  # resolve the RF cycle
    sp = str(cfg["beam"].get("species", "positrons"))
    kq = (1.0 if sp.startswith("positron") else -1.0) * E_CHARGE / (M_E * C_LIGHT)

    machine = [int(a) for a in args] if args else list(range(first, first + n_sec))
    indices = sorted({m - first for m in machine})
    if any(idx < 0 or idx >= n_sec for idx in indices):
        sys.exit(f"usage: python sim/autophase_impact.py "
                 f"[{' '.join(str(first + j) for j in range(n_sec))}] [--dry-run] [--verify]")

    # ── Pencil-ise the converter positron core (longitudinal scan; matches sim/autophase.py) ──────
    # The bare core's ~600 mrad divergence would scrape the bore before any section ends (no capture
    # optic here), so a divergent bunch yields zero survivors at every phase. Redirect each particle's
    # momentum onto +z (preserving energy), centre on-axis: the scan then measures pure energy gain.
    P_in, info = drv.load_converter_core(cfg)
    pmag = np.sqrt(np.asarray(P_in.px) ** 2 + np.asarray(P_in.py) ** 2 + np.asarray(P_in.pz) ** 2)
    z = np.asarray(P_in.z, float)
    u = pmag / MC2_EV                                        # u = gamma*beta = pc/mc^2 (px=py=0 now)
    w = np.asarray(P_in.weight, float)
    if z.size > SUBSAMPLE:                                   # weighted resample -> equal weight
        rng = np.random.default_rng(cfg["beam"]["rng_seed"])
        sel = rng.choice(z.size, SUBSAMPLE, replace=False, p=w / w.sum())
        z, u, w = z[sel], u[sel], np.ones(SUBSAMPLE)
    print(f"\nPositron core: {z.size} parts (pencil-ised for the longitudinal scan), "
          f"<KE>_in {info['ke_in_mev']:.3f} MeV, beta_min {info['beta_min_core']:.5f}\n", flush=True)

    secf = _build_section_fields(drv, cfg)

    # Persistent bunch state (z, u, t) advances section by section: each section's crest is found in
    # the field of all EARLIER sections (pinned to their found/frozen crest), entering with the right
    # energy and absolute arrival time -- the latter is what the absolute theta0 crest depends on.
    found = [float(s["crest_phase_deg"]) for s in cfg["sections"]]
    zc, uc, tc = z, u, 0.0
    for i in range(max(indices) + 1):
        sec = cfg["sections"][i]
        sf = secf[i]
        z_stop = sf["z_exit"] + 0.02                         # clear the field, then coast
        if i in indices:
            n_c = min(COARSE_SUB, zc.size)                   # zc is a random resample -> first n_c is a subset

            def scan(phi, sub, sf=sf, zc=zc, uc=uc, tc=tc, z_stop=z_stop):
                _, ue, _ = _push(zc[sub], uc[sub], tc, sf["lines"], math.radians(phi),
                                 z_stop, omega, kq, dt)
                return _mean_ke_mev(ue, w[sub])
            old = found[i]
            crest, ke_out = find_crest(lambda b: scan(b, slice(0, n_c)),
                                       lambda b: scan(b, slice(None)))
            found[i] = crest
            print(f"── Section {i + first} ({sec['name']}) ──────────────────────────")
            print(f"  scale={sec['field_scale']:.4g}, deck z [{sf['z_entry']:.2f}, {sf['z_exit']:.2f}] m")
            print(f"  crest crest_phase_deg = {crest:.4f}°  (was {old:.4f}°, Δ {crest - old:+.4f}°)")
            print(f"  model ⟨KE⟩ {ke_out:.3f} MeV (relative; ~1.7× Impact-T magnitude, crest-only)\n",
                  flush=True)
        # Advance the persistent state through this section at its (found or frozen) crest.
        zc, uc, tc = _push(zc, uc, tc, sf["lines"], math.radians(found[i]), z_stop, omega, kq, dt)

    print(f"Final model ⟨KE⟩ after section {max(indices) + first}: {_mean_ke_mev(uc, w):.3f} MeV "
          f"(relative model energy; not the stage energy)", flush=True)

    if not dry:
        crests = [found[idx] for idx in indices]
        new_txt, changes = _phase_yaml_lines(CONFIG, crests, indices)
        with open(CONFIG, "w") as fh:
            fh.write(new_txt)
        for idx, old, new in changes:
            print(f"Wrote {CONFIG}: sec {idx + first} crest_phase_deg {old} -> {new}")
        if changes:
            print("\nSetpoints updated. Re-run sim/linac5-8.py to propagate the new crests.")
    else:
        print("--dry-run: YAML unchanged.")

    if verify:
        _verify_impact(drv, cfg, found, max(indices))


if __name__ == "__main__":
    main()
