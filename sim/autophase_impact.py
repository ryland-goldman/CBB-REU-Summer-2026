"""
Auto-phase the linac 4-8 (Impact-T) RF sections for the POSITRON beam and rewrite the YAML.

A standalone tool — NOT wired into the chain. The frozen per-section crest_phase_deg in
config/linac4-8.yaml were derived for ELECTRONS; positrons (q=+e) crest ~180 deg away and their
lower injection energy/velocity shifts the ABSOLUTE Impact-T theta0 (referenced to t=0), so the
crest must be re-found on the real deck with the positron core injected. For each section it builds
the deck TRUNCATED to that section (earlier sections pinned to their already-found crest), scans the
section base phase, and takes the bunch-averaged exit-energy maximum (coarse -> fine -> parabolic).

Runtime caveat: this drives Impact-T O(sections x scan-points) times — minutes to tens of minutes.
DO NOT wire it into main(). See docs/linac4-8.md (Positron mode).

  python sim/autophase_impact.py            # phase sections 4 5 6 7 8, rewrite the YAML
  python sim/autophase_impact.py 4 5        # only sections 4-5
  python sim/autophase_impact.py --dry-run  # scan + report, write nothing
"""

import importlib.util
import math
import os
import re
import shutil
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np

from sim.helpers.tools import C_LIGHT, prepare_env

CONFIG = "config/linac4-8.yaml"
WORKDIR = "logs/diags/autophase_impact"

COARSE_STEP_DEG = 30.0       # wide 0..360 scan resolution
FINE_HALF_DEG = 15.0         # +-window about the coarse max
FINE_STEP_DEG = 5.0          # fine scan resolution
FEW_SURVIVORS_FRAC = 0.10    # WARN below this surviving fraction (crest maximises the surviving core)


def _load_driver():
    """Import sim/linac4-8.py (hyphenated, not a normal module name) for its deck builders so the
    deck, header, and positron-core handoff are byte-identical to what the driver runs."""
    path = os.path.join(os.path.dirname(__file__), "linac4-8.py")
    spec = importlib.util.spec_from_file_location("linac48_driver", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _exit_ke_mev(drv, I):
    """Run the configured deck and return (bunch-mean exit KE [MeV], surviving macro count).

    Re-uses the SAME Impact object across phases: only theta0 changed + reconfigure(), so the deck
    geometry and initial particles are untouched. (If lume-impact 0.11.0 ever proves unreliable
    re-running one Impact object, rebuild inside the scan loop instead.)"""
    I.run()
    if not I.finished or I.error:
        raise RuntimeError(f"Impact-T did not finish cleanly (finished={I.finished}, "
                           f"error={I.error})")
    P = I.particles.get("final_particles") if hasattr(I.particles, "get") else (
        I.particles["final_particles"] if "final_particles" in I.particles else None)
    if P is None or P.n_particle == 0:                        # whole bunch lost -> deprioritise phase
        return float("-inf"), 0
    return float(P["mean_energy"]) / 1e6 - drv.MC2_MEV, int(P.n_particle)


def find_crest(scan):
    """Crest base phase [deg, in [0,360)] from a phase->KE callable `scan(phi)->ke_mev`.

    Coarse 0..360 scan, +-FINE_HALF window fine scan about the max, parabolic refine of the top 3
    (mirrors sim/autophase.find_crest). Returns (crest_deg, exit_ke_at_crest_mev, survivors_at_crest).
    """
    coarse = np.arange(0.0, 360.0, COARSE_STEP_DEG)
    cg = np.array([scan(b)[0] for b in coarse])
    c0 = coarse[int(np.argmax(cg))]

    fine = c0 + np.arange(-FINE_HALF_DEG, FINE_HALF_DEG + 1e-9, FINE_STEP_DEG)
    fg = np.array([scan(b)[0] for b in fine])
    i = int(np.argmax(fg))
    if 0 < i < len(fine) - 1:                                 # parabolic vertex of the top 3
        y0, y1, y2 = fg[i - 1], fg[i], fg[i + 1]
        denom = y0 - 2 * y1 + y2
        d = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
        crest = fine[i] + d * (fine[1] - fine[0])
    else:
        crest = fine[i]
    crest = float(crest % 360.0)
    ke, surv = scan(crest)                                    # exact value at the refined crest
    return crest, ke, surv


def _truncated_ntstep(drv, cfg, sections):
    """A trimmed Ntstep for a deck truncated to `sections`: Impact-T runs a FIXED Ntstep, and the
    production 200000 is sized for the full ~29 m deck — wasteful (slow) for a short truncated deck.
    Size it from the truncated length so the bunch clears the deck with margin."""
    _I, total_len, _bounds = drv.build_impact({**cfg, "sections": sections}, workdir=WORKDIR)
    dt = float(cfg["deck"]["dt"])
    n = math.ceil(total_len / (0.6 * C_LIGHT * dt)) * 1.3     # 0.6c floor + 30% margin
    return int(n)


def _phase_yaml_lines(path, found, indices):
    """Section-aware writeback: the sections are inline-flow mappings
    `  - {name: "CU 5", ..., crest_phase_deg: 68.7787}    # sec 4`, so autophase.set_yaml_param's
    top-level `key:` regex does NOT match. Walk the file; for each `^\\s*- \\{` line carrying a
    `crest_phase_deg`, if its running section index is in `indices`, substitute the number (keeping
    every other field + the trailing `# sec N` comment). Returns the rewritten text and a list of
    (index, old_str, new_str) for the report."""
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
            crest = new_to[sec_i]
            m = pat.search(line)
            old = m.group(2) if m else "?"
            new = f"{crest:.4f}"
            lines[li] = pat.sub(lambda mm: mm.group(1) + new + mm.group(3), line, count=1)
            changes.append((sec_i, old, new))
        sec_i += 1
    return "".join(lines), changes


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    dry = "--dry-run" in sys.argv or "-n" in sys.argv

    prepare_env()
    drv = _load_driver()
    cfg = drv.load_config(CONFIG)
    first = cfg["lattice"]["first_section"]
    n_sec = len(cfg["sections"])

    machine = [int(a) for a in args] if args else list(range(first, first + n_sec))
    indices = sorted({m - first for m in machine})
    if any(idx < 0 or idx >= n_sec for idx in indices):
        sys.exit(f"usage: python sim/autophase_impact.py "
                 f"[{' '.join(str(first + j) for j in range(n_sec))}] [--dry-run]")

    if os.path.isdir(WORKDIR):
        shutil.rmtree(WORKDIR)
    os.makedirs(WORKDIR, exist_ok=True)

    P_in, info = drv.load_sec3_core(cfg)
    ke_in = info["ke_in_mev"]
    n_in = int(P_in.n_particle)

    # Crest is a LONGITUDINAL quantity. The bare positron core's ~600 mrad divergence scrapes on the
    # bore before any section ends (no capture optic modelled here), so a divergent bunch yields zero
    # survivors at every phase. Pencil-ise it -- redirect each particle's momentum onto +z (preserving
    # its energy) and centre it on-axis -- so the scan measures pure energy gain, exactly as the WarpX
    # autophase.py 1D longitudinal model does. The full divergent beam is what sim/linac4-8.py runs.
    pmag = np.sqrt(np.asarray(P_in.px) ** 2 + np.asarray(P_in.py) ** 2 + np.asarray(P_in.pz) ** 2)
    P_in.x = np.zeros(n_in)
    P_in.y = np.zeros(n_in)
    P_in.px = np.zeros(n_in)
    P_in.py = np.zeros(n_in)
    P_in.pz = pmag
    print(f"\nPositron core: {n_in} parts (pencil-ised for the longitudinal scan), "
          f"<KE>_in {ke_in:.3f} MeV, beta_min {info['beta_min_core']:.5f}\n", flush=True)

    found = []                                                # crest per section, machine order
    found_idx = []
    max_target = max(indices)
    for i in range(max_target + 1):
        sec = cfg["sections"][i]
        if i not in indices:                                  # not requested: keep the YAML crest
            found.append(float(sec["crest_phase_deg"]))
            continue

        # Deck truncated from the start through section i; earlier sections (whether requested or
        # not) are pinned to found[j] so this section's crest is found in their actual field.
        cfg_i = dict(cfg)
        cfg_i["sections"] = cfg["sections"][:i + 1]
        cfg_i["deck"] = {**cfg["deck"], "ntstep": _truncated_ntstep(drv, cfg, cfg_i["sections"])}

        I, total_len_i, _bounds = drv.build_impact(cfg_i, workdir=WORKDIR)
        for j in range(i):                                    # pin earlier sections to their crest
            gname = drv._ensure_section_group(I, cfg, j)
            drv._set_group_scale(I, gname, cfg["sections"][j]["field_scale"])
            drv._set_section_phase(I, cfg, j, found[j])
        gname_i = drv._ensure_section_group(I, cfg, i)
        drv._set_group_scale(I, gname_i, sec["field_scale"])
        I.initial_particles = P_in

        def scan(phi):
            drv._set_section_phase(I, cfg, i, phi)
            I.configure()
            return _exit_ke_mev(drv, I)

        old = float(sec["crest_phase_deg"])
        print(f"-- Section {i + first} ({sec['name']}) -------------------------------------")
        print(f"  scale={sec['field_scale']:.4g}, truncated deck {total_len_i:.2f} m, "
              f"Ntstep={cfg_i['deck']['ntstep']}", flush=True)
        crest, ke_out, surv = find_crest(scan)
        found.append(crest)
        found_idx.append(i)

        surv_frac = surv / n_in if n_in else 0.0
        print(f"  crest crest_phase_deg = {crest:.4f} deg  (was {old:.4f} deg, "
              f"Delta {crest - old:+.4f} deg)")
        print(f"  exit <KE> {ke_out:.3f} MeV (Delta from inject {ke_out - ke_in:+.3f} MeV); "
              f"survivors {surv}/{n_in} ({surv_frac * 100:.1f}%)", flush=True)
        if surv_frac < FEW_SURVIVORS_FRAC:
            print(f"  WARN: only {surv_frac * 100:.1f}% survive — the degraded positron beam may "
                  f"lose most particles; the crest then maximises the SURVIVING core.", flush=True)
        print(flush=True)

    if dry:
        print("--dry-run: YAML unchanged.")
        return

    crests = [found[idx] for idx in found_idx]
    new_txt, changes = _phase_yaml_lines(CONFIG, crests, found_idx)
    with open(CONFIG, "w") as fh:
        fh.write(new_txt)
    for idx, old, new in changes:
        print(f"Wrote {CONFIG}: sec {idx + first} crest_phase_deg {old} -> {new}")
    if changes:
        print("\nSetpoints updated. Re-run sim/linac4-8.py to propagate the new crests.")


if __name__ == "__main__":
    main()
