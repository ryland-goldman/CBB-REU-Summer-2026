"""Cornell Linac e+/e- converter target (G4beamline). main(): read the linac3 exit electron beam
(the 3->4 boundary) -> resample to the incident-event count + write a BLTrackFile -> generate the
.g4bl deck -> run g4bl (7 mm W target, brems->pair production) -> sample the exit plane -> keep the
forward positron core -> openPMD handoff OUT (`positrons`, +e) + injection_summary.json.

Run `python sim/converter.py [n_events]` (the optional arg overrides config n_events for a quick
test). sim/plot/converter.py makes the figures. See docs/converter.md.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

os.environ.setdefault("OMP_NUM_THREADS", os.environ.get("OMP_THREADS", "1"))

import json
import math
import shutil
import subprocess

import numpy as np
import yaml

from sim.helpers.tools import MC2_EV, M_E, E_CHARGE as Q_E, prepare_env
from sim.helpers.loadparticles import read_warpx_dump, write_openpmd_particles, upstream_exit_lab_z
from sim.helpers.tqdmwrapper import g4bl_progress
from sim.helpers import g4bl

CONFIG = "config/converter.yaml"
MC2_MEV = MC2_EV / 1e6
_G4BL_FALLBACK = "/Applications/G4beamline-3.08.app/Contents/MacOS/g4bl"


def load_config(path=CONFIG):
    """Load the converter YAML (sci-notation numeric fields coerced as in linac4-8.py)."""
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    cfg["physics"]["n_events"] = int(cfg["physics"]["n_events"])
    return cfg


def _resolve_g4bl(exe):
    """Resolve the g4bl binary: PATH first, then the literal path, then the macOS app bundle."""
    return shutil.which(exe) or (exe if os.path.exists(exe) else None) or (
        _G4BL_FALLBACK if os.path.exists(_G4BL_FALLBACK) else None) or _raise_no_g4bl(exe)


def _raise_no_g4bl(exe):
    raise RuntimeError(
        f"g4bl not found (config g4bl.exe={exe!r}); install G4beamline 3.08 and add it to PATH "
        f"or set g4bl.exe to its absolute path")


def load_incident_beam(cfg):
    """Read the linac3 exit electrons and resample (with replacement, charge-weighted) to n_events
    incident tracks placed just upstream of the target front face. Returns
    (ParticleGroup incident, q_incident_C, w_evt_C, z_inject_lab_m).
    """
    diag = cfg["io"]["sec3_particles"]
    n_events = cfg["physics"]["n_events"]
    rng = np.random.default_rng(cfg["physics"]["rng_seed"])

    P = read_warpx_dump(diag)                                  # full sec3 exit (last dump), electrons
    q_incident = float(P["charge"])                            # honest incident charge denominator
    z_inject_lab = upstream_exit_lab_z(cfg["io"]["sec3_summary"], float(P["mean_z"]))

    w = np.asarray(P.weight, dtype=float)
    idx = rng.choice(P.n_particle, size=n_events, replace=True, p=w / w.sum())
    Pin = P[idx]                                               # bootstrap (each event seeds its own RNG)

    # Place the bunch head `front_clearance_mm` upstream of the target front face (g4bl injects
    # tracks at their file z; relativistic e- drift straight into the target).
    g = cfg["geometry"]
    z_start_m = (g["target_front_z_mm"] - g["front_clearance_mm"]) * 1e-3
    Pin.z = Pin.z - float(np.max(Pin.z)) + z_start_m

    w_evt = q_incident / n_events                             # uniform charge per incident event
    return Pin, q_incident, w_evt, z_inject_lab


def coil_current_density(b_tesla, inner_mm, outer_mm, length_mm):
    """Conductor current density [A/mm^2] giving a central field of `b_tesla` for a uniform thick
    solenoid (inner/outer radius, length [mm]). Exact thick-solenoid on-axis centre formula
    B = mu0*J*b*ln[(a2+sqrt(a2^2+b^2))/(a1+sqrt(a1^2+b^2))], b=L/2; matches g4bl's coil to <1e-5."""
    mu0 = 4e-7 * math.pi
    a1, a2, b = inner_mm * 1e-3, outer_mm * 1e-3, length_mm * 1e-3 / 2.0
    fac = mu0 * b * math.log((a2 + math.hypot(a2, b)) / (a1 + math.hypot(a1, b)))
    return b_tesla / fac / 1e6                                # A/m^2 -> A/mm^2


def _solenoid_deck(cfg, front):
    """The capture-solenoid lines (a real `coil` + `solenoid`, so g4bl computes the Maxwellian field
    with end fringe) for the deck, or '' when disabled. Coil spans `start_z_mm`..+`length_mm` rel. to
    the target front; current density is auto-solved so the central field equals `b_tesla`."""
    sol = cfg.get("solenoid", {})
    if not sol.get("enabled", False):
        return ""
    L = sol["length_mm"]
    center = front + sol["start_z_mm"] + L / 2.0
    j = coil_current_density(sol["b_tesla"], sol["inner_radius_mm"], sol["outer_radius_mm"], L)
    return (
        f"coil Capture innerRadius={sol['inner_radius_mm']:g} outerRadius={sol['outer_radius_mm']:g} "
        f"length={L:g} material=Cu\n"
        f"solenoid CaptureSol coilName=Capture current={j:.6g}\n"
        f"place CaptureSol z={center:g}\n")


def write_g4bl_deck(cfg, path, beam_in, out_file):
    """Generate the .g4bl deck (7 mm W target + capture solenoid + forward sampling plane). Returns
    sample_z_mm."""
    g, ph, det = cfg["geometry"], cfg["physics"], cfg["detector"]
    front, L = g["target_front_z_mm"], g["target_length_mm"]
    center = front + L / 2.0
    sol = cfg.get("solenoid", {})
    # Sample plane: when the capture coil is on, drift `exit_drift_mm` past its exit so the handoff is
    # field-free (the exit fringe has died -> the fringe focusing kick is fully applied). Else sample
    # `back_clearance_mm` past the target back face.
    if sol.get("enabled", False):
        sample_z = front + sol["start_z_mm"] + sol["length_mm"] + sol["exit_drift_mm"]
    else:
        sample_z = front + L + det["back_clearance_mm"]
    deck = (
        "#  converter.g4bl -- generated by sim/converter.py -- DO NOT EDIT\n"
        f"#  e+/e- converter: linac3 exit electrons -> {L:g} mm {g['target_material']} target "
        f"-> positron sampling plane\n"
        f"physics {ph['list']} minRangeCut={ph['min_range_cut_mm']:g}\n"
        f"beam ascii file={beam_in} nEvents={ph['n_events']}\n"
        f"cylinder Target outerRadius={g['target_radius_mm']:g} length={L:g} "
        f"material={g['target_material']} maxStep={ph['max_step_mm']:g} color=1,0,0\n"
        f"place Target z={center:g}\n"
        f"{_solenoid_deck(cfg, front)}"
        f"virtualdetector Sample radius={det['radius_mm']:g} length=1 format=ASCII "
        f"file={out_file} color=0,1,0\n"
        f"place Sample z={sample_z:g}\n")
    with open(path, "w") as fh:
        fh.write(deck)
    return sample_z


def run_g4bl(cfg, workdir):
    """Run g4bl on the generated deck (cwd=workdir, so its relative paths + g4beamline.root land
    there), streaming a tqdm bar from its `Event N Completed` stdout. Raises on failure."""
    exe = _resolve_g4bl(cfg["g4bl"]["exe"])
    deck = cfg["io"]["deck_name"]
    print(f"Running g4bl ({cfg['physics']['n_events']} events, {cfg['physics']['list']}) ...",
          flush=True)
    proc = subprocess.Popen([exe, deck], cwd=workdir, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, bufsize=1)
    tail = []
    for line in g4bl_progress(proc.stdout, cfg["physics"]["n_events"], desc="converter"):
        tail.append(line)
        if len(tail) > 40:
            tail.pop(0)
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f"g4bl exited {proc.returncode}:\n{''.join(tail)}")


def _lepton_cuts(df, cfg):
    """Forward (Pz>0) + min-KE mask on a single-species BLTrackFile DataFrame, plus its KE [MeV]."""
    ke = g4bl.bltrack_ke_mev(df)
    keep = ke >= cfg["detector"]["min_ke_mev"]
    if cfg["detector"].get("forward_only", True):
        keep = keep & (df["Pz"].values > 0.0)
    return keep, ke


def main():
    prepare_env()
    cfg = load_config()
    if len(sys.argv) > 1:                                     # quick-test override: python ... N
        cfg["physics"]["n_events"] = int(sys.argv[1])

    outdir, workdir = cfg["io"]["outdir"], cfg["io"]["workdir"]
    if os.path.isdir(workdir):
        shutil.rmtree(workdir)
    os.makedirs(outdir, exist_ok=True)                       # makes workdir too (outdir is under it)

    # ── Incident e- beam -> BLTrackFile ────────────────────────────────────────────────────
    Pin, q_incident, w_evt, z_inject_lab = load_incident_beam(cfg)
    n_events = cfg["physics"]["n_events"]
    beam_in = os.path.join(workdir, cfg["io"]["beam_in_file"])
    df_in = g4bl.particlegroup_to_bltrack_df(
        Pin, pdgid=g4bl.PDG_ELECTRON, event_ids=np.arange(1, n_events + 1, dtype=np.int64))
    g4bl.write_bltrackfile(df_in, beam_in)
    print(f"Incident: {q_incident*1e12:.2f} pC sec3-exit e- resampled to {n_events} events "
          f"(inject lab-z {z_inject_lab:.3f} m)", flush=True)

    # ── Generate + run the g4bl deck ───────────────────────────────────────────────────────
    sample_z_mm = write_g4bl_deck(cfg, os.path.join(workdir, cfg["io"]["deck_name"]),
                                  cfg["io"]["beam_in_file"], cfg["io"]["out_file"])
    run_g4bl(cfg, workdir)

    out_path = os.path.join(workdir, cfg["io"]["out_file"])
    if not os.path.isfile(out_path) or os.path.getsize(out_path) == 0:
        raise RuntimeError(f"g4bl produced no sampling output at {out_path}")
    df = g4bl.read_bltrackfile(out_path)

    # ── Split species, apply the forward/KE cut to the positron core ───────────────────────
    df_pos = df[df["PDGid"].values == g4bl.PDG_POSITRON].reset_index(drop=True)
    df_ele = df[df["PDGid"].values == g4bl.PDG_ELECTRON].reset_index(drop=True)
    n_gamma = int(np.count_nonzero(df["PDGid"].values == g4bl.PDG_GAMMA))
    if len(df_pos) == 0:
        raise RuntimeError("no positrons reached the sampling plane -- raise n_events or check the deck")
    keep, ke_pos_all = _lepton_cuts(df_pos, cfg)
    df_keep = df_pos[keep].reset_index(drop=True)
    ke_pos = ke_pos_all[keep]
    if len(df_keep) == 0:
        raise RuntimeError(
            f"no positrons survived the cut (min_ke={cfg['detector']['min_ke_mev']} MeV, "
            f"forward_only); raise n_events or relax the cut")

    # ── Yield bookkeeping: each kept e+ inherits its event's charge -> q_pos = q_in * yield ──
    q_pos = float(len(df_keep) * w_evt)
    yield_pos = len(df_keep) / n_events
    yield_pos_raw = len(df_pos) / n_events

    pg_pos = g4bl.bltrack_df_to_particlegroup(df_keep, species="positron",
                                              weight_C=np.full(len(df_keep), w_evt))
    pg_pos.z = pg_pos.z - float(pg_pos["mean_z"])             # converter-local frame, centroid at 0

    # ── Handoff OUT ────────────────────────────────────────────────────────────────────────
    part_dir = os.path.join(outdir, "particles")
    os.makedirs(part_dir, exist_ok=True)
    write_openpmd_particles(pg_pos, part_dir, iteration=0, species="positrons",
                            charge=+Q_E, mass=M_E)

    pt = np.hypot(df_keep["Px"].values, df_keep["Py"].values)
    div_mrad = 1e3 * pt / np.abs(df_keep["Pz"].values)
    rr_mm = np.hypot(df_keep["x"].values, df_keep["y"].values)
    _sol = cfg.get("solenoid", {})
    sol_on = bool(_sol.get("enabled", False))
    inj = dict(
        q_injected_C=q_incident,                             # honest denominator (full sec3 exit e-)
        q_incident_C=q_incident, n_events=n_events,
        z_inject_lab_m=float(z_inject_lab + sample_z_mm * 1e-3),
        z_inject_mean_m=0.0,
        target_material=cfg["geometry"]["target_material"],
        target_length_mm=cfg["geometry"]["target_length_mm"],
        target_radius_mm=cfg["geometry"]["target_radius_mm"],
        physics_list=cfg["physics"]["list"],
        min_range_cut_mm=cfg["physics"]["min_range_cut_mm"],
        max_step_mm=cfg["physics"]["max_step_mm"],
        min_ke_mev_cut=cfg["detector"]["min_ke_mev"],
        forward_only=bool(cfg["detector"].get("forward_only", True)),
        # Capture coil (the focusing optic; peak field, length, auto-solved current density)
        capture_solenoid_on=sol_on,
        capture_solenoid_b_tesla=float(_sol.get("b_tesla", 0.0)) if sol_on else 0.0,
        capture_solenoid_length_mm=float(_sol.get("length_mm", 0.0)) if sol_on else 0.0,
        capture_solenoid_current_a_per_mm2=(
            coil_current_density(_sol["b_tesla"], _sol["inner_radius_mm"],
                                 _sol["outer_radius_mm"], _sol["length_mm"]) if sol_on else 0.0),
        # Yields (per incident electron)
        yield_positron=yield_pos, yield_positron_raw=yield_pos_raw,
        yield_electron=len(df_ele) / n_events, yield_gamma=n_gamma / n_events,
        n_positron_out=int(len(df_keep)), n_positron_raw=int(len(df_pos)),
        n_electron_out=int(len(df_ele)), n_gamma_out=n_gamma,
        q_positron_out_C=q_pos,
        # Positron core spectra (kept set)
        ke_pos_mean_mev=float(np.mean(ke_pos)), ke_pos_sigma_mev=float(np.std(ke_pos)),
        ke_pos_min_mev=float(np.min(ke_pos)), ke_pos_max_mev=float(np.max(ke_pos)),
        sigma_r_pos_mm=float(np.sqrt(np.mean(rr_mm ** 2))),
        div_pos_rms_mrad=float(np.sqrt(np.mean(div_mrad ** 2))),
        sample_z_mm=float(sample_z_mm),
    )
    with open(os.path.join(outdir, "injection_summary.json"), "w") as fh:
        json.dump(inj, fh, indent=2)

    print(f"\nDone. {len(df_pos)} e+ at plane, {len(df_keep)} forward >= "
          f"{cfg['detector']['min_ke_mev']} MeV ({q_pos*1e15:.2f} fC); "
          f"yield {yield_pos:.4f} e+/e-; <KE>_e+ {inj['ke_pos_mean_mev']:.2f} MeV "
          f"(sigma {inj['ke_pos_sigma_mev']:.2f}); div_rms {inj['div_pos_rms_mrad']:.1f} mrad "
          f"-> {outdir}/", flush=True)


if __name__ == "__main__":
    main()
