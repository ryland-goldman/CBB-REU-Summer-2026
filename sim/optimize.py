"""
Multi-objective optimization of the Cornell linac chain with Xopt/CNSGA.

Reads config/xopt.yaml (VOCS + CNSGA + executor), then for each population member: builds a unique
LINACSIM_OUT_DIR sandbox (isolation plan), seeds the frozen upstream dump (every stage upstream of
the earliest varied one -- the freeze boundary),
writes the variables into the sandbox config COPY (NEVER the canonical config/ -- phase variables
go to the crest_offset_deg / PHASE_OFFSET_DEG knobs autophase never rewrites, NOT crest_phase_deg/
PHASE_DEG), runs the chain as a subprocess, and reads the objectives/constraints from the (extended)
linac5-8 injection_summary.json. CNSGA returns the charge<->emittance Pareto front. See
docs/optimization-plan.md and docs/per-eval-isolation-plan.md.

  python sim/optimize.py            # run the optimization (config/xopt.yaml)

Objectives: maximize q_out_C, minimize eps_n=sqrt(eps_n_x*eps_n_y). Constraints: ke_out_mev,
sigma_E_rel, sigma_x_mm, transmission_core. A failed/empty-beam eval returns NaN objectives (CNSGA's
discard mechanism), so a diverged member never kills the worker.
"""

import hashlib
import math
import os
import re
import shutil
import signal
import subprocess
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np
import yaml

from sim.helpers.tools import REPO_ROOT
from sim.helpers.sandbox import make_out_dir

XOPT_CONFIG = os.path.join(REPO_ROOT, "config/xopt.yaml")

# Objective + constraint output names (the keys Xopt reads from each eval). Failure => all NaN.
OUTPUT_KEYS = ("q_out_C", "eps_n", "ke_out_mev", "sigma_E_rel", "sigma_x_mm", "transmission_core")

# Variable -> sandbox-config target. Each entry is (relative yaml path, addressing). Addressing is
# ("block", key) for a top-level-ish scalar `key: value` (keys are unique within the file), or
# ("section", index, key) for a key inside linac5-8.yaml's i-th inline-flow `- {...}` section dict.
# Phase variables map to the OFFSET keys (crest_offset_deg / PHASE_OFFSET_DEG) -- autophase owns the
# crest itself and would overwrite it. field_scale/quad/solenoid are amplitude knobs autophase never
# touches. PHASE_OFFSET_DEG is absent from the canonical linac1.yaml and is inserted under `params:`.
OVERRIDES = {
    "l5_field_scale_0":   ("config/linac5-8.yaml", ("section", 0, "field_scale")),
    "l5_field_scale_1":   ("config/linac5-8.yaml", ("section", 1, "field_scale")),
    "l5_field_scale_2":   ("config/linac5-8.yaml", ("section", 2, "field_scale")),
    "l5_field_scale_3":   ("config/linac5-8.yaml", ("section", 3, "field_scale")),
    "l5_phase_off_0":     ("config/linac5-8.yaml", ("section", 0, "crest_offset_deg")),
    "l5_phase_off_1":     ("config/linac5-8.yaml", ("section", 1, "crest_offset_deg")),
    "l5_phase_off_2":     ("config/linac5-8.yaml", ("section", 2, "crest_offset_deg")),
    "l5_phase_off_3":     ("config/linac5-8.yaml", ("section", 3, "crest_offset_deg")),
    "l5_quad_k1_1":       ("config/linac5-8.yaml", ("section", 1, "quad_k1")),
    "l5_quad_k1_2":       ("config/linac5-8.yaml", ("section", 2, "quad_k1")),
    "l5_sol_b_5":         ("config/linac5-8.yaml", ("section", 0, "solenoid_b_tesla")),
    "l5_sol_b_6":         ("config/linac5-8.yaml", ("section", 1, "solenoid_b_tesla")),
    "conv_target_len_mm": ("config/converter.yaml", ("block", "target_length_mm")),
    "conv_sol_b_tesla":   ("config/converter.yaml", ("block", "b_tesla")),
    "conv_exit_drift_mm": ("config/converter.yaml", ("block", "exit_drift_mm")),
    "l1_power_mw":        ("config/linac1.yaml", ("block", "POWER_MW")),
    "l1_phase_off":       ("config/linac1.yaml", ("block", "PHASE_OFFSET_DEG", "params")),
    "inj_i_sol0":         ("config/injector.yaml", ("block", "I_SOL0")),
    "inj_i_lens0e":       ("config/injector.yaml", ("block", "I_LENS0E")),
    "inj_preb1_kw":       ("config/injector.yaml", ("block", "PREB1_KW")),
    "inj_preb2_kw":       ("config/injector.yaml", ("block", "PREB2_KW")),
    "inj_preb1_phi":      ("config/injector.yaml", ("block", "PREB1_PHI_OFF")),
    "inj_preb2_phi":      ("config/injector.yaml", ("block", "PREB2_PHI_OFF")),
}

# Pipeline stage order (mirrors sim/main.py STAGES) and which stage each OVERRIDES config feeds.
# Every stage UPSTREAM of the earliest varied one is parameter-independent across the whole
# population: it runs ONCE and each eval symlinks its frozen dump instead of recomputing it. The
# scope=full VOCS varies injector onward, so cathode+gun (~16 min/eval) are frozen.
STAGE_ORDER = ["cathode", "gun", "injector", "linac1", "linac2", "linac3", "linac4",
               "converter", "linac5-8"]
CONFIG_STAGE = {
    "config/cathode.yaml": "cathode", "config/gun.yaml": "gun",
    "config/injector.yaml": "injector", "config/linac1.yaml": "linac1",
    "config/converter.yaml": "converter", "config/linac5-8.yaml": "linac5-8",
}
# Boundary `--from` stage -> (frozen-upstream dump symlinked into each sandbox, the `--to` stage of
# the one-time prefix prerun that produces it). Keyed by every stage that reads a previous stage's
# openPMD dump; a boundary not listed (cathode/gun) means nothing upstream is frozen.
FREEZE_SEED = {
    "injector":  ("logs/diags/gun",           "gun"),
    "linac1":    ("logs/diags/injector",      "injector"),
    "converter": ("logs/diags/linac1-4/sec4", "linac4"),
    "linac5-8":  ("logs/diags/converter",     "converter"),
}


def freeze_boundary(cfg):
    """The earliest pipeline stage any ACTIVE variable targets (= the first `--from` stage each eval
    runs); everything before it is frozen. None if a variable targets the very first stages
    (cathode/gun) -- then the chain runs from the top with nothing frozen."""
    stages = {CONFIG_STAGE[OVERRIDES[n][0]] for n in cfg["vocs"]["variables"]}
    if not stages:
        return None
    first = min(stages, key=STAGE_ORDER.index)
    return first if first in FREEZE_SEED else None


# ── Config ────────────────────────────────────────────────────────────────────────
def load_xopt_config(path=XOPT_CONFIG):
    """Load config/xopt.yaml, coercing variable bounds to float (PyYAML 1.1 leaves unsigned-exponent
    forms like `2.4e7` as strings, which VOCS would reject)."""
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    cfg["vocs"]["variables"] = {
        name: [float(lo), float(hi)] for name, (lo, hi) in cfg["vocs"]["variables"].items()}
    # scope=downstream runs only converter + linac5-8 (`--from converter`), so variables targeting
    # upstream stage configs (linac1/injector) never reach a running stage -- inert CNSGA dimensions
    # that also defeat sandbox dedup. Drop them here so the YAML can document all knobs and `full`
    # auto-includes them.
    if cfg["run"]["scope"] == "downstream":
        active = {"config/converter.yaml", "config/linac5-8.yaml"}
        inert = [n for n in cfg["vocs"]["variables"] if OVERRIDES[n][0] not in active]
        for n in inert:
            del cfg["vocs"]["variables"][n]
        if inert:
            print(f"scope=downstream: dropped {len(inert)} upstream variables ({', '.join(inert)})",
                  flush=True)
    return cfg


# ── Sandbox + override writing ──────────────────────────────────────────────────────
def eval_sandbox(inputs):
    """The ONE unique-per-eval sandbox dir, on NODE-LOCAL scratch. Each stage writes hundreds of
    openPMD diagnostic dumps; concurrent evals sharing /nfs saturate the NFS server (the whole job
    stalls in D/disk-sleep). Prefer $LINACSIM_RUNS_DIR, else $TMPDIR (SGE per-job, auto-cleaned), else
    /tmp -- all node-local. Deterministic in the inputs so a resumed run reuses it."""
    h = hashlib.sha1(repr(sorted(inputs.items())).encode()).hexdigest()[:16]
    base = os.environ.get("LINACSIM_RUNS_DIR") or os.environ.get("TMPDIR") or "/tmp"
    return os.path.join(base, "linac_runs", h)


def _fmt(value):
    """Unambiguous YAML scalar for an override value (repr keeps full float precision)."""
    return repr(float(value))


def _set_block_key(text, key, value, parent=None):
    """Set/insert a block scalar `  key: value` (keys are file-unique). If absent and `parent` is
    given, insert `key: value` indented two under the `parent:` line; else raise."""
    pat = re.compile(rf"^(\s*{re.escape(key)}:\s*)(\S+)(.*)$", re.M)
    if pat.search(text):
        return pat.sub(lambda m: m.group(1) + value + m.group(3), text, count=1)
    if parent is None:
        raise KeyError(f"{key} not found and no parent to insert under")
    ppat = re.compile(rf"^(\s*){re.escape(parent)}:\s*$", re.M)
    pm = ppat.search(text)
    if not pm:
        raise KeyError(f"parent {parent} not found to insert {key} under")
    eol = text.index("\n", pm.end()) + 1
    return text[:eol] + f"{pm.group(1)}  {key}: {value}\n" + text[eol:]


def _set_section_key(text, idx, key, value):
    """Set/insert `key: value` inside linac5-8.yaml's idx-th inline-flow `- {...}` section dict,
    preserving every other field and the trailing `# sec N` comment (so autophase_impact's
    section-aware crest regex still matches)."""
    lines = text.splitlines(keepends=True)
    sec_i = -1
    for li, line in enumerate(lines):
        if not re.match(r"^\s*- \{", line):
            continue
        sec_i += 1
        if sec_i != idx:
            continue
        if re.search(rf"\b{re.escape(key)}:", line):
            lines[li] = re.sub(rf"({re.escape(key)}:\s*)(\S+?)(\s*[}},])",
                               lambda m: m.group(1) + value + m.group(3), line, count=1)
        else:                                       # insert before the closing brace
            lines[li] = re.sub(r"\}", f", {key}: {value}}}", line, count=1)
        return "".join(lines)
    raise KeyError(f"section index {idx} not found")


def apply_overrides(out_dir, inputs):
    """Write each variable into the SANDBOX config copy (never the canonical config/). Text edits,
    not yaml.dump, so the files keep their structure (autophase rewrites the same copies in place,
    and linac5-8.yaml's inline-flow sections must survive for autophase_impact's regex)."""
    edits = {}                                      # path -> list of addressing tuples + value
    for name, value in inputs.items():
        spec = OVERRIDES.get(name)
        if spec is None:
            continue
        rel, addr = spec
        edits.setdefault(rel, []).append((addr, value))
    for rel, ops in edits.items():
        path = os.path.join(out_dir, rel)
        with open(path) as fh:
            text = fh.read()
        for addr, value in ops:
            v = _fmt(value)
            if addr[0] == "section":
                text = _set_section_key(text, addr[1], addr[2], v)
            else:                                   # ("block", key[, parent])
                text = _set_block_key(text, addr[1], v, parent=addr[2] if len(addr) > 2 else None)
        with open(path, "w") as fh:
            fh.write(text)


def seed_upstream(out_dir, from_stage):
    """Symlink the frozen upstream dump the `--from` boundary stage reads into the sandbox (read-only,
    identical across the population): converter reads linac1-4/sec4, injector reads gun, etc.
    from_stage is None => the chain runs from the top and regenerates everything (no-op)."""
    if from_stage is None:
        return
    rel = FREEZE_SEED[from_stage][0]
    src = os.path.join(REPO_ROOT, rel)
    if not os.path.isdir(src):
        raise FileNotFoundError(
            f"freeze boundary '{from_stage}' needs the frozen upstream dump at {src} -- run the chain "
            f"through it once (e.g. `python sim/main.py`) before optimizing.")
    dst = os.path.join(out_dir, rel)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if not (os.path.islink(dst) or os.path.exists(dst)):
        os.symlink(src, dst)


# ── Run the chain + read objectives ─────────────────────────────────────────────────
def run_chain(out_dir, from_stage, timeout_s):
    """Run the (partial) chain as a subprocess with LINACSIM_OUT_DIR=out_dir. `--from <boundary>`
    starts at the earliest varied stage off the seeded frozen dump (full -> injector, downstream ->
    converter); `--no-plots` drops the figures the optimizer never reads. cwd=out_dir mirrors how
    sim/main.py launches its own stages; the script path is absolutized against REPO_ROOT."""
    argv = [sys.executable, os.path.join(REPO_ROOT, "sim/main.py"), "--no-plots"]
    if from_stage is not None:
        argv += ["--from", from_stage]
    env = dict(os.environ)
    env["LINACSIM_OUT_DIR"] = out_dir
    env["PYTHONPATH"] = REPO_ROOT + os.pathsep + env.get("PYTHONPATH", "")
    env.setdefault("OMP_NUM_THREADS", "1")
    # start_new_session=True makes main.py its own process-group leader; on timeout kill the WHOLE
    # group. subprocess.run(timeout=) would SIGKILL only main.py and ORPHAN its WarpX/g4bl
    # grandchildren (reparented to init), which keep burning the node's cores long after the eval is
    # abandoned -- the failure mode that collapsed the first 128-worker run.
    proc = subprocess.Popen(argv, cwd=out_dir, env=env, start_new_session=True)
    try:
        rc = proc.wait(timeout=timeout_s)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait()
        raise
    if rc != 0:
        raise subprocess.CalledProcessError(rc, argv)


def read_summary(out_dir):
    """The (richer) linac5-8 injection_summary.json from this sandbox."""
    import json
    path = os.path.join(out_dir, "logs/diags/linac5-8/main/injection_summary.json")
    with open(path) as fh:
        return json.load(fh)


def _finite(x):
    return x is not None and isinstance(x, (int, float)) and math.isfinite(x)


def penalty_outputs():
    """CNSGA discards rows with NaN OBJECTIVES (not a flag); NaN every objective + constraint name."""
    return {k: float("nan") for k in OUTPUT_KEYS}


def evaluate(inputs):
    """One CNSGA evaluation: sandbox -> seed -> overrides -> run -> read objectives. Any failure
    (timeout, run error, OR a successful-but-empty beam) routes to penalty NaN so the worker survives."""
    cfg = load_xopt_config()
    from_stage = freeze_boundary(cfg)               # freeze every stage upstream of the first varied one
    timeout_s = cfg["run"].get("eval_timeout_s", 1800)
    keep_failed = cfg["run"].get("keep_failed_sandboxes", True)
    keep_success = cfg["run"].get("keep_successful_sandboxes", False)

    out_dir = eval_sandbox(inputs)
    make_out_dir(out_dir)                            # config/ copy + empty logs/ + fieldmaps symlink
    failed = False
    try:
        seed_upstream(out_dir, from_stage)
        apply_overrides(out_dir, inputs)            # into the COPY
        run_chain(out_dir, from_stage, timeout_s)
        m = read_summary(out_dir)
        if not (_finite(m.get("eps_n_x_m")) and _finite(m.get("eps_n_y_m"))
                and _finite(m.get("ke_out_mev"))):
            raise RuntimeError("empty beam / no survivors")     # NaN keys would propagate otherwise
        out = {k: m.get(k) for k in OUTPUT_KEYS if k != "eps_n"}
        out["eps_n"] = (m["eps_n_x_m"] * m["eps_n_y_m"]) ** 0.5
    except Exception:
        failed = True
        out = penalty_outputs()
    finally:
        # The summary is already read into `out`; the sandbox tree (a full chain's diags) is otherwise
        # dead weight. Prune unless retained for triage (failures) / inspection (successes).
        if (failed and not keep_failed) or (not failed and not keep_success):
            shutil.rmtree(out_dir, ignore_errors=True)
    return out


# ── Xopt wiring ─────────────────────────────────────────────────────────────────────
def make_executor(evcfg):
    """Build the Xopt Evaluator executor. `process` (DEFAULT, local) is a ProcessPoolExecutor sized
    to max_workers (else nproc-2 on the 16-core interactive node). `dask-sge` targets the CLASSE
    cluster -- it runs SGE/qsub (NOT HTCondor, NOT SLURM), so the cluster path is
    dask_jobqueue.SGECluster (the import is deferred into this branch so the process path needs no
    dask-jobqueue). SGECluster typically needs a site `queue`/`resource_spec`; pass them via the
    evaluator config when wiring the real CLASSE submit (cores=1 is one chain per slot)."""
    from concurrent.futures import ProcessPoolExecutor
    kind = evcfg.get("executor", "process")
    if kind == "process":
        n = evcfg.get("max_workers") or max(1, (os.cpu_count() or 2) - 2)
        return ProcessPoolExecutor(max_workers=n)
    if kind == "dask-sge":
        # Each dask worker is a separate 1-core qsub job that MUST reproduce submit.sge's environment:
        # activate the scratch conda env and re-set the env vars, else ~/.local shadows the env
        # (PYTHONNOUSERSITE) and WarpX/HDF5 misbehave. Site params come from the evaluator config.
        # `job_script_prologue` is dask-jobqueue >= 0.8 (older: `env_extra`).
        from dask_jobqueue import SGECluster
        from dask.distributed import Client
        # Prefer sourcing cluster_env.sh (the ONE site env: conda + OMP/HDF5/PYTHONNOUSERSITE + g4bl
        # PATH + Geant4 data vars) so workers match the controller exactly (the converter needs g4bl).
        env_setup = evcfg.get("env_setup")
        if env_setup:
            prologue = [f"source {env_setup}"]
        else:                                        # minimal fallback (NO g4bl/G4 data -> converter fails)
            conda_env = evcfg.get("conda_env")
            prologue = ["source ~/miniforge3/etc/profile.d/conda.sh",
                        f"conda activate {conda_env}" if conda_env else "",
                        "export OMP_NUM_THREADS=1",
                        "export HDF5_USE_FILE_LOCKING=FALSE",
                        "export PYTHONNOUSERSITE=1"]  # ~/.local lume/openpmd shadow the env otherwise
            prologue = [ln for ln in prologue if ln]
        cluster = SGECluster(
            cores=1, processes=1, memory=evcfg.get("memory", "4GB"),
            queue=evcfg.get("queue"),
            resource_spec=evcfg.get("resource_spec", "mem_free=4G"),
            walltime=evcfg.get("walltime", "48:00:00"),
            job_script_prologue=[ln for ln in prologue if ln])
        cluster.scale(jobs=int(evcfg.get("max_workers", 64)))
        return Client(cluster).get_executor()
    raise ValueError(f"unknown executor {kind!r} (process | dask-sge)")


def save_pareto(X, outdir):
    """Checkpoint the full data + the non-dominated FEASIBLE Pareto front to outdir. The front is
    computed over ALL configured objectives (q_out_C, eps_n, and any demoted beam-quality objectives),
    each normalised to a minimise sense via its MAXIMIZE/MINIMIZE direction, so it stays correct as the
    objective set changes."""
    os.makedirs(outdir, exist_ok=True)
    if X.data is None or len(X.data) == 0:
        return
    X.data.to_csv(os.path.join(outdir, "all_evaluations.csv"))
    from xopt.vocs import get_feasibility_data    # a module function in xopt 3.1.1, NOT a VOCS method
    try:
        feas = X.data[get_feasibility_data(X.vocs, X.data)["feasible"]]
    except Exception:
        feas = X.data
    objs = dict(X.vocs.objectives)                  # name -> Maximize/MinimizeObjective instance
    cols = list(objs)
    sub = feas.dropna(subset=cols)
    if len(sub) == 0:
        return
    # Direction from the objective TYPE name: str(MaximizeObjective()) is "dtype=None" (no "MAX"), so
    # key off type(...).__name__ ("MaximizeObjective"). MAXIMIZE -> sign -1 to recast as a minimise.
    signs = np.array([-1.0 if type(objs[c]).__name__.upper().startswith("MAX") else 1.0 for c in cols])
    pts = sub[cols].to_numpy() * signs              # now everything is minimise
    keep = []
    for i, a in enumerate(pts):                     # a is non-dominated if no b is <= in all, < in one
        dominated = any(np.all(b <= a) and np.any(b < a) for j, b in enumerate(pts) if j != i)
        if not dominated:
            keep.append(i)
    sub.iloc[keep].to_csv(os.path.join(outdir, "pareto_front.csv"))


def prebuild_fieldmaps():
    """Build the shared WarpX field maps ONCE before launching the population. The maps are
    parameter-independent for the optimized variables (gun voltage is not optimized), so every
    sandbox reads them via its fieldmaps symlink; the per-output skip-guard makes a re-run a no-op.
    Without this, the first concurrent wave of evals would race building the shared fieldmaps/h5.
    Only the upstream WarpX stages (gun/injector/linac1-4) build maps, so this is a no-op need for
    scope=downstream (converter is g4bl, linac5-8 reads the committed rfdata template)."""
    from sim.helpers.buildfields import (
        build_gun_field, build_injector_fields, build_linac_slac_fields)
    os.chdir(REPO_ROOT)                              # maps write to the shared REPO_ROOT/fieldmaps/h5
    gun_v = yaml.safe_load(open(os.path.join(REPO_ROOT, "config/gun.yaml")))["params"]["GUN_VOLTAGE"]
    build_gun_field(float(gun_v))
    build_injector_fields()
    build_linac_slac_fields()


def prefreeze_upstream(from_stage):
    """Run the frozen prefix (cathode..the stage before `from_stage`) ONCE at REPO_ROOT so every eval
    symlinks its dump instead of recomputing a parameter-independent prefix (scope=full freezes
    cathode+gun, ~16 min/eval). Idempotent: a no-op once the seed dump exists (resume / prior run)."""
    if from_stage is None:
        return
    rel, to_stage = FREEZE_SEED[from_stage]
    if os.path.isdir(os.path.join(REPO_ROOT, rel)):
        print(f"frozen upstream prefix present ({rel}) -- skipping prerun", flush=True)
        return
    print(f"freezing upstream prefix (cathode..{to_stage}) once -> {rel} ...", flush=True)
    env = dict(os.environ)
    env.pop("LINACSIM_OUT_DIR", None)               # write the shared REPO_ROOT dump, not a sandbox
    env["PYTHONPATH"] = REPO_ROOT + os.pathsep + env.get("PYTHONPATH", "")
    env.setdefault("OMP_NUM_THREADS", "1")
    subprocess.run([sys.executable, os.path.join(REPO_ROOT, "sim/main.py"),
                    "--to", to_stage, "--no-plots"], cwd=REPO_ROOT, env=env, check=True)


def main():
    from xopt import Xopt, VOCS, Evaluator
    from xopt.generators.ga.cnsga import CNSGAGenerator

    cfg = load_xopt_config()
    if cfg["run"]["scope"] == "full":               # upstream stages build the shared maps -> prebuild
        print("prebuilding shared field maps (fieldmaps/h5) ...", flush=True)
        prebuild_fieldmaps()
        prefreeze_upstream(freeze_boundary(cfg))    # run cathode..gun once; evals symlink the dump
    opt_dir = os.path.join(REPO_ROOT, "logs", "opt")
    os.makedirs(opt_dir, exist_ok=True)             # holds the per-step data.csv checkpoint + Pareto output

    vocs = VOCS(**cfg["vocs"])
    gen = CNSGAGenerator(vocs=vocs, population_size=cfg["generator"]["population_size"])
    ev = Evaluator(function=evaluate, executor=make_executor(cfg["evaluator"]),
                   max_workers=cfg["evaluator"].get("max_workers", 1))
    X = Xopt(generator=gen, evaluator=ev)           # Xopt 3.x derives vocs from the generator
    # Do NOT set X.dump_file: its per-step auto-dump serializes the CNSGA generator's population
    # DataFrame through pandas to_json(orient="columns"), which raises "index must be unique" under
    # this pandas/py3.14 once a generation rolls over. Checkpoint X.data to CSV each step instead.
    data_csv = os.path.join(opt_dir, "data.csv")

    max_eval = cfg["xopt"]["max_evaluations"]
    print(f"CNSGA: scope={cfg['run']['scope']}, population={cfg['generator']['population_size']}, "
          f"max_evaluations={max_eval}, {len(vocs.variables)} variables -> {opt_dir}/", flush=True)
    while (len(X.data) if X.data is not None else 0) < max_eval:
        X.step()                                    # CNSGA is async -- keeps the executor saturated
        if X.data is not None and len(X.data):
            X.data.reset_index(drop=True).to_csv(data_csv)   # per-step checkpoint (unique index)
    try:
        save_pareto(X, opt_dir)
    except Exception as e:
        print(f"save_pareto skipped ({type(e).__name__}); data.csv holds all evaluations", flush=True)
    print(f"Done. {len(X.data)} evaluations. Pareto + checkpoints in {opt_dir}/", flush=True)


if __name__ == "__main__":
    main()
