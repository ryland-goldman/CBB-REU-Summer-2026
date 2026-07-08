# Plan: `sim/optimize.py` + `config/xopt.yaml` — multi-objective beam optimization (Xopt/CNSGA)

Implementation plan for a multi-objective optimizer over the full simulation chain, maximizing beam
quality (brightness / low emittance / low energy spread / high energy / small spot / high charge)
using **Xopt** with the **CNSGA** generator. Builds directly on the per-eval isolation mechanism in
`docs/per-eval-isolation-plan.md` (the `LINACSIM_OUT_DIR` sandbox) — that is a **prerequisite**.

## Deliverables

- `config/xopt.yaml` — declarative VOCS (variables/objectives/constraints) + CNSGA + executor config.
- `sim/optimize.py` — loads `config/xopt.yaml`, builds the Xopt run, defines `evaluate()`
  (sandbox → apply overrides → run chain → read objectives from the summary), runs CNSGA,
  checkpoints the Pareto front.
- **Extend `sim/linac5-8.py`'s `injection_summary.json` writer** to emit the exit beam-quality
  metrics it currently omits (emittance, energy spread, spot size). The optimizer then just reads
  the summary — no separate computation — and the metrics are available for plotting/reporting too.
- **Add a per-section `crest_offset_deg` knob to `sim/linac5-8.py`** (and a `PHASE_OFFSET_DEG` to the
  section-1 path of `sim/linac1-4.py`) — a phase offset the driver **adds on top of the autophased
  crest** at apply time. This is mandatory: the optimizer's phase variables must target a key that
  autophase never rewrites (see "Crest-offset wiring" below — autophase owns `crest_phase_deg` /
  sec-1 `PHASE_DEG`, so the optimizer cannot write those directly).
- **Add a `--from <stage>` (partial-chain) flag to `sim/main.py`** to support `scope: downstream`
  (see "Execution").
- `xopt` added to `requirements.txt` (not in the CBB env today — `import xopt` fails).

## Prerequisite & dependency

- **Per-eval isolation must land first.** `evaluate()` relies on `LINACSIM_OUT_DIR` to give every
  population member its own `config/` copy + `logs/` tree (see the isolation plan). Without it,
  concurrent evals collide.
- **Install:** `conda install -c conda-forge xopt` (or `pip install xopt`); pin in
  `requirements.txt`. CNSGA lives at `xopt.generators.ga.cnsga.CNSGAGenerator`.

## What "a good beam" means here → objectives & constraints

CNSGA returns a **Pareto front**, but high-dimensional objective space degrades fast (a front over
6 objectives barely converges). So keep **2–3 true objectives** and push the rest to **constraints**
or fold into a derived **brightness** figure of merit. All evaluated at the **final exit**
(`linac5-8` particles dump):

**Recommended objectives (minimize form):**
1. **Maximize transmitted charge** `q_out_C` → objective `-q_out_C`.
2. **Minimize transverse emittance** `eps_n = sqrt(eps_n_x * eps_n_y)` [m·rad].
   (Equivalent single-axis brightness leverage; energy spread handled as a constraint.)

**Optional 3rd objective (brightness FoM, if a 3-axis front is wanted):**
- **Maximize 5-D brightness** `B = q_out_C / (eps_n_x * eps_n_y * sigma_E_rel)` → objective `-B`.
  (Captures charge, transverse, and energy-spread quality in one scalar; use *either* this as a
  single headline objective with charge+emittance as constraints, *or* the 2-objective set above —
  not both, to keep the front legible.)

**Constraints (Xopt `constraints:` — feasible when satisfied):**
- `ke_out_mev >= KE_FLOOR` (e.g. ≥ 200 MeV) — high energy.
- `sigma_E_rel <= DE_TOL` (e.g. ≤ 0.01) — low energy spread (σ_E/⟨E⟩).
- `sigma_x_mm <= SPOT_TOL` (e.g. ≤ 2 mm) — small beam size.
- `transmission_core >= TRANS_FLOOR` (linac5-8 `n_out/n_in`, already in the summary) — don't win on
  quality by throwing the beam away (the `q_out_C` objective already pushes absolute charge).

This is the standard accelerator pattern: optimize the genuine trade pair (charge ↔ emittance),
constrain the rest. The exact objective/constraint split is **declared in `config/xopt.yaml`** so it
can be retuned without touching code.

## Extend the final `injection_summary.json` (single source of truth)

**Gap to close:** `logs/diags/linac5-8/main/injection_summary.json` reports `ke_out_mev`,
`transmission_core`, `q_out_C`, `core_charge_frac`, `beta_min_core`, `n_out` — but **no emittance,
energy spread, or spot size at the exit** (emittance is computed only in `sim/plot/linac1-4.py` via
`metrics.screen_profile`).

**Decision (per the request): add these to the summary writer in `sim/linac5-8.py`, rather than
recomputing them in a separate optimizer-only extractor.** `linac5-8.py` already holds the final
beam in memory when it writes the summary, so this is a few lines, and the metrics become available
for plotting/reporting too — useful independent of the optimizer. New keys to emit:

| key | meaning | how |
|---|---|---|
| `eps_n_x_m`, `eps_n_y_m` | normalized transverse emittance [m·rad] | **`ParticleGroup.norm_emit_x/y`** (P_out is already a ParticleGroup — the centered, charge-weighted normalized emittance, no hand math). If computing manually, use the **centered** form `sqrt(var_q·var_u − cov²)` with mean-subtracted moments and `u=γβ` (matching `metrics.py:91-96`); the uncentered `<q²><u²>−<qu>²` is wrong for any off-axis/tilted beam |
| `sigma_E_mev`, `sigma_E_rel` | abs / relative energy spread | `KE=(sqrt(1+ux²+uy²+uz²)-1)·MC2`; `sigma_E_rel = std(KE)/mean(KE)` |
| `sigma_x_mm`, `sigma_y_mm` | RMS spot | `sqrt(<x²>-<x>²)` |

Implement once as a small helper (e.g. `metrics.beam_quality(pg)`) so the same numbers are reusable
by plotters; `linac5-8.py` calls it and merges the dict into the summary it already writes. The
charge & transmission objectives (`q_out_C`, `transmission_core`, `ke_out_mev`) are **already** in
the summary. The optimizer's `evaluate()` then just **reads the (now richer) summary JSON** — no
separate computation.

**No-survivor guard (required):** when the bunch fully scrapes, `sim/linac5-8.py:542-547` sets
`P_out = None` and `ke_out = None` and still writes a valid summary (no exception). So
`metrics.beam_quality` must be called as `beam_quality(P_out) if P_out is not None else {…NaN…}` and
emit `NaN` for every new key (and `ke_out_mev`). Otherwise `evaluate()`'s
`(m["eps_n_x_m"]*m["eps_n_y_m"])**0.5` raises `KeyError` on an empty beam — which (see `evaluate()`
below) must route to the penalty path, not kill the worker.

(Same pattern could later extend the converter / linac1-4 summaries for consistency; the converter
already reports `ke_pos_sigma_mev`, `sigma_r_pos_mm`, `div_pos_rms_mrad`. Out of scope here.)

## `config/xopt.yaml` (declarative VOCS + run config)

```yaml
# Xopt CNSGA configuration for the Cornell linac chain. sim/optimize.py reads this.
xopt:
  max_evaluations: 4000          # or drive by generations * population in code
generator:                       # CNSGAGenerator is bound in code; only population_size is read here
  population_size: 64            # match to available cores/slots (see isolation plan)
evaluator:
  executor: process              # process (local, default) | dask-sge (CLASSE cluster); see "Execution"
  max_workers: 64
vocs:
  variables:                     # name: [lo, hi]  — bounds in physical units
    # --- final-stage beam-quality knobs (highest leverage, cheapest if scope=downstream) ---
    l5_field_scale_0: [2.4e7, 2.9e7]   # linac5-8 sections[0].field_scale  (CEA 4)
    l5_field_scale_1: [2.4e7, 2.9e7]   # sections[1].field_scale           (CEA 5)
    l5_field_scale_2: [2.2e7, 2.7e7]   # sections[2].field_scale           (CU 3)
    l5_field_scale_3: [2.2e7, 2.7e7]   # sections[3].field_scale           (CU 4)
    l5_phase_off_0:   [-30, 30]        # crest OFFSET [deg] on sections[0] (added to autophased crest)
    l5_phase_off_1:   [-30, 30]
    l5_phase_off_2:   [-30, 30]
    l5_phase_off_3:   [-30, 30]
    l5_quad_scale_1:  [0.2, 2.5]       # sections[1].quad_scale (× the sec6 CESR exit-triplet gradients)
    l5_quad_scale_2:  [0.2, 2.5]       # sections[2].quad_scale (× the sec7 CESR exit-triplet gradients)
    # applied sec-5/6 field = this × (conv_sol_b_tesla/0.7022) via solenoid_tracking — range bounds the pre-scale value
    l5_sol_b_5:       [0.05, 0.6]      # sections[0].solenoid_b_tesla (264 A machine flat-top = 0.243 T)
    l5_sol_b_6:       [0.05, 0.6]      # sections[1].solenoid_b_tesla
    # --- converter (positron yield + emittance into linac5-8) ---
    conv_target_len_mm: [4.0, 9.0]     # converter geometry.target_length_mm
    conv_sol_b_tesla:   [0.35, 1.4]    # converter solenoid.b_tesla (con_sol peak; 0.7022 = 3300 A)
    conv_exit_drift_mm: [125, 250]     # converter solenoid.exit_drift_mm (back face -> plane; >121 = past the map)
    # --- capture (linac1) — include only when scope=full ---
    l1_power_mw:  [14, 24]             # linac1 POWER_MW
    l1_phase_off: [-30, 30]            # crest OFFSET [deg] on linac1 sec-1 (added to autophased PHASE_DEG)
    # --- injector match — include only when scope=full (expensive upstream) ---
    inj_i_sol0:   [20, 60]             # injector I_SOL0
    inj_i_lens0e: [0, 20]             # injector I_LENS0E
    inj_preb1_kw: [4, 14]             # injector PREB1_KW
    inj_preb2_kw: [4, 16]             # injector PREB2_KW
  objectives:
    q_out_C: MAXIMIZE
    eps_n: MINIMIZE
    sigma_E_rel: MINIMIZE   # demoted from constraint: the positron exit sits ~9% / ~8 mm, so a
    sigma_x_mm: MINIMIZE    # <1% / <2 mm hard cut admits no member (empty front). Trade on the front.
  constraints:
    ke_out_mev:        [GREATER_THAN, 200.0]
    transmission_core: [GREATER_THAN, 0.005]   # linac5-8 internal n_out/n_in (already emitted); the
                                               # q_out_C MAXIMIZE objective already pushes end-to-end charge
run:
  scope: downstream              # downstream (converter+linac5-8) | full (cathode→linac5-8)
  eval_timeout_s: 1800
  keep_failed_sandboxes: true
```

A `scope` switch keeps the variable set / chain length tractable: **start `downstream`** (re-run
only converter + linac5-8 off a frozen `linac1-4/sec4` electron dump — minutes per eval, the
final-beam knobs dominate quality), then widen to `full` once the front looks right.

## `sim/optimize.py` (structure)

```python
def eval_sandbox(inputs):
    """The ONE unique-per-eval dir helper (both plans reference this single name)."""
    return f"{REPO_ROOT}/logs/runs/{hashlib.sha1(repr(sorted(inputs.items())).encode()).hexdigest()[:16]}"

def evaluate(inputs):
    out_dir = eval_sandbox(inputs)                 # unique LINACSIM_OUT_DIR
    make_out_dir(out_dir)                           # isolation-plan helper: config/ copy + logs/ + fieldmaps symlink
    seed_upstream(out_dir, scope)                   # downstream scope only: copy/symlink the frozen sec4 dump in
    apply_overrides(out_dir, inputs)                # write vars into the COPY (see below)
    try:
        run_chain(out_dir, scope)                  # subprocess sim/main.py [--from converter], env=LINACSIM_OUT_DIR,
                                                   #   timeout=eval_timeout_s
        m = read_summary(out_dir)                  # the now-richer linac5-8 injection_summary.json
        if m.get("eps_n_x_m") is None or m.get("ke_out_mev") is None:
            raise RunFailed("empty beam")          # no survivors → penalty (NaN keys would KeyError otherwise)
        out = {**m, "eps_n": (m["eps_n_x_m"]*m["eps_n_y_m"])**0.5}
    except Exception:                              # timeout, run failure, OR successful-but-empty beam
        out = penalty_outputs()                    # NaN for every objective → CNSGA discards the row
    finally:
        if not keep_failed: shutil.rmtree(out_dir, ignore_errors=True)
    return out

def penalty_outputs():
    # CNSGA discards rows by NaN OBJECTIVES (not a flag). Return NaN for each objective/constraint name.
    return {"q_out_C": float("nan"), "eps_n": float("nan"), "ke_out_mev": float("nan"),
            "sigma_E_rel": float("nan"), "sigma_x_mm": float("nan"), "transmission_core": float("nan")}

def main():
    os.makedirs(f"{REPO_ROOT}/logs/opt", exist_ok=True)   # dump dir must exist before the CSV checkpoint
    cfg = yaml.safe_load(open("config/xopt.yaml"))
    vocs = VOCS(**cfg["vocs"])
    gen  = CNSGAGenerator(vocs=vocs, population_size=cfg["generator"]["population_size"])
    ev   = Evaluator(function=evaluate, executor=make_executor(cfg["evaluator"]))
    X = Xopt(vocs=vocs, generator=gen, evaluator=ev)      # X.dump_file left unset — see "Execution"
    while len(X.data) < cfg["xopt"]["max_evaluations"]:   # len(X.data) is the eval count (verify vs X.n_evaluations)
        X.step()                                   # CNSGA is async — keeps the executor saturated
        X.data.to_csv("logs/opt/xopt_data.csv")    # per-step progress record (no reload path — a restart starts fresh)
    save_pareto(X, "logs/opt/")
```

`penalty_outputs` encodes failure as **NaN objectives** (CNSGA's actual discard mechanism) — there is
no `xopt_error` return flag (Xopt sets that itself if the function *raises*; here we catch and return
NaN so the worker survives). Named helpers still to implement: `apply_overrides`, `seed_upstream`,
`run_chain`, `read_summary`, `make_executor`, `save_pareto` (specified in the sections around this).

**`save_pareto`'s direction handling:** `str(MaximizeObjective())` does **not** contain `"MAX"` (it
prints as `dtype=None`), so the minimize/maximize direction for each objective must be read off
`type(objective).__name__` rather than the object's string form.

### `apply_overrides(out_dir, inputs)` — how variables reach the sim (no canonical-YAML edit)

Each variable maps to a key in the **sandbox `config/` copy** (never the repo's canonical config —
same pattern autophase already uses to write crests). Concretely:

| variable | sandbox file · key |
|---|---|
| `l5_field_scale_{i}` | `config/linac5-8.yaml` · `sections[i].field_scale` |
| `l5_phase_off_{i}`   | `config/linac5-8.yaml` · `sections[i].crest_offset_deg` *(NEW key; driver adds it to the autophased `crest_phase_deg`. Do NOT write `crest_phase_deg` — autophase overwrites it. See note.)* |
| `l5_quad_scale_{i}`  | `config/linac5-8.yaml` · `sections[i].quad_scale` |
| `l5_sol_b_{5,6}`     | `config/linac5-8.yaml` · `sections[0|1].solenoid_b_tesla` |
| `conv_target_len_mm` | `config/converter.yaml` · `geometry.target_length_mm` |
| `conv_sol_b_tesla`   | `config/converter.yaml` · `solenoid.b_tesla` |
| `conv_exit_drift_mm` | `config/converter.yaml` · `solenoid.exit_drift_mm` |
| `l1_power_mw` | `config/linac1.yaml` · `params.POWER_MW` |
| `l1_phase_off` | `config/linac1.yaml` · `params.PHASE_OFFSET_DEG` *(NEW key; sec-1 driver adds it to the autophased `PHASE_DEG`. Do NOT write `PHASE_DEG` — `autophase.py` overwrites it for section 1. This key does not exist in the canonical `linac1.yaml`; `apply_overrides` inserts it under the sandbox copy's `params:` block on first use.)* |
| `inj_*` | `config/injector.yaml` · `params.I_SOL0` / `I_LENS0E` / `PREB1_KW` / `PREB2_KW` |

**Crest-offset vs autophase ordering (the critical integration point).** `sim/main.py` runs
autophase **before** each linac stage, and autophase **overwrites** the very keys the optimizer would
otherwise set: `sim/autophase_impact.py` rewrites `sections[i].crest_phase_deg`, and `sim/autophase.py`
rewrites section-1's `PHASE_DEG` (`PHASE_KEY[1]`). Since `evaluate()` runs `apply_overrides` *then*
launches `main.py` (which autophases), writing those keys directly would be silently discarded.
**Resolution (the deliverable above):** the optimizer writes phase **offsets** to NEW keys autophase
never touches — `sections[i].crest_offset_deg` (linac5-8) and `params.PHASE_OFFSET_DEG` (linac1
sec-1) — and the drivers apply `effective = autophased_crest + offset`. This keeps the crest physical
and the offset meaningful across the whole population.

Sections 2/3/4 need **no** new key: autophase writes `CREST_PHASE_DEG` while the detune knob is the
separate `PHASE_DEG` (`linac2/3/4.yaml`) — so if those phases are ever optimized they map straight to
`params.PHASE_DEG`. (The current VOCS does not optimize them.) `field_scale`/`FIELD_SCALE` are always
safe — autophase never touches amplitude.

## Partial chain (`scope: downstream`) — what it requires

`sim/main.py` today iterates **all** `STAGES` unconditionally (main.py:145) — there is no partial
run. `scope: downstream` (re-run only converter + linac5-8 off a frozen `linac1-4/sec4` electron
dump) therefore needs two pieces, both owned by this plan:

1. **`run_chain(out_dir, scope)`** → for `downstream`, invoke `sim/main.py --from converter`; for
   `full`, plain `sim/main.py`. The `--from <stage>` flag (deliverable above) slices `STAGES` to
   start at the named stage and skips autophase for skipped stages.
2. **`seed_upstream(out_dir, scope)`** → for `downstream`, the converter reads
   `logs/diags/linac1-4/sec4/main/particles`, which does **not** exist in a fresh sandbox
   (`make_out_dir` creates an empty `logs/`). So seed it: symlink (or copy) the frozen repo
   `logs/diags/linac1-4/sec4/` into `<out_dir>/logs/diags/linac1-4/sec4/` before running. The
   upstream dump is fixed across the population (only converter+downstream vars move), so a read-only
   symlink is safe. For `full`, `seed_upstream` is a no-op.

## Freeze boundary — running upstream stages once per population

Every stage upstream of the earliest actively-varied stage is parameter-independent across the
whole population (none of its inputs move), so it is run **once** (`prefreeze_upstream`) and each
eval symlinks the frozen dump instead of recomputing it:

- **`scope: full`** varies from the injector onward, so cathode+gun are prefrozen — avoiding their
  combined ~16 min/eval on every population member.
- **`scope: downstream`** varies only converter/linac5-8 configs; `prefreeze_upstream` also drops
  (as inert) any VOCS variable that targets an upstream config, since those configs never reach a
  running stage under `--from converter` anyway.
- The shared WarpX field maps (`fieldmaps/h5`) are likewise parameter-independent for every
  currently-optimized variable (e.g. gun voltage is not a VOCS knob), so `sim/optimize.py` builds
  them once before the population starts — otherwise the first wave of concurrent evals would race
  to build the same maps.

## Execution (local → cluster)

- **Local smoke test:** `executor: process`, small `population_size` (8), `scope: downstream`,
  `max_evaluations: ~40` — confirms the loop, sandboxing, and objective extraction end to end.
- **Cluster:** `executor: dask-sge` → a `dask_jobqueue.SGECluster` client; one worker = one slot
  running a sandboxed chain. **CLASSE runs SGE/qsub, not HTCondor or SLURM** (lnx201.classe.cornell.edu;
  submit with `qsub`/`qstat`). `SGECluster` typically needs a site `queue`/`resource_spec` — pass them
  via the evaluator config. The `SGECluster` import lives inside the `dask-sge` branch of
  `make_executor` so the default `process` path has no dask dependency. CNSGA is async, so workers stay saturated (no per-generation
  barrier). Wall-clock ≈ generations × slowest eval; total ≈ population × generations evals (see the
  feasibility discussion). The manual `X.data` CSV dump (see below) is a per-step checkpoint-of-record:
  it preserves every evaluation for analysis and potential manual re-seeding, but `main()` has no
  reload path (and the CSV omits the CNSGA population state), so a restarted run re-optimizes from
  scratch.
- `run_chain` launches the chain with `start_new_session=True` so an eval timeout can `SIGKILL` the
  whole process group (`main.py` plus its WarpX/g4bl grandchildren) — killing only the `main.py` PID
  orphans those grandchildren (reparented to init), which then keep burning cluster cores
  indefinitely.
- **`X.dump_file` is intentionally left unset**, despite the pseudocode above: Xopt's per-step
  auto-dump serializes the CNSGA generator's population `DataFrame` via pandas
  `to_json(orient="columns")`, which raises `"index must be unique"` under this project's
  pandas/py3.14 combination once a generation rolls over. `X.data` is checkpointed to CSV manually
  each step instead.
- `RNG_SEED=0` / `OMP_NUM_THREADS=1` stay per eval (deterministic, reproducible, restart-safe).

## Validation

1. **Summary-metrics check:** re-run `sim/linac5-8.py` on the existing upstream dump; confirm the
   new `eps_n_x_m/eps_n_y_m`, `sigma_E_rel`, `sigma_x_mm` keys appear in the summary, are positive
   and finite, and `sigma_E_mev/ke_out_mev` is physically sane.
2. **Single `evaluate()`:** call with the baseline (current config values) under a temp
   `LINACSIM_OUT_DIR`; confirm it returns finite objectives and leaves the canonical `config/`
   unchanged (`git diff` clean).
3. **Failure path:** force a diverging input (e.g. `conv_target_len_mm` huge); confirm the timeout →
   penalty path returns rather than hanging the population.
4. **Tiny CNSGA run:** 8×5 downstream; confirm a Pareto front (charge ↔ emittance) and checkpoint
   files under `logs/opt/`.

## Out of scope / follow-ups

- The isolation mechanism itself (separate plan — prerequisite).
- Surrogate/Bayesian acceleration (Xopt also offers MOBO; could warm-start CNSGA later).
- Extending the converter / linac1-4 summaries with the same beam-quality keys (only `linac5-8` is
  in scope here).

## Reference: exact knobs (grounding the VOCS)

| stage | file | key(s) | role |
|---|---|---|---|
| linac5-8 | `config/linac5-8.yaml` | `sections[i].{field_scale, crest_phase_deg, quad_scale, solenoid_b_tesla}`, `rf.phase_deg` | final accel + capture optics |
| converter | `config/converter.yaml` | `geometry.target_length_mm/target_radius_mm`, `solenoid.b_tesla`, `solenoid.exit_drift_mm` | e⁺ yield / emittance / size |
| linac1 | `config/linac1.yaml` | `params.POWER_MW`, `params.PHASE_DEG` | capture into the relativistic chain |
| linac2-4 | `config/linac{2,3,4}.yaml` | `params.FIELD_SCALE`, `params.PHASE_DEG` (detune-from-crest) | back-half energy/phase |
| injector | `config/injector.yaml` | `params.{I_SOL0,I_LENS0A..E,PREB1_KW,PREB2_KW,PREB1_PHI_OFF,PREB2_PHI_OFF}` | bunching + transverse match |
| gun/cathode | `config/{gun,cathode}.yaml` | `params.GUN_VOLTAGE`, cathode emission | low priority (fixed operating point) |

Objective sources: **all** read from `logs/diags/linac5-8/main/injection_summary.json` —
`q_out_C`, `transmission_core`, `ke_out_mev` already present; `eps_n_x_m/eps_n_y_m`, `sigma_E_rel`,
`sigma_x_mm` added by extending `sim/linac5-8.py`'s summary writer (`metrics.beam_quality`).
