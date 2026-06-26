# Plan: per-eval isolation via the `LINACSIM_OUT_DIR` environment variable

Implementation plan for isolating **configs and outputs per run** so an Xopt/CNSGA population can
evaluate many full-chain runs concurrently without racing on shared files. Isolation is driven by an
**environment variable**, **not** by editing any `config/*.yaml` in place.

## Goal & constraints

- One evaluation = one full chain (`cathode → … → linac5-8`) over a parameter set.
- N evaluations run **concurrently** (cluster). They must not collide on: the canonical
  `config/*.yaml`, the `logs/diags/**` outputs, `injection_summary.json`, WarpX `warpx_used_inputs`,
  Impact-T `fort.*`.
- **Hard constraint (from the request):** isolation is achieved via an env var, never by mutating
  the checked-in YAML (e.g. do **not** rewrite `io.outdir` / `io.workdir` in the YAML to redirect a
  run). The canonical `config/*.yaml` stays byte-identical across all evals.

## The contract

```
LINACSIM_OUT_DIR=<DIR>     # if set: read config from <DIR>/config, write output to <DIR>/logs
                           # unset: use <repo>/config and <repo>/logs (today's behavior, unchanged)
```

- `fieldmaps/` is **not** redirected — it stays shared at the repo root (committed inputs + built
  maps are identical across evals; no reason to copy them per run).
- Unset ⇒ resolves to `REPO_ROOT` ⇒ **fully backward-compatible**: a plain `python sim/main.py` or
  `python sim/linac1-4.py 2` behaves exactly as today.
- An env var (vs a CLI flag) is the right mechanism here: `sim/main.py` spawns every stage and
  autophase step as a subprocess, and env vars **inherit automatically** — the optimizer sets
  `LINACSIM_OUT_DIR` once and every stage of that eval picks it up, with **no argv threading** and
  **no per-stage flag parsing** (which would otherwise collide with `linac1-4.py`'s positional
  section number and `converter.py`'s positional `n_events`).
- One tradeoff to note: env vars are "sticky" — a stray exported `LINACSIM_OUT_DIR` silently
  redirects a manual run. The `REPO_ROOT` default keeps normal runs safe; document the var.

## Core idea: redirect the working directory, not the path literals

`prepare_env()` (`sim/helpers/tools.py:58`) ends with `os.chdir(REPO_ROOT)`, and every stage reads
config from `"config/<stage>.yaml"` and writes to `"logs/diags/<stage>/…"` — all **relative**
literals resolved against the cwd (see the inventory at the end). Therefore:

> Point `prepare_env()`'s chdir at `LINACSIM_OUT_DIR` (default `REPO_ROOT`). Every `config/…` and
> `logs/…` literal then resolves under that directory — with **zero edits to the hundreds of
> hardcoded path strings**. Keep `fieldmaps/` shared by symlinking `<DIR>/fieldmaps → <repo>/fieldmaps`
> so the relative `fieldmaps/…` literals still hit the shared maps.

The code still lives in the repo (imported via `sys.path`); only the **cwd** (config/output) moves.
These two roles — *where the code is* (`REPO_ROOT`, on `sys.path`) vs *where the run reads/writes*
(`LINACSIM_OUT_DIR`, the cwd) — must stay separate.

This is preferred over rewriting every literal to `f"{OUT_DIR}/config/…"` / `f"{OUT_DIR}/logs/…"`:
those would touch ~12 path constants in `linac1-4.py` alone, plus the WarpX `path=` kwargs, the
Impact-T `cfg["io"]` reads, and every `loadparticles` upstream-read constant. The chdir approach
changes ~3 functions total. With cwd at `REPO_ROOT` instead, `fieldmaps/…` resolves for free but
`config`/`logs` would need the mass rewrite — so chdir-into-`OUT_DIR` + one fieldmaps symlink is the
strictly smaller change.

## Code changes (small, enumerated)

### 1. `sim/helpers/tools.py` — resolver + chdir target

```python
def out_root():
    """This run's config/output root: $LINACSIM_OUT_DIR if set, else REPO_ROOT. Always absolute."""
    return os.path.abspath(os.environ.get("LINACSIM_OUT_DIR", REPO_ROOT))
```

Change `prepare_env()` so the **chdir target is `out_root()`** while **`sys.path` stays
`REPO_ROOT`** (imports must always find the code):

```python
    root = out_root()
    if os.getcwd() != root:
        os.chdir(root)
    if REPO_ROOT not in sys.path:          # unchanged: code lives in the repo, not the sandbox
        sys.path.insert(0, REPO_ROOT)
```

`prepare_env()` must run **before** any stage touches a relative path — it already does (every
driver calls it before `WarpX(...)` / `load_config`). No other stage code changes.

### 2. New helper `sim/helpers/sandbox.py` — build the sandbox

```python
def _link(target, link):                           # idempotent symlink
    if os.path.islink(link) or os.path.exists(link):
        return
    os.symlink(target, link)

def make_out_dir(out_dir, src_root=REPO_ROOT):
    """Populate a LINACSIM_OUT_DIR sandbox: own config/ copy + empty logs/, fieldmaps shared.
    NO-OP when out_dir == src_root (the unset / REPO_ROOT case — never sandbox the repo itself)."""
    if os.path.abspath(out_dir) == os.path.abspath(src_root):
        return                                     # plain run: leave repo config/logs/fieldmaps as-is
    os.makedirs(out_dir, exist_ok=True)
    shutil.copytree(f"{src_root}/config", f"{out_dir}/config", dirs_exist_ok=True)  # isolated copy
    os.makedirs(f"{out_dir}/logs", exist_ok=True)
    _link(f"{src_root}/fieldmaps", f"{out_dir}/fieldmaps")     # shared maps, NOT redirected
```

`make_out_dir(REPO_ROOT)` is an explicit **no-op** so a plain `python sim/main.py` (var unset →
`out_root()==REPO_ROOT`) neither re-copies config nor tries `_link(repo/fieldmaps, repo/fieldmaps)`
(which would otherwise be a self-link / `FileExistsError`). Backward compatibility is preserved.

- The `config/` **copy** is what makes isolation env-driven and YAML-edit-free: autophase rewrites
  the **sandbox copy's** crest (next section), never the canonical `config/`.
- `logs/` starts empty so every stage's `logs/diags/**` and `injection_summary.json` are private.
- The single `fieldmaps` symlink realizes "fieldmaps stays at repo root" while the cwd is the
  sandbox.

### 3. Autophase scripts — already compatible, just honor the cwd

`sim/autophase.py` (`CONFIG = {N: f"config/linac{N}.yaml"}`, `set_yaml_param` rewrites it) and
`sim/autophase_impact.py` (`CONFIG = "config/linac5-8.yaml"`) open **relative** `config/...`. Once
`prepare_env()` chdirs into `LINACSIM_OUT_DIR`, those reads/writes hit the **sandbox copy**
automatically. The only change: ensure they call `prepare_env()` (or `os.chdir(out_root())`) at
startup. This keeps the crest-then-offset design intact — autophase writes the crest into the
sandbox config; Xopt's offset/field-scale overrides apply on top per the optimization plan (out of
scope here).

### 4. `sim/main.py` — set the var once; subprocesses inherit; **redirect main.py itself**

- Resolve `out_dir = out_root()` (i.e. `LINACSIM_OUT_DIR` or `REPO_ROOT`). **Do not** invent a
  `logs/runs/<id>` when the var is unset — a plain run must stay at `REPO_ROOT` (earlier drafts'
  "or generate" branch is wrong; it would redirect normal runs and break backward-compat).
- Call `make_out_dir(out_dir)` before the first stage **only to build a not-yet-populated sandbox**:
  guard it as `if out_dir != REPO_ROOT and not os.path.isdir(f"{out_dir}/config"): make_out_dir(out_dir)`.
  This matters because the optimizer's `evaluate()` **pre-builds the sandbox and writes overrides
  into `<out_dir>/config` before launching `main.py`** — an unconditional `make_out_dir` here would
  `copytree` the canonical `config/` back over the sandbox and **wipe the optimizer's overrides**
  (`copytree(..., dirs_exist_ok=True)` overwrites file-by-file). The `config/` existence check makes
  `main.py` build the sandbox for a standalone `LINACSIM_OUT_DIR=… python sim/main.py`, but defer to
  a caller that already prepared one. (The `make_out_dir` no-op on `REPO_ROOT` is a second guard.)
- In the subprocess builder (`subprocess.run([sys.executable, *argv], cwd=REPO_ROOT, env=env)`,
  ~line 78): set `env["LINACSIM_OUT_DIR"] = out_dir`, change `cwd=REPO_ROOT` → `cwd=out_dir`, **and
  absolutize the script path against `REPO_ROOT`**. `STAGES` stores *relative* script paths
  (`"sim/cathode.py"`, `["sim/autophase.py","1"]`, `"sim/plot/*.py"` — main.py:33-43); Python resolves
  a command-line script path against **cwd**, not `PYTHONPATH`, so `cwd=out_dir` alone makes every
  subprocess look for `<out_dir>/sim/cathode.py` (absent in the sandbox) and fail before any import.
  Build the command as `cmd = [sys.executable, os.path.join(REPO_ROOT, argv[0]), *argv[1:]]` (argv[0]
  is the script, rest are its args), keep `cwd=out_dir`, keep `PYTHONPATH=REPO_ROOT` (main.py:73).
  Then the absolute script path loads, `__file__`-based `sys.path` inserts still point at `REPO_ROOT`,
  and relative `config/`/`logs/` I/O resolves in the sandbox. `run_subprocess` builds
  `env = dict(os.environ)` (main.py:70), so the env line propagates to every stage and autophase
  subprocess. The `cwd=out_dir` change is **required**, not cosmetic: two entry points touch relative
  paths **before** `prepare_env()` runs, so relying on the in-process chdir is insufficient —
  - `sim/injector.py:37` resolves `GUN_DIAG` with an `os.path.isdir("logs/diags/gun/handoff")` check
    at **module import**, before `main()`/`prepare_env()`. With `cwd=REPO_ROOT` it would test the
    shared repo and pick the wrong upstream beam; with `cwd=out_dir` it tests the sandbox.
  - The plotters `sim/plot/{cathode,converter,linac1-4,linac5-8}.py` never call `prepare_env()`/chdir
    at all — they read `logs/diags/…` and `savefig` to `logs/plots/…` relative to cwd. With
    `cwd=out_dir` their I/O lands in the sandbox; otherwise they read/write the shared repo and race.
  `PYTHONPATH=REPO_ROOT` is already set (main.py:73) so imports still resolve with cwd moved. The
  in-process `prepare_env()` chdir stays (it covers standalone stage runs), now redundantly aligned.
- **`sim/main.py` does NOT go through `prepare_env()`** — it hardcodes `os.chdir(REPO_ROOT)` at
  **main.py:134**, and its own outputs are written relative to cwd: the run-log
  `logs/pipeline/log_<ts>.log` (main.py:136) and `beam_summary()` reading `logs/diags/…`
  (main.py:154-157). Under a sandboxed run these would escape to the **shared** repo `logs/` —
  concurrent evals collide on the pipeline log (timestamped only to the second) and the end-of-run
  summary reads foreign/stale data. **Fix:** replace the line-134 `os.chdir(REPO_ROOT)` with
  `os.chdir(out_root())` (import `out_root` from `helpers.tools`), placed after `make_out_dir`, so
  main.py's own log + summary resolve into the sandbox. (`sys.path` already carries `REPO_ROOT` via
  the stages' `__file__`-based insert + `PYTHONPATH`, so imports still resolve.)

## Field maps (shared ⇒ prebuild AND add a build skip-guard)

`buildfields.py` writes built maps to the relative `fieldmaps/h5/`, which now resolves through the
symlink to the **shared** repo `fieldmaps/h5/`. Two problems the drivers create, both **must** be
fixed or concurrent evals corrupt each other's maps:

1. **No idempotence — every run truncates the maps.** `build_gun_field` / `build_injector_fields` /
   `build_linac_slac_fields` have no existence guard: `write_thetamode_series` always opens with
   `io.Access.create` (recreate/truncate), and the drivers call them **unconditionally every run**
   (`gun.py:189`, `injector.py:115`, `linac1-4.py:134`). So even after a one-time prebuild, each
   concurrent eval **rewrites** the shared `fieldmaps/h5/*.h5` while siblings read → partial-read /
   create-race corruption. Prebuilding alone does **not** make them read-only. **Fix (required), two
   parts:**
   - **Per-output-file skip-guard.** A single `os.path.exists()` per `build_*` is insufficient:
     `build_injector_fields` emits **8** maps (`preb1_EB`, `preb2_EB`, 6 solenoids) and
     `build_linac_slac_fields` emits 2. Guard **each** output file (skip the individual
     `write_thetamode_series` when its target `.h5` already exists), or gate the whole `build_*` on a
     sentinel that is written last and implies all of its maps are complete.
   - **An explicit one-time prebuild BEFORE the population launches.** A skip-guard alone does not
     prevent the *first wave*: if N evals start before any `h5` exists, all see `exists()==False` and
     build concurrently into the shared dir. So the optimizer must build the maps once up front
     (e.g. a single `python sim/main.py` on the baseline, or a dedicated `python -c "buildfields…"`
     step in `sim/optimize.py`'s startup, before submitting any eval). After that, every eval's guard
     short-circuits to pure reads.
   The alternative to both is a **per-`OUT_DIR` real `h5/`** (symlink only `gdf/` + `rfdata/`, make
   `h5/` a real dir) — no shared state, no guard needed, at the cost of rebuilding maps per eval.
2. **The gun map is parameter-DEPENDENT.** `gun_E.h5` is scaled by `-GUN_VOLTAGE/GUN_MAP_VOLTAGE` at
   build time (`buildfields.py:113`, `gun.py:185-189`) — unlike the RF/solenoid maps, which are
   normalized (per-kW / per-Ampere) with amplitudes applied at **runtime**, so those are genuinely
   parameter-independent and safe to prebuild-and-share. Consequence: **if `GUN_VOLTAGE` is ever an
   Xopt variable, the gun must use a private per-`OUT_DIR` `fieldmaps/h5/`** — a shared `gun_E.h5`
   would cross-contaminate evals, and a skip-guard would be *worse* (it freezes the first writer's
   voltage for all). With the optimization plan's current variable set `GUN_VOLTAGE` is **not**
   optimized, so prebuild-and-share + skip-guard is safe today; revisit the moment it is added.

## Edge cases / things to verify

- **`io.outdir` / `io.workdir` must be relative.** `linac5-8.py:484` and `converter.py:175` read
  these from YAML. If relative (e.g. `logs/diags/linac5-8/main`), the chdir redirects them for free.
  **Audit:** confirm none is absolute in the committed YAML; if any is, make it relative (a config
  fix — the value stays identical for every eval, so this is not isolation-by-edit).
- **`sys.path` vs `cwd`.** Keep `REPO_ROOT` on `sys.path` (imports) while cwd is `LINACSIM_OUT_DIR`.
  Do not collapse the two.
- **`warpx_used_inputs`** and Impact-T `fort.*` write to cwd ⇒ land in the sandbox ⇒ isolated.
- **Robustness for the fleet:** a diverged eval must still return a value, and its sandbox should be
  removable. Have the optimizer wrapper `make_out_dir` → run → parse → optionally
  `shutil.rmtree(out_dir)` (keep on failure for triage). Per-eval timeout + penalty is a separate
  optimizer concern.
- **HDF5 file locking** already disabled (`HDF5_USE_FILE_LOCKING=FALSE`); separate sandboxes remove
  the remaining cross-eval contention on `logs/`.

## Xopt `evaluate()` integration (how it gets used)

```python
def evaluate(inputs):
    out_dir = eval_sandbox(inputs)                              # unique per eval (see optimization-plan)
    make_out_dir(out_dir)
    apply_overrides(out_dir, inputs)          # write offsets/field_scales into the COPY (opt. plan)
    subprocess.run([sys.executable, "sim/main.py"],
                   env={**os.environ, "LINACSIM_OUT_DIR": out_dir}, timeout=...)
    return parse_objectives(f"{out_dir}/logs/diags")             # the isolated summaries
```

Concurrency falls out of distinct `LINACSIM_OUT_DIR` values — no shared mutable state between
population members (fieldmaps is read-only-shared after prebuild).

## Validation

1. **Backward-compat:** `python sim/main.py` (var unset) reproduces today's `logs/diags/**`
   byte-for-byte (determinism: `RNG_SEED=0`, `OMP=1`).
2. **Single sandboxed run:** `LINACSIM_OUT_DIR=/tmp/run0 python sim/main.py` writes nothing under the
   repo `logs/` or `config/`; all output under `/tmp/run0/logs`; canonical `config/*.yaml` unchanged
   (`git diff` clean).
3. **Concurrency:** launch two sandboxed full chains at once with different params; confirm distinct
   outputs, no file-lock errors, each `injection_summary.json` matches its own params.
4. **Autophase isolation:** confirm a sandboxed run's autophase rewrites only
   `<OUT_DIR>/config/linac*.yaml`, never the repo copy.

## Out of scope (separate plans)

- The Xopt VOCS, crest-**offset** + `field_scale` variable wiring, objectives/constraints.
- Per-eval timeout/penalty + diverged-run handling.
- The SGE/Dask executor and CNSGA generator config.

## Reference: current path resolution (why chdir suffices)

| Concern | Where set | Form |
|---|---|---|
| cwd / chdir | `tools.py:72` `prepare_env` | `os.chdir(REPO_ROOT)` ← **the one redirect point → `out_root()`** |
| config (WarpX) | `cathode.py:40`, `gun.py:36`, `injector.py:34`, `linac1-4.py:137` | relative `"config/<stage>.yaml"` |
| config (Impact/conv) | `linac5-8.py:41`, `converter.py:30` | relative, `yaml.safe_load` |
| config (autophase) | `autophase.py:43`, `autophase_impact.py:51` | relative; rewritten in place ⇒ rewrites the copy |
| WarpX output | `…path=DIAG_DIR/"logs/diags/…"` (cathode:41, gun:39, injector:35, linac1-4:138) | relative |
| Impact/conv output | `cfg["io"]["outdir"]/["workdir"]` (linac5-8:484, converter:175) | relative (verify) |
| upstream reads | `loadparticles` + module consts (linac1-4.py:50-57, etc.) | relative `"logs/diags/…"` |
| field maps | `buildfields.py:22-23` | relative `"fieldmaps/{gdf,h5}"` ⇒ **shared via symlink** |
| subprocess | `main.py:78` | `env=…` ⇒ inject `LINACSIM_OUT_DIR` (inherited by all stages) |

`config/…` and `logs/…` are relative ⇒ one chdir into `LINACSIM_OUT_DIR` isolates them;
`fieldmaps/…` is kept shared by the symlink.
