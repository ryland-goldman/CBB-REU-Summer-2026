---
name: physics-review
description: Full physics review of the CHESS beam-simulation chain. Spawns an agent team of accelerator-physics, WarpX/PIC, and Python experts to find bugs, quality, and performance issues and cross-examine each other, then runs a workflow to confirm/fix/review, commits, and gives a final merge recommendation.
disable-model-invocation: true
argument-hint: [pr-number]
allowed-tools: Bash(git *), Bash(gh *), Bash(conda *), Task, TaskCreate, TaskUpdate, TaskGet, TaskList, SendMessage, TeamCreate, TeamDelete, Agent, Workflow, Read, Write, Edit, Grep, Glob
---

# Full Physics Review

End-to-end review of the front-end Linac chain. Drive all four phases below in order. Do **not** skip the cross-examination or the merge gate.

## Review target (auto-detected)

- Requested PR override (arg): `$ARGUMENTS`
- Current branch: !`git rev-parse --abbrev-ref HEAD`
- Active PR for branch: !`gh pr view --json number,title,state,headRefName,baseRefName 2>/dev/null || echo "none"`
- Open PRs in repo: !`gh pr list --state open --json number,title,headRefName 2>/dev/null || echo "none"`
- Uncommitted changes: !`git status --short`
- Diff stat vs main: !`git diff main...HEAD --stat 2>/dev/null | tail -40`

**Scope rule:**
- If `$ARGUMENTS` names a PR number, review that PR (`gh pr checkout <n>` first if not already on it).
- Else if an active PR exists for the current branch, review **that PR's diff** (`git diff <baseRefName>...HEAD`).
- Else review the **full codebase** (all four stages).
State the chosen scope in one line before proceeding.

## Prerequisites

- Agent teams must be enabled (`CLAUDE_CODE_EXPERIMENTAL_AGENT_TEAMS=1`). If team creation fails with a "not enabled" error, stop and tell the user to add it to `~/.claude/settings.json` `env` and restart.
- Read `CLAUDE.md` and the `README.md` of every in-scope stage (plus immediate neighbors) before spawning the team, per the repo's read-aggressively rule. Teammates load `CLAUDE.md` automatically but **not** your conversation history — put stage-specific context in each spawn prompt.
- Only one team can exist at a time; if a prior team is still running, clean it up first.

---

## Phase 1 — Agent-team review with cross-examination

Create an agent team and spawn these teammates. Give each a spawn prompt naming the exact stage(s)/files in scope, the relevant `README.md`, and the WarpX gotchas from `CLAUDE.md`. Tell each to record findings as `file:line` with a one-line rationale.

| Teammate | Name it | Domain focus |
|----------|---------|--------------|
| Accelerator physicist | `accel` | probe the overall start-to-end integration of the simulation and verify the accuracy of output figures |
| Emission & space-charge physicist | `emission` | cathode — SCL diode, Child–Langmuir current, PPC/flux injection, space-charge at low energy |
| RF & bunching physicist | `rf` | prebuncher velocity bunching (power/phase, ~1 mm bunch) + linac TW section energy gain (~37 MeV), RF field-map scaling/summing |
| Beam-optics & transport physicist | `optics` | electrostatics (~148 keV), solenoid focusing, emittance/envelope, **inter-stage handoff contract** (charge renorm, openPMD beam I/O) |
| WarpX / PIC expert | `warpx` | geometry binding (2D vs RZ), MLMG/Poisson solver convergence, CFL & step counts, ≈3:1 cell-aspect requirement, thetaMode openPMD axis order, external-field load conventions, `LoadAppliedField` ordering |
| Python / scientific-software engineer | `pyeng` | bugs, unit/constant correctness, dead code, API/`config()` consistency, performance knobs, doc drift vs `CLAUDE.md`/READMEs |

Add a numerical-methods, diagnostics, or other specialists only if scope clearly warrants it.

Run it in two rounds:

1. **Independent review.** Each teammate reviews its domain and appends findings to a shared findings doc at `pipeline/logs/physics_review_${CLAUDE_SESSION_ID}_{$AGENT_NAME}_findings.md`
2. **Cross-examination (mandatory).** Have teammates message each other to challenge findings like a scientific debate — each tries to *disprove* the others' claims and flag anything missed at a stage boundary (physics that looks fine in one stage but breaks the downstream contract). A finding survives only if it withstands challenge. Resolve duplicates and downgrade/withdraw refuted items.

Then synthesize a single **deduplicated, severity-ranked issue list** from the findings docs. Each issue: id, title, `file:line`, severity, why-it's-real (post-debate), suggested fix direction. Shut down teammates and clean up the team before Phase 2.

---

## Phase 2 — Confirm / fix / review workflow

Launch a **dynamic `Workflow`** over the surviving issue list (this skill is your opt-in to call `Workflow`). Recommended shape — adapt to the issue count:

- **Confirm** (parallel, read-only): one skeptic agent per issue re-verifies it against the current tree; drop any it cannot reproduce. `log()` anything dropped.
- **Fix** (per confirmed issue): implement the change. If fixes touch **disjoint** files, run them in parallel with `isolation: 'worktree'`; if any two fixes share a file, serialize those. Behavior-changing physics fixes must note what re-validation is needed (energy gain, Child–Langmuir current, bunching).
- **Review** (per fix): an independent reviewer checks correctness, that it didn't break the inter-stage contract or the ≈3:1 aspect/MLMG constraints, and that no docs drifted. Reject → loop back to Fix once.

Use `pipeline()` so each issue flows confirm→fix→review independently. Apply verified fixes to the working tree (consolidate from worktrees if used). Keep `pipeline/`, READMEs, `CLAUDE.md`, `FIGURES.md`, and `requirements.txt` in sync per the repo's doc-sync rule — a fix isn't done until its docs match.

---

## Phase 3 — Commit

Commit the verified fixes. If on `main`, create a branch first (`physics-review-fixes-<short-desc>`). Follow the repo's commit convention (stage `*.py` + `README.md`; `git add -f <stage>/results/*.png` only if figures were regenerated; never commit `diags/`, `.h5`, or logs). Group logically; clear messages.

---

## Phase 4 — Final review & merge gate

Run **one** `Agent` (general-purpose) as a fresh reviewer over the committed diff. Give it the original issue list and the commits. It must return an explicit **MERGE** or **DO NOT MERGE** with reasons, plus any residual risks.

- If **MERGE** and a PR is in scope: merge it (`gh pr merge <n> --squash` unless the user prefers otherwise; confirm the method if ambiguous). If no PR, merge the fixes branch into `main`.
- If **DO NOT MERGE**: do not merge. Summarize the blockers and stop.

Report the final recommendation, what was merged (or why not), and the path to the findings doc.
