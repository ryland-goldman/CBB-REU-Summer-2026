"""Cornell Linac — Sections 2–8 stage (Impact-T generic constant-gradient TW).

The rest of the straight electron line to CHESS, after `linac_sec1`: seven S-band
(2856 MHz) traveling-wave sections (CEA 2/3/4/5 + CU 3/4/5) chained into ONE
Impact-T deck and integrated as one time-ordered beam. Reads `linac_sec1`'s
captured ~25 MeV exit beam (its last/exit openPMD dump) and accelerates it on-crest
to ≈308 MeV at the default 11 MW klystron point (the real tapered-bore survivors;
the full-beam ⟨KE⟩ is ≈309 MeV, the bore scrapes lower-energy off-axis particles).

Unlike the WarpX stages, this stage is an external serial Impact-T run
(`ImpactTexe`) driven through lume-impact, so it runs IN-PROCESS via `ImpactStage`
(no pywarpx global-geometry binding ⇒ no per-stage subprocess), while still reusing
the pipeline's repo-root chdir + fd-limit raise + shared log.

No field maps: sections 2–8 have no GPT/CST maps, so the field shape reuses the
shipped lume-impact `rfdata4–7` traveling-wave template (shape only), rescaled per
section to the calibrated gradient. Space charge is OFF by default (γ>49 at entry;
`SPACE_CHARGE=True` opts into an exploratory single-bunch Impact-T SC run). Quads are
present at real lengths but OFF (K1=0) for the headline beam — the A→T calibration
is undocumented, so the FODO line is a separate, clearly-labeled exploratory figure.

Usage:
    import linac_rest
    linac_rest.config(POWER_MW=11.0)           # optional overrides (one power convention)
    linac_rest.run()                           # build deck + Impact-T + plots
    linac_rest.run(plots=False)                # build + sim only
    linac_rest.plot()                          # plots only

Parameter names match the module-level constants in
`linac_rest/build_linac_rest_lattice.py` and `linac_rest/linac_rest_sim.py`.
"""

from pipeline._impact_runner import ImpactStage

# Diags dir for run() (the single operating point); a scan can override via OUTDIR.
DEFAULT_OUTDIR = "linac_rest/diags/main"

_stage = ImpactStage(
    name="linac_rest",
    build_module="linac_rest.build_linac_rest_lattice",
    sim_module="linac_rest.linac_rest_sim",
    plot_module="linac_rest.plot_linac_rest",
)
config = _stage.config
run = _stage.run
plot = _stage.plot


def resolve_outdir():
    """Return the diags dir the next run() will write to (OUTDIR override or default).

    Used by `pipeline/run_pipeline.py` so the final-beam summary reads the same
    directory the sim wrote, without importing the lume-impact-laden sim module.
    """
    return _stage._params.get("OUTDIR") or DEFAULT_OUTDIR
