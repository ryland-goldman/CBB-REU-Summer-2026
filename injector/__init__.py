"""Cornell Linac injector stage facade: config()/run()/plot().

See injector/README.md for physics, parameters, and gotchas.
"""

from pipeline._runner import Stage

# Single source of truth for the default diags dir; imported by injector_sim.py so the
# parent can resolve the path without importing the pywarpx-laden sim module.
DEFAULT_OUTDIR = "injector/diags/main"

_stage = Stage(
    name="injector",
    build_module="injector.build_injector_field",
    sim_module="injector.injector_sim",
    plot_module="injector.plot_injector",
)
config = _stage.config
run = _stage.run
plot = _stage.plot


def faithful_gpt_deck():
    """Print how to switch to the GPT master-deck GUI working point.

    The operating point now lives in injector/injector.yaml; to reproduce the deck's GUI-saved
    phases (crest-referenced Preb-1 −70° / Preb-2 −45°, reversed-install +π absorbed into the
    crest reference) set these `params:` keys and re-run. CAVEAT: this re-introduces Preb-1's
    non-zero net kick, which desyncs the analytic Preb-2 inter-cavity timing — a hardened study
    needs the two-pass timing fix (see injector/README.md).
    """
    print("Edit injector/injector.yaml params: PHASE: crest, PREB1_PHI_OFF: -70.0, "
          "PREB2_PHI_OFF: -45.0, PREB2_REV_PHASE: 0.0")


def resolve_outdir():
    """Return the diags dir the next run() will write to (OUTDIR override or default)."""
    return _stage._params.get("OUTDIR") or DEFAULT_OUTDIR
