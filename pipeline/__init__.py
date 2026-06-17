"""Cornell Linac pipeline package; exposes the cross-stage figure generator.

See pipeline/README.md for the pipeline overview and outputs.
"""


def plot_chain():
    """Generate the cross-stage beam-evolution figures into the repo-root results/."""
    # importlib by full dotted path: this function shadows the plot_chain submodule.
    import importlib
    _mod = importlib.import_module("pipeline.plot_chain")
    _mod.main()
