"""Cornell Linac WarpX pipeline package.

Exposes the cross-stage figure generator so it can be run standalone:

    python -c "import pipeline; pipeline.plot_chain()"

(It is also called automatically at the end of pipeline/run_pipeline.py:main().)
"""


def plot_chain():
    """Generate the cross-stage beam-evolution figures into the repo-root results/."""
    # Import the submodule by its full dotted path: a bare `from pipeline import
    # plot_chain` would resolve to THIS function (it shadows the submodule in the package
    # namespace), so use importlib to get the module object.
    import importlib
    _mod = importlib.import_module("pipeline.plot_chain")
    _mod.main()
