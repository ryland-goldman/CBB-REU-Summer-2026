# openPMD-viewer documentation

`openPMD-viewer` contains a set of tools to load and visualize the contents of a set
of [openPMD](http://www.openpmd.org/#/start) files (typically, a timeseries).

The routines of `openPMD-viewer` can be used in two ways:

- Using the **Python API**, in order to write a script that loads the data and
  produces a set of pre-defined plots.
- Using the **interactive GUI inside a Jupyter Notebook**, in order to interactively
  visualize the data.

## Installation

You can install openPMD-viewer with `pip`:

```
pip install openpmd-viewer
```

or alternatively with `conda`:

```
conda install -c conda-forge openpmd-viewer
```

## Usage

The notebooks in [tutorials/tutorials.md](tutorials/tutorials.md) demonstrate how to
use both the API and the interactive GUI.

If you wish to use the **interactive GUI**, the installation of `openPMD-viewer`
provides a convenient executable which automatically **creates a new pre-filled
notebook** and **opens it in a browser**. To use it, type in a regular terminal:

```
openPMD_notebook
```

## Contents

- [Tutorials](tutorials/tutorials.md)
- [API reference](api_reference/api_reference.md)
