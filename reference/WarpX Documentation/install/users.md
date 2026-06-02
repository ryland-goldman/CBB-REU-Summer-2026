<a id="install-methods"></a>

# Installation Methods

<style>
.rst-content section>img {
    width: 30px;
    margin-bottom: 0;
    margin-top: 0;
    margin-right: 15px;
    margin-left: 15px;
    float: left;
}
</style>

Our community is here to help – please [report installation issues](https://github.com/BLAST-WarpX/warpx/issues/new) if you encounter any.

Please choose **one** of the installation methods below to get started.

## HPC Systems

If you want to use WarpX on a specific high-performance computing (HPC) system, please go directly to our [HPC system-specific documentation](hpc.md#install-hpc).

<a id="install-methods-conda"></a>

## Using the conda-forge Package

A package for WarpX is available via [conda-forge](https://conda-forge.org/download/).

```bash
mamba create -n warpx -c conda-forge warpx
mamba activate warpx
```

#### NOTE
The `warpx` package on conda-forge does not yet provide [GPU support](https://github.com/conda-forge/warpx-feedstock/issues/89).

<a id="install-methods-spack"></a>

## Using the Spack Package

Packages for WarpX are available via the [Spack](https://spack.readthedocs.io) package manager.
The `warpx` package installs executables. The `warpx +python` variant also builds Python bindings, which can be used with [PICMI](https://github.com/picmi-standard/picmi).

```bash
# optional: activate Spack binary caches
spack mirror add rolling https://binaries.spack.io/develop
spack buildcache keys --install --trust

# see `spack info py-warpx` for build options.
# optional arguments:       -mpi compute=cuda
spack install warpx +python
spack load warpx +python
```

See `spack info warpx` and [the official Spack tutorial](https://spack-tutorial.readthedocs.io) for more information.

<a id="install-methods-pypi"></a>

## Using the PyPI Package

If you have the [WarpX dependencies](dependencies.md#install-dependencies) installed, you can use `pip` to install WarpX (with PICMI) [from source](cmake.md#install-build-cmake):

```bash
python3 -m pip install -U pip
python3 -m pip install -U build packaging setuptools[core] wheel
python3 -m pip install -U cmake

python3 -m pip wheel -v git+https://github.com/BLAST-WarpX/warpx.git
python3 -m pip install *whl
```

Pre-compiled binary packages will be published on [PyPI](https://pypi.org/) in the future for faster installs.
Please consider using [conda](#install-methods-conda) in the meantime.

<a id="install-methods-brew"></a>

## Using the Brew Package

#### NOTE
Coming soon.

<a id="install-methods-cmake"></a>

## From Source with CMake

After installing the [WarpX dependencies](dependencies.md#install-dependencies), you can also install WarpX from source with [CMake](https://cmake.org/):

```bash
# get the source code
git clone https://github.com/BLAST-WarpX/warpx.git $HOME/src/warpx
cd $HOME/src/warpx

# configure
cmake -S . -B build

# optional: change configuration
ccmake build

# compile
#   on Windows:          --config RelWithDebInfo
cmake --build build -j 4

# executables for WarpX are now in build/bin/
```

For more details on how to configure WarpX from source, please see the section [Build from Source](cmake.md#install-build-cmake).

<a id="install-tips-macos"></a>

## Tips for macOS Users

See also: A. Huebl, [Working With Multiple Package Managers](https://collegeville.github.io/CW20/WorkshopResources/WhitePapers/huebl-working-with-multiple-pkg-mgrs.pdf), [Collegeville Workshop (CW20)](https://collegeville.github.io/CW20/), 2020
