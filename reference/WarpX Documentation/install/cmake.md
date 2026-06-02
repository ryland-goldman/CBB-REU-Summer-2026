<a id="install-build-cmake"></a>

# Build from Source

[CMake](https://cmake.org) is our primary build system.
If you are new to CMake, we recommend starting with [this concise tutorial](https://hsf-training.github.io/hsf-training-cmake-webpage/) from the HEP Software Foundation.
For those primarily interested in building the project, focus on these key sections: [1. Introduction](https://hsf-training.github.io/hsf-training-cmake-webpage/01-intro/index.html), [2. Building with CMake](https://hsf-training.github.io/hsf-training-cmake-webpage/02-building/index.html), and [9. Finding Packages](https://hsf-training.github.io/hsf-training-cmake-webpage/09-findingpackages/index.html).

Please refer to the [HPC Systems](hpc.md#install-hpc) section.

Software dependencies of WarpX.

Configuration, compilation and install.

All build configuration options.

<a id="install-build-dependencies"></a>

## Install Dependencies

To begin, obtain a copy of the WarpX source code:

```bash
git clone https://github.com/BLAST-WarpX/warpx.git $HOME/src/warpx
cd $HOME/src/warpx
```

WarpX relies on [widely-used third-party software](dependencies.md#install-dependencies).
Below, you’ll find instructions for installing these dependencies using various package managers.
To ensure compatibility, pick **one** package manager for your development workflows.

### Install with conda-forge

Conda provides a convenient way to install dependencies across Linux, macOS, and Windows platforms.

[conda-forge](https://conda-forge.org/download/) is a community-led collection of recipes, build infrastructure and distributions for the conda package manager, offering cross-platform compatibility at the user level.

### With MPI (only Linux/macOS)

```bash
conda create -n warpx-cpu-mpich-dev -c conda-forge blaspp boost ccache cmake compilers git lapackpp "openpmd-api=*=mpi_mpich*" openpmd-viewer packaging pytest python python-build make numpy pandas scipy setuptools yt "fftw=*=mpi_mpich*" pkg-config matplotlib mamba mpich mpi4py ninja pip virtualenv wheel
conda activate warpx-cpu-mpich-dev

# compile WarpX with -DWarpX_MPI=ON
# for pip, use: export WARPX_MPI=ON
```

### Without MPI

```bash
conda create -n warpx-cpu-dev -c conda-forge blaspp boost ccache cmake compilers git lapackpp openpmd-api openpmd-viewer packaging pytest python python-build make numpy pandas scipy setuptools yt fftw pkg-config matplotlib mamba ninja pip virtualenv wheel
conda activate warpx-cpu-dev

# compile WarpX with -DWarpX_MPI=OFF
# for pip, use: export WARPX_MPI=OFF
```

For OpenMP support, you will further need:

### Linux

```bash
mamba install -c conda-forge libgomp
```

### macOS or Windows

```bash
mamba install -c conda-forge llvm-openmp
```

For Nvidia CUDA GPU support, you will need to have [a recent CUDA driver installed](https://developer.nvidia.com/cuda-downloads) or you can lower the CUDA version of [the Nvidia cuda package](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#conda-installation) and [conda-forge to match your drivers](https://docs.cupy.dev/en/stable/install.html#install-cupy-from-conda-forge) and then add these packages:

```bash
mamba install -c nvidia -c conda-forge cuda cuda-nvtx-dev cupy
```

More info for [CUDA-enabled ML packages](https://twitter.com/jeremyphoward/status/1697435241152127369).

### Install with Spack

Spack provides another option for installing dependencies on Linux and macOS systems.

[Spack](https://spack.readthedocs.io) is a flexible, user-level package manager designed primarily for Linux, with growing support for macOS and planned future support for Windows.

To begin, download a [WarpX Spack desktop development environment](https://github.com/BLAST-WarpX/warpx/blob/development/Tools/machines/desktop) configuration.
For most desktop development work, we recommend using the OpenMP environment for CPUs, unless you have a supported GPU device.

* **Debian/Ubuntu** Linux:
  * OpenMP: `system=ubuntu; compute=openmp` (CPUs)
  * CUDA: `system=ubuntu; compute=cuda` (Nvidia GPUs)
  * ROCm: `system=ubuntu; compute=rocm` (AMD GPUs)
  * SYCL: *todo* (Intel GPUs)
* **macOS**: first, prepare with `brew install gpg2; brew install gcc`
  * OpenMP: `system=macos; compute=openmp`

If you already [installed Spack](https://spack.io), we recommend to activate its [binary caches](https://spack.io/spack-binary-packages/) for faster builds:

```bash
spack mirror add rolling https://binaries.spack.io/develop
spack buildcache keys --install --trust
```

Now install the WarpX dependencies in a new WarpX development environment:

```bash
# download environment file
curl -sLO https://raw.githubusercontent.com/BLAST-WarpX/warpx/development/Tools/machines/desktop/spack-${system}-${compute}.yaml

# create new development environment
spack env create warpx-${compute}-dev spack-${system}-${compute}.yaml
spack env activate warpx-${compute}-dev

# installation
spack install
python3 -m pip install jupyter matplotlib numpy openpmd-api openpmd-viewer pandas scipy virtualenv yt
```

In new terminal sessions, re-activate the environment with

```bash
spack env activate warpx-openmp-dev
```

again.
Replace `openmp` with the equivalent you chose.

Compile WarpX with `-DWarpX_MPI=ON`.
For `pip`, use `export WARPX_MPI=ON`.

### Install with Brew

Brew can be used to install dependencies on Linux and macOS.

[Homebrew (Brew)](https://brew.sh) is a user-level package manager primarily for [Apple macOS](https://en.wikipedia.org/wiki/MacOS), but also supports Linux.

```bash
brew update
brew tap openpmd/openpmd
brew install adios2      # for openPMD
brew install ccache
brew install cmake
brew install fftw        # for PSATD
brew install git
brew install hdf5-mpi    # for openPMD
brew install libomp
brew unlink gcc
brew link --force libomp
brew install pkg-config  # for fftw
brew install open-mpi
brew install openblas    # for PSATD in RZ
brew install openpmd-api # for openPMD
```

If you also want to compile with PSATD in RZ, you need to manually install BLAS++ and LAPACK++:

```bash
sudo mkdir -p /usr/local/bin/
sudo curl -L -o /usr/local/bin/cmake-easyinstall https://raw.githubusercontent.com/ax3l/cmake-easyinstall/main/cmake-easyinstall
sudo chmod a+x /usr/local/bin/cmake-easyinstall

cmake-easyinstall --prefix=/usr/local git+https://github.com/icl-utk-edu/blaspp.git \
    -Duse_openmp=OFF -Dbuild_tests=OFF -DCMAKE_VERBOSE_MAKEFILE=ON
cmake-easyinstall --prefix=/usr/local git+https://github.com/icl-utk-edu/lapackpp.git \
    -Duse_cmake_find_lapack=ON -Dbuild_tests=OFF -DCMAKE_VERBOSE_MAKEFILE=ON
```

Compile WarpX with `-DWarpX_MPI=ON`.
For `pip`, use `export WARPX_MPI=ON`.

### Install with APT

The [Advanced Package Tool (APT)](https://en.wikipedia.org/wiki/APT_(software)) is a system-level package manager on Debian-based Linux distributions, including Ubuntu.

### With MPI (only Linux/macOS)

```bash
sudo apt update
sudo apt install build-essential ccache cmake g++ git libfftw3-mpi-dev libfftw3-dev libhdf5-openmpi-dev libopenmpi-dev pkg-config python3 python3-dev python3-matplotlib python3-mpi4py python3-numpy python3-pandas python3-pip python3-scipy python3-venv

# optional:
# for CUDA, either install
#   https://developer.nvidia.com/cuda-downloads (preferred)
# or, if your Debian/Ubuntu is new enough, use the packages
#   sudo apt install nvidia-cuda-dev libcub-dev

# compile WarpX with -DWarpX_MPI=ON
# for pip, use: export WARPX_MPI=ON
```

### Without MPI

```bash
sudo apt update
sudo apt install build-essential ccache cmake g++ git libfftw3-dev libfftw3-dev libhdf5-dev pkg-config python3 python3-dev python3-matplotlib python3-numpy python3-pandas python3-pip python3-scipy python3-venv

# optional:
# for CUDA, either install
#   https://developer.nvidia.com/cuda-downloads (preferred)
# or, if your Debian/Ubuntu is new enough, use the packages
#   sudo apt install nvidia-cuda-dev libcub-dev

# compile WarpX with -DWarpX_MPI=OFF
# for pip, use: export WARPX_MPI=OFF
```

<a id="install-build-code"></a>

## Build the Code

<a id="build-the-executable-with-cmake"></a>

### Build the Executable with CMake

To build WarpX from the source directory, execute these commands:

```bash
# Configure the build system
# Additional options available, such as:
#   -DWarpX_PYTHON=ON
#   -DCMAKE_INSTALL_PREFIX=$HOME/sw/warpx
cmake -S . -B build

# Compile using four parallel threads
cmake --build build -j 4
```

**That’s it!**
The 3D WarpX binary is now available in `build/bin/` and is ready to [run](../usage/how_to_run.md#usage-run) with any [3D example input file](../usage/examples.md#usage-examples).
You can either run the binary directly from this location or copy it to another directory.

For a system-wide installation, use the following command:

```bash
# for default install paths, you will need administrator rights, e.g. with sudo:
cmake --build build --target install
```

You can inspect and modify build options after running `cmake -S . -B build` with either

```bash
ccmake build
```

or by adding arguments with `-D<OPTION>=<VALUE>` to the first CMake call.
For example, this builds WarpX in all geometries, enables Python bindings and Nvidia GPU (CUDA) support:

```bash
cmake -S . -B build -DWarpX_DIMS="1;2;3;RZ;RCYLINDER;RSPHERE" -DWarpX_COMPUTE=CUDA
```

An executable WarpX binary with the current compile-time options encoded in its file name will be created in `build/bin/`.
Note that you need separate binaries to run 1D, 2D, 3D, RZ, RCYLINDER, RSPHERE geometry inputs scripts.
Additionally, a [symbolic link](https://en.wikipedia.org/wiki/Symbolic_link) named `warpx` can be found in that directory, which points to the last built WarpX executable.

More details on running simulations are in the section [Run WarpX](../usage/how_to_run.md#usage-run).
Alternatively, read on and also build our PICMI Python interface.

<a id="install-build-python-cmake"></a>

### Build the Python Interface with CMake

#### NOTE
First, ensure your Python development environment is up-to-date:

```bash
python3 -m pip install -U pip
python3 -m pip install -U build packaging setuptools[core] wheel
python3 -m pip install -U cmake
python3 -m pip install -r requirements.txt
```

To build the PICMI Python bindings, configure WarpX to generate a library and install it using our `pip_install` *CMake target*:

```bash
# Configure with all WarpX dimensionalities and Python support enabled
cmake -S . -B build_py -DWarpX_DIMS="1;2;3;RZ;RCYLINDER;RSPHERE" -DWarpX_PYTHON=ON

# Build and install the Python package
cmake --build build_py --target pip_install -j 4
```

**That’s it!**
You can now [run a first 3D PICMI script](../usage/python.md#usage-picmi) from our [examples](../usage/examples.md#usage-examples).

Developers could now change the WarpX source code and then call the build line again to refresh the Python installation.

<a id="install-build-python-pip"></a>

### Build the Python Interface with pip

This section is relevant for Python package management, mainly for maintainers or people that rather like to interact only with `pip`.

One can build and install `pywarpx` from the root of the WarpX source tree:

```bash
python3 -m pip wheel -v .
python3 -m pip install pywarpx*whl
```

This will call the CMake logic above implicitly.
Using this workflow has the advantage that it can build and package up multiple libraries with varying `WarpX_DIMS` into one `pywarpx` package.

<a id="install-build-options"></a>

## Build Options

### Configure your Compiler

To use a specific compiler instead of the system default, set the appropriate environment variables.
For instance, to use Clang/LLVM:

```bash
export CC=$(which clang)
export CXX=$(which clang++)
```

For CUDA development, specify both the CUDA compiler and the host C++ compiler:

```bash
export CUDACXX=$(which nvcc)
export CUDAHOSTCXX=$(which clang++)
```

We also support adding [additional compiler flags via environment variables](https://cmake.org/cmake/help/latest/manual/cmake-language.7.html#cmake-language-environment-variables) such as [CXXFLAGS](https://cmake.org/cmake/help/latest/envvar/CXXFLAGS.html)/[LDFLAGS](https://cmake.org/cmake/help/latest/envvar/LDFLAGS.html):

```bash
# example: treat all compiler warnings as errors
export CXXFLAGS="-Werror"
```

#### NOTE
Please clean your build directory with `rm -rf build/` after changing the compiler.
Now call `cmake -S . -B build` (+ further options) again to re-initialize the build configuration.

### CMake

| CMake Option                | Default & Values                 | Description                                                                                                                          |
|-----------------------------|----------------------------------|--------------------------------------------------------------------------------------------------------------------------------------|
| `CMAKE_BUILD_TYPE`          | RelWithDebInfo/**Release**/Debug | [Type of build, symbols & optimizations](https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html)                         |
| `CMAKE_INSTALL_PREFIX`      | system-dependent path            | [Install path prefix](https://cmake.org/cmake/help/latest/variable/CMAKE_INSTALL_PREFIX.html)                                        |
| `CMAKE_VERBOSE_MAKEFILE`    | ON/**OFF**                       | [Print all compiler commands to the terminal during build](https://cmake.org/cmake/help/latest/variable/CMAKE_VERBOSE_MAKEFILE.html) |
| `WarpX_APP`                 | **ON**/OFF                       | Build the WarpX executable application                                                                                               |
| `WarpX_ASCENT`              | ON/**OFF**                       | Ascent in situ visualization                                                                                                         |
| `WarpX_CATALYST`            | ON/**OFF**                       | Catalyst in situ visualization                                                                                                       |
| `WarpX_COMPUTE`             | NOACC/**OMP**/CUDA/SYCL/HIP      | On-node, accelerated computing backend                                                                                               |
| `WarpX_DIMS`                | **3**/2/1/RZ/RCYLINDER/RSPHERE   | Simulation dimensionality. Use `"1;2;3;RZ;RCYLINDER;RSPHERE"` for all.                                                               |
| `WarpX_EB`                  | **ON**/OFF                       | Embedded boundary support (not supported in RZ, RCYLINDER, and RSPHERE  yet)                                                         |
| `WarpX_PETSC`               | ON/**OFF**                       | PETSc linear/nonlinear solvers via AMReX                                                                                             |
| `WarpX_IPO`                 | ON/**OFF**                       | Compile WarpX with interprocedural optimization (aka LTO)                                                                            |
| `WarpX_LIB`                 | ON/**OFF**                       | Build WarpX as a library, e.g., for PICMI Python                                                                                     |
| `WarpX_MPI`                 | **ON**/OFF                       | Multi-node support (message-passing)                                                                                                 |
| `WarpX_MPI_THREAD_MULTIPLE` | **ON**/OFF                       | MPI thread-multiple support, i.e. for `async_io`                                                                                     |
| `WarpX_OPENPMD`             | **ON**/OFF                       | openPMD I/O (HDF5, ADIOS)                                                                                                            |
| `WarpX_PRECISION`           | SINGLE/**DOUBLE**                | Floating point precision (single/double)                                                                                             |
| `WarpX_PARTICLE_PRECISION`  | SINGLE/**DOUBLE**                | Particle floating point precision (single/double), defaults to WarpX_PRECISION value if not set                                      |
| `WarpX_FASTMATH`            | ON/**OFF**                       | Enable fast-math optimizations                                                                                                       |
| `WarpX_FFT`                 | ON/**OFF**                       | FFT-based solvers                                                                                                                    |
| `WarpX_PYTHON`              | ON/**OFF**                       | Python bindings                                                                                                                      |
| `WarpX_QED`                 | **ON**/OFF                       | QED support (requires PICSAR)                                                                                                        |
| `WarpX_QED_TABLE_GEN`       | ON/**OFF**                       | QED table generation support (requires PICSAR and Boost)                                                                             |
| `WarpX_QED_TOOLS`           | ON/**OFF**                       | Build external tool to generate QED lookup tables (requires PICSAR and Boost)                                                        |
| `WarpX_QED_TABLES_GEN_OMP`  | **AUTO**/ON/OFF                  | Enables OpenMP support for QED lookup tables generation                                                                              |
| `WarpX_SENSEI`              | ON/**OFF**                       | SENSEI in situ visualization                                                                                                         |
| `Python_EXECUTABLE`         | (newest found)                   | Path to Python executable                                                                                                            |
| `PY_PIP_OPTIONS`            | `-v`                             | Additional options for `pip`, e.g., `-vvv;-q`                                                                                        |
| `PY_PIP_INSTALL_OPTIONS`    |                                  | Additional options for `pip install`, e.g., `--user;-q`                                                                              |

WarpX can be configured in further detail with options from AMReX, which are documented in the AMReX manual:

* [general AMReX build options](https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#customization-options)
* [GPU-specific options](https://amrex-codes.github.io/amrex/docs_html/GPU.html#building-gpu-support).

**Developers** might be interested in additional options that control dependencies of WarpX.
By default, the most important dependencies of WarpX are automatically downloaded for convenience:

| CMake Option              | Default & Values                             | Description                                                                                                    |
|---------------------------|----------------------------------------------|----------------------------------------------------------------------------------------------------------------|
| `BUILD_SHARED_LIBS`       | ON/**OFF**                                   | [Build shared libraries for dependencies](https://cmake.org/cmake/help/latest/variable/BUILD_SHARED_LIBS.html) |
| `WarpX_CCACHE`            | **ON**/OFF                                   | Search and use CCache to speed up rebuilds.                                                                    |
| `WarpX_UNITY_BUILD`       | ON/**OFF**                                   | WarpX library as unity build (single TU)                                                                       |
| `AMReX_CUDA_PTX_VERBOSE`  | ON/**OFF**                                   | Print CUDA code generation statistics from `ptxas`.                                                            |
| `WarpX_amrex_src`         | *None*                                       | Path to AMReX source directory (preferred if set)                                                              |
| `WarpX_amrex_repo`        | `https://github.com/AMReX-Codes/amrex.git`   | Repository URI to pull and build AMReX from                                                                    |
| `WarpX_amrex_branch`      | *we set and maintain a compatible commit*    | Repository branch for `WarpX_amrex_repo`                                                                       |
| `WarpX_amrex_internal`    | **ON**/OFF                                   | Needs a pre-installed AMReX library if set to `OFF`                                                            |
| `WarpX_openpmd_src`       | *None*                                       | Path to openPMD-api source directory (preferred if set)                                                        |
| `WarpX_openpmd_repo`      | `https://github.com/openPMD/openPMD-api.git` | Repository URI to pull and build openPMD-api from                                                              |
| `WarpX_openpmd_branch`    | `0.17.0`                                     | Repository branch for `WarpX_openpmd_repo`                                                                     |
| `WarpX_openpmd_internal`  | **ON**/OFF                                   | Needs a pre-installed openPMD-api library if set to `OFF`                                                      |
| `WarpX_picsar_src`        | *None*                                       | Path to PICSAR source directory (preferred if set)                                                             |
| `WarpX_picsar_repo`       | `https://github.com/ECP-WarpX/picsar.git`    | Repository URI to pull and build PICSAR from                                                                   |
| `WarpX_picsar_branch`     | *we set and maintain a compatible commit*    | Repository branch for `WarpX_picsar_repo`                                                                      |
| `WarpX_picsar_internal`   | **ON**/OFF                                   | Needs a pre-installed PICSAR library if set to `OFF`                                                           |
| `WarpX_pyamrex_src`       | *None*                                       | Path to PICSAR source directory (preferred if set)                                                             |
| `WarpX_pyamrex_repo`      | `https://github.com/AMReX-Codes/pyamrex.git` | Repository URI to pull and build pyAMReX from                                                                  |
| `WarpX_pyamrex_branch`    | *we set and maintain a compatible commit*    | Repository branch for `WarpX_pyamrex_repo`                                                                     |
| `WarpX_pyamrex_internal`  | **ON**/OFF                                   | Needs a pre-installed pyAMReX library if set to `OFF`                                                          |
| `WarpX_PYTHON_IPO`        | **ON**/OFF                                   | Build Python w/ interprocedural/link optimization (IPO/LTO)                                                    |
| `WarpX_pybind11_src`      | *None*                                       | Path to pybind11 source directory (preferred if set)                                                           |
| `WarpX_pybind11_repo`     | `https://github.com/pybind/pybind11.git`     | Repository URI to pull and build pybind11 from                                                                 |
| `WarpX_pybind11_branch`   | *we set and maintain a compatible commit*    | Repository branch for `WarpX_pybind11_repo`                                                                    |
| `WarpX_pybind11_internal` | **ON**/OFF                                   | Needs a pre-installed pybind11 library if set to `OFF`                                                         |
| `WarpX_TEST_CLEANUP`      | ON/**OFF**                                   | Clean up automated test directories                                                                            |
| `WarpX_TEST_DEBUGGER`     | ON/**OFF**                                   | Run automated tests without AMReX signal handling (to attach debuggers)                                        |
| `WarpX_TEST_FPETRAP`      | ON/**OFF**                                   | Run automated tests with FPE-trapping runtime parameters                                                       |
| `WarpX_BACKTRACE_INFO`    | ON/**OFF**                                   | Compile with -g1 for minimal debug symbols (currently used in CI tests)                                        |

For example, one can also build against a local AMReX copy.
Assuming AMReX’ source is located in `$HOME/src/amrex`, add the `cmake` argument `-DWarpX_amrex_src=$HOME/src/amrex`.
Relative paths are also supported, e.g. `-DWarpX_amrex_src=../amrex`.

Or build against an AMReX feature branch of a colleague.
Assuming your colleague pushed AMReX to `https://github.com/WeiqunZhang/amrex/` in a branch `new-feature` then pass to `cmake` the arguments: `-DWarpX_amrex_repo=https://github.com/WeiqunZhang/amrex.git -DWarpX_amrex_branch=new-feature`.
More details on this [workflow are described here](../developers/how_to_compile_locally.md#developers-local-compile-src).

You can speed up the install further if you pre-install these dependencies, e.g. with a package manager.
Set `-DWarpX_<dependency-name>_internal=OFF` and add installation prefix of the dependency to the environment variable [CMAKE_PREFIX_PATH](https://cmake.org/cmake/help/latest/envvar/CMAKE_PREFIX_PATH.html).
Please see the [introduction to CMake](#install-build-cmake) if this sounds new to you.
More details on this [workflow are described here](../developers/how_to_compile_locally.md#developers-local-compile-findpackage).

If you re-compile often, consider installing the [Ninja](https://github.com/ninja-build/ninja/wiki/Pre-built-Ninja-packages) build system.
Pass `-G Ninja` to the CMake configuration call to speed up parallel compiles.

### pip

Environment variables can be used to control the build step:

| Environment Variable        | Default & Values                 | Description                                                                                     |
|-----------------------------|----------------------------------|-------------------------------------------------------------------------------------------------|
| `WARPX_COMPUTE`             | NOACC/**OMP**/CUDA/SYCL/HIP      | On-node, accelerated computing backend                                                          |
| `WARPX_DIMS`                | `"1;2;3;RZ;RCYLINDER;RSPHERE"`   | Simulation dimensionalities (semicolon-separated list)                                          |
| `WARPX_EB`                  | **ON**/OFF                       | Embedded boundary support (not supported in RZ, RCYLINDER, and RSPHERE yet)                     |
| `WARPX_PETSC`               | ON/**OFF**                       | PETSc linear/nonlinear solvers via AMReX                                                        |
| `WARPX_MPI`                 | ON/**OFF**                       | Multi-node support (message-passing)                                                            |
| `WARPX_OPENPMD`             | **ON**/OFF                       | openPMD I/O (HDF5, ADIOS)                                                                       |
| `WARPX_PRECISION`           | SINGLE/**DOUBLE**                | Floating point precision (single/double)                                                        |
| `WARPX_PARTICLE_PRECISION`  | SINGLE/**DOUBLE**                | Particle floating point precision (single/double), defaults to WarpX_PRECISION value if not set |
| `WARPX_FFT`                 | ON/**OFF**                       | FFT-based solvers                                                                               |
| `WARPX_QED`                 | **ON**/OFF                       | PICSAR QED (requires PICSAR)                                                                    |
| `WARPX_QED_TABLE_GEN`       | ON/**OFF**                       | QED table generation (requires PICSAR and Boost)                                                |
| `BUILD_PARALLEL`            | `2`                              | Number of threads to use for parallel builds                                                    |
| `BUILD_SHARED_LIBS`         | ON/**OFF**                       | Build shared libraries for dependencies                                                         |
| `HDF5_USE_STATIC_LIBRARIES` | ON/**OFF**                       | Prefer static libraries for HDF5 dependency (openPMD)                                           |
| `ADIOS_USE_STATIC_LIBS`     | ON/**OFF**                       | Prefer static libraries for ADIOS1 dependency (openPMD)                                         |
| `WARPX_AMREX_SRC`           | *None*                           | Absolute path to AMReX source directory (preferred if set)                                      |
| `WARPX_AMREX_REPO`          | *None (uses cmake default)*      | Repository URI to pull and build AMReX from                                                     |
| `WARPX_AMREX_BRANCH`        | *None (uses cmake default)*      | Repository branch for `WARPX_AMREX_REPO`                                                        |
| `WARPX_AMREX_INTERNAL`      | **ON**/OFF                       | Needs a pre-installed AMReX library if set to `OFF`                                             |
| `WARPX_OPENPMD_SRC`         | *None*                           | Absolute path to openPMD-api source directory (preferred if set)                                |
| `WARPX_OPENPMD_INTERNAL`    | **ON**/OFF                       | Needs a pre-installed openPMD-api library if set to `OFF`                                       |
| `WARPX_PICSAR_SRC`          | *None*                           | Absolute path to PICSAR source directory (preferred if set)                                     |
| `WARPX_PICSAR_INTERNAL`     | **ON**/OFF                       | Needs a pre-installed PICSAR library if set to `OFF`                                            |
| `WARPX_PYAMREX_SRC`         | *None*                           | Absolute path to pyAMReX source directory (preferred if set)                                    |
| `WARPX_PYAMREX_INTERNAL`    | **ON**/OFF                       | Needs a pre-installed pyAMReX library if set to `OFF`                                           |
| `WARPX_PYTHON_IPO`          | **ON**/OFF                       | Build Python w/ interprocedural/link optimization (IPO/LTO)                                     |
| `WARPX_PYBIND11_SRC`        | *None*                           | Absolute path to pybind11 source directory (preferred if set)                                   |
| `WARPX_PYBIND11_INTERNAL`   | **ON**/OFF                       | Needs a pre-installed pybind11 library if set to `OFF`                                          |
| `WARPX_CCACHE_PROGRAM`      | First found `ccache` executable. | Set to `NO` to disable CCache.                                                                  |
| `PYWARPX_LIB_DIR`           | *None*                           | If set, search for pre-built WarpX C++ libraries (see below)                                    |

Note that we currently change the `WARPX_MPI` default intentionally to `OFF`, to simplify a first install from source.

Some hints and workflows follow.
Developers, that want to test a change of the source code but did not change the `pywarpx` version number, can force a reinstall via:

```bash
python3 -m pip install --force-reinstall --no-deps -v .
```

Some Developers like to code directly against a local copy of AMReX, changing both code-bases at a time:

```bash
WARPX_AMREX_SRC=$PWD/../amrex python3 -m pip install --force-reinstall --no-deps -v .
```

Additional environment control as common for CMake ([see above](#install-build-cmake)) can be set as well, e.g. `CC`, CXX\`, and `CMAKE_PREFIX_PATH` hints.
So another sophisticated example might be: use Clang as the compiler, build with local source copies of PICSAR and AMReX, support the FFT-based solvers, MPI and openPMD, hint a parallel HDF5 installation in `$HOME/sw/hdf5-parallel-1.10.4`, and only build 2D and 3D geometry:

```bash
CC=$(which clang) CXX=$(which clang++) WARPX_AMREX_SRC=$PWD/../amrex WARPX_PICSAR_SRC=$PWD/../picsar WARPX_FFT=ON WARPX_MPI=ON WARPX_DIMS="2;3" CMAKE_PREFIX_PATH=$HOME/sw/hdf5-parallel-1.10.4:$CMAKE_PREFIX_PATH python3 -m pip install --force-reinstall --no-deps -v .
```

Here we wrote this all in one line, but one can also set all environment variables in a development environment and keep the pip call nice and short as in the beginning.
Note that you need to use absolute paths for external source trees, because pip builds in a temporary directory, e.g. `export WARPX_AMREX_SRC=$HOME/src/amrex`.

All of this can also be run from CMake.
This is the workflow most developers will prefer as it allows rapid re-compiles:

```bash
# build WarpX executables and libraries
cmake -S . -B build_py -DWarpX_DIMS="1;2;3;RZ;RCYLINDER;RSPHERE" -DWarpX_PYTHON=ON

# build & install Python only
cmake --build build_py -j 4 --target pip_install
```

There is also a `--target pip_install_nodeps` option that [skips pip-based dependency checks](../developers/how_to_compile_locally.md#developers-local-compile-pylto).

WarpX release managers might also want to generate a self-contained source package that can be distributed to exotic architectures:

```bash
python setup.py sdist --dist-dir .
python3 -m pip wheel -v pywarpx-*.tar.gz
python3 -m pip install *whl
```

The above steps can also be executed in one go to build from source on a machine:

```bash
python3 setup.py sdist --dist-dir .
python3 -m pip install -v pywarpx-*.tar.gz
```

Last but not least, you can uninstall `pywarpx` as usual with:

```bash
python3 -m pip uninstall pywarpx
```
