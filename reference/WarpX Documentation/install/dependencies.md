<a id="install-dependencies"></a>

# List of Dependencies

WarpX depends on the following popular third party software.
Please see installation instructions below.

- a mature [C++20](https://en.wikipedia.org/wiki/C%2B%2B20) compiler, e.g., GCC 12+, Clang 14, NVCC 12.4, MSVC 19.39 or newer
- [CMake 3.24.0+](https://cmake.org)
- [Git 2.18+](https://git-scm.com)
- [AMReX](https://amrex-codes.github.io): we automatically download and compile a copy of AMReX
- [PICSAR](https://github.com/ECP-WarpX/picsar): we automatically download and compile a copy of PICSAR

and for Python bindings:

- [pyAMReX](https://github.com/AMReX-Codes/pyamrex): we automatically download and compile a copy of pyAMReX
- [pybind11](https://github.com/pybind/pybind11): we automatically download and compile a copy of pybind11

Optional dependencies include:

- [MPI 3.0+](https://www.mpi-forum.org/docs/): for multi-node and/or multi-GPU execution
- for on-node accelerated compute *one of either*:
  - [OpenMP 3.1+](https://www.openmp.org): for threaded CPU execution or
  - [CUDA Toolkit 12.2+](https://developer.nvidia.com/cuda-downloads): for Nvidia GPU support (see [matching host-compilers](https://gist.github.com/ax3l/9489132)) or
  - [ROCm 6.0+](https://gpuopen.com/learn/amd-lab-notes/amd-lab-notes-rocm-installation-readme/): for AMD GPU support
  - [oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html): for Intel GPU support
- [FFTW3](http://www.fftw.org): for spectral solver (PSATD or IGF) support when running on CPU or SYCL
  - also needs the `pkg-config` tool on Unix
- [BLAS++](https://github.com/icl-utk-edu/blaspp) and [LAPACK++](https://github.com/icl-utk-edu/lapackpp): for spectral solver (PSATD) support in RZ geometry
- [Boost 1.66.0+](https://www.boost.org/): for QED lookup tables generation support
- [openPMD-api 0.17.0+](https://github.com/openPMD/openPMD-api): we automatically download and compile a copy of openPMD-api for openPMD I/O support
  - see [optional I/O backends](https://github.com/openPMD/openPMD-api#dependencies), i.e., ADIOS2 and/or HDF5
- [Ascent 0.8.0+](https://ascent.readthedocs.io): for in situ 3D visualization
- [CCache](https://ccache.dev): to speed up rebuilds (For CUDA support, needs version 3.7.9+ and 4.2+ is recommended)
- [Ninja](https://ninja-build.org): for faster parallel compiles
- [Python 3.11+](https://www.python.org)
  - [mpi4py](https://mpi4py.readthedocs.io)
  - [numpy](https://numpy.org)
  - [periodictable](https://periodictable.readthedocs.io)
  - [picmistandard](https://picmi-standard.github.io)
  - [lasy](https://lasydoc.readthedocs.io)
  - see our `requirements.txt` file for compatible versions

If you are on a high-performance computing (HPC) system, then [please see our separate HPC documentation](hpc.md#install-hpc).

For all other systems, we recommend to use a **package dependency manager**:
Pick *one* of the installation methods below to install all dependencies for WarpX development in a consistent manner.
