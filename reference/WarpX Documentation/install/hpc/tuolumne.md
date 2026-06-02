<a id="building-tuolumne"></a>

# Tuolumne (LLNL)

The [Tuolumne AMD GPU cluster](https://hpc.llnl.gov/hardware/compute-platforms/tuolumne) (short name *tuo*) is located at LLNL.
Tuolumne is an unclassified sibling system of [El Capitan](https://hpc.llnl.gov/hardware/compute-platforms/el-capitan), sharing the same architecture.

El Capitan & Tuolumne provide four AMD MI300A APUs per compute node.

## Introduction

If you are new to this system, **please see the following resources**:

* [Tuolumne overview](https://hpc.llnl.gov/hardware/compute-platforms/tuolumne)
* [LLNL user account](https://lc.llnl.gov/lorenz/mylc/mylc.cgi) (login required)
* [Tuolumne user guide](https://hpc.llnl.gov/documentation/user-guides/using-el-capitan-systems)
* Batch system: [Flux](https://lc.llnl.gov/confluence/display/ELCAPEA/Running+Jobs)
* [Jupyter service](https://lc.llnl.gov/jupyter) ([documentation](https://lc.llnl.gov/confluence/display/LC/JupyterHub+and+Jupyter+Notebook), login required)
* [Production directories](https://lc.llnl.gov/confluence/display/ELCAPEA/File+Systems) (login required):
  * `/p/lustre5/${USER}`: personal directory on the parallel filesystem (also: `lustre2`)
  * Note that the `$HOME` directory and the `/usr/workspace/${USER}` space are NFS mounted and *not* suitable for production quality data generation.

## Login

```bash
ssh tuolumne.llnl.gov
```

<a id="building-tuolumne-preparation"></a>

## Preparation

Use the following commands to download the WarpX source code:

```bash
git clone https://github.com/BLAST-WarpX/warpx.git /p/lustre5/${USER}/tuolumne/src/warpx
```

On Tuolumne, we usually accelerate all computations with the GPU cores of the MI300A APU.
For development purposes, you can also limit yourself to the CPU cores of the MI300A.

### GPU

We use system software modules, add environment hints and further dependencies via the file `$HOME/tuolumne_mi300a_warpx.profile`.
Create it now:

```bash
cp /p/lustre5/${USER}/tuolumne/src/warpx/Tools/machines/tuolumne-llnl/tuolumne_mi300a_warpx.profile.example $HOME/tuolumne_mi300a_warpx.profile
```

### Script Details

```bash
# please set your project account
export proj=""  # change me!

# remembers the location of this script
export MY_PROFILE=$(cd $(dirname $BASH_SOURCE) && pwd)"/"$(basename $BASH_SOURCE)
# early access: not yet used
# if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your $MY_PROFILE file! Please edit its line 2 to continue!"; return; fi

# required dependencies
module load rocm/6.4.1
module load cmake/3.29.2

# optional: faster builds
# ccache is system provided
module load ninja/1.10.2

# optional: for QED support with detailed tables
# TODO: no Boost module found

# optional: for openPMD and PSATD+RZ support
SW_DIR="/p/lustre5/${USER}/tuolumne/warpx/mi300a"
export CMAKE_PREFIX_PATH=${SW_DIR}/hdf5-1.14.1.2:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/c-blosc-2.15.1:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/adios2-2.10.2:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/blaspp-2024.05.31:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/lapackpp-2024.05.31:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/petsc-3.24.0:$CMAKE_PREFIX_PATH

export LD_LIBRARY_PATH=${SW_DIR}/hdf5-1.14.1.2/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/c-blosc-2.15.1/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/adios2-2.10.2/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/blaspp-2024.05.31/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/lapackpp-2024.05.31/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/petsc-3.24.0/lib:$LD_LIBRARY_PATH

export PATH=${SW_DIR}/adios2-2.10.2/bin:${PATH}
export PATH=${SW_DIR}/hdf5-1.14.1.2/bin:${PATH}

# python
module load cray-python/3.11.7

if [ -d "${SW_DIR}/venvs/warpx-tuolumne-mi300a" ]
then
  source ${SW_DIR}/venvs/warpx-tuolumne-mi300a/bin/activate
fi

# an alias to request an interactive batch node for one hour
#   for parallel execution, start on the batch node: srun <command>
alias getNode="salloc -N 1 -t 1:00:00"
# an alias to run a command on a batch node for up to 30min
#   usage: runNode <command>
alias runNode="srun -N 1 --ntasks-per-node=4 -t 0:30:00"

# GPU-aware MPI
export MPICH_GPU_SUPPORT_ENABLED=1
export LDFLAGS="${PE_MPICH_GTL_DIR_amd_gfx942} ${PE_MPICH_GTL_LIBS_amd_gfx942} -Wl,-rpath,${CRAYLIBS_X86_64}"

# Transparent huge pages on CPU
# https://hpc.llnl.gov/documentation/user-guides/using-el-capitan-systems/introduction-and-quickstart/pro-tips
#export LDFLAGS="${LDFLAGS} -lhugetlbfs"
#export HSA_XNACK=1
#export HUGETLB_MORECORE=yes

# optimize ROCm/HIP compilation for MI300A
export AMREX_AMD_ARCH=gfx942

# compiler environment hints
export CC=$(which amdclang)
export CXX=$(which amdclang++)
export FC=$(which amdflang)
export HIPCXX=${CXX}
```

Edit the 2nd line of this script, which sets the `export proj=""` variable.
**Currently, this is unused and can be kept empty.**
Once project allocation becomes required, e.g., if you are member of the project `abcde`, then run `vi $HOME/tuolumne_mi300a_warpx.profile`.
Enter the edit mode by typing `i` and edit line 2 to read:

```bash
export proj="abcde"
```

Exit the `vi` editor with `Esc` and then type `:wq` (write & quit).

#### IMPORTANT
Now, and as the first step on future logins to Tuolumne, activate these environment settings:

```bash
source $HOME/tuolumne_mi300a_warpx.profile
```

Finally, since Tuolumne does not yet provide software modules for some of our dependencies, install them once:

> ```bash
> bash /p/lustre5/${USER}/tuolumne/src/warpx/Tools/machines/tuolumne-llnl/install_mi300a_dependencies.sh
> source /p/lustre5/${USER}/tuolumne/warpx/mi300a/venvs/warpx-tuolumne-mi300a/bin/activate
> ```

> ### Script Details

> ```bash
> #!/bin/bash
> #
> # Copyright 2024 The WarpX Community
> #
> # This file is part of WarpX.
> #
> # Author: Axel Huebl
> # License: BSD-3-Clause-LBNL

> # Exit on first error encountered #############################################
> #
> set -eu -o pipefail


> # Check: ######################################################################
> #
> #   Was tuolumne_mi300a_warpx.profile sourced and configured correctly?
> #   early access: not yet used!
> #if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your tuolumne_mi300a_warpx.profile file! Please edit its line 2 to continue!"; exit 1; fi


> # Remove old dependencies #####################################################
> #
> SRC_DIR="/p/lustre5/${USER}/tuolumne/src"
> SW_DIR="/p/lustre5/${USER}/tuolumne/warpx/mi300a"
> rm -rf ${SW_DIR}
> mkdir -p ${SW_DIR}

> # remove common user mistakes in python, located in .local instead of a venv
> python3 -m pip uninstall -qq -y pywarpx
> python3 -m pip uninstall -qq -y warpx
> python3 -m pip uninstall -qqq -y mpi4py 2>/dev/null || true


> # General extra dependencies ##################################################
> #

> # tmpfs build directory: avoids issues often seen with $HOME and is faster
> build_dir=$(mktemp -d)
> build_procs=24

> # C-Blosc2 (I/O compression)
> if [ -d ${SRC_DIR}/c-blosc2 ]
> then
>   cd ${SRC_DIR}/c-blosc2
>   git fetch --prune
>   git checkout v2.15.1
>   cd -
> else
>   git clone -b v2.15.1 https://github.com/Blosc/c-blosc2.git ${SRC_DIR}/c-blosc2
> fi
> cmake \
>     --fresh                        \
>     -S ${SRC_DIR}/c-blosc2         \
>     -B ${build_dir}/c-blosc2-build \
>     -DBUILD_SHARED_LIBS=OFF        \
>     -DBUILD_TESTS=OFF              \
>     -DBUILD_BENCHMARKS=OFF         \
>     -DBUILD_EXAMPLES=OFF           \
>     -DBUILD_FUZZERS=OFF            \
>     -DBUILD_STATIC=OFF             \
>     -DDEACTIVATE_AVX2=OFF          \
>     -DDEACTIVATE_AVX512=OFF        \
>     -DWITH_SANITIZER=OFF           \
>     -DCMAKE_INSTALL_PREFIX=${SW_DIR}/c-blosc-2.15.1
> cmake \
>     --build ${build_dir}/c-blosc2-build \
>     --target install                    \
>     --parallel ${build_procs}
> rm -rf ${build_dir}/c-blosc2-build

> # HDF5
> if [ -d ${SRC_DIR}/hdf5 ]
> then
>   cd ${SRC_DIR}/hdf5
>   git fetch --prune
>   git checkout hdf5-1_14_1-2
>   cd -
> else
>   git clone -b hdf5-1_14_1-2 https://github.com/HDFGroup/hdf5.git ${SRC_DIR}/hdf5
> fi
> cmake \
>     --fresh                      \
>     -S ${SRC_DIR}/hdf5           \
>     -B ${build_dir}/hdf5-build   \
>     -DBUILD_SHARED_LIBS=OFF        \
>     -DBUILD_TESTING=OFF          \
>     -DHDF5_ENABLE_PARALLEL=ON    \
>     -DCMAKE_INSTALL_PREFIX=${SW_DIR}/hdf5-1.14.1.2
> cmake \
>     --build ${build_dir}/hdf5-build \
>     --target install                \
>     --parallel ${build_procs}
> rm -rf ${build_dir}/hdf5-build

> # ADIOS2
> if [ -d ${SRC_DIR}/adios2 ]
> then
>   cd ${SRC_DIR}/adios2
>   git fetch --prune
>   git checkout v2.10.2
>   cd -
> else
>   git clone -b v2.10.2 https://github.com/ornladios/ADIOS2.git ${SRC_DIR}/adios2
> fi
> cmake \
>     --fresh                      \
>     -S ${SRC_DIR}/adios2         \
>     -B ${build_dir}/adios2-build \
>     -DADIOS2_USE_Blosc2=ON       \
>     -DADIOS2_USE_Campaign=OFF    \
>     -DADIOS2_USE_Fortran=OFF     \
>     -DADIOS2_USE_Python=OFF      \
>     -DADIOS2_USE_ZeroMQ=OFF      \
>     -DBUILD_SHARED_LIBS=OFF      \
>     -DCMAKE_INSTALL_PREFIX=${SW_DIR}/adios2-2.10.2
> cmake \
>     --build ${build_dir}/adios2-build \
>     --target install                  \
>     --parallel ${build_procs}
> rm -rf ${build_dir}/adios2-build

> # BLAS++ (for PSATD+RZ)
> if [ -d ${SRC_DIR}/blaspp ]
> then
>   cd ${SRC_DIR}/blaspp
>   git fetch --prune
>   git checkout v2024.05.31
>   cd -
> else
>   git clone -b v2024.05.31 https://github.com/icl-utk-edu/blaspp.git ${SRC_DIR}/blaspp
> fi
> cmake \
>     --fresh \
>     -S ${SRC_DIR}/blaspp                      \
>     -B ${build_dir}/blaspp-tuolumne-mi300a-build \
>     -Duse_openmp=OFF                          \
>     -Dgpu_backend=hip                         \
>     -DBUILD_SHARED_LIBS=OFF                   \
>     -DCMAKE_CXX_STANDARD=20                   \
>     -DCMAKE_INSTALL_PREFIX=${SW_DIR}/blaspp-2024.05.31
> cmake \
>     --build ${build_dir}/blaspp-tuolumne-mi300a-build \
>     --target install                               \
>     --parallel ${build_procs}
> rm -rf ${build_dir}/blaspp-tuolumne-mi300a-build

> # LAPACK++ (for PSATD+RZ)
> if [ -d ${SRC_DIR}/lapackpp ]
> then
>   cd ${SRC_DIR}/lapackpp
>   git fetch --prune
>   git checkout v2024.05.31
>   cd -
> else
>   git clone -b v2024.05.31 https://github.com/icl-utk-edu/lapackpp.git ${SRC_DIR}/lapackpp
> fi
> cmake \
>     --fresh                                     \
>     -S ${SRC_DIR}/lapackpp                      \
>     -B ${build_dir}/lapackpp-tuolumne-mi300a-build \
>     -DCMAKE_CXX_STANDARD=20                     \
>     -Dgpu_backend=hip                           \
>     -Dbuild_tests=OFF                           \
>     -DBUILD_SHARED_LIBS=OFF                     \
>     -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON      \
>     -DCMAKE_INSTALL_PREFIX=${SW_DIR}/lapackpp-2024.05.31
> cmake \
>     --build ${build_dir}/lapackpp-tuolumne-mi300a-build \
>     --target install                                 \
>     --parallel ${build_procs}
> rm -rf ${build_dir}/lapackpp-tuolumne-mi300a-build

> # PETSC
> if [ -d ${SRC_DIR}/petsc ]
> then
>   cd ${SRC_DIR}/petsc
>   git fetch --prune
>   git checkout v3.24.0
>   cd -
> else
>   git clone -b v3.24.0 https://gitlab.com/petsc/petsc.git ${SRC_DIR}/petsc
> fi
> cd ${SRC_DIR}/petsc
> ./configure               \
>     COPTFLAGS="-g -O3"    \
>     FOPTFLAGS="-g -O3"    \
>     CXXOPTFLAGS="-g -O2"  \
>     HIPOPTFLAGS="-g -O3"  \
>     LDFLAGS+="${LDFLAGS}" \
>     --prefix=${SW_DIR}/petsc-3.24.0  \
>     --with-batch                     \
>     --with-cmake=1                   \
>     --with-cuda=0                    \
>     --with-hip=1                     \
>     --with-hip-dir=${ROCM_PATH}      \
>     --with-fortran-bindings=0        \
>     --with-fftw=0                    \
>     --download-kokkos                \
>     --download-kokkos-kernels        \
>     --with-make-np=${build_procs}    \
>     --with-mpi-dir=${MPICH_DIR}      \
>     --with-clean=1                   \
>     --with-debugging=0               \
>     --with-x=0                       \
>     --with-zlib=1
> make all
> make install
> cd -

> # Python ######################################################################
> #
> # sometimes, the Tuolumne PIP Index is down
> export PIP_EXTRA_INDEX_URL="https://pypi.org/simple"

> python3 -m pip install --upgrade pip
> # python3 -m pip cache purge || true  # Cache disabled on system
> rm -rf ${SW_DIR}/venvs/warpx-tuolumne-mi300a
> python3 -m venv ${SW_DIR}/venvs/warpx-tuolumne-mi300a
> source ${SW_DIR}/venvs/warpx-tuolumne-mi300a/bin/activate
> python3 -m pip install --upgrade pip
> python3 -m pip install --upgrade build
> python3 -m pip install --upgrade packaging
> python3 -m pip install --upgrade wheel
> python3 -m pip install --upgrade setuptools[core]
> python3 -m pip install --upgrade "cython>=3.0"
> python3 -m pip install --upgrade numpy
> python3 -m pip install --upgrade pandas
> python3 -m pip install --upgrade scipy
> python3 -m pip install --upgrade mpi4py --no-cache-dir --no-build-isolation --no-binary mpi4py
> python3 -m pip install --upgrade openpmd-api
> python3 -m pip install --upgrade openpmd-viewer
> python3 -m pip install --upgrade matplotlib
> python3 -m pip install --upgrade yt
> # install or update WarpX dependencies such as picmistandard
> python3 -m pip install --upgrade -r ${SRC_DIR}/warpx/requirements.txt
> # cupy for ROCm
> #   https://docs.cupy.dev/en/stable/install.html#building-cupy-for-rocm-from-source
> #   https://docs.cupy.dev/en/stable/install.html#using-cupy-on-amd-gpu-experimental
> #   https://github.com/cupy/cupy/issues/7830
> #   https://github.com/cupy/cupy/pull/8457
> #   https://github.com/cupy/cupy/pull/8319
> #python3 -m pip install --upgrade "cython<3"
> #HIPCC=${CXX} \
> #CXXFLAGS="-I${ROCM_PATH}/include/hipblas -I${ROCM_PATH}/include/hipsparse -I${ROCM_PATH}/include/hipfft -I${ROCM_PATH}/include/rocsolver -I${ROCM_PATH}/include/rccl -I${ROCM_PATH}/include/thrust" \
> #CUPY_INSTALL_USE_HIP=1  \
> #ROCM_HOME=${ROCM_PATH}  \
> #HCC_AMDGPU_TARGET=${AMREX_AMD_ARCH}  \
> #  python3 -m pip install -v cupy
> #python3 -m pip install --upgrade "cython>=3"


> # for ML dependencies, see install_mi300a_ml.sh

> # remove build temporary directory
> rm -rf ${build_dir}
> ```

> ### AI/ML Dependencies (Optional)

> If you plan to run AI/ML workflows depending on PyTorch et al., run the next step as well.
> This will take a while and should be skipped if not needed.

> ```bash
> bash /p/lustre5/${USER}/tuolumne/src/warpx/Tools/machines/tuolumne-llnl/install_mi300a_ml.sh
> ```

> ### Script Details

> ```bash
> #!/bin/bash
> #
> # Copyright 2024 The WarpX Community
> #
> # This file is part of WarpX.
> #
> # Author: Axel Huebl
> # License: BSD-3-Clause-LBNL

> # Exit on first error encountered #############################################
> #
> set -eu -o pipefail


> # Check: ######################################################################
> #
> #   Was tuolumne_mi300a_warpx.profile sourced and configured correctly?
> #   early access: not yet used!
> #if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your tuolumne_mi300a_warpx.profile file! Please edit its line 2 to continue!"; exit 1; fi


> # Remove old dependencies #####################################################
> #
> SRC_DIR="/p/lustre5/${USER}/tuolumne/src"
> SW_DIR="/p/lustre5/${USER}/tuolumne/warpx/mi300a"

> # remove common user mistakes in python, located in .local instead of a venv
> python3 -m pip uninstall -qqq -y torch 2>/dev/null || true


> # Python ML ###################################################################
> #
> # for basic python dependencies, see install_mi300a_dependencies.sh

> # sometimes, the Lassen PIP Index is down
> export PIP_EXTRA_INDEX_URL="https://pypi.org/simple"

> source ${SW_DIR}/venvs/warpx-tuolumne-mi300a/bin/activate

> python3 -m pip install --upgrade torch torchvision --index-url https://download.pytorch.org/whl/rocm6.1
> python3 -m pip install --upgrade scikit-learn
> python3 -m pip install --upgrade "optimas[all]"
> ```

### CPU

We use system software modules, add environment hints and further dependencies via the file `$HOME/tuolumne_cpu_warpx.profile`.
Create it now:

```bash
cp /p/lustre5/${USER}/tuolumne/src/warpx/Tools/machines/tuolumne-llnl/tuolumne_cpu_warpx.profile.example $HOME/tuolumne_cpu_warpx.profile
```

### Script Details

```bash
# please set your project account
export proj=""  # change me!

# remembers the location of this script
export MY_PROFILE=$(cd $(dirname $BASH_SOURCE) && pwd)"/"$(basename $BASH_SOURCE)
# early access: not yet used
# if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your $MY_PROFILE file! Please edit its line 2 to continue!"; return; fi

# required dependencies
module unload cray-libsci
module load cmake/3.29.2
module load cray-fftw/3.3.10.11

# optional: faster builds
# ccache is system provided
module load ninja/1.10.2

# optional: for QED support with detailed tables
# TODO: no Boost module found

# optional: for openPMD and PSATD+RZ support
SW_DIR="/p/lustre5/${USER}/tuolumne/warpx/cpu"
export CMAKE_PREFIX_PATH=${SW_DIR}/hdf5-1.14.1.2:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/c-blosc-2.15.1:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/adios2-2.10.2:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/blaspp-2024.05.31:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/lapackpp-2024.05.31:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/petsc-3.24.0:$CMAKE_PREFIX_PATH

export LD_LIBRARY_PATH=${SW_DIR}/hdf5-1.14.1.2/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/c-blosc-2.15.1/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/adios2-2.10.2/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/blaspp-2024.05.31/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/lapackpp-2024.05.31/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/petsc-3.24.0/lib:$LD_LIBRARY_PATH

export PATH=${SW_DIR}/adios2-2.10.2/bin:${PATH}
export PATH=${SW_DIR}/hdf5-1.14.1.2/bin:${PATH}

# python
module load cray-python/3.11.7

if [ -d "${SW_DIR}/venvs/warpx-tuolumne-cpu" ]
then
  source ${SW_DIR}/venvs/warpx-tuolumne-cpu/bin/activate
fi

# an alias to request an interactive batch node for one hour
#   for parallel execution, start on the batch node: srun <command>
alias getNode="salloc -N 1 -t 1:00:00"
# an alias to run a command on a batch node for up to 30min
#   usage: runNode <command>
alias runNode="srun -N 1 --ntasks-per-node=4 -t 0:30:00"

# Transparent huge pages on CPU
# https://hpc.llnl.gov/documentation/user-guides/using-el-capitan-systems/introduction-and-quickstart/pro-tips
#export LDFLAGS="${LDFLAGS} -lhugetlbfs"
#export HSA_XNACK=1
#export HUGETLB_MORECORE=yes

# optimize compilation for MI300A CPU part (Zen4)
export CXXFLAGS="-march=znver4"
export CFLAGS="-march=znver4"

# compiler environment hints
export CC=$(which cc)
export CXX=$(which CC)
export FC=$(which ftn)
```

Edit the 2nd line of this script, which sets the `export proj=""` variable.
**Currently, this is unused and can be kept empty.**
Once project allocation becomes required, e.g., if you are member of the project `abcde`, then run `vi $HOME/tuolumne_cpu_warpx.profile`.
Enter the edit mode by typing `i` and edit line 2 to read:

```bash
export proj="abcde"
```

Exit the `vi` editor with `Esc` and then type `:wq` (write & quit).

#### IMPORTANT
Now, and as the first step on future logins to Tuolumne, activate these environment settings:

```bash
source $HOME/tuolumne_cpu_warpx.profile
```

Finally, since Tuolumne does not yet provide software modules for some of our dependencies, install them once:

> ```bash
> bash /p/lustre5/${USER}/tuolumne/src/warpx/Tools/machines/tuolumne-llnl/install_cpu_dependencies.sh
> source /p/lustre5/${USER}/tuolumne/warpx/cpu/venvs/warpx-tuolumne-cpu/bin/activate
> ```

> ### Script Details

> ```bash
> #!/bin/bash
> #
> # Copyright 2024 The WarpX Community
> #
> # This file is part of WarpX.
> #
> # Author: Axel Huebl
> # License: BSD-3-Clause-LBNL

> # Exit on first error encountered #############################################
> #
> set -eu -o pipefail


> # Check: ######################################################################
> #
> #   Was tuolumne_cpu_warpx.profile sourced and configured correctly?
> #   early access: not yet used!
> #if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your tuolumne_cpu_warpx.profile file! Please edit its line 2 to continue!"; exit 1; fi


> # Remove old dependencies #####################################################
> #
> SRC_DIR="/p/lustre5/${USER}/tuolumne/src"
> SW_DIR="/p/lustre5/${USER}/tuolumne/warpx/cpu"
> rm -rf ${SW_DIR}
> mkdir -p ${SW_DIR}

> # remove common user mistakes in python, located in .local instead of a venv
> python3 -m pip uninstall -qq -y pywarpx
> python3 -m pip uninstall -qq -y warpx
> python3 -m pip uninstall -qqq -y mpi4py 2>/dev/null || true


> # General extra dependencies ##################################################
> #

> # tmpfs build directory: avoids issues often seen with $HOME and is faster
> build_dir=$(mktemp -d)
> build_procs=24

> # C-Blosc2 (I/O compression)
> if [ -d ${SRC_DIR}/c-blosc2 ]
> then
>   cd ${SRC_DIR}/c-blosc2
>   git fetch --prune
>   git checkout v2.15.1
>   cd -
> else
>   git clone -b v2.15.1 https://github.com/Blosc/c-blosc2.git ${SRC_DIR}/c-blosc2
> fi
> cmake \
>     --fresh                        \
>     -S ${SRC_DIR}/c-blosc2         \
>     -B ${build_dir}/c-blosc2-build \
>     -DBUILD_SHARED_LIBS=OFF        \
>     -DBUILD_TESTS=OFF              \
>     -DBUILD_BENCHMARKS=OFF         \
>     -DBUILD_EXAMPLES=OFF           \
>     -DBUILD_FUZZERS=OFF            \
>     -DBUILD_STATIC=OFF             \
>     -DDEACTIVATE_AVX2=OFF          \
>     -DDEACTIVATE_AVX512=OFF        \
>     -DWITH_SANITIZER=OFF           \
>     -DCMAKE_INSTALL_PREFIX=${SW_DIR}/c-blosc-2.15.1
> cmake \
>     --build ${build_dir}/c-blosc2-build \
>     --target install                    \
>     --parallel ${build_procs}
> rm -rf ${build_dir}/c-blosc2-build

> # HDF5
> if [ -d ${SRC_DIR}/hdf5 ]
> then
>   cd ${SRC_DIR}/hdf5
>   git fetch --prune
>   git checkout hdf5-1_14_1-2
>   cd -
> else
>   git clone -b hdf5-1_14_1-2 https://github.com/HDFGroup/hdf5.git ${SRC_DIR}/hdf5
> fi
> cmake \
>     --fresh                      \
>     -S ${SRC_DIR}/hdf5           \
>     -B ${build_dir}/hdf5-build   \
>     -DBUILD_SHARED_LIBS=OFF        \
>     -DBUILD_TESTING=OFF          \
>     -DHDF5_ENABLE_PARALLEL=ON    \
>     -DCMAKE_INSTALL_PREFIX=${SW_DIR}/hdf5-1.14.1.2
> cmake \
>     --build ${build_dir}/hdf5-build \
>     --target install                \
>     --parallel ${build_procs}
> rm -rf ${build_dir}/hdf5-build

> # ADIOS2
> if [ -d ${SRC_DIR}/adios2 ]
> then
>   cd ${SRC_DIR}/adios2
>   git fetch --prune
>   git checkout v2.10.2
>   cd -
> else
>   git clone -b v2.10.2 https://github.com/ornladios/ADIOS2.git ${SRC_DIR}/adios2
> fi
> cmake \
>     --fresh                      \
>     -S ${SRC_DIR}/adios2         \
>     -B ${build_dir}/adios2-build \
>     -DADIOS2_USE_Blosc2=ON       \
>     -DADIOS2_USE_Campaign=OFF    \
>     -DADIOS2_USE_Fortran=OFF     \
>     -DADIOS2_USE_Python=OFF      \
>     -DADIOS2_USE_ZeroMQ=OFF      \
>     -DBUILD_SHARED_LIBS=OFF      \
>     -DCMAKE_INSTALL_PREFIX=${SW_DIR}/adios2-2.10.2
> cmake \
>     --build ${build_dir}/adios2-build \
>     --target install                  \
>     --parallel ${build_procs}
> rm -rf ${build_dir}/adios2-build

> # BLAS++ (for PSATD+RZ)
> if [ -d ${SRC_DIR}/blaspp ]
> then
>   cd ${SRC_DIR}/blaspp
>   git fetch --prune
>   git checkout v2024.05.31
>   cd -
> else
>   git clone -b v2024.05.31 https://github.com/icl-utk-edu/blaspp.git ${SRC_DIR}/blaspp
> fi
> cmake \
>     --fresh \
>     -S ${SRC_DIR}/blaspp                      \
>     -B ${build_dir}/blaspp-tuolumne-cpu-build \
>     -Duse_openmp=ON                           \
>     -Dgpu_backend=OFF                         \
>     -DBUILD_SHARED_LIBS=OFF                   \
>     -DCMAKE_CXX_STANDARD=20                   \
>     -DCMAKE_INSTALL_PREFIX=${SW_DIR}/blaspp-2024.05.31
> cmake \
>     --build ${build_dir}/blaspp-tuolumne-cpu-build \
>     --target install                               \
>     --parallel ${build_procs}
> rm -rf ${build_dir}/blaspp-tuolumne-cpu-build

> # LAPACK++ (for PSATD+RZ)
> if [ -d ${SRC_DIR}/lapackpp ]
> then
>   cd ${SRC_DIR}/lapackpp
>   git fetch --prune
>   git checkout v2024.05.31
>   cd -
> else
>   git clone -b v2024.05.31 https://github.com/icl-utk-edu/lapackpp.git ${SRC_DIR}/lapackpp
> fi
> cmake \
>     --fresh                                     \
>     -S ${SRC_DIR}/lapackpp                      \
>     -B ${build_dir}/lapackpp-tuolumne-cpu-build \
>     -DCMAKE_CXX_STANDARD=20                     \
>     -Dgpu_backend=OFF                           \
>     -Dbuild_tests=OFF                           \
>     -DBUILD_SHARED_LIBS=OFF                     \
>     -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON      \
>     -DCMAKE_INSTALL_PREFIX=${SW_DIR}/lapackpp-2024.05.31
> cmake \
>     --build ${build_dir}/lapackpp-tuolumne-cpu-build \
>     --target install                                 \
>     --parallel ${build_procs}
> rm -rf ${build_dir}/lapackpp-tuolumne-cpu-build

> # PETSC
> if [ -d ${SRC_DIR}/petsc ]
> then
>   cd ${SRC_DIR}/petsc
>   git fetch --prune
>   git checkout v3.24.0
>   cd -
> else
>   git clone -b v3.24.0 https://gitlab.com/petsc/petsc.git ${SRC_DIR}/petsc
> fi
> cd ${SRC_DIR}/petsc
> ./configure    \
>     CC=${CC}   \
>     CXX=${CXX} \
>     FC=${FC}   \
>     COPTFLAGS="-g -O3"   \
>     FOPTFLAGS="-g -O3"   \
>     CXXOPTFLAGS="-g -O2" \
>     --prefix=${SW_DIR}/petsc-3.24.0  \
>     --with-batch                     \
>     --with-cmake=1                   \
>     --with-cuda=0                    \
>     --with-hip=0                     \
>     --with-fortran-bindings=0        \
>     --with-fftw=1                    \
>     --with-fftw-dir=${FFTW_ROOT}     \
>     --with-make-np=${build_procs}    \
>     ---with-openmp-kernels=1         \
>     --with-clean=1                   \
>     --with-debugging=0               \
>     --with-x=0                       \
>     --with-zlib=1
> make all
> make install
> cd -

> # Python ######################################################################
> #
> # sometimes, the Tuolumne PIP Index is down
> export PIP_EXTRA_INDEX_URL="https://pypi.org/simple"

> python3 -m pip install --upgrade pip
> # python3 -m pip cache purge || true  # Cache disabled on system
> rm -rf ${SW_DIR}/venvs/warpx-tuolumne-cpu
> python3 -m venv ${SW_DIR}/venvs/warpx-tuolumne-cpu
> source ${SW_DIR}/venvs/warpx-tuolumne-cpu/bin/activate
> python3 -m pip install --upgrade pip
> python3 -m pip install --upgrade build
> python3 -m pip install --upgrade packaging
> python3 -m pip install --upgrade wheel
> python3 -m pip install --upgrade setuptools[core]
> python3 -m pip install --upgrade "cython>=3.0"
> python3 -m pip install --upgrade numpy
> python3 -m pip install --upgrade pandas
> python3 -m pip install --upgrade scipy
> python3 -m pip install --upgrade mpi4py --no-cache-dir --no-build-isolation --no-binary mpi4py
> python3 -m pip install --upgrade openpmd-api
> python3 -m pip install --upgrade openpmd-viewer
> python3 -m pip install --upgrade matplotlib
> python3 -m pip install --upgrade yt
> # install or update WarpX dependencies such as picmistandard
> python3 -m pip install --upgrade -r ${SRC_DIR}/warpx/requirements.txt

> # for ML dependencies, see install_cpu_ml.sh (TODO)

> # remove build temporary directory
> rm -rf ${build_dir}
> ```

<a id="building-tuolumne-compilation"></a>

## Compilation

Use the following [cmake commands](../cmake.md#install-build-cmake) to compile the application executable:

### GPU

```bash
cd /p/lustre5/${USER}/tuolumne/src/warpx

cmake --fresh -S . -B build_tuolumne -DWarpX_COMPUTE=HIP -DWarpX_FFT=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_tuolumne -j 24
```

The WarpX application executables are now in `/p/lustre5/${USER}/tuolumne/src/warpx/build_tuolumne/bin/`.
Additionally, the following commands will install WarpX as a Python module:

```bash
cmake --fresh -S . -B build_tuolumne_py -DWarpX_COMPUTE=HIP -DWarpX_FFT=ON -DWarpX_APP=OFF -DWarpX_PYTHON=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_tuolumne_py -j 24 --target pip_install
```

### CPU

```bash
cd /p/lustre5/${USER}/tuolumne/src/warpx

cmake --fresh -S . -B build_tuolumne_cpu -DWarpX_COMPUTE=OMP -DWarpX_FFT=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_tuolumne_cpu -j 24
```

The WarpX application executables are now in `/p/lustre5/${USER}/tuolumne/src/warpx/build_tuolumne_cpu/bin/`.
Additionally, the following commands will install WarpX as a Python module:

```bash
cmake --fresh -S . -B build_tuolumne_cpu_py -DWarpX_COMPUTE=OMP -DWarpX_FFT=ON -DWarpX_APP=OFF -DWarpX_PYTHON=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_tuolumne_cpu_py -j 24 --target pip_install
```

Now, you can [submit tuolumne compute jobs](#running-cpp-tuolumne) for WarpX [Python (PICMI) scripts](../../usage/python.md#usage-picmi) ([example scripts](../../usage/examples.md#usage-examples)).
Or, you can use the WarpX executables to submit tuolumne jobs ([example inputs](../../usage/examples.md#usage-examples)).
For executables, you can reference their location in your [job script](#running-cpp-tuolumne) or copy them to a location in `$PROJWORK/$proj/`.

<a id="building-tuolumne-update"></a>

## Update WarpX & Dependencies

If you already installed WarpX in the past and want to update it, start by getting the latest source code:

```bash
cd /p/lustre5/${USER}/tuolumne/src/warpx

# read the output of this command - does it look ok?
git status

# get the latest WarpX source code
git fetch
git pull

# read the output of these commands - do they look ok?
git status
git log     # press q to exit
```

And, if needed,

- [update the tuolumne_mi300a_warpx.profile file](#building-tuolumne-preparation),
- log out and into the system, activate the now updated environment profile as usual,
- [execute the dependency install scripts](#building-tuolumne-preparation).

As a last step [rebuild WarpX](#building-tuolumne-compilation).

<a id="running-cpp-tuolumne"></a>

## Running

<a id="running-cpp-tuolumne-mi300a-apus"></a>

### MI300A APUs (128GB)

[Each compute node](https://hpc.llnl.gov/documentation/user-guides/using-el-capitan-systems/introduction-and-quickstart/pro-tips) is divided into 4 sockets, each with:

* 1 MI300A APU (incl. 1 GPU),
* 21 available user CPU cores, with 3 cores reserved for the OS (2 hardware threads per core)
* 128GB HBM3 memory (a single NUMA domain)

The batch script below can be used to run a WarpX simulation on 1 node with 4 APUs on the supercomputer Tuolumne at LLNL.
Replace descriptions between chevrons `<>` by relevant values, for instance `<input file>` could be `plasma_mirror_inputs`.
WarpX runs with one MPI rank per GPU and uses 21 (of 24) CPU cores (3 are reserved for the system).

The batch script below also [sends WarpX a signal](../../usage/parameters.md#running-cpp-parameters-signal) when the simulations gets close to the walltime of the job, to shut down cleanly.
Adjust the `FLUX_WT_SIG` and `WARPX_WT` to modify or disable this behavior as needed.

### GPU

```bash
#!/bin/bash

# Copyright 2025 The WarpX Community
#
# This file is part of WarpX.
#
# Authors: Axel Huebl, Andreas Kemp
# License: BSD-3-Clause-LBNL

### Flux directives ###
#flux: --setattr=bank=mstargt
#flux: --job-name=hemi
#flux: --nodes=16
#flux: --time-limit=360s
#flux: --queue=pbatch
#              pdebug
#flux: --exclusive
#flux: --error=WarpX.e{{id}}
#flux: --output=WarpX.o{{id}}

# Not yet tested: Transparent huge pages on CPU
# https://hpc.llnl.gov/documentation/user-guides/using-el-capitan-systems/introduction-and-quickstart/pro-tips
#      --setattr=thp=always

# executable & inputs file or python interpreter & PICMI script here
EXE="./warpx.2d"
INPUTS="./inputs_hist_10.input"

# clean shutdown close to walltime (or checkpoint)
# https://warpx.readthedocs.io/en/latest/usage/parameters.html#signal-handling
FLUX_WT_SIG="--signal=SIGUSR1@120s"
WARPX_WT="warpx.break_signals=USR1"

# enviroment setup
if [[ -z "${MY_PROFILE}" ]]; then
    echo "WARNING: FORGOT TO"
    echo "   source $HOME/tuolumne_mi300a_warpx.profile"
    echo "before submission. Doing that now."

    source $HOME/tuolumne_mi300a_warpx.profile
fi

# pin to closest NIC to GPU
export MPICH_OFI_NIC_POLICY=GPU

# Not yet tested: Transparent huge pages on CPU
# https://hpc.llnl.gov/documentation/user-guides/using-el-capitan-systems/introduction-and-quickstart/pro-tips
#export HSA_XNACK=1
#export HUGETLB_MORECORE=yes

# threads for OpenMP and threaded compressors per MPI rank
#   note: 21 physical cores per socket maximum (system reserves 3)
export OMP_NUM_THREADS=21

# GPU-aware MPI optimizations
GPU_AWARE_MPI="amrex.use_gpu_aware_mpi=1"

# start MPI parallel processes
NNODES=$(flux resource list -s up -no {nnodes})
flux run ${FLUX_WT_SIG} --exclusive --nodes=${NNODES} \
  --tasks-per-node=4 \
  ${EXE} ${INPUTS} \
  ${GPU_AWARE_MPI} ${WARPX_WT} \
  > output.txt
```

To run a simulation, copy the lines above to a file `tuolumne_mi300a.flux` and run

```bash
flux batch tuolumne_mi300a.flux
```

### CPU

```bash
#!/bin/bash

# Copyright 2025 The WarpX Community
#
# This file is part of WarpX.
#
# Authors: Axel Huebl, Andreas Kemp
# License: BSD-3-Clause-LBNL

### Flux directives ###
#flux: --setattr=bank=mstargt
#flux: --job-name=hemi
#flux: --nodes=16
#flux: --time-limit=360s
#flux: --queue=pbatch
#              pdebug
#flux: --exclusive
#flux: --error=WarpX.e{{id}}
#flux: --output=WarpX.o{{id}}

# Not yet tested: Transparent huge pages on CPU
# https://hpc.llnl.gov/documentation/user-guides/using-el-capitan-systems/introduction-and-quickstart/pro-tips
#      --setattr=thp=always

# executable & inputs file or python interpreter & PICMI script here
EXE="./warpx.2d"
INPUTS="./inputs_hist_10.input"

# clean shutdown close to walltime (or checkpoint)
# https://warpx.readthedocs.io/en/latest/usage/parameters.html#signal-handling
FLUX_WT_SIG="--signal=SIGUSR1@120s"
WARPX_WT="warpx.break_signals=USR1"

# enviroment setup
if [[ -z "${MY_PROFILE}" ]]; then
    echo "WARNING: FORGOT TO"
    echo "   source $HOME/tuolumne_cpu_warpx.profile"
    echo "before submission. Doing that now."

    source $HOME/tuolumne_cpu_warpx.profile
fi

# pin to closest NIC to APU
#export MPICH_OFI_NIC_POLICY=APU

# Not yet tested: Transparent huge pages on CPU
# https://hpc.llnl.gov/documentation/user-guides/using-el-capitan-systems/introduction-and-quickstart/pro-tips
#export HSA_XNACK=1
#export HUGETLB_MORECORE=yes

# threads for OpenMP and threaded compressors per MPI rank
#   note: 21 physical cores per socket maximum (system reserves 3)
export OMP_NUM_THREADS=21

# start MPI parallel processes
NNODES=$(flux resource list -s up -no {nnodes})
flux run ${FLUX_WT_SIG} --exclusive --nodes=${NNODES} \
  --tasks-per-node=4           \
  ${EXE} ${INPUTS} ${WARPX_WT} \
  > output.txt
```

To run a simulation, copy the lines above to a file `tuolumne_cpu.flux` and run

```bash
flux batch tuolumne_cpu.flux
```

to submit the job.
