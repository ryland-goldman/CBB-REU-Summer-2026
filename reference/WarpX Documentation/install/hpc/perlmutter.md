<a id="building-perlmutter"></a>

# Perlmutter (NERSC)

The [Perlmutter cluster](https://docs.nersc.gov/systems/perlmutter/) is located at NERSC.

## Introduction

If you are new to this system, **please see the following resources**:

* [NERSC user guide](https://docs.nersc.gov/)
* Batch system: [Slurm](https://docs.nersc.gov/systems/perlmutter/#running-jobs)
* [Jupyter service](https://jupyter.nersc.gov) ([documentation](https://docs.nersc.gov/services/jupyter/))
* [Filesystems](https://docs.nersc.gov/filesystems/):
  * `$HOME`: per-user directory, use only for inputs, source and scripts; backed up (40GB)
  * `${CFS}/m3239/`: [community file system](https://docs.nersc.gov/filesystems/community/) for users in the project `m3239` (or equivalent); moderate performance (20TB default)
  * `$PSCRATCH`: per-user [production directory](https://docs.nersc.gov/filesystems/perlmutter-scratch/); very fast for parallel jobs; purged every 8 weeks (20TB default)

<a id="building-perlmutter-preparation"></a>

## Preparation

Use the following commands to download the WarpX source code:

```bash
git clone https://github.com/BLAST-WarpX/warpx.git $HOME/src/warpx
```

On Perlmutter, you can run either on GPU nodes with fast A100 GPUs (recommended) or CPU nodes.

### A100 GPUs

We use system software modules, add environment hints and further dependencies via the file `$HOME/perlmutter_gpu_warpx.profile`.
Create it now:

```bash
cp $HOME/src/warpx/Tools/machines/perlmutter-nersc/perlmutter_gpu_warpx.profile.example $HOME/perlmutter_gpu_warpx.profile
```

### Script Details

```bash
# please set your project account
export proj=""  # change me! GPU projects must end in "..._g"

# remembers the location of this script
export MY_PROFILE=$(cd $(dirname $BASH_SOURCE) && pwd)"/"$(basename $BASH_SOURCE)
if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your $MY_PROFILE file! Please edit its line 2 to continue!"; return; fi

# required dependencies
module load gpu
module load PrgEnv-gnu
module load craype
module load craype-x86-milan
module load craype-accel-nvidia80
module load cudatoolkit
module load gcc-native/13.2  # default gcc-native/14 breaks pybind11 builds with NVCC 12.9.41
module load cmake/3.30.2

# missing modules installed here
export SW_DIR=${PSCRATCH}/storage/sw/warpx/perlmutter/gpu

# optional: for QED support with detailed tables
export CMAKE_PREFIX_PATH=${SW_DIR}/boost-1.82.0:${CMAKE_PREFIX_PATH}
export LD_LIBRARY_PATH=${SW_DIR}/boost-1.82.0/lib:${LD_LIBRARY_PATH}

# optional: for openPMD and PSATD+RZ support
module load cray-hdf5-parallel/1.12.2.9
export CMAKE_PREFIX_PATH=${SW_DIR}/c-blosc-1.21.1:${CMAKE_PREFIX_PATH}
export CMAKE_PREFIX_PATH=${SW_DIR}/adios2-2.10.2:${CMAKE_PREFIX_PATH}
export CMAKE_PREFIX_PATH=${SW_DIR}/blaspp-2024.05.31:${CMAKE_PREFIX_PATH}
export CMAKE_PREFIX_PATH=${SW_DIR}/lapackpp-2024.05.31:${CMAKE_PREFIX_PATH}

export LD_LIBRARY_PATH=${SW_DIR}/c-blosc-1.21.1/lib64:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${SW_DIR}/adios2-2.10.2/lib64:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${SW_DIR}/blaspp-2024.05.31/lib64:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${SW_DIR}/lapackpp-2024.05.31/lib64:${LD_LIBRARY_PATH}

export PATH=${SW_DIR}/adios2-2.10.2/bin:${PATH}

# optional: CCache
export PATH=${SW_DIR}/ccache-4.10.2:$PATH

# optional: for Python bindings or libEnsemble
module load cray-python/3.11.5

if [ -d "${SW_DIR}/venvs/warpx-gpu" ]
then
  source ${SW_DIR}/venvs/warpx-gpu/bin/activate
fi

# an alias to request an interactive batch node for one hour
#   for parallel execution, start on the batch node: srun <command>
alias getNode="salloc -N 1 --ntasks-per-node=4 -t 1:00:00 -q interactive -C gpu --gpu-bind=none -c 32 -G 4 -A $proj"
# an alias to run a command on a batch node for up to 30min
#   usage: runNode <command>
alias runNode="srun -N 1 --ntasks-per-node=4 -t 0:30:00 -q interactive -C gpu --gpu-bind=none -c 32 -G 4 -A $proj"

# necessary to use CUDA-Aware MPI and run a job
export CRAY_ACCEL_TARGET=nvidia80

# optimize CUDA compilation for A100
export AMREX_CUDA_ARCH=8.0

# optimize CPU microarchitecture for AMD EPYC 3rd Gen (Milan/Zen3)
# note: the cc/CC/ftn wrappers below add those
export CXXFLAGS="-march=znver3"
export CFLAGS="-march=znver3"

# compiler environment hints
export CC=cc
export CXX=CC
export FC=ftn
export CUDACXX=$(which nvcc)
export CUDAHOSTCXX=CC
```

Edit the 2nd line of this script, which sets the `export proj=""` variable.
Perlmutter GPU projects must end in `..._g`.
For example, if you are member of the project `m3239`, then run `nano $HOME/perlmutter_gpu_warpx.profile` and edit line 2 to read:

```bash
export proj="m3239_g"
```

Exit the `nano` editor with `Ctrl` + `O` (save) and then `Ctrl` + `X` (exit).

#### IMPORTANT
Now, and as the first step on future logins to Perlmutter, activate these environment settings:

```bash
source $HOME/perlmutter_gpu_warpx.profile
```

Finally, since Perlmutter does not yet provide software modules for some of our dependencies, install them once:

```bash
bash $HOME/src/warpx/Tools/machines/perlmutter-nersc/install_gpu_dependencies.sh
source ${PSCRATCH}/storage/sw/warpx/perlmutter/gpu/venvs/warpx-gpu/bin/activate
```

### Script Details

```bash
#!/bin/bash
#
# Copyright 2023 The WarpX Community
#
# This file is part of WarpX.
#
# Author: Axel Huebl
# License: BSD-3-Clause-LBNL

# Exit on first error encountered #############################################
#
set -eu -o pipefail


# Check: ######################################################################
#
#   Was perlmutter_gpu_warpx.profile sourced and configured correctly?
if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your perlmutter_gpu_warpx.profile file! Please edit its line 2 to continue!"; exit 1; fi


# Check $proj variable is correct and has a corresponding CFS directory #######
#
if [ ! -d "${CFS}/${proj%_g}/" ]
then
    echo "WARNING: The directory ${CFS}/${proj%_g}/ does not exist!"
    echo "Is the \$proj environment variable of value \"$proj\" correctly set? "
    echo "Please edit line 2 of your perlmutter_gpu_warpx.profile file to continue!"
    exit
fi


# Remove old dependencies #####################################################
#
SW_DIR="${PSCRATCH}/storage/sw/warpx/perlmutter/gpu"
rm -rf ${SW_DIR}
mkdir -p ${SW_DIR}

# remove common user mistakes in python, located in .local instead of a venv
python3 -m pip uninstall -qq -y pywarpx
python3 -m pip uninstall -qq -y warpx
python3 -m pip uninstall -qqq -y mpi4py 2>/dev/null || true


# General extra dependencies ##################################################
#

# build parallelism
PARALLEL=16

# tmpfs build directory: avoids issues often seen with $HOME and is faster
build_dir=$(mktemp -d)

# CCache
curl -Lo ccache.tar.xz https://github.com/ccache/ccache/releases/download/v4.10.2/ccache-4.10.2-linux-x86_64.tar.xz
tar -xf ccache.tar.xz
mv ccache-4.10.2-linux-x86_64 ${SW_DIR}/ccache-4.10.2
rm -rf ccache.tar.xz

# Boost (QED tables)
rm -rf $HOME/src/boost-temp
mkdir -p $HOME/src/boost-temp
curl -Lo $HOME/src/boost-temp/boost.tar.gz https://archives.boost.io/release/1.82.0/source/boost_1_82_0.tar.gz
tar -xzf $HOME/src/boost-temp/boost.tar.gz -C $HOME/src/boost-temp
cd $HOME/src/boost-temp/boost_1_82_0
./bootstrap.sh --with-libraries=math --prefix=${SW_DIR}/boost-1.82.0
./b2 cxxflags="-std=c++20" install -j ${PARALLEL}
cd -
rm -rf $HOME/src/boost-temp

# c-blosc (I/O compression)
if [ -d $HOME/src/c-blosc ]
then
  cd $HOME/src/c-blosc
  git fetch --prune
  git checkout v1.21.1
  cd -
else
  git clone -b v1.21.1 https://github.com/Blosc/c-blosc.git $HOME/src/c-blosc
fi
rm -rf $HOME/src/c-blosc-pm-gpu-build
cmake -S $HOME/src/c-blosc -B ${build_dir}/c-blosc-pm-gpu-build -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF -DDEACTIVATE_AVX2=OFF -DCMAKE_INSTALL_PREFIX=${SW_DIR}/c-blosc-1.21.1
cmake --build ${build_dir}/c-blosc-pm-gpu-build --target install --parallel ${PARALLEL}
rm -rf ${build_dir}/c-blosc-pm-gpu-build

# ADIOS2
if [ -d $HOME/src/adios2 ]
then
  cd $HOME/src/adios2
  git fetch --prune
  git checkout v2.10.2
  cd -
else
  git clone -b v2.10.2 https://github.com/ornladios/ADIOS2.git $HOME/src/adios2
fi
rm -rf $HOME/src/adios2-pm-gpu-build
cmake -S $HOME/src/adios2 -B ${build_dir}/adios2-pm-gpu-build -DADIOS2_USE_Blosc=ON -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DADIOS2_USE_ZeroMQ=OFF -DCMAKE_INSTALL_PREFIX=${SW_DIR}/adios2-2.10.2
cmake --build ${build_dir}/adios2-pm-gpu-build --target install -j ${PARALLEL}
rm -rf ${build_dir}/adios2-pm-gpu-build

# BLAS++ (for PSATD+RZ)
if [ -d $HOME/src/blaspp ]
then
  cd $HOME/src/blaspp
  git fetch --prune
  git checkout v2024.05.31
  cd -
else
  git clone -b v2024.05.31 https://github.com/icl-utk-edu/blaspp.git $HOME/src/blaspp
fi
rm -rf $HOME/src/blaspp-pm-gpu-build
CXX=$(which CC) cmake -S $HOME/src/blaspp -B ${build_dir}/blaspp-pm-gpu-build -Duse_openmp=OFF -Dgpu_backend=cuda -DCMAKE_CXX_STANDARD=20 -DCMAKE_INSTALL_PREFIX=${SW_DIR}/blaspp-2024.05.31
cmake --build ${build_dir}/blaspp-pm-gpu-build --target install --parallel ${PARALLEL}
rm -rf ${build_dir}/blaspp-pm-gpu-build

# LAPACK++ (for PSATD+RZ)
if [ -d $HOME/src/lapackpp ]
then
  cd $HOME/src/lapackpp
  git fetch --prune
  git checkout v2024.05.31
  cd -
else
  git clone -b v2024.05.31 https://github.com/icl-utk-edu/lapackpp.git $HOME/src/lapackpp
fi
rm -rf $HOME/src/lapackpp-pm-gpu-build
CXX=$(which CC) CXXFLAGS="-DLAPACK_FORTRAN_ADD_" cmake -S $HOME/src/lapackpp -B ${build_dir}/lapackpp-pm-gpu-build -DCMAKE_CXX_STANDARD=20 -Dbuild_tests=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX=${SW_DIR}/lapackpp-2024.05.31
cmake --build ${build_dir}/lapackpp-pm-gpu-build --target install --parallel ${PARALLEL}
rm -rf ${build_dir}/lapackpp-pm-gpu-build

# Python ######################################################################
#
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade virtualenv
python3 -m pip cache purge
rm -rf ${SW_DIR}/venvs/warpx-gpu
python3 -m venv ${SW_DIR}/venvs/warpx-gpu
source ${SW_DIR}/venvs/warpx-gpu/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade build
python3 -m pip install --upgrade packaging
python3 -m pip install --upgrade wheel
python3 -m pip install --upgrade setuptools[core]
python3 -m pip install --upgrade cython
python3 -m pip install --upgrade numpy
python3 -m pip install --upgrade pandas
python3 -m pip install --upgrade scipy
MPICC="cc -target-accel=nvidia80 -shared" python3 -m pip install --upgrade mpi4py --no-cache-dir --no-build-isolation --no-binary mpi4py
python3 -m pip install --upgrade openpmd-api
python3 -m pip install --upgrade matplotlib
python3 -m pip install --upgrade yt
# install or update WarpX dependencies
python3 -m pip install --upgrade -r $HOME/src/warpx/requirements.txt
python3 -m pip install --upgrade cupy-cuda12x  # CUDA 12 compatible wheel
# optimas (based on libEnsemble & ax->botorch->gpytorch->pytorch)
python3 -m pip install --upgrade torch  # CUDA 12 compatible wheel
python3 -m pip install --upgrade optimas[all]
python3 -m pip install --upgrade lasy

# remove build temporary directory
rm -rf ${build_dir}
```

### CPU Nodes

We use system software modules, add environment hints and further dependencies via the file `$HOME/perlmutter_cpu_warpx.profile`.
Create it now:

```bash
cp $HOME/src/warpx/Tools/machines/perlmutter-nersc/perlmutter_cpu_warpx.profile.example $HOME/perlmutter_cpu_warpx.profile
```

### Script Details

```bash
# please set your project account
export proj=""  # change me!

# remembers the location of this script
export MY_PROFILE=$(cd $(dirname $BASH_SOURCE) && pwd)"/"$(basename $BASH_SOURCE)
if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your $MY_PROFILE file! Please edit its line 2 to continue!"; return; fi

# required dependencies
module load cpu
module load cmake/3.30.2
module load cray-fftw/3.3.10.6

# missing modules installed here
export SW_DIR=${PSCRATCH}/storage/sw/warpx/perlmutter/cpu

# optional: for QED support with detailed tables
export CMAKE_PREFIX_PATH=${SW_DIR}/boost-1.82.0:${CMAKE_PREFIX_PATH}
export LD_LIBRARY_PATH=${SW_DIR}/boost-1.82.0/lib:${LD_LIBRARY_PATH}

# optional: for openPMD and PSATD+RZ support
module load cray-hdf5-parallel/1.12.2.9
export CMAKE_PREFIX_PATH=${SW_DIR}/c-blosc-1.21.1:${CMAKE_PREFIX_PATH}
export CMAKE_PREFIX_PATH=${SW_DIR}/adios2-2.10.2:${CMAKE_PREFIX_PATH}
export CMAKE_PREFIX_PATH=${SW_DIR}/blaspp-2024.05.31:${CMAKE_PREFIX_PATH}
export CMAKE_PREFIX_PATH=${SW_DIR}/lapackpp-2024.05.31:${CMAKE_PREFIX_PATH}

export LD_LIBRARY_PATH=${SW_DIR}/c-blosc-1.21.1/lib64:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${SW_DIR}/adios2-2.10.2/lib64:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${SW_DIR}/blaspp-2024.05.31/lib64:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${SW_DIR}/lapackpp-2024.05.31/lib64:${LD_LIBRARY_PATH}

export PATH=${SW_DIR}/adios2-2.10.2/bin:${PATH}

# optional: CCache
export PATH=${SW_DIR}/ccache-4.10.2:$PATH

# optional: for Python bindings or libEnsemble
module load cray-python/3.11.5

if [ -d "${SW_DIR}/venvs/warpx-cpu" ]
then
  source ${SW_DIR}/venvs/warpx-cpu/bin/activate
fi

# an alias to request an interactive batch node for one hour
#   for parallel execution, start on the batch node: srun <command>
alias getNode="salloc --nodes 1 --qos interactive --time 01:00:00 --constraint cpu --account=$proj"
# an alias to run a command on a batch node for up to 30min
#   usage: runNode <command>
alias runNode="srun --nodes 1 --qos interactive --time 01:00:00 --constraint cpu $proj"

# optimize CPU microarchitecture for AMD EPYC 3rd Gen (Milan/Zen3)
# note: the cc/CC/ftn wrappers below add those
export CXXFLAGS="-march=znver3"
export CFLAGS="-march=znver3"

# compiler environment hints
export CC=cc
export CXX=CC
export FC=ftn
```

Edit the 2nd line of this script, which sets the `export proj=""` variable.
For example, if you are member of the project `m3239`, then run `nano $HOME/perlmutter_cpu_warpx.profile` and edit line 2 to read:

```bash
export proj="m3239"
```

Exit the `nano` editor with `Ctrl` + `O` (save) and then `Ctrl` + `X` (exit).

#### IMPORTANT
Now, and as the first step on future logins to Perlmutter, activate these environment settings:

```bash
source $HOME/perlmutter_cpu_warpx.profile
```

Finally, since Perlmutter does not yet provide software modules for some of our dependencies, install them once:

```bash
bash $HOME/src/warpx/Tools/machines/perlmutter-nersc/install_cpu_dependencies.sh
source ${PSCRATCH}/storage/sw/warpx/perlmutter/cpu/venvs/warpx-cpu/bin/activate
```

### Script Details

```bash
#!/bin/bash
#
# Copyright 2023 The WarpX Community
#
# This file is part of WarpX.
#
# Author: Axel Huebl
# License: BSD-3-Clause-LBNL

# Exit on first error encountered #############################################
#
set -eu -o pipefail


# Check: ######################################################################
#
#   Was perlmutter_cpu_warpx.profile sourced and configured correctly?
if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your perlmutter_cpu_warpx.profile file! Please edit its line 2 to continue!"; exit 1; fi


# Check $proj variable is correct and has a corresponding CFS directory #######
#
if [ ! -d "${CFS}/${proj}/" ]
then
    echo "WARNING: The directory ${CFS}/${proj}/ does not exist!"
    echo "Is the \$proj environment variable of value \"$proj\" correctly set? "
    echo "Please edit line 2 of your perlmutter_cpu_warpx.profile file to continue!"
    exit
fi


# Remove old dependencies #####################################################
#
SW_DIR="${PSCRATCH}/storage/sw/warpx/perlmutter/cpu"
rm -rf ${SW_DIR}
mkdir -p ${SW_DIR}

# remove common user mistakes in python, located in .local instead of a venv
python3 -m pip uninstall -qq -y pywarpx
python3 -m pip uninstall -qq -y warpx
python3 -m pip uninstall -qqq -y mpi4py 2>/dev/null || true


# General extra dependencies ##################################################
#

# build parallelism
PARALLEL=16

# tmpfs build directory: avoids issues often seen with $HOME and is faster
build_dir=$(mktemp -d)

# CCache
curl -Lo ccache.tar.xz https://github.com/ccache/ccache/releases/download/v4.10.2/ccache-4.10.2-linux-x86_64.tar.xz
tar -xf ccache.tar.xz
mv ccache-4.10.2-linux-x86_64 ${SW_DIR}/ccache-4.10.2
rm -rf ccache.tar.xz

# Boost (QED tables)
rm -rf $HOME/src/boost-temp
mkdir -p $HOME/src/boost-temp
curl -Lo $HOME/src/boost-temp/boost.tar.gz https://archives.boost.io/release/1.82.0/source/boost_1_82_0.tar.gz
tar -xzf $HOME/src/boost-temp/boost.tar.gz -C $HOME/src/boost-temp
cd $HOME/src/boost-temp/boost_1_82_0
./bootstrap.sh --with-libraries=math --prefix=${SW_DIR}/boost-1.82.0
./b2 cxxflags="-std=c++20" install -j ${PARALLEL}
cd -
rm -rf $HOME/src/boost-temp

# c-blosc (I/O compression)
if [ -d $HOME/src/c-blosc ]
then
  cd $HOME/src/c-blosc
  git fetch --prune
  git checkout v1.21.1
  cd -
else
  git clone -b v1.21.1 https://github.com/Blosc/c-blosc.git $HOME/src/c-blosc
fi
rm -rf $HOME/src/c-blosc-pm-cpu-build
cmake -S $HOME/src/c-blosc -B ${build_dir}/c-blosc-pm-cpu-build -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF -DDEACTIVATE_AVX2=OFF -DCMAKE_INSTALL_PREFIX=${SW_DIR}/c-blosc-1.21.1
cmake --build ${build_dir}/c-blosc-pm-cpu-build --target install --parallel ${PARALLEL}
rm -rf ${build_dir}/c-blosc-pm-cpu-build

# ADIOS2
if [ -d $HOME/src/adios2 ]
then
  cd $HOME/src/adios2
  git fetch --prune
  git checkout v2.10.2
  cd -
else
  git clone -b v2.10.2 https://github.com/ornladios/ADIOS2.git $HOME/src/adios2
fi
rm -rf $HOME/src/adios2-pm-cpu-build
cmake -S $HOME/src/adios2 -B ${build_dir}/adios2-pm-cpu-build -DADIOS2_USE_Blosc=ON -DADIOS2_USE_CUDA=OFF -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DADIOS2_USE_ZeroMQ=OFF -DCMAKE_INSTALL_PREFIX=${SW_DIR}/adios2-2.10.2
cmake --build ${build_dir}/adios2-pm-cpu-build --target install -j ${PARALLEL}
rm -rf ${build_dir}/adios2-pm-cpu-build

# BLAS++ (for PSATD+RZ)
if [ -d $HOME/src/blaspp ]
then
  cd $HOME/src/blaspp
  git fetch --prune
  git checkout v2024.05.31
  cd -
else
  git clone -b v2024.05.31 https://github.com/icl-utk-edu/blaspp.git $HOME/src/blaspp
fi
rm -rf $HOME/src/blaspp-pm-cpu-build
CXX=$(which CC) cmake -S $HOME/src/blaspp -B ${build_dir}/blaspp-pm-cpu-build -Duse_openmp=ON -Dgpu_backend=OFF -DCMAKE_CXX_STANDARD=20 -DCMAKE_INSTALL_PREFIX=${SW_DIR}/blaspp-2024.05.31
cmake --build ${build_dir}/blaspp-pm-cpu-build --target install --parallel ${PARALLEL}
rm -rf ${build_dir}/blaspp-pm-cpu-build

# LAPACK++ (for PSATD+RZ)
if [ -d $HOME/src/lapackpp ]
then
  cd $HOME/src/lapackpp
  git fetch --prune
  git checkout v2024.05.31
  cd -
else
  git clone -b v2024.05.31 https://github.com/icl-utk-edu/lapackpp.git $HOME/src/lapackpp
fi
rm -rf $HOME/src/lapackpp-pm-cpu-build
CXX=$(which CC) CXXFLAGS="-DLAPACK_FORTRAN_ADD_" cmake -S $HOME/src/lapackpp -B ${build_dir}/lapackpp-pm-cpu-build -DCMAKE_CXX_STANDARD=20 -Dbuild_tests=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX=${SW_DIR}/lapackpp-2024.05.31
cmake --build ${build_dir}/lapackpp-pm-cpu-build --target install --parallel ${PARALLEL}
rm -rf ${build_dir}/lapackpp-pm-cpu-build

# Python ######################################################################
#
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade virtualenv
python3 -m pip cache purge
rm -rf ${SW_DIR}/venvs/warpx-cpu
python3 -m venv ${SW_DIR}/venvs/warpx-cpu
source ${SW_DIR}/venvs/warpx-cpu/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade build
python3 -m pip install --upgrade packaging
python3 -m pip install --upgrade wheel
python3 -m pip install --upgrade setuptools[core]
python3 -m pip install --upgrade cython
python3 -m pip install --upgrade numpy
python3 -m pip install --upgrade pandas
python3 -m pip install --upgrade scipy
MPICC="cc -shared" python3 -m pip install --upgrade mpi4py --no-cache-dir --no-build-isolation --no-binary mpi4py
python3 -m pip install --upgrade openpmd-api
python3 -m pip install --upgrade matplotlib
python3 -m pip install --upgrade yt
# install or update WarpX dependencies
python3 -m pip install --upgrade -r $HOME/src/warpx/requirements.txt
# optimas (based on libEnsemble & ax->botorch->gpytorch->pytorch)
python3 -m pip install --upgrade torch --index-url https://download.pytorch.org/whl/cpu
python3 -m pip install --upgrade optimas[all]


# remove build temporary directory
rm -rf ${build_dir}
```

<a id="building-perlmutter-compilation"></a>

## Compilation

Use the following [cmake commands](../cmake.md#install-build-cmake) to compile the application executable:

### A100 GPUs

```bash
cd $HOME/src/warpx
rm -rf build_pm_gpu

cmake -S . -B build_pm_gpu -DWarpX_COMPUTE=CUDA -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_pm_gpu -j 16
```

The WarpX application executables are now in `$HOME/src/warpx/build_pm_gpu/bin/`.
Additionally, the following commands will install WarpX as a Python module:

```bash
cd $HOME/src/warpx
rm -rf build_pm_gpu_py

cmake -S . -B build_pm_gpu_py -DWarpX_COMPUTE=CUDA -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_APP=OFF -DWarpX_PYTHON=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_pm_gpu_py -j 16 --target pip_install
```

### CPU Nodes

```bash
cd $HOME/src/warpx
rm -rf build_pm_cpu

cmake -S . -B build_pm_cpu -DWarpX_COMPUTE=OMP -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_pm_cpu -j 16
```

The WarpX application executables are now in `$HOME/src/warpx/build_pm_cpu/bin/`.
Additionally, the following commands will install WarpX as a Python module:

```bash
rm -rf build_pm_cpu_py

cmake -S . -B build_pm_cpu_py -DWarpX_COMPUTE=OMP -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_APP=OFF -DWarpX_PYTHON=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_pm_cpu_py -j 16 --target pip_install
```

Now, you can [submit Perlmutter compute jobs](#running-cpp-perlmutter) for WarpX [Python (PICMI) scripts](../../usage/python.md#usage-picmi) ([example scripts](../../usage/examples.md#usage-examples)).
Or, you can use the WarpX executables to submit Perlmutter jobs ([example inputs](../../usage/examples.md#usage-examples)).
For executables, you can reference their location in your [job script](#running-cpp-perlmutter) or copy them to a location in `$PSCRATCH`.

<a id="building-perlmutter-update"></a>

## Update WarpX & Dependencies

If you already installed WarpX in the past and want to update it, start by getting the latest source code:

```bash
cd $HOME/src/warpx

# read the output of this command - does it look ok?
git status

# get the latest WarpX source code
git fetch
git pull

# read the output of these commands - do they look ok?
git status
git log # press q to exit
```

And, if needed,

- [update the perlmutter_gpu_warpx.profile or perlmutter_cpu_warpx files](#building-perlmutter-preparation),
- log out and into the system, activate the now updated environment profile as usual,
- [execute the dependency install scripts](#building-perlmutter-preparation).

As a last step, clean the build directory `rm -rf $HOME/src/warpx/build_pm_*` and rebuild WarpX.

<a id="running-cpp-perlmutter"></a>

## Running

### A100 (40GB) GPUs

The batch script below can be used to run a WarpX simulation on multiple nodes (change `-N` accordingly) on the supercomputer Perlmutter at NERSC.
This partition has up to [1536 nodes](https://docs.nersc.gov/systems/perlmutter/architecture/).

The batch script is set up for Python (PICMI) by default (`EXE=python3`, `INPUTS=run_script.py`).
To use the WarpX executable instead, edit the script and comment the Python lines and uncomment the `EXE`/`INPUTS` lines.
Source your `perlmutter_gpu_warpx.profile` before submitting (the script will try to source it if you forgot).
Replace descriptions between chevrons `<>` by relevant values, for instance `<input file>` could be `plasma_mirror_inputs`.
Note that we run one MPI rank per GPU.

```bash
#!/bin/bash -l

# Copyright 2021-2023 Axel Huebl, Kevin Gott
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

#SBATCH -t 00:10:00
#SBATCH -N 2
#SBATCH -J WarpX
#    note: <proj> must end on _g
#SBATCH -A <proj>
#SBATCH -q regular
# A100 40GB (most nodes)
#SBATCH -C gpu
# A100 80GB (256 nodes)
#S BATCH -C gpu&hbm80g
#SBATCH --exclusive
#SBATCH --cpus-per-task=32
# ideally single:1, but NERSC cgroups issue
#SBATCH --gpu-bind=none
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH -o WarpX.o%j
#SBATCH -e WarpX.e%j

# python interpreter & script here
EXE=python3
INPUTS=run_script.py
# or executable & inputs file
#EXE=./warpx
#INPUTS=inputs

# environment setup
if [[ -z "${MY_PROFILE}" ]]; then
    echo "WARNING: FORGOT TO"
    echo "   source $HOME/perlmutter_gpu_warpx.profile"
    echo "before submission. Doing that now."

    source $HOME/perlmutter_gpu_warpx.profile
fi

# pin to closest NIC to GPU
export MPICH_OFI_NIC_POLICY=GPU

# threads for OpenMP and threaded compressors per MPI rank
#   note: 16 avoids hyperthreading (32 virtual cores, 16 physical)
export OMP_NUM_THREADS=16

# GPU-aware MPI optimizations
export AMREX_DEFAULT_INIT="amrex.use_gpu_aware_mpi=1"

# CUDA visible devices are ordered inverse to local task IDs
#   Reference: nvidia-smi topo -m
srun --cpu-bind=cores bash -c "
    export CUDA_VISIBLE_DEVICES=\$((3-SLURM_LOCALID));
    ${EXE} ${INPUTS}" \
  > output.txt
```

To run a simulation, copy the lines above to a file `perlmutter_gpu.sbatch` and run

```bash
sbatch perlmutter_gpu.sbatch
```

to submit the job.

### A100 (80GB) GPUs

Perlmutter has [256 nodes](https://docs.nersc.gov/systems/perlmutter/architecture/) that provide 80 GB HBM per A100 GPU.
In the A100 (40GB) batch script, replace `-C gpu` with `-C gpu&hbm80g` to use these large-memory GPUs.

### CPU Nodes

The Perlmutter CPU partition has up to [3072 nodes](https://docs.nersc.gov/systems/perlmutter/architecture/), each with 2x AMD EPYC 7763 CPUs.
The batch script is set up for Python (PICMI) by default; to use the executable instead, edit `EXE`/`INPUTS`. Source your `perlmutter_cpu_warpx.profile` before submitting (the script will try to source it if you forgot).

```bash
#!/bin/bash -l

# Copyright 2021-2023 WarpX
#
# This file is part of WarpX.
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL

#SBATCH -t 00:10:00
#SBATCH -N 2
#SBATCH -J WarpX
#SBATCH -A <proj>
#SBATCH -q regular
#SBATCH -C cpu
# 8 cores per chiplet, 2x SMP
#SBATCH --cpus-per-task=16
#SBATCH --ntasks-per-node=16
#SBATCH --exclusive
#SBATCH -o WarpX.o%j
#SBATCH -e WarpX.e%j

# python interpreter & script here
EXE=python3
INPUTS=run_script.py
# or executable & inputs file
#EXE=./warpx
#INPUTS=inputs_small

# environment setup
if [[ -z "${MY_PROFILE}" ]]; then
    echo "WARNING: FORGOT TO"
    echo "   source $HOME/perlmutter_cpu_warpx.profile"
    echo "before submission. Doing that now."

    source $HOME/perlmutter_cpu_warpx.profile
fi

# each CPU node on Perlmutter (NERSC) has 64 hardware cores with
# 2x Hyperthreading/SMP
# https://en.wikichip.org/wiki/amd/epyc/7763
# https://www.amd.com/en/products/cpu/amd-epyc-7763
# Each CPU is made up of 8 chiplets, each sharing 32MB L3 cache.
# This will be our MPI rank assignment (2x8 is 16 ranks/node).

# threads for OpenMP and threaded compressors per MPI rank
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

srun --cpu-bind=cores \
  ${EXE} ${INPUTS} \
  > output.txt
```

<a id="post-processing-perlmutter"></a>

## Post-Processing

For post-processing, most users use Python via NERSC’s [Jupyter service](https://jupyter.nersc.gov) ([documentation](https://docs.nersc.gov/services/jupyter/)).

As a one-time preparatory setup, log into Perlmutter via SSH and do *not* source the WarpX profile script above.
Create your own Conda environment and [Jupyter kernel](https://docs.nersc.gov/services/jupyter/how-to-guides/#how-to-use-a-conda-environment-as-a-python-kernel) for post-processing:

```bash
module load python

conda config --set auto_activate_base false

# create conda environment
rm -rf $HOME/.conda/envs/warpx-pm-postproc
conda create --yes -n warpx-pm-postproc -c conda-forge mamba conda-libmamba-solver
conda activate warpx-pm-postproc
conda config --set solver libmamba
mamba install --yes -c conda-forge python ipykernel ipympl matplotlib numpy pandas yt openpmd-viewer openpmd-api h5py fast-histogram dask dask-jobqueue pyarrow

# create Jupyter kernel
rm -rf $HOME/.local/share/jupyter/kernels/warpx-pm-postproc/
python -m ipykernel install --user --name warpx-pm-postproc --display-name WarpX-PM-PostProcessing
echo -e '#!/bin/bash\nmodule load python\nconda activate warpx-pm-postproc\nexec "$@"' > $HOME/.local/share/jupyter/kernels/warpx-pm-postproc/kernel-helper.sh
chmod a+rx $HOME/.local/share/jupyter/kernels/warpx-pm-postproc/kernel-helper.sh
KERNEL_STR=$(jq '.argv |= ["{resource_dir}/kernel-helper.sh"] + .' $HOME/.local/share/jupyter/kernels/warpx-pm-postproc/kernel.json | jq '.argv[1] = "python"')
echo ${KERNEL_STR} | jq > $HOME/.local/share/jupyter/kernels/warpx-pm-postproc/kernel.json

exit
```

When opening a Jupyter notebook on [https://jupyter.nersc.gov](https://jupyter.nersc.gov), just select `WarpX-PM-PostProcessing` from the list of available kernels on the top right of the notebook.

Additional software can be installed later on, e.g., in a Jupyter cell using `!mamba install -y -c conda-forge ...`.
Software that is not available via conda can be installed via `!python -m pip install ...`.
