<a id="building-pitzer"></a>

# Pitzer (OSC)

The [Pitzer cluster](https://www.osc.edu/supercomputing/computing/pitzer) is located at the Ohio Supercomputer Center (OSC). It is currently the main CPU/GPU cluster at OSC. However, the [Cardinal cluster](https://www.osc.edu/resources/technical_support/supercomputers/cardinal) is soon going to take over Pitzer to become the next major CPU/GPU cluster at OSC in the second half of 2024. A list of all OSC clusters can be found [here](https://www.osc.edu/services/cluster_computing).

The Pitzer cluster offers a variety of partitions suitable for different computational needs, including GPU nodes, CPU nodes, and nodes with large memory capacities. For more information on the specifications and capabilities of these partitions, visit the [Ohio Supercomputer Center’s Pitzer page](https://www.osc.edu/supercomputing/computing/pitzer).

## Introduction

If you are new to this system, **please see the following resources**:

* [Pitzer user guide](https://www.osc.edu/resources/getting_started/new_user_resource_guide)
* Batch system: [Slurm](https://www.osc.edu/supercomputing/batch-processing-at-osc)
* [Jupyter service](https://www.osc.edu/vocabulary/documentation/jupyter)
* [Filesystems](https://www.osc.edu/supercomputing/storage-environment-at-osc/storage-hardware/overview_of_file_systems):
  * `$HOME`: per-user directory, use only for inputs, source, and scripts; backed up (500GB)
  * `/fs/ess`: per-project storage directory, use for long-term storage of data and analysis; backed up (1-5TB)
  * `/fs/scratch`: per-project production directory; fast I/O for parallel jobs; not backed up (100TB)

<a id="building-pitzer-preparation"></a>

## Preparation

Use the following commands to download the WarpX source code:

```bash
git clone https://github.com/BLAST-WarpX/warpx.git $HOME/src/warpx
```

On Pitzer, you can run either on GPU nodes with V100 GPUs or CPU nodes.

### V100 GPUs

We use system software modules, add environment hints and further dependencies via the file `$HOME/pitzer_v100_warpx.profile`.
Create it now:

```bash
cp $HOME/src/warpx/Tools/machines/pitzer-osc/pitzer_v100_warpx.profile.example $HOME/pitzer_v100_warpx.profile
```

### Script Details

```bash
# please set your project account
export proj="" # change me!

# remembers the location of this script
export MY_V100_PROFILE=$(cd $(dirname $BASH_SOURCE) && pwd)"/"$(basename $BASH_SOURCE)
if [ -z ${proj-} ]; then
  echo "WARNING: The 'proj' variable is not yet set in your $MY_V100_PROFILE file! Please edit its line 2 to continue!"
  return
fi

export SW_DIR="${HOME}/sw/osc/pitzer/v100"

module purge
module load cmake/3.25.2
module load intel/19.0.5
module load cuda/11.8.0
module load openmpi-cuda/4.1.5-hpcx
module load gcc-compatibility/11.2.0

# optional: for python binding support
module load miniconda3/24.1.2-py310
export VENV_NAME="warpx-pitzer-v100"
if [ -d "${SW_DIR}/venvs/${VENV_NAME}" ]; then
  source ${SW_DIR}/venvs/${VENV_NAME}/bin/activate
fi

# an alias to request an interactive batch node for one hour
#   for parallel execution, start on the batch node: srun <command>
alias getNode="salloc -N 1 --ntasks-per-node=2 --cpus-per-task=20 --gpus-per-task=v100:1 -t 1:00:00 -A $proj"
# an alias to run a command on a batch node for up to 30min
#   usage: runNode <command>
alias runNode="srun -N 1 --ntasks-per-node=2 --cpus-per-task=20 --gpus-per-task=v100:1 -t 1:00:00 -A $proj"

# optional: for PSATD in RZ geometry support
export CMAKE_PREFIX_PATH=${SW_DIR}/blaspp-2024.05.31:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/lapackpp-2024.05.31:$CMAKE_PREFIX_PATH
export LD_LIBRARY_PATH=${SW_DIR}/blaspp-2024.05.31/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/lapackpp-2024.05.31/lib64:$LD_LIBRARY_PATH

# optional: for QED lookup table generation support
# use self-installed boost
export CMAKE_PREFIX_PATH=${SW_DIR}/boost-1.82.0:$CMAKE_PREFIX_PATH
export LD_LIBRARY_PATH=${SW_DIR}/boost-1.82.0/lib:$LD_LIBRARY_PATH

# optional: for openPMD support (hdf5 and adios2)
# use self-installed hdf5
module load hdf5/1.12.0

export CMAKE_PREFIX_PATH=${SW_DIR}/c-blosc-1.21.6:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/adios2-2.10.2:$CMAKE_PREFIX_PATH
export LD_LIBRARY_PATH=${SW_DIR}/c-blosc-1.21.6/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/adios2-2.10.2/lib64:$LD_LIBRARY_PATH
export PATH=${SW_DIR}/adios2-2.10.2/bin:${PATH}

# avoid relocation truncation error which result from large executable size
export CUDAFLAGS="--host-linker-script=use-lcs" # https://github.com/BLAST-WarpX/warpx/pull/3673
export AMREX_CUDA_ARCH=7.0 # 7.0: V100, 8.0: V100, 9.0: H100 https://github.com/BLAST-WarpX/warpx/issues/3214

# compiler environment hints
export CC=$(which gcc)
export CXX=$(which g++)
export FC=$(which gfortran)
export CUDACXX=$(which nvcc)
export CUDAHOSTCXX=${CXX}
```

Edit the 2nd line of this script, which sets the `export proj=""` variable.
For example, if you are member of the project `pas2024`, then run `nano $HOME/pitzer_v100_warpx.profile` and edit line 2 to read:

```bash
export proj="pas2024"
```

Exit the `nano` editor with `Ctrl` + `O` (save) and then `Ctrl` + `X` (exit).

#### IMPORTANT
Now, and as the first step on future logins to pitzer, activate these environment settings:

```bash
source $HOME/pitzer_v100_warpx.profile
```

Finally, since pitzer does not yet provide software modules for some of our dependencies, install them once:

```bash
bash $HOME/src/warpx/Tools/machines/pitzer-osc/install_v100_dependencies.sh
```

### Script Details

```bash
#!/bin/bash
#
# Copyright 2024 The WarpX Community
#
# This file is part of WarpX.
#
# Author: Zhongwei Wang
# License: BSD-3-Clause-LBNL

# Exit on first error encountered #############################################
#
set -eu -o pipefail

# Check: ######################################################################
#
#   Was pitzer_v100_warpx.profile sourced and configured correctly?
if [ -z ${proj-} ]; then
  echo "WARNING: The 'proj' variable is not yet set in your pitzer_v100_warpx.profile file! Please edit its line 2 to continue!"
  exit 1
fi

# Remove old dependencies #####################################################
#
rm -rf ${SW_DIR}
mkdir -p ${SW_DIR}

# remove common user mistakes in python, located in .local instead of a venv
python3 -m pip uninstall -qq -y pywarpx
python3 -m pip uninstall -qq -y warpx
python3 -m pip uninstall -qqq -y mpi4py 2>/dev/null || true

# General extra dependencies ##################################################
#
SRC_DIR="${HOME}/src"
build_dir=$(mktemp -d)

# boost (for QED table generation support)
cd ${SRC_DIR}
wget https://archives.boost.io/release/1.82.0/source/boost_1_82_0.tar.gz
tar -xzvf boost_1_82_0.tar.gz
rm -rf boost_1_82_0.tar.gz
cd -

cd ${SRC_DIR}/boost_1_82_0
./bootstrap.sh --prefix=${SW_DIR}/boost-1.82.0
./b2 install
cd -

# BLAS++ (for PSATD+RZ)
if [ -d ${SRC_DIR}/blaspp ]; then
  cd ${SRC_DIR}/blaspp
  git fetch
  git checkout v2024.05.31
  cd -
else
  git clone -b v2024.05.31 https://github.com/icl-utk-edu/blaspp.git ${SRC_DIR}/blaspp
fi
rm -rf ${build_dir}/blaspp-pitzer-v100-build
CXX=$(which CC) cmake -S ${SRC_DIR}/blaspp \
  -B ${build_dir}/blaspp-pitzer-v100-build \
  -Duse_openmp=ON \
  -Dgpu_backend=cuda \
  -DCMAKE_CXX_STANDARD=20 \
  -DCMAKE_INSTALL_PREFIX=${SW_DIR}/blaspp-2024.05.31
cmake --build ${build_dir}/blaspp-pitzer-v100-build --target install --parallel 16
rm -rf ${build_dir}/blaspp-pitzer-v100-build

# LAPACK++ (for PSATD+RZ)
if [ -d ${SRC_DIR}/lapackpp ]; then
  cd ${SRC_DIR}/lapackpp
  git fetch
  git checkout v2024.05.31
  cd -
else
  git clone -b v2024.05.31 https://github.com/icl-utk-edu/lapackpp.git ${SRC_DIR}/lapackpp
fi
rm -rf ${build_dir}/lapackpp-pitzer-v100-build
CXX=$(which CC) CXXFLAGS="-DLAPACK_FORTRAN_ADD_" cmake -S ${SRC_DIR}/lapackpp \
  -B ${build_dir}/lapackpp-pitzer-v100-build \
  -DCMAKE_CXX_STANDARD=20 \
  -Dbuild_tests=OFF \
  -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON \
  -DCMAKE_INSTALL_PREFIX=${SW_DIR}/lapackpp-2024.05.31
cmake --build ${build_dir}/lapackpp-pitzer-v100-build --target install --parallel 16
rm -rf ${build_dir}/lapackpp-pitzer-v100-build

# c-blosc (I/O compression, for openPMD)
if [ -d ${SRC_DIR}/c-blosc ]; then
  cd ${SRC_DIR}/c-blosc
  git fetch --prune
  git checkout v1.21.6
  cd -
else
  git clone -b v1.21.6 https://github.com/Blosc/c-blosc.git ${SRC_DIR}/c-blosc
fi
rm -rf ${build_dir}/c-blosc-pitzer-build
cmake -S ${SRC_DIR}/c-blosc \
  -B ${build_dir}/c-blosc-pitzer-build \
  -DBUILD_TESTS=OFF \
  -DBUILD_BENCHMARKS=OFF \
  -DDEACTIVATE_AVX2=OFF \
  -DCMAKE_INSTALL_PREFIX=${SW_DIR}/c-blosc-1.21.6
cmake --build ${build_dir}/c-blosc-pitzer-build --target install --parallel 16
rm -rf ${build_dir}/c-blosc-pitzer-build

# ADIOS2 (for openPMD)
if [ -d ${SRC_DIR}/adios2 ]; then
  cd ${SRC_DIR}/adios2
  git fetch --prune
  git checkout v2.10.2
  cd -
else
  git clone -b v2.10.2 https://github.com/ornladios/ADIOS2.git ${SRC_DIR}/adios2
fi
rm -rf ${build_dir}/adios2-pitzer-build
cmake -S ${SRC_DIR}/adios2 \
  -B ${build_dir}/adios2-pitzer-build \
  -DBUILD_TESTING=OFF \
  -DADIOS2_BUILD_EXAMPLES=OFF \
  -DADIOS2_USE_Blosc=ON \
  -DADIOS2_USE_Fortran=OFF \
  -DADIOS2_USE_Python=OFF \
  -DADIOS2_USE_SST=OFF \
  -DADIOS2_USE_ZeroMQ=OFF \
  -DCMAKE_INSTALL_PREFIX=${SW_DIR}/adios2-2.10.2
cmake --build ${build_dir}/adios2-pitzer-build --target install -j 16
rm -rf ${build_dir}/adios2-pitzer-build

rm -rf ${build_dir}

# Python ######################################################################
#
python3 -m pip install --upgrade --user virtualenv
rm -rf ${SW_DIR}/venvs/${VENV_NAME}
python3 -m venv ${SW_DIR}/venvs/${VENV_NAME}
source ${SW_DIR}/venvs/${VENV_NAME}/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip cache purge
python3 -m pip install --upgrade build
python3 -m pip install --upgrade packaging
python3 -m pip install --upgrade wheel
python3 -m pip install --upgrade setuptools[core]
python3 -m pip install --upgrade cython
python3 -m pip install --upgrade numpy
python3 -m pip install --upgrade pandas
python3 -m pip install --upgrade scipy
python3 -m pip install --upgrade mpi4py --no-cache-dir --no-build-isolation --no-binary mpi4py
python3 -m pip install --upgrade openpmd-api
python3 -m pip install --upgrade matplotlib
python3 -m pip install --upgrade yt

# install or update WarpX dependencies such as picmistandard
python3 -m pip install --upgrade -r ${SRC_DIR}/warpx/requirements.txt

# ML dependencies
python3 -m pip install --upgrade torch
```

### CPU Nodes

We use system software modules, add environment hints and further dependencies via the file `$HOME/pitzer_cpu_warpx.profile`.
Create it now:

```bash
cp $HOME/src/warpx/Tools/machines/pitzer-osc/pitzer_cpu_warpx.profile.example $HOME/pitzer_cpu_warpx.profile
```

### Script Details

```bash
# please set your project account
export proj="" # change me!

# remembers the location of this script
export MY_CPU_PROFILE=$(cd $(dirname $BASH_SOURCE) && pwd)"/"$(basename $BASH_SOURCE)
if [ -z ${proj-} ]; then
  echo "WARNING: The 'proj' variable is not yet set in your $MY_CPU_PROFILE file! Please edit its line 2 to continue!"
  return
fi

export SW_DIR="${HOME}/sw/osc/pitzer/cpu"

module purge
module load cmake/3.25.2
module load gnu/12.3.0
module load openmpi/4.1.5-hpcx

# optional: for python binding support
module load miniconda3/24.1.2-py310
export VENV_NAME="warpx-pitzer-cpu"
if [ -d "${SW_DIR}/venvs/${VENV_NAME}" ]; then
  source ${SW_DIR}/venvs/${VENV_NAME}/bin/activate
fi

# an alias to request an interactive batch node for one hour
#   for parallel execution, start on the batch node: srun <command>
alias getNode="salloc -N 1 --ntasks-per-node=2 --cpus-per-task=20 -t 1:00:00 -A $proj"
# an alias to run a command on a batch node for up to 30min
#   usage: runNode <command>
alias runNode="srun -N 1 --ntasks-per-node=2 --cpus-per-task=20 -t 1:00:00 -A $proj"

# optional: for PSATD in RZ geometry support
export CMAKE_PREFIX_PATH=${SW_DIR}/blaspp-2024.05.31:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/lapackpp-2024.05.31:$CMAKE_PREFIX_PATH
export LD_LIBRARY_PATH=${SW_DIR}/blaspp-2024.05.31/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/lapackpp-2024.05.31/lib64:$LD_LIBRARY_PATH

# optional: for QED lookup table generation support
# use self-installed boost
export CMAKE_PREFIX_PATH=${SW_DIR}/boost-1.82.0:$CMAKE_PREFIX_PATH
export LD_LIBRARY_PATH=${SW_DIR}/boost-1.82.0/lib:$LD_LIBRARY_PATH

# optional: for openPMD support (hdf5 and adios2)
module load hdf5/1.12.2
export CMAKE_PREFIX_PATH=${SW_DIR}/c-blosc-1.21.6:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/adios2-2.10.2:$CMAKE_PREFIX_PATH
export LD_LIBRARY_PATH=${SW_DIR}/c-blosc-1.21.6/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/adios2-2.10.2/lib64:$LD_LIBRARY_PATH
export PATH=${SW_DIR}/adios2-2.10.2/bin:${PATH}

# compiler environment hints
export CC=$(which gcc)
export CXX=$(which g++)
export FC=$(which gfortran)
```

Edit the 2nd line of this script, which sets the `export proj=""` variable.
For example, if you are member of the project `pas2024`, then run `nano $HOME/pitzer_cpu_warpx.profile` and edit line 2 to read:

```bash
export proj="pas2024"
```

Exit the `nano` editor with `Ctrl` + `O` (save) and then `Ctrl` + `X` (exit).

#### IMPORTANT
Now, and as the first step on future logins to pitzer, activate these environment settings:

```bash
source $HOME/pitzer_cpu_warpx.profile
```

Finally, since pitzer does not yet provide software modules for some of our dependencies, install them once:

```bash
bash $HOME/src/warpx/Tools/machines/pitzer-osc/install_cpu_dependencies.sh
```

### Script Details

```bash
#!/bin/bash
#
# Copyright 2024 The WarpX Community
#
# This file is part of WarpX.
#
# Author: Zhongwei Wang
# License: BSD-3-Clause-LBNL

# Exit on first error encountered #############################################
#
set -eu -o pipefail

# Check: ######################################################################
#
#   Was pitzer_cpu_warpx.profile sourced and configured correctly?
if [ -z ${proj-} ]; then
  echo "WARNING: The 'proj' variable is not yet set in your pitzer_cpu_warpx.profile file! Please edit its line 2 to continue!"
  exit 1
fi

# Remove old dependencies #####################################################
#
rm -rf ${SW_DIR}
mkdir -p ${SW_DIR}

# remove common user mistakes in python, located in .local instead of a venv
python3 -m pip uninstall -qq -y pywarpx
python3 -m pip uninstall -qq -y warpx
python3 -m pip uninstall -qqq -y mpi4py 2>/dev/null || true

# General extra dependencies ##################################################
#
SRC_DIR="${HOME}/src"
build_dir=$(mktemp -d)

# boost (for QED table generation support)
cd ${SRC_DIR}
wget https://archives.boost.io/release/1.82.0/source/boost_1_82_0.tar.gz
tar -xzvf boost_1_82_0.tar.gz
rm -rf boost_1_82_0.tar.gz
cd -

cd ${SRC_DIR}/boost_1_82_0
./bootstrap.sh --prefix=$SW_DIR/boost-1.82.0
./b2 install
cd -

# BLAS++ (for PSATD+RZ)
if [ -d ${SRC_DIR}/blaspp ]; then
  cd ${SRC_DIR}/blaspp
  git fetch --prune
  git checkout v2024.05.31
  cd -
else
  git clone -b v2024.05.31 https://github.com/icl-utk-edu/blaspp.git ${SRC_DIR}/blaspp
fi
rm -rf ${build_dir}/blaspp-pitzer-cpu-build
CXX=$(which CC) cmake -S ${SRC_DIR}/blaspp \
  -B ${build_dir}/blaspp-pitzer-cpu-build \
  -Duse_openmp=ON \
  -Dgpu_backend=OFF \
  -DCMAKE_CXX_STANDARD=20 \
  -DCMAKE_INSTALL_PREFIX=${SW_DIR}/blaspp-2024.05.31
cmake --build ${build_dir}/blaspp-pitzer-cpu-build --target install --parallel 16
rm -rf ${build_dir}/blaspp-pitzer-cpu-build

# LAPACK++ (for PSATD+RZ)
if [ -d ${SRC_DIR}/lapackpp ]; then
  cd ${SRC_DIR}/lapackpp
  git fetch --prune
  git checkout v2024.05.31
  cd -
else
  git clone -b v2024.05.31 https://github.com/icl-utk-edu/lapackpp.git ${SRC_DIR}/lapackpp
fi
rm -rf ${build_dir}/lapackpp-pitzer-cpu-build
CXX=$(which CC) CXXFLAGS="-DLAPACK_FORTRAN_ADD_" cmake -S ${SRC_DIR}/lapackpp \
  -B ${build_dir}/lapackpp-pitzer-cpu-build \
  -DCMAKE_CXX_STANDARD=20 \
  -Dbuild_tests=OFF \
  -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON \
  -DCMAKE_INSTALL_PREFIX=${SW_DIR}/lapackpp-2024.05.31
cmake --build ${build_dir}/lapackpp-pitzer-cpu-build --target install --parallel 16
rm -rf ${build_dir}/lapackpp-pitzer-cpu-build

# c-blosc (I/O compression, for openPMD)
if [ -d ${SRC_DIR}/c-blosc ]; then
  cd ${SRC_DIR}/c-blosc
  git fetch --prune
  git checkout v1.21.6
  cd -
else
  git clone -b v1.21.6 https://github.com/Blosc/c-blosc.git ${SRC_DIR}/c-blosc
fi
rm -rf ${build_dir}/c-blosc-pitzer-build
cmake -S ${SRC_DIR}/c-blosc \
  -B ${build_dir}/c-blosc-pitzer-build \
  -DBUILD_TESTS=OFF \
  -DBUILD_BENCHMARKS=OFF \
  -DDEACTIVATE_AVX2=OFF \
  -DCMAKE_INSTALL_PREFIX=${SW_DIR}/c-blosc-1.21.6
cmake --build ${build_dir}/c-blosc-pitzer-build --target install --parallel 16
rm -rf ${build_dir}/c-blosc-pitzer-build

# ADIOS2 (for openPMD)
if [ -d ${SRC_DIR}/adios2 ]; then
  cd ${SRC_DIR}/adios2
  git fetch --prune
  git checkout v2.10.2
  cd -
else
  git clone -b v2.10.2 https://github.com/ornladios/ADIOS2.git ${SRC_DIR}/adios2
fi
rm -rf ${build_dir}/adios2-pitzer-build
cmake -S ${SRC_DIR}/adios2 \
  -B ${build_dir}/adios2-pitzer-build \
  -DBUILD_TESTING=OFF \
  -DADIOS2_BUILD_EXAMPLES=OFF \
  -DADIOS2_USE_Blosc=ON \
  -DADIOS2_USE_Fortran=OFF \
  -DADIOS2_USE_Python=OFF \
  -DADIOS2_USE_SST=OFF \
  -DADIOS2_USE_ZeroMQ=OFF \
  -DCMAKE_INSTALL_PREFIX=${SW_DIR}/adios2-2.10.2
cmake --build ${build_dir}/adios2-pitzer-build --target install -j 16
rm -rf ${build_dir}/adios2-pitzer-build

rm -rf ${build_dir}

# Python ######################################################################
#
python3 -m pip install --upgrade --user virtualenv
rm -rf ${SW_DIR}/venvs/${VENV_NAME}
python3 -m venv ${SW_DIR}/venvs/${VENV_NAME}
source ${SW_DIR}/venvs/${VENV_NAME}/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip cache purge
python3 -m pip install --upgrade build
python3 -m pip install --upgrade packaging
python3 -m pip install --upgrade wheel
python3 -m pip install --upgrade setuptools[core]
python3 -m pip install --upgrade cython
python3 -m pip install --upgrade numpy
python3 -m pip install --upgrade pandas
python3 -m pip install --upgrade scipy
python3 -m pip install --upgrade mpi4py --no-cache-dir --no-build-isolation --no-binary mpi4py
python3 -m pip install --upgrade openpmd-api
python3 -m pip install --upgrade matplotlib
python3 -m pip install --upgrade yt

# install or update WarpX dependencies such as picmistandard
python3 -m pip install --upgrade -r ${SRC_DIR}/warpx/requirements.txt

# ML dependencies
python3 -m pip install --upgrade torch
```

<a id="building-pitzer-compilation"></a>

## Compilation

Use the following [cmake commands](../cmake.md#install-build-cmake) to compile the application executable:

### V100 GPUs

```bash
cd $HOME/src/warpx
rm -rf build_v100

cmake -S . -B build_v100 -DWarpX_COMPUTE=CUDA -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_v100 -j 48
```

The WarpX application executables are now in `$HOME/src/warpx/build_v100/bin/`. Additionally, the following commands will install WarpX as a Python module:

```bash
cd $HOME/src/warpx
rm -rf build_v100_py

cmake -S . -B build_v100_py -DWarpX_COMPUTE=CUDA -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_APP=OFF -DWarpX_PYTHON=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_v100_py -j 48 --target pip_install
```

### CPU Nodes

```bash
cd $HOME/src/warpx
rm -rf build

cmake -S . -B build -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build -j 48
```

The WarpX application executables are now in `$HOME/src/warpx/build/bin/`. Additionally, the following commands will install WarpX as a Python module:

```bash
cd $HOME/src/warpx
rm -rf build_py

cmake -S . -B build_py -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_APP=OFF -DWarpX_PYTHON=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_py -j 48 --target pip_install
```

Now, you can [submit Pitzer compute jobs](#running-pitzer) for WarpX [Python (PICMI) scripts](../../usage/python.md#usage-picmi) ([example scripts](../../usage/examples.md#usage-examples)). Or, you can use the WarpX executables to submit Pitzer jobs ([example inputs](../../usage/examples.md#usage-examples)). For executables, you can reference their location in your [job script](#running-pitzer) or copy them to a location in `/scratch`.

<a id="building-pitzer-update"></a>

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

- [update the pitzer_cpu_warpx.profile file](#building-pitzer-preparation),
- log out and into the system, activate the now updated environment profile as usual,
- [execute the dependency install scripts](#building-pitzer-preparation).

As a last step, clean the build directory `rm -rf $HOME/src/warpx/build_*` and rebuild WarpX.

<a id="running-pitzer"></a>

## Running

### V100 GPUs

Pitzer’s GPU partition includes:

- 32 nodes, each equipped with two V100 (16GB) GPUs.
- 42 nodes, each with two V100 (32GB) GPUs.
- 4 large memory nodes, each with quad V100 (32GB) GPUs.

To run a WarpX simulation on the GPU nodes, use the batch script provided below. Adjust the `-N` parameter in the script to match the number of nodes you intend to use. Each node in this partition supports running one MPI rank per GPU.

```bash
#!/bin/bash
#SBATCH --time=0:20:00
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --cpus-per-task=24
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=closest
#SBATCH --job-name=<job_name>
#SBATCH --account=<project_id>
#SBATCH --output=./logs/%x_%j.out
#SBATCH --error=./logs/%x_%j.err

# Pitzer cluster has 32 GPU nodes with dual Intel Xeon 6148 and dual V100 (16GB) GPUs and 42 nodes with dual Intel Xeon 8268 and dual V100 (32GB) GPUs. https://www.osc.edu/resources/technical_support/supercomputers/pitzer

source ${HOME}/pitzer_v100_warpx.profile
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "GPU Information:"
nvidia-smi

GPU_AWARE_MPI="amrex.use_gpu_aware_mpi=1"

# executable & inputs file or python interpreter & PICMI script here
EXE=${HOME}/src/warpx/build_v100/bin/warpx.2d
INPUTS=inputs

srun --cpu-bind=cores ${EXE} ${INPUTS} ${GPU_AWARE_MPI} >./logs/${SLURM_JOB_NAME}_${SLURM_JOBID}.log 2>&1
```

After preparing your script, submit your job with the following command:

```bash
sbatch pitzer_v100.sbatch
```

### CPU Nodes

For CPU-based computations, Pitzer offers:

- 224 nodes, each with dual Intel Xeon Gold 6148 CPUs and 192 GB RAM.
- 340 nodes, each with dual Intel Xeon Platinum 8268 CPUs and 192 GB RAM.
- 16 large memory nodes.

To submit a job to the CPU partition, use the provided batch script. Ensure you have copied the script to your working directory.

```bash
#!/bin/bash
#SBATCH --time=0:20:00
#SBATCH --nodes=1 --ntasks-per-node=6
#SBATCH --cpus-per-task=8
#SBATCH --job-name=<job_name>
#SBATCH --account=<project_id>
#SBATCH --output=./logs/%x_%j.out
#SBATCH --error=./logs/%x_%j.err

# Pitzer cluster has 224 CPU nodes equipped with dual Intel Xeon 6148 (40 cores per node) and 340 CPU nodes with dual Intel Xeon 8268 (48 cores per node). https://www.osc.edu/resources/technical_support/supercomputers/pitzer

source ${HOME}/pitzer_cpu_warpx.profile
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# executable & inputs file or python interpreter & PICMI script here
EXE=${HOME}/src/warpx/build/bin/warpx.2d
INPUTS=inputs

srun --cpu-bind=cores ${EXE} ${INPUTS} >./logs/${SLURM_JOB_NAME}_${SLURM_JOBID}.log 2>&1
```

Submit your job with:

```bash
sbatch pitzer_cpu.sbatch
```

<a id="post-processing-osc"></a>

## Post-Processing

For post-processing, many users prefer to use the online [Jupyter service](https://ondemand.osc.edu/pun/sys/dashboard/batch_connect/sessions) ([documentation](https://www.osc.edu/vocabulary/documentation/jupyter)) that is directly connected to the cluster’s fast filesystem.

#### NOTE
This section is a stub and contributions are welcome.
We can document further details, e.g., which recommended post-processing Python software to install or how to customize Jupyter kernels here.
