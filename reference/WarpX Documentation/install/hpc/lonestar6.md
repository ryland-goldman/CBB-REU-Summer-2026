<a id="building-lonestar6"></a>

# Lonestar6 (TACC)

The [Lonestar6 cluster](https://portal.tacc.utexas.edu/user-guides/lonestar6) is located at [TACC](https://www.tacc.utexas.edu).

## Introduction

If you are new to this system, **please see the following resources**:

* [TACC user guide](https://portal.tacc.utexas.edu/user-guides/)
* Batch system: [Slurm](https://portal.tacc.utexas.edu/user-guides/lonestar6#job-management)
* [Jupyter service](https://tacc.github.io/ctls2017/docs/intro_to_python/intro_to_python_011_jupyter.html)
* [Filesystem directories](https://portal.tacc.utexas.edu/user-guides/lonestar6#managing-files-on-lonestar6):
  * `$HOME`: per-user home directory, backed up (10 GB)
  * `$WORK`: per-user production directory, not backed up, not purged, Lustre (1 TB)
  * `$SCRATCH`: per-user production directory, not backed up, purged every 10 days, Lustre (no limits, 8PByte total)

## Installation

Use the following commands to download the WarpX source code and switch to the correct branch:

```bash
git clone https://github.com/BLAST-WarpX/warpx.git $WORK/src/warpx
```

We use system software modules, add environment hints and further dependencies via the file `$HOME/lonestar6_warpx_a100.profile`.
Create it now:

```bash
cp $HOME/src/warpx/Tools/machines/lonestar6-tacc/lonestar6_warpx_a100.profile.example $HOME/lonestar6_warpx_a100.profile
```

### Script Details

```bash
# please set your project account
#export proj="<yourProject>"  # change me

# required dependencies
module purge
module load TACC
module load gcc/11.2.0
module load cuda/12.2
module load cmake
module load mvapich2

# optional: for QED support with detailed tables
module load boost/1.84

# optional: for openPMD and PSATD+RZ support
module load phdf5/1.10.4

SW_DIR="${WORK}/sw/lonestar6/sw/lonestar6/a100"
export CMAKE_PREFIX_PATH=${SW_DIR}/c-blosc-1.21.1:${CMAKE_PREFIX_PATH}
export CMAKE_PREFIX_PATH=${SW_DIR}/adios2-2.10.2:${CMAKE_PREFIX_PATH}
export CMAKE_PREFIX_PATH=${SW_DIR}/blaspp-2024.05.31:${CMAKE_PREFIX_PATH}
export CMAKE_PREFIX_PATH=${SW_DIR}/lapackpp-2024.05.31:${CMAKE_PREFIX_PATH}

export LD_LIBRARY_PATH=${SW_DIR}/c-blosc-1.21.1/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/adios2-2.10.2/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/blaspp-2024.05.31/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/lapackpp-2024.05.31/lib64:$LD_LIBRARY_PATH

export PATH=${SW_DIR}/adios2-2.10.2/bin:${PATH}

# optional: CCache
#module load ccache  # TODO: request from support

# optional: for Python bindings or libEnsemble
module load python3/3.9.7

if [ -d "$WORK/sw/lonestar6/a100/venvs/warpx-a100" ]
then
  source $WORK/sw/lonestar6/a100/venvs/warpx-a100/bin/activate
fi

# an alias to request an interactive batch node for one hour
#   for parallel execution, start on the batch node: srun <command>
alias getNode="salloc -N 1 --ntasks-per-node=2 -t 1:00:00 -p gpu-100 --gpu-bind=single:1 -c 32 -G 2 -A $proj"
# an alias to run a command on a batch node for up to 30min
#   usage: runNode <command>
alias runNode="srun -N 1 --ntasks-per-node=2 -t 0:30:00 -p gpu-100 --gpu-bind=single:1 -c 32 -G 2 -A $proj"

# optimize CUDA compilation for A100
export AMREX_CUDA_ARCH=8.0

# compiler environment hints
export CC=$(which gcc)
export CXX=$(which g++)
export FC=$(which gfortran)
export CUDACXX=$(which nvcc)
export CUDAHOSTCXX=${CXX}
```

Edit the 2nd line of this script, which sets the `export proj=""` variable.
For example, if you are member of the project `abcde`, then run `nano $HOME/lonestar6_warpx_a100.profile` and edit line 2 to read:

```bash
export proj="abcde"
```

Exit the `nano` editor with `Ctrl` + `O` (save) and then `Ctrl` + `X` (exit).

#### IMPORTANT
Now, and as the first step on future logins to Lonestar6, activate these environment settings:

```bash
source $HOME/lonestar6_warpx_a100.profile
```

Finally, since Lonestar6 does not yet provide software modules for some of our dependencies, install them once:

```bash
bash $HOME/src/warpx/Tools/machines/lonestar6-tacc/install_a100_dependencies.sh
source ${SW_DIR}/venvs/warpx-a100/bin/activate
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
#   Was lonestar6_warpx_a100.profile sourced and configured correctly?
if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your lonestar6_warpx_a100.profile file! Please edit its line 2 to continue!"; exit 1; fi


# Remove old dependencies #####################################################
#
SW_DIR="${WORK}/sw/lonestar6/sw/lonestar6/a100"
rm -rf ${SW_DIR}
mkdir -p ${SW_DIR}

# remove common user mistakes in python, located in .local instead of a venv
python3 -m pip uninstall -qq -y pywarpx
python3 -m pip uninstall -qq -y warpx
python3 -m pip uninstall -qqq -y mpi4py 2>/dev/null || true


# General extra dependencies ##################################################
#

# tmpfs build directory: avoids issues often seen with $HOME and is faster
build_dir=$(mktemp -d)

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
rm -rf $HOME/src/c-blosc-a100-build
cmake -S $HOME/src/c-blosc -B ${build_dir}/c-blosc-a100-build -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF -DDEACTIVATE_AVX2=OFF -DCMAKE_INSTALL_PREFIX=${SW_DIR}/c-blosc-1.21.1
cmake --build ${build_dir}/c-blosc-a100-build --target install --parallel 16
rm -rf ${build_dir}/c-blosc-a100-build

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
rm -rf $HOME/src/adios2-a100-build
cmake -S $HOME/src/adios2 -B ${build_dir}/adios2-a100-build -DADIOS2_USE_Blosc=ON -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DADIOS2_USE_ZeroMQ=OFF -DCMAKE_INSTALL_PREFIX=${SW_DIR}/adios2-2.10.2
cmake --build ${build_dir}/adios2-a100-build --target install -j 16
rm -rf ${build_dir}/adios2-a100-build

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
rm -rf $HOME/src/blaspp-a100-build
cmake -S $HOME/src/blaspp -B ${build_dir}/blaspp-a100-build -Duse_openmp=OFF -Dgpu_backend=cuda -DCMAKE_CXX_STANDARD=20 -DCMAKE_INSTALL_PREFIX=${SW_DIR}/blaspp-2024.05.31
cmake --build ${build_dir}/blaspp-a100-build --target install --parallel 16
rm -rf ${build_dir}/blaspp-a100-build

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
rm -rf $HOME/src/lapackpp-a100-build
CXXFLAGS="-DLAPACK_FORTRAN_ADD_" cmake -S $HOME/src/lapackpp -B ${build_dir}/lapackpp-a100-build -DCMAKE_CXX_STANDARD=20 -Dbuild_tests=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX=${SW_DIR}/lapackpp-2024.05.31
cmake --build ${build_dir}/lapackpp-a100-build --target install --parallel 16
rm -rf ${build_dir}/lapackpp-a100-build

# Python ######################################################################
#
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade virtualenv
python3 -m pip cache purge
rm -rf ${SW_DIR}/venvs/warpx-a100
python3 -m venv ${SW_DIR}/venvs/warpx-a100
source ${SW_DIR}/venvs/warpx-a100/bin/activate
python3 -m pip install --upgrade pip
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
# install or update WarpX dependencies
python3 -m pip install --upgrade -r $HOME/src/warpx/requirements.txt
#python3 -m pip install --upgrade cupy-cuda12x  # CUDA 12 compatible wheel
# optimas (based on libEnsemble & ax->botorch->gpytorch->pytorch)
#python3 -m pip install --upgrade torch  # CUDA 12 compatible wheel
#python3 -m pip install --upgrade optimas[all]


# remove build temporary directory
rm -rf ${build_dir}
```

<a id="building-lonestar6-compilation"></a>

## Compilation

Use the following [cmake commands](../cmake.md#install-build-cmake) to compile the application executable:

```bash
cd $HOME/src/warpx
rm -rf build_pm_gpu

cmake -S . -B build_gpu -DWarpX_COMPUTE=CUDA -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_gpu -j 16
```

The WarpX application executables are now in `$HOME/src/warpx/build_gpu/bin/`.
Additionally, the following commands will install WarpX as a Python module:

```bash
cd $HOME/src/warpx
rm -rf build_pm_gpu_py

cmake -S . -B build_gpu_py -DWarpX_COMPUTE=CUDA -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_APP=OFF -DWarpX_PYTHON=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_gpu_py -j 16 --target pip_install
```

Now, you can [submit Lonestar6 compute jobs](#running-cpp-lonestar6) for WarpX [Python (PICMI) scripts](../../usage/python.md#usage-picmi) ([example scripts](../../usage/examples.md#usage-examples)).
Or, you can use the WarpX executables to submit Lonestar6 jobs ([example inputs](../../usage/examples.md#usage-examples)).
For executables, you can reference their location in your [job script](#running-cpp-lonestar6) or copy them to a location in `$WORK` or `$SCRATCH`.

<a id="running-cpp-lonestar6"></a>

## Running

<a id="running-cpp-lonestar6-a100-gpus"></a>

### A100 GPUs (40 GB)

[84 GPU nodes, each with 2 A100 GPUs (40 GB)](https://portal.tacc.utexas.edu/user-guides/lonestar6#system-gpu).

The batch script below can be used to run a WarpX simulation on multiple nodes (change `-N` accordingly) on the supercomputer lonestar6 at tacc.
Replace descriptions between chevrons `<>` by relevant values, for instance `<input file>` could be `plasma_mirror_inputs`.
Note that we run one MPI rank per GPU.

```bash
#!/bin/bash -l

# Copyright 2021-2022 Axel Huebl, Kevin Gott
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
#SBATCH -C gpu
#SBATCH --exclusive
#SBATCH --cpus-per-task=32
#SBATCH --gpu-bind=none
#SBATCH --gpus-per-node=4
#SBATCH -o WarpX.o%j
#SBATCH -e WarpX.e%j

# executable & inputs file or python interpreter & PICMI script here
EXE=./warpx
INPUTS=inputs_small

# pin to closest NIC to GPU
export MPICH_OFI_NIC_POLICY=GPU

# threads for OpenMP and threaded compressors per MPI rank
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# depends on https://github.com/BLAST-WarpX/warpx/issues/2009
#GPU_AWARE_MPI="amrex.the_arena_is_managed=0 amrex.use_gpu_aware_mpi=1"
GPU_AWARE_MPI=""

# CUDA visible devices are ordered inverse to local task IDs
#   Reference: nvidia-smi topo -m
srun --cpu-bind=cores bash -c "
    export CUDA_VISIBLE_DEVICES=\$((3-SLURM_LOCALID));
    ${EXE} ${INPUTS} ${GPU_AWARE_MPI}" \
  > output.txt
```

To run a simulation, copy the lines above to a file `lonestar6.sbatch` and run

```bash
sbatch lonestar6_a100.sbatch
```

to submit the job.
