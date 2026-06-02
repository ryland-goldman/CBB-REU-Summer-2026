<a id="building-polaris"></a>

# Polaris (ALCF)

The [Polaris cluster](https://docs.alcf.anl.gov/polaris/getting-started/) is located at ALCF.

## Introduction

If you are new to this system, **please see the following resources**:

* [ALCF user guide](https://docs.alcf.anl.gov/)
* Batch system: [PBS](https://docs.alcf.anl.gov/running-jobs/job-and-queue-scheduling/)
* [Filesystems](https://docs.alcf.anl.gov/data-management/filesystem-and-storage/file-systems/)

<a id="building-polaris-preparation"></a>

## Preparation

Use the following commands to download the WarpX source code:

```bash
git clone https://github.com/BLAST-WarpX/warpx.git $HOME/src/warpx
```

On Polaris, you can run either on GPU nodes with fast A100 GPUs (recommended) or CPU nodes.

### A100 GPUs

We use system software modules, add environment hints and further dependencies via the file `$HOME/polaris_gpu_warpx.profile`.
Create it now:

```bash
cp $HOME/src/warpx/Tools/machines/polaris-alcf/polaris_gpu_warpx.profile.example $HOME/polaris_gpu_warpx.profile
```

### Script Details

```bash
# Set the project name
export proj=""  # change me!

# swap to GNU programming environment (with gcc 13.2)
module load PrgEnv-gnu

# swap to the Milan cray package
module load craype-x86-milan

# extra modules
module use /soft/modulefiles
module load spack-pe-gnu

# add cuda
module load cuda/12.6
module load cudatoolkit-standalone/12.6
module load craype-accel-nvidia80

# required dependencies
module load cmake

# optional: for QED support with detailed tables
module load boost

# optional: for openPMD and PSATD+RZ support
module load cray-hdf5-parallel
module load cray-libsci/25.03.0
export CMAKE_PREFIX_PATH=/home/${USER}/sw/polaris/gpu/c-blosc-1.21.1:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=/home/${USER}/sw/polaris/gpu/adios2-2.10.2:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=/home/${USER}/sw/polaris/gpu/blaspp-2024.05.31:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=/home/${USER}/sw/polaris/gpu/lapackpp-2024.05.31:$CMAKE_PREFIX_PATH

export LD_LIBRARY_PATH=/home/${USER}/sw/polaris/gpu/c-blosc-1.21.1/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/${USER}/sw/polaris/gpu/adios2-2.10.2/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/${USER}/sw/polaris/gpu/blaspp-2024.05.31/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/${USER}/sw/polaris/gpu/lapackpp-2024.05.31/lib64:$LD_LIBRARY_PATH

export PATH=/home/${USER}/sw/polaris/gpu/adios2-2.10.2/bin:${PATH}

# optional: for Python bindings or libEnsemble
module load cray-python/3.11.7

if [ -d "/home/${USER}/sw/polaris/gpu/venvs/warpx" ]
then
  source /home/${USER}/sw/polaris/gpu/venvs/warpx/bin/activate
fi

# necessary to use CUDA-Aware MPI and run a job
export CRAY_ACCEL_TARGET=nvidia80

# optimize CUDA compilation for A100
export AMREX_CUDA_ARCH=8.0

# optimize CPU microarchitecture for AMD EPYC 3rd Gen (Milan/Zen3)
# note: the cc/CC/ftn wrappers below add those
export CXXFLAGS="-march=znver3"
export CFLAGS="-march=znver3"

# compiler environment hints
export CC=$(which gcc-12)
export CXX=$(which g++-12)
export CUDACXX=$(which nvcc)
export CUDAHOSTCXX=${CXX}
```

Edit the 2nd line of this script, which sets the `export proj=""` variable.
For example, if you are member of the project `proj_name`, then run `nano $HOME/polaris_gpu_warpx.profile` and edit line 2 to read:

```bash
export proj="proj_name"
```

Exit the `nano` editor with `Ctrl` + `O` (save) and then `Ctrl` + `X` (exit).

#### IMPORTANT
Now, and as the first step on future logins to Polaris, activate these environment settings:

```bash
source $HOME/polaris_gpu_warpx.profile
```

Finally, since Polaris does not yet provide software modules for some of our dependencies, install them once:

```bash
bash $HOME/src/warpx/Tools/machines/polaris-alcf/install_gpu_dependencies.sh
source $HOME/sw/polaris/gpu/venvs/warpx/bin/activate
```

### Script Details

```bash
#!/bin/bash
#
# Copyright 2024-2025 The WarpX Community
#
# This file is part of WarpX.
#
# Authors: Axel Huebl, Roelof Groenewald
# License: BSD-3-Clause-LBNL

# Exit on first error encountered #############################################
#
set -eu -o pipefail

# Check: ######################################################################
#
#   Was polaris_gpu_warpx.profile sourced and configured correctly?
if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your polaris_gpu_warpx.profile file! Please edit its line 2 to continue!"; exit 1; fi

# Remove old dependencies #####################################################
#
SW_DIR="/home/${USER}/sw/polaris/gpu"
rm -rf ${SW_DIR}
mkdir -p ${SW_DIR}

# remove common user mistakes in python, located in .local instead of a venv
python3 -m pip uninstall -qq -y pywarpx
python3 -m pip uninstall -qq -y warpx
python3 -m pip uninstall -qqq -y mpi4py 2>/dev/null || true

# General extra dependencies ##################################################
#

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
cmake -S $HOME/src/c-blosc -B $HOME/src/c-blosc-pm-gpu-build -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF -DDEACTIVATE_AVX2=OFF -DCMAKE_INSTALL_PREFIX=${SW_DIR}/c-blosc-1.21.1
cmake --build $HOME/src/c-blosc-pm-gpu-build --target install --parallel 16
rm -rf $HOME/src/c-blosc-pm-gpu-build

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
cmake -S $HOME/src/adios2 -B $HOME/src/adios2-pm-gpu-build -DADIOS2_USE_Blosc=ON -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DADIOS2_USE_ZeroMQ=OFF -DCMAKE_INSTALL_PREFIX=${SW_DIR}/adios2-2.10.2
cmake --build $HOME/src/adios2-pm-gpu-build --target install -j 16
rm -rf $HOME/src/adios2-pm-gpu-build

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
CXX=$(which CC) cmake -S $HOME/src/blaspp -B $HOME/src/blaspp-pm-gpu-build -Duse_openmp=OFF -Dgpu_backend=cuda -DCMAKE_CXX_STANDARD=20 -DCMAKE_INSTALL_PREFIX=${SW_DIR}/blaspp-2024.05.31
cmake --build $HOME/src/blaspp-pm-gpu-build --target install --parallel 16
rm -rf $HOME/src/blaspp-pm-gpu-build

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
CXX=$(which CC) CXXFLAGS="-DLAPACK_FORTRAN_ADD_" cmake -S $HOME/src/lapackpp -B $HOME/src/lapackpp-pm-gpu-build -DCMAKE_CXX_STANDARD=20 -Dbuild_tests=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX=${SW_DIR}/lapackpp-2024.05.31
cmake --build $HOME/src/lapackpp-pm-gpu-build --target install --parallel 16
rm -rf $HOME/src/lapackpp-pm-gpu-build

# Python ######################################################################
#
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade virtualenv
python3 -m pip cache purge
rm -rf ${SW_DIR}/venvs/warpx
python3 -m venv ${SW_DIR}/venvs/warpx
source ${SW_DIR}/venvs/warpx/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade build
python3 -m pip install --upgrade packaging
python3 -m pip install --upgrade wheel
python3 -m pip install --upgrade setuptools[core]
python3 -m pip install --upgrade cython
python3 -m pip install --upgrade numpy
python3 -m pip install --upgrade pandas
python3 -m pip install --upgrade scipy
MPICC="CC -target-accel=nvidia80 -shared" python3 -m pip install --upgrade mpi4py --no-cache-dir --no-build-isolation --no-binary mpi4py
python3 -m pip install --upgrade openpmd-api
python3 -m pip install --upgrade matplotlib
python3 -m pip install --upgrade yt
# install or update WarpX dependencies such as picmistandard
python3 -m pip install --upgrade -r $HOME/src/warpx/requirements.txt
python3 -m pip install cupy-cuda12x  # CUDA 12.6 compatible wheel
# optimas (based on libEnsemble & ax->botorch->gpytorch->pytorch)
python3 -m pip install --upgrade torch  # CUDA 12.6 compatible wheel
python3 -m pip install --upgrade optimas[all]
python3 -m pip install --upgrade lasy
```

### CPU Nodes

*Under construction*

<a id="building-polaris-compilation"></a>

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

*Under construction*

Now, you can [submit Polaris compute jobs](#running-cpp-polaris) for WarpX [Python (PICMI) scripts](../../usage/python.md#usage-picmi) ([example scripts](../../usage/examples.md#usage-examples)).
Or, you can use the WarpX executables to submit Polaris jobs ([example inputs](../../usage/examples.md#usage-examples)).
For executables, you can reference their location in your [job script](#running-cpp-polaris) or copy them to a location in `$PSCRATCH`.

<a id="building-polaris-update"></a>

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

- [update the polaris_gpu_warpx.profile or polaris_cpu_warpx files](#building-polaris-preparation),
- log out and into the system, activate the now updated environment profile as usual,
- [execute the dependency install scripts](#building-polaris-preparation).

As a last step, clean the build directory `rm -rf $HOME/src/warpx/build_pm_*` and rebuild WarpX.

<a id="running-cpp-polaris"></a>

## Running

### A100 (40GB) GPUs

The batch script below can be used to run a WarpX simulation on multiple nodes (change `<NODES>` accordingly) on the supercomputer Polaris at ALCF.

Replace descriptions between chevrons `<>` by relevant values, for instance `<input file>` could be `plasma_mirror_inputs`.
Note that we run one MPI rank per GPU.

```bash
#!/bin/bash -l

#PBS -A <proj>
#PBS -l select=<NODES>:system=polaris
#PBS -l place=scatter
#PBS -l walltime=0:10:00
#PBS -l filesystems=home:eagle
#PBS -q debug
#PBS -N test_warpx

# Set required environment variables
# support gpu-aware-mpi
# export MPICH_GPU_SUPPORT_ENABLED=1

# Change to working directory
echo Working directory is $PBS_O_WORKDIR
cd ${PBS_O_WORKDIR}

echo Jobid: $PBS_JOBID
echo Running on host `hostname`
echo Running on nodes `cat $PBS_NODEFILE`

# executable & inputs file or python interpreter & PICMI script here
EXE=./warpx
INPUTS=input1d

# MPI and OpenMP settings
NNODES=`wc -l < $PBS_NODEFILE`
NRANKS_PER_NODE=4
NDEPTH=1
NTHREADS=1

NTOTRANKS=$(( NNODES * NRANKS_PER_NODE ))
echo "NUM_OF_NODES= ${NNODES} TOTAL_NUM_RANKS= ${NTOTRANKS} RANKS_PER_NODE= ${NRANKS_PER_NODE} THREADS_PER_RANK= ${NTHREADS}"

mpiexec -np ${NTOTRANKS} ${EXE} ${INPUTS} > output.txt
```

To run a simulation, copy the lines above to a file `polaris_gpu.pbs` and run

```bash
qsub polaris_gpu.pbs
```

to submit the job.

### CPU Nodes

*Under construction*
