<a id="building-aurora"></a>

# Aurora (ALCF)

The [Aurora cluster](https://docs.alcf.anl.gov/aurora/) is located at ALCF.

## Introduction

If you are new to this system, **please see the following resources**:

* [ALCF user guide](https://docs.alcf.anl.gov/)
* Batch system: [PBS](https://docs.alcf.anl.gov/running-jobs/)
* [Filesystems](https://docs.alcf.anl.gov/data-management/filesystem-and-storage/):
  * `/lus/flare/projects/$proj/`: shared with all members of a project, Lustre

<a id="building-aurora-preparation"></a>

## Preparation

Use the following commands to download the WarpX source code:

```bash
git clone https://github.com/BLAST-WarpX/warpx.git $HOME/src/warpx
```

We use system software modules, add environment hints and further dependencies via the file `$HOME/aurora_warpx.profile`.
Create it now:

```bash
cp $HOME/src/warpx/Tools/machines/aurora-alcf/aurora_warpx.profile.example $HOME/aurora_warpx.profile
```

### Script Details

```bash
# Set the project name
export proj=""  # change me!

# required dependencies
module load cmake
module load adios2/2.10.2-cpu

# optional: for QED support with detailed tables
module load boost

# optional: for openPMD and PSATD+RZ support
module load hdf5/1.14.5
export CMAKE_PREFIX_PATH=/home/${USER}/sw/aurora/gpu/blaspp-2024.05.31:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=/home/${USER}/sw/aurora/gpu/lapackpp-2024.05.31:$CMAKE_PREFIX_PATH

export LD_LIBRARY_PATH=/home/${USER}/sw/aurora/gpu/blaspp-2024.05.31/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/${USER}/sw/aurora/gpu/lapackpp-2024.05.31/lib64:$LD_LIBRARY_PATH

# necessary for Python and AI/ML on Aurora:
module load frameworks

if [ -d "/home/${USER}/sw/aurora/gpu/venvs/warpx-aurora" ]
then
  source /home/${USER}/sw/aurora/gpu/venvs/warpx-aurora/bin/activate
fi

# necessary to use build or run with GPU-aware MPICH
export MPIR_CVAR_ENABLE_GPU=1
```

Edit the 2nd line of this script, which sets the `export proj=""` variable.
For example, if you are member of the project `proj_name`, then run `nano $HOME/aurora_warpx.profile` and edit line 2 to read:

```bash
export proj="proj_name"
```

Exit the `nano` editor with `Ctrl` + `O` (save) and then `Ctrl` + `X` (exit).

#### IMPORTANT
Now, and as the first step on future logins to Aurora, activate these environment settings:

```bash
source $HOME/aurora_warpx.profile
```

Finally, since Aurora does not yet provide software modules for some of our dependencies, install them once:

```bash
bash $HOME/src/warpx/Tools/machines/aurora-alcf/install_dependencies.sh
source /home/${USER}/sw/aurora/gpu/venvs/warpx-aurora/bin/activate
```

### Script Details

```bash
#!/bin/bash
#
# Copyright 2025 The WarpX Community
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
#   Was aurora_warpx.profile sourced and configured correctly?
if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your aurora_warpx.profile file! Please edit its line 2 to continue!"; exit 1; fi

# Remove old dependencies #####################################################
#
SW_DIR="/home/${USER}/sw/aurora/gpu"
rm -rf ${SW_DIR}
mkdir -p ${SW_DIR}

# remove common user mistakes in python, located in .local instead of a venv
python3 -m pip uninstall -qq -y pywarpx
python3 -m pip uninstall -qq -y warpx
python3 -m pip uninstall -qqq -y mpi4py 2>/dev/null || true

# General extra dependencies ##################################################
#

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
rm -rf $HOME/src/blaspp-aurora-gpu-build
CXX=icpx CXXFLAGS="-qmkl" cmake -S $HOME/src/blaspp -B $HOME/src/blaspp-aurora-gpu-build -Duse_openmp=OFF -Dgpu_backend=sycl -DCMAKE_CXX_STANDARD=20 -DCMAKE_INSTALL_PREFIX=${SW_DIR}/blaspp-2024.05.31 -DCMAKE_EXE_LINKER_FLAGS="-qmkl"
cmake --build $HOME/src/blaspp-aurora-gpu-build --target install --parallel 16
rm -rf $HOME/src/blaspp-aurora-gpu-build

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
rm -rf $HOME/src/lapackpp-aurora-gpu-build
CXX=icpx CXXFLAGS="-DLAPACK_FORTRAN_ADD_ -qmkl" cmake -S $HOME/src/lapackpp -B $HOME/src/lapackpp-aurora-gpu-build -DCMAKE_CXX_STANDARD=20 -Dbuild_tests=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX=${SW_DIR}/lapackpp-2024.05.31 -DCMAKE_EXE_LINKER_FLAGS="-qmkl"
cmake --build $HOME/src/lapackpp-aurora-gpu-build --target install --parallel 16
rm -rf $HOME/src/lapackpp-aurora-gpu-build

# Python ######################################################################
#
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade virtualenv
python3 -m pip cache purge
rm -rf ${SW_DIR}/venvs/warpx-aurora
python3 -m venv ${SW_DIR}/venvs/warpx-aurora
source ${SW_DIR}/venvs/warpx-aurora/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade build
python3 -m pip install --upgrade packaging
python3 -m pip install --upgrade wheel
python3 -m pip install --upgrade setuptools[core]
python3 -m pip install --upgrade "cython>=3.0"
python3 -m pip install --upgrade numpy
python3 -m pip install --upgrade h5py
python3 -m pip install --upgrade pandas
python3 -m pip install --upgrade scipy
# Is next line OK/needed for Aurora?
MPICC="mpicxx -fsycl -fsycl-targets=spir64_gen -Xsycl-target-backend \\\"-device pvc\\\" -shared" python3 -m pip install --upgrade mpi4py --no-cache-dir --no-build-isolation --no-binary mpi4py
python3 -m pip install --upgrade openpmd-api
python3 -m pip install --upgrade matplotlib
python3 -m pip install --upgrade yt
# install or update WarpX dependencies such as picmistandard
python3 -m pip install --upgrade -r $HOME/src/warpx/requirements.txt
# optional: for libEnsemble
#python3 -m pip install -r $HOME/src/warpx-aurora/Tools/LibEnsemble/requirements.txt
# optional: for optimas (based on libEnsemble & ax->botorch->gpytorch->pytorch)
#python3 -m pip install --upgrade torch  # should get from frameworks module on Aurora
#python3 -m pip install -r $HOME/src/warpx/Tools/optimas/requirements.txt
```

<a id="building-aurora-compilation"></a>

## Compilation

Use the following [cmake commands](../cmake.md#install-build-cmake) to compile the application executable:

```bash
cd $HOME/src/warpx
rm -rf build_aurora

cmake -S . -B build_aurora -DWarpX_COMPUTE=SYCL -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_aurora -j 16
```

The WarpX application executables are now in `$HOME/src/warpx/build_aurora/bin/`.
Additionally, the following commands will install WarpX as a Python module:

```bash
cd $HOME/src/warpx
rm -rf build_aurora_py

cmake -S . -B build_aurora_py -DWarpX_COMPUTE=SYCL -DWarpX_FFT=OFF -DWarpX_QED_TABLE_GEN=ON -DWarpX_APP=OFF -DWarpX_PYTHON=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_aurora_py -j 16 --target pip_install
```

Now, you can [submit Aurora compute jobs](#running-cpp-aurora) for WarpX [Python (PICMI) scripts](../../usage/python.md#usage-picmi) ([example scripts](../../usage/examples.md#usage-examples)).
Or, you can use the WarpX executables to submit Aurora jobs ([example inputs](../../usage/examples.md#usage-examples)).
For executables, you can reference their location in your [job script](#running-cpp-aurora) or copy them to a location in `/lus/flare/projects/$proj/`.

<a id="building-aurora-update"></a>

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

- [update the aurora_warpx.profile or aurora_cpu_warpx files](#building-aurora-preparation),
- log out and into the system, activate the now updated environment profile as usual,
- [execute the dependency install scripts](#building-aurora-preparation).

As a last step, clean the build directory `rm -rf $HOME/src/warpx/build_aurora*` and rebuild WarpX.

<a id="running-cpp-aurora"></a>

## Running

The batch script below can be used to run a WarpX simulation on multiple nodes (change `<NODES>` accordingly) on the supercomputer Aurora at ALCF.

Replace descriptions between chevrons `<>` by relevant values, for instance `<input file>` could be `plasma_mirror_inputs`.
Note that we run one MPI rank per GPU.

```bash
#!/bin/bash -l

#PBS -A <proj>
#PBS -l select=<NODES>:system=aurora
#PBS -W run_count=17
#PBS -l walltime=0:10:00
#PBS -l filesystems=home:flare
#PBS -q debug
#PBS -N test_warpx

# Set required environment variables
# support gpu-aware-mpi
# export MPIR_CVAR_ENABLE_GPU=1

# Change to working directory
echo Working directory is $PBS_O_WORKDIR
cd ${PBS_O_WORKDIR}

echo Jobid: $PBS_JOBID
echo Running on host `hostname`
echo Running on nodes `cat $PBS_NODEFILE`

# On Aurora, must load module environment in job script:
source $HOME/aurora_warpx.profile

# executable & inputs file or python interpreter & PICMI script here
EXE=./warpx
INPUTS=input1d

# MPI and OpenMP settings
NNODES=`wc -l < $PBS_NODEFILE`
NRANKS_PER_NODE=12 # 1 rank per PVC tile (a.k.a. stack; like a gcd on an AMD GPU)
NTHREADS=1
# Avoid core 0 on socket 0 and core 52 on socket 1, where OS threads run:
CPU_RANK_BIND="--cpu-bind=list:1-8:9-16:17-24:25-32:33-40:41-48:53-60:61-68:69-76:77-84:85-92:93-100"
# Convention: gpuNumber.tileNumber
GPU_RANK_BIND="--gpu-bind=list:0.0:0.1:1.0:1.1:2.0:2.1:3.0:3.1:4.0:4.1:5.0:5.1"

NTOTRANKS=$(( NNODES * NRANKS_PER_NODE ))
echo "NUM_OF_NODES= ${NNODES} TOTAL_NUM_RANKS= ${NTOTRANKS} RANKS_PER_NODE= ${NRANKS_PER_NODE} THREADS_PER_RANK= ${NTHREADS}"

mpiexec --np ${NTOTRANKS} -ppn ${NRANKS_PER_NODE} ${CPU_RANK_BIND} ${GPU_RANK_BIND} -envall ${EXE} ${INPUTS} > output.txt
```

To run a simulation, copy the lines above to a file `aurora.pbs` and run

```bash
qsub aurora.pbs
```

to submit the job.
