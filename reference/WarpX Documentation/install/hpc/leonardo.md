<a id="building-leonardo"></a>

# Leonardo (CINECA)

The [Leonardo cluster](https://leonardo-supercomputer.cineca.eu/) is hosted at [CINECA](https://www.cineca.it/en).

On Leonardo, each one of the 3456 compute nodes features a custom Atos Bull Sequana XH21355 “Da Vinci” blade, composed of:

* 1 x CPU Intel Ice Lake Xeon 8358 32 cores 2.60 GHz
* 512 (8 x 64) GB RAM DDR4 3200 MHz
* 4 x NVidia custom Ampere A100 GPU 64GB HBM2
* 2 x NVidia HDR 2×100 GB/s cards

## Introduction

If you are new to this system, **please see the following resources**:

* [Leonardo website](https://leonardo-supercomputer.cineca.eu/)
* [Leonardo user guide](https://docs.hpc.cineca.it/hpc/leonardo.html)

Storage organization:

* `$HOME`: permanent, backed up, user specific (50 GB quota)
* `$CINECA_SCRATCH`: temporary, user specific, no backup, a large disk for the storage of run time data and files, automatic cleaning procedure of data older than 40 days
* `$PUBLIC`: permanent, no backup (50 GB quota)
* `$WORK`: permanent, project specific, no backup

<a id="building-leonardo-preparation"></a>

## Preparation

Use the following commands to download the WarpX source code:

```bash
git clone https://github.com/BLAST-WarpX/warpx.git $HOME/src/warpx
```

We use system software modules, add environment hints and further dependencies via the file `$HOME/leonardo_gpu_warpx.profile`.
Create it now:

```bash
cp $HOME/src/warpx/Tools/machines/leonardo-cineca/leonardo_gpu_warpx.profile.example $HOME/leonardo_gpu_warpx.profile
```

### Script Details

```bash
# required dependencies
module load profile/base
module load cmake/3.24.3
module load gmp/6.2.1
module load mpfr/4.1.0
module load mpc/1.2.1
module load gcc/11.3.0
module load cuda/11.8
module load zlib/1.2.13--gcc--11.3.0
module load openmpi/4.1.4--gcc--11.3.0-cuda-11.8

# optional: for QED support with detailed tables
module load boost/1.80.0--openmpi--4.1.4--gcc--11.3.0

# optional: for openPMD and PSATD+RZ support
module load openblas/0.3.21--gcc--11.3.0
export CMAKE_PREFIX_PATH=/leonardo/prod/spack/03/install/0.19/linux-rhel8-icelake/gcc-11.3.0/c-blosc-1.21.1-aifmix6v5lwxgt7rigwoebalrgbcnv26:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=$HOME/sw/adios2-2.10.2:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=$HOME/sw/blaspp-2024.05.31:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=$HOME/sw/lapackpp-2024.05.31:$CMAKE_PREFIX_PATH

export LD_LIBRARY_PATH=/leonardo/prod/spack/03/install/0.19/linux-rhel8-icelake/gcc-11.3.0/c-blosc-1.21.1-aifmix6v5lwxgt7rigwoebalrgbcnv26/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/sw/adios2-2.10.2/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/sw/blaspp-2024.05.31/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/sw/lapackpp-2024.05.31/lib64:$LD_LIBRARY_PATH

export PATH=$HOME/sw/adios2-2.10.2/bin:$PATH

# optional: for Python bindings or libEnsemble
module load python/3.10.8--gcc--11.3.0

if [ -d "$HOME/sw/venvs/warpx" ]
then
  source $HOME/sw/venvs/warpx/bin/activate
fi

# optimize CUDA compilation for A100
export AMREX_CUDA_ARCH=8.0

# compiler environment hints
export CXX=$(which g++)
export CC=$(which gcc)
export FC=$(which gfortran)
export CUDACXX=$(which nvcc)
export CUDAHOSTCXX=${CXX}
```

#### IMPORTANT
Now, and as the first step on future logins to Leonardo, activate these environment settings:

```bash
source $HOME/leonardo_gpu_warpx.profile
```

Finally, since Leonardo does not yet provide software modules for some of our dependencies, install them once:

```bash
bash $HOME/src/warpx/Tools/machines/leonardo-cineca/install_gpu_dependencies.sh
source $HOME/sw/venvs/warpx/bin/activate
```

### Script Details

```bash
#!/bin/bash
#
# Copyright 2023 The WarpX Community
#
# This file is part of WarpX.
#
# Author: Axel Huebl, Marta Galbiati
# License: BSD-3-Clause-LBNL

set -eu -o pipefail


# Check: ######################################################################
#
#   Was leonardo_gpu_warpx.profile sourced and configured correctly?
#


# Remove old dependencies #####################################################
#
SW_DIR="$HOME/sw"
rm -rf ${SW_DIR}
mkdir -p ${SW_DIR}


# General extra dependencies ##################################################
#

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
rm -rf $HOME/src/adios2-gpu-build
cmake -S $HOME/src/adios2 -B $HOME/src/adios2-gpu-build -DADIOS2_USE_Blosc=ON -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DADIOS2_USE_ZeroMQ=OFF -DCMAKE_INSTALL_PREFIX=${SW_DIR}/adios2-2.10.2
cmake --build $HOME/src/adios2-gpu-build --target install -j 16
rm -rf $HOME/src/adios2-gpu-build


# BLAS++ (for PSATD+RZ)
if [ -d $HOME/src/blaspp ]
then
  cd $HOME/src/blaspp
  git fetch
  git checkout v2024.05.31
  cd -
else
  git clone -b v2024.05.31 https://github.com/icl-utk-edu/blaspp.git $HOME/src/blaspp
fi
rm -rf $HOME/src/blaspp-gpu-build
CXX=$(which g++) cmake -S $HOME/src/blaspp -B $HOME/src/blaspp-gpu-build -Duse_openmp=OFF -Dgpu_backend=cuda -DCMAKE_CXX_STANDARD=20 -DCMAKE_INSTALL_PREFIX=${SW_DIR}/blaspp-2024.05.31
cmake --build $HOME/src/blaspp-gpu-build --target install --parallel 16
rm -rf $HOME/src/blaspp-gpu-build


# LAPACK++ (for PSATD+RZ)
if [ -d $HOME/src/lapackpp ]
then
  cd $HOME/src/lapackpp
  git fetch
  git checkout v2024.05.31
  cd -
else
  git clone -b v2024.05.31 https://github.com/icl-utk-edu/lapackpp.git $HOME/src/lapackpp
fi
rm -rf $HOME/src/lapackpp-gpu-build
CXX=$(which CC) CXXFLAGS="-DLAPACK_FORTRAN_ADD_" cmake -S $HOME/src/lapackpp -B $HOME/src/lapackpp-gpu-build -DCMAKE_CXX_STANDARD=20 -Dbuild_tests=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX=${SW_DIR}/lapackpp-2024.05.31
cmake --build $HOME/src/lapackpp-gpu-build --target install --parallel 16
rm -rf $HOME/src/lapackpp-gpu-build


# Python ######################################################################
#
rm -rf ${SW_DIR}/venvs/warpx
python3 -m venv ${SW_DIR}/venvs/warpx
source ${SW_DIR}/venvs/warpx/bin/activate
python3 -m ensurepip --upgrade
python3 -m pip cache purge
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade build
python3 -m pip install --upgrade packaging
python3 -m pip install --upgrade wheel
python3 -m pip install --upgrade setuptools[core]
python3 -m pip install --upgrade cython
python3 -m pip install --upgrade numpy
python3 -m pip install --upgrade pandas
python3 -m pip install --upgrade scipy
MPICC="gcc -shared" python3 -m pip install --upgrade mpi4py --no-cache-dir --no-build-isolation --no-binary mpi4py
python3 -m pip install --upgrade openpmd-api
python3 -m pip install --upgrade matplotlib
python3 -m pip install --upgrade yt
# install or update WarpX dependencies such as picmistandard
python3 -m pip install --upgrade -r $HOME/src/warpx/requirements.txt
# optional: for libEnsemble
python3 -m pip install -r $HOME/src/warpx/Tools/LibEnsemble/requirements.txt
# optional: for optimas (based on libEnsemble & ax->botorch->gpytorch->pytorch)
python3 -m pip install --upgrade torch  # CUDA 11.8 compatible wheel
python3 -m pip install -r $HOME/src/warpx/Tools/optimas/requirements.txt
```

<a id="building-leonardo-compilation"></a>

## Compilation

Use the following [cmake commands](../cmake.md#install-build-cmake) to compile the application executable:

```bash
cd $HOME/src/warpx
rm -rf build_gpu

cmake -S . -B build_gpu -DWarpX_COMPUTE=CUDA -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_gpu -j 16
```

The WarpX application executables are now in `$HOME/src/warpx/build_gpu/bin/`.
Additionally, the following commands will install WarpX as a Python module:

```bash
cd $HOME/src/warpx
rm -rf build_gpu_py

cmake -S . -B build_gpu_py -DWarpX_COMPUTE=CUDA -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_PYTHON=ON -DWarpX_APP=OFF -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_gpu_py -j 16 --target pip_install
```

Now, you can [submit Leonardo compute jobs](#running-leonardo) for WarpX [Python (PICMI) scripts](../../usage/python.md#usage-picmi) ([example scripts](../../usage/examples.md#usage-examples)).
Or, you can use the WarpX executables to submit Leonardo jobs ([example inputs](../../usage/examples.md#usage-examples)).
For executables, you can reference their location in your [job script](#running-leonardo) or copy them to a location in `$CINECA_SCRATCH`.

<a id="building-leonardo-update"></a>

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
git log     # press q to exit
```

And, if needed,

- [update the leonardo_gpu_warpx.profile file](#building-leonardo-preparation),
- log out and into the system, activate the now updated environment profile as usual,
- [execute the dependency install scripts](#building-leonardo-preparation).

As a last step, clean the build directories `rm -rf $HOME/src/warpx/build_gpu*` and rebuild WarpX.

<a id="running-leonardo"></a>

## Running

The batch script below can be used to run a WarpX simulation on multiple nodes on Leonardo.
Replace descriptions between chevrons `<>` by relevant values.
Note that we run one MPI rank per GPU.

```bash
#!/usr/bin/bash
#SBATCH --time=02:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --ntasks-per-socket=4
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=4
#SBATCH --gpus-per-task=1
#SBATCH --mem=494000
#SBATCH --partition=boost_usr_prod
#SBATCH --job-name=<job name>
#SBATCH --gres=gpu:4
#SBATCH --err=job.err
#SBATCH --out=job.out
#SBATCH --account=<project id>
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<mail>

cd /leonardo_scratch/large/userexternal/<username>/<directory>
srun /leonardo/home/userexternal/<username>/src/warpx/build_gpu/bin/warpx.2d <input file> > output.txt
```

To run a simulation, copy the lines above to a file `job.sh` and run

```bash
sbatch job.sh
```

to submit the job.

<a id="post-processing-leonardo"></a>

## Post-Processing

For post-processing, activate the environment settings:

```bash
source $HOME/leonardo_gpu_warpx.profile
```

and run python scripts.
