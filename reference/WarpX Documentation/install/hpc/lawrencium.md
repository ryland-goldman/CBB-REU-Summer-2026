<a id="building-lawrencium"></a>

# Lawrencium (LBNL)

The [Lawrencium cluster](http://scs.lbl.gov/Systems) is located at LBNL.

## Introduction

If you are new to this system, **please see the following resources**:

* [Lawrencium user guide](https://sites.google.com/a/lbl.gov/high-performance-computing-services-group/lbnl-supercluster/lawrencium)
* Batch system: [Slurm](https://sites.google.com/a/lbl.gov/high-performance-computing-services-group/scheduler/slurm-usage-instructions)
* [Production directories](https://sites.google.com/a/lbl.gov/high-performance-computing-services-group/lbnl-supercluster/lawrencium#backup):
  * `/global/scratch/users/$USER/`: production directory
  * `/global/home/groups/$GROUP/`: group production directory
  * `/global/home/users/$USER`: home directory (10 GB)

## Installation

Use the following commands to download the WarpX source code and switch to the correct branch:

```bash
git clone https://github.com/BLAST-WarpX/warpx.git $HOME/src/warpx
```

We use the following modules and environments on the system (`$HOME/lawrencium_warpx.profile`).

```bash
# please set your project account
#export proj="<yourProject>"  # change me, e.g., ac_blast

# required dependencies
module load cmake/3.27.7
module load gcc/11.4.0
module load cuda/12.2.1
module load openmpi/4.1.6

# optional: for QED support with detailed tables
module load boost/1.83.0

# optional: for openPMD and PSATD+RZ support
module load hdf5/1.14.3

export CMAKE_PREFIX_PATH=$HOME/sw/v100/c-blosc-1.21.1:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=$HOME/sw/v100/adios2-2.10.2:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=$HOME/sw/v100/blaspp-2024.05.31:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=$HOME/sw/v100/lapackpp-2024.05.31:$CMAKE_PREFIX_PATH

export PATH=$HOME/sw/v100/adios2-2.10.2/bin:$PATH

# optional: CCache
#module load ccache  # missing

# optional: for Python bindings or libEnsemble
module load python/3.11.6-gcc-11.4.0

if [ -d "$HOME/sw/v100/venvs/warpx" ]
then
  source $HOME/sw/v100/venvs/warpx/bin/activate
fi

# an alias to request an interactive batch node for one hour
#   for parallel execution, start on the batch node: srun <command>
alias getNode="salloc -N 1 -t 1:00:00 --qos=es_debug --partition=es1 --constraint=es1_v100 --gres=gpu:1 --cpus-per-task=4 -A $proj"
# an alias to run a command on a batch node for up to 30min
#   usage: runNode <command>
alias runNode="srun -N 1 -t 1:00:00 --qos=es_debug --partition=es1 --constraint=es1_v100 --gres=gpu:1 --cpus-per-task=4 -A $proj"

# optimize CUDA compilation for 1080 Ti (deprecated)
#export AMREX_CUDA_ARCH=6.1
# optimize CUDA compilation for V100
export AMREX_CUDA_ARCH=7.0
# optimize CUDA compilation for 2080 Ti
#export AMREX_CUDA_ARCH=7.5

# compiler environment hints
export CXX=$(which g++)
export CC=$(which gcc)
export FC=$(which gfortran)
export CUDACXX=$(which nvcc)
export CUDAHOSTCXX=${CXX}
```

We recommend to store the above lines in a file, such as `$HOME/lawrencium_warpx.profile`, and load it into your shell after a login:

```bash
source $HOME/lawrencium_warpx.profile
```

And since Lawrencium does not yet provide a module for them, install ADIOS2, BLAS++ and LAPACK++:

```bash
# c-blosc (I/O compression)
git clone -b v1.21.1 https://github.com/Blosc/c-blosc.git src/c-blosc
rm -rf src/c-blosc-v100-build
cmake -S src/c-blosc -B src/c-blosc-v100-build -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF -DDEACTIVATE_AVX2=OFF -DCMAKE_INSTALL_PREFIX=$HOME/sw/v100/c-blosc-1.21.1
cmake --build src/c-blosc-v100-build --target install --parallel 12

# ADIOS2
git clone -b v2.8.3 https://github.com/ornladios/ADIOS2.git src/adios2
rm -rf src/adios2-v100-build
cmake -S src/adios2 -B src/adios2-v100-build -DADIOS2_USE_Blosc=ON -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DADIOS2_USE_ZeroMQ=OFF -DCMAKE_INSTALL_PREFIX=$HOME/sw/v100/adios2-2.8.3
cmake --build src/adios2-v100-build --target install -j 12

# BLAS++ (for PSATD+RZ)
git clone https://github.com/icl-utk-edu/blaspp.git src/blaspp
rm -rf src/blaspp-v100-build
cmake -S src/blaspp -B src/blaspp-v100-build -Duse_openmp=OFF -Dgpu_backend=cuda -DCMAKE_CXX_STANDARD=20 -DCMAKE_INSTALL_PREFIX=$HOME/sw/v100/blaspp-master
cmake --build src/blaspp-v100-build --target install --parallel 12

# LAPACK++ (for PSATD+RZ)
git clone https://github.com/icl-utk-edu/lapackpp.git src/lapackpp
rm -rf src/lapackpp-v100-build
cmake -S src/lapackpp -B src/lapackpp-v100-build -DCMAKE_CXX_STANDARD=20 -Dgpu_backend=cuda -Dbuild_tests=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX=$HOME/sw/v100/lapackpp-master -Duse_cmake_find_lapack=ON -DBLAS_LIBRARIES=${LAPACK_DIR}/lib/libblas.a -DLAPACK_LIBRARIES=${LAPACK_DIR}/lib/liblapack.a
cmake --build src/lapackpp-v100-build --target install --parallel 12
```

Optionally, download and install Python packages for [PICMI](../../usage/python.md#usage-picmi) or dynamic ensemble optimizations ([libEnsemble](https://libensemble.readthedocs.io/en/main/)):

```bash
python3 -m pip install --user --upgrade pip
python3 -m pip install --user virtualenv
python3 -m pip cache purge
rm -rf $HOME/sw/v100/venvs/warpx
python3 -m venv $HOME/sw/v100/venvs/warpx
source $HOME/sw/v100/venvs/warpx/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade build
python3 -m pip install --upgrade packaging
python3 -m pip install --upgrade wheel
python3 -m pip install --upgrade setuptools[core]
python3 -m pip install --upgrade cython
python3 -m pip install --upgrade numpy
python3 -m pip install --upgrade pandas
python3 -m pip install --upgrade scipy
python3 -m pip install --upgrade mpi4py --no-build-isolation --no-binary mpi4py
python3 -m pip install --upgrade openpmd-api
python3 -m pip install --upgrade matplotlib
python3 -m pip install --upgrade yt
# optional: for libEnsemble
python3 -m pip install -r $HOME/src/warpx/Tools/LibEnsemble/requirements.txt
```

Then, `cd` into the directory `$HOME/src/warpx` and use the following commands to compile the application executable:

```bash
cd $HOME/src/warpx
rm -rf build

cmake -S . -B build -DWarpX_DIMS="1;2;RZ;3" -DWarpX_COMPUTE=CUDA -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON
cmake --build build -j 12
```

The general [cmake compile-time options](../cmake.md#install-build-cmake) apply as usual.

**That’s it!**
A 3D WarpX executable is now in `build/bin/` and [can be run](#running-cpp-lawrencium) with a [3D example inputs file](../../usage/examples.md#usage-examples).
Most people execute the binary directly or copy it out to a location in `/global/scratch/users/$USER/`.

For a *full PICMI install*, follow the [instructions for Python (PICMI) bindings](../cmake.md#install-build-python-cmake):

```bash
# PICMI build
cd $HOME/src/warpx

# install or update dependencies
python3 -m pip install -r requirements.txt

# compile parallel PICMI interfaces in 3D, 2D, 1D and RZ
WARPX_MPI=ON WARPX_COMPUTE=CUDA WARPX_FFT=ON BUILD_PARALLEL=12 python3 -m pip install --force-reinstall --no-deps -v .
```

Or, if you are *developing*, do a quick PICMI install of a *single geometry* (see: [WarpX_DIMS](../cmake.md#install-build-options)) using:

```bash
# find dependencies & configure
cmake -S . -B build -DWarpX_COMPUTE=CUDA -DWarpX_FFT=ON -DWarpX_APP=OFF -DWarpX_PYTHON=ON -DWarpX_DIMS=RZ

# build and then call "python3 -m pip install ..."
cmake --build build --target pip_install -j 12
```

<a id="running-cpp-lawrencium"></a>

## Running

<a id="running-cpp-lawrencium-v100-gpus"></a>

### V100 GPUs (16 GB)

12 nodes with each two NVIDIA V100 GPUs.

```bash
#!/bin/bash -l

# Copyright 2023 The WarpX Community
#
# Author: Axel Huebl
# License: BSD-3-Clause-LBNL

#SBATCH -t 00:10:00
#SBATCH -N 2
#SBATCH --job-name=WarpX
#SBATCH --account=<proj>
#SBATCH --qos=es_normal
# 2xV100 nodes
#SBATCH --partition=es1
#SBATCH --constraint=es1_v100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH -o WarpX.o%j
#SBATCH -e WarpX.e%j
#S BATCH --mail-type=all
#S BATCH --mail-user=yourmail@lbl.gov

# executable & inputs file or python interpreter & PICMI script here
EXE=./warpx
INPUTS=inputs_3d

srun ${EXE} ${INPUTS} \
  > output_${SLURM_JOB_ID}.txt
```

To run a simulation, copy the lines above to a file `v100.sbatch` and run

```bash
sbatch lawrencium_v100.sbatch
```

<a id="running-cpp-lawrencium-2080ti-gpus"></a>

### 2080 Ti GPUs (10 GB)

18 nodes with each four NVIDIA 2080 TI GPUs.
These are most interesting if you run in single precision.

Use `--constraint=es1_2080ti --cpus-per-task=2` in the above template to run on those nodes.
