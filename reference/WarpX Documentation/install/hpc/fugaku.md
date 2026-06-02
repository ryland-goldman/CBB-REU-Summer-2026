<a id="building-fugaku"></a>

# Fugaku (Riken)

The [Fugaku cluster](https://docs.nersc.gov/systems/perlmutter/) is located at the Riken Center for Computational Science (Japan).

## Introduction

If you are new to this system, **please see the following resources**:

* [Fugaku user guide](https://www.r-ccs.riken.jp/en/fugaku/user-guide/)

<a id="building-fugaku-preparation"></a>

## Preparation

Use the following commands to download the WarpX source code and switch to the correct branch:

```bash
git clone https://github.com/BLAST-WarpX/warpx.git $HOME/src/warpx
```

Compiling WarpX on Fugaku is more practical on a compute node. Use the following commands to acquire a compute node for one hour:

```bash
pjsub --interact -L "elapse=02:00:00" -L "node=1" --sparam "wait-time=300" --mpi "max-proc-per-node=48" --all-mount-gfscache
```

We use system software modules, add environment hints and further dependencies via the file `$HOME/fugaku_warpx.profile`.
Create it now, modify it if needed, and source it (it will take few minutes):

```bash
cp $HOME/src/warpx/Tools/machines/fugaku-riken/fugaku_warpx.profile.example $HOME/fugaku_warpx.profile
source $HOME/fugaku_warpx.profile
```

### Script Details

```bash
. /vol0004/apps/oss/spack/share/spack/setup-env.sh

# required dependencies
spack load cmake@3.24.3%fj@4.10.0 arch=linux-rhel8-a64fx

# avoid harmless warning messages "[WARN] xos LPG [...]"
export LD_LIBRARY_PATH=/lib64:$LD_LIBRARY_PATH

# optional: faster builds
spack load ninja@1.11.1%fj@4.10.0

# optional: for PSATD
spack load fujitsu-fftw@1.1.0%fj@4.10.0

# optional: for QED lookup table generation support
spack load boost@1.80.0%fj@4.8.1/zc5pwgc

# optional: for openPMD support
spack load hdf5@1.12.2%fj@4.8.1/im6lxev
export CMAKE_PREFIX_PATH=${HOME}/sw/fugaku/a64fx/c-blosc-1.21.1-install:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${HOME}/sw/fugaku/a64fx/adios2-2.10.2-install:$CMAKE_PREFIX_PATH

# compiler environment hints
export CC=$(which mpifcc)
export CXX=$(which mpiFCC)
export FC=$(which mpifrt)
export CFLAGS="-O3 -Nclang -Nlibomp -Klib -g -DNDEBUG"
export CXXFLAGS="-O3 -Nclang -Nlibomp -Klib -g -DNDEBUG"

# avoid harmless warning messages "[WARN] xos LPG [...]"
export LD_LIBRARY_PATH=/lib64:$LD_LIBRARY_PATH
```

Finally, since Fugaku does not yet provide software modules for some of our dependencies, install them once:

```bash
bash $HOME/src/warpx/Tools/machines/fugaku-riken/install_dependencies.sh
```

### Script Details

```bash
#!/bin/bash
#
# Copyright 2023 The WarpX Community
#
# This file is part of WarpX.
#
# Author: Axel Huebl, Luca Fedeli
# License: BSD-3-Clause-LBNL

# Exit on first error encountered #############################################
#
set -eu -o pipefail

# Remove old dependencies #####################################################
#
SRC_DIR="${HOME}/src/"
SW_DIR="${HOME}/sw/fugaku/a64fx/"
rm -rf ${SW_DIR}
mkdir -p ${SW_DIR}
mkdir -p ${SRC_DIR}

# General extra dependencies ##################################################
#

# c-blosc (I/O compression)
if [ -d ${SRC_DIR}/c-blosc ]
then
  cd ${SRC_DIR}/c-blosc
  git fetch --prune
  git checkout v1.21.1
  cd -
else
  git clone -b v1.21.1 https://github.com/Blosc/c-blosc.git ${SRC_DIR}/c-blosc
fi
  rm -rf ${SRC_DIR}/c-blosc-fugaku-build
  cmake -S ${SRC_DIR}/c-blosc -B ${SRC_DIR}/c-blosc-fugaku-build -DBUILD_SHARED_LIBS=OFF -DBUILD_SHARED=OFF -DBUILD_STATIC=ON -DBUILD_TESTS=OFF -DBUILD_FUZZERS=OFF -DBUILD_BENCHMARKS=OFF -DDEACTIVATE_AVX2=OFF -DCMAKE_INSTALL_PREFIX=${SW_DIR}/c-blosc-1.21.1-install
  cmake --build ${SRC_DIR}/c-blosc-fugaku-build --target install --parallel 48
  rm -rf ${SRC_DIR}/c-blosc-fugaku-build

# ADIOS2 (I/O)
if [ -d ${SRC_DIR}/c-blosc ]
then
  cd ${SRC_DIR}/adios2
  git fetch --prune
  git checkout v2.10.2
  cd -
else
  git clone -b v2.10.2 https://github.com/ornladios/ADIOS2.git ${SRC_DIR}/adios2
fi
rm -rf ${SRC_DIR}/adios2-fugaku-build
cmake -S ${SRC_DIR}/adios2 -B ${SRC_DIR}/adios2-fugaku-build -DBUILD_SHARED_LIBS=OFF -DADIOS2_USE_Blosc=ON -DBUILD_TESTING=OFF -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DADIOS2_USE_ZeroMQ=OFF -DCMAKE_INSTALL_PREFIX=${SW_DIR}/adios2-2.10.2-install
cmake --build ${SRC_DIR}/adios2-fugaku-build --target install -j 48
rm -rf ${SRC_DIR}/adios2-fugaku-build
```

<a id="building-fugaku-compilation"></a>

## Compilation

Use the following [cmake commands](../cmake.md#install-build-cmake) to compile the application executable:

```bash
cd $HOME/src/warpx
rm -rf build

export CC=$(which mpifcc)
export CXX=$(which mpiFCC)
export CFLAGS="-Nclang"
export CXXFLAGS="-Nclang"

cmake -S . -B build -DWarpX_COMPUTE=OMP \
    -DWarpX_DIMS="1;2;3" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS_RELEASE="-Ofast" \
    -DAMReX_DIFFERENT_COMPILER=ON \
    -DWarpX_MPI_THREAD_MULTIPLE=OFF

cmake --build build -j 48
```

**That’s it!**
A 3D WarpX executable is now in `build/bin/` and [can be run](#running-cpp-fugaku) with a [3D example inputs file](../../usage/examples.md#usage-examples).

<a id="running-cpp-fugaku"></a>

## Running

<a id="running-cpp-fugaku-a64fx-cpus"></a>

### A64FX CPUs

In non-interactive runs, you can use pjsub submit.sh where submit.sh can be adapted from:

```bash
#!/bin/bash
#PJM -L "node=48"
#PJM -L "rscgrp=small"
#PJM -L "elapse=0:30:00"
#PJM -s
#PJM -L "freq=2200,eco_state=2"
#PJM --mpi "max-proc-per-node=12"
#PJM -x PJM_LLIO_GFSCACHE=/vol0004:/vol0003
#PJM --llio localtmp-size=10Gi
#PJM --llio sharedtmp-size=10Gi

export NODES=48
export MPI_RANKS=$((NODES * 12))
export OMP_NUM_THREADS=4

export EXE="./warpx"
export INPUT="i.3d"

export XOS_MMM_L_PAGING_POLICY=demand:demand:demand

# Add HDF5 library path to LD_LIBRARY_PATH
# This is done manually to avoid calling spack during the run,
# since this would take a significant amount of time.
export LD_LIBRARY_PATH=/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/hdf5-1.12.2-im6lxevf76cu6cbzspi4itgz3l4gncjj/lib:$LD_LIBRARY_PATH

# Broadcast WarpX executable to all the nodes
llio_transfer ${EXE}

mpiexec -stdout-proc ./output.%j/%/1000r/stdout -stderr-proc ./output.%j/%/1000r/stderr -n ${MPI_RANKS} ${EXE} ${INPUT}

llio_transfer --purge ${EXE}
```

Note: the `Boost Eco Mode` mode that is set in this example increases the default frequency of the A64FX
from 2 GHz to 2.2 GHz, while at the same time switching off one of the two floating-point arithmetic
pipelines. Some preliminary tests with WarpX show that this mode achieves performances similar to those of
the normal mode but with a reduction of the energy consumption of approximately 20%.
