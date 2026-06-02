<a id="building-adastra"></a>

# Adastra (CINES)

The [Adastra cluster](https://www.cines.fr/calcul/adastra/) is located at CINES (France).
Each node contains 4 AMD MI250X GPUs, each with 2 Graphics Compute Dies (GCDs) for a total of 8 GCDs per node.
You can think of the 8 GCDs as 8 separate GPUs, each having 64 GB of high-bandwidth memory (HBM2E).

## Introduction

If you are new to this system, **please see the following resources**:

* [Adastra user guide](https://dci.dci-gitlab.cines.fr/webextranet/user_support/)
* Batch system: [Slurm](https://dci.dci-gitlab.cines.fr/webextranet/user_support/index.html?highlight=sbatch#running-jobs)
* [Production directories](https://dci.dci-gitlab.cines.fr/webextranet/data_and_storage/index.html#data-and-storage):
  * `$SHAREDSCRATCHDIR`: meant for short-term data storage, shared with all members of a project, purged every 30 days (17.6 TB default quota)
  * `$SCRATCHDIR`: meant for short-term data storage, single user, purged every 30 days
  * `$SHAREDWORKDIR`: meant for mid-term data storage, shared with all members of a project, never purged (4.76 TB default quota)
  * `$WORKDIR`: meant for mid-term data storage, single user, never purged
  * `$STORE` : meant for long term storage, single user, never purged, backed up
  * `$SHAREDHOMEDIR` : meant for scripts and tools, shared with all members of a project, never purged, backed up
  * `$HOME` : meant for scripts and tools, single user, never purged, backed up

<a id="building-adastra-preparation"></a>

## Preparation

The following instructions will install WarpX in the `$WORKDIR` directory.
On Adastra the Home folder has relatively stringent space quota and inode quota shared
among all the members of a given project. Therefore, installing WarpX in a different location
is advisable.

Use the following commands to download the WarpX source code:

```bash
# If you have multiple projects, activate the project that you want to use with:
#
# myproject -a YOUR_PROJECT_NAME
#
git clone https://github.com/BLAST-WarpX/warpx.git $WORKDIR/src/warpx
```

We use system software modules, add environment hints and further dependencies via the file `$HOME/adastra_warpx.profile`.
Create it now:

```bash
cp $WORKDIR/src/warpx/Tools/machines/adastra-cines/adastra_warpx.profile.example $HOME/adastra_warpx.profile
```

### Script Details

```bash
# please set your project account and uncomment the following two lines
#export proj=your_project_id
#myproject -a $proj

#path of the directory where software is installed (by default equal to $WORKDIR/sw)
export SW_DIR=${WORKDIR}/sw/

# required dependencies
module purge
module load cpe/25.09
module load craype-accel-amd-gfx90a craype-x86-trento
module load rocm/6.4.3
module load PrgEnv-amd
module load cray-mpich/9.0.1
module load develop
module load cmake/4.0.3

# optional: faster builds
module load ninja

# optional: for PSATD in RZ geometry support
export CMAKE_PREFIX_PATH=${SW_DIR}/adastra/gpu/blaspp-2024.05.31:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/adastra/gpu/lapackpp-2024.05.31:$CMAKE_PREFIX_PATH
export LD_LIBRARY_PATH=${SW_DIR}/adastra/gpu/blaspp-2024.05.31/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/adastra/gpu/lapackpp-2024.05.31/lib64:$LD_LIBRARY_PATH

# optional: for QED lookup table generation support
export CMAKE_PREFIX_PATH=${SW_DIR}/adastra/gpu/boost-1.88.0:${CMAKE_PREFIX_PATH}
export LD_LIBRARY_PATH=${SW_DIR}/adastra/gpu/boost-1.88.0/lib:${LD_LIBRARY_PATH}

# optional: for openPMD support
module load cray-hdf5-parallel
export CMAKE_PREFIX_PATH=${SW_DIR}/adastra/gpu/c-blosc-2.23.0:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/adastra/gpu/adios2-2.10.2:$CMAKE_PREFIX_PATH
export PATH=${SW_DIR}/adastra/gpu/adios2-2.10.2/bin:${PATH}

# optional: for Python bindings or libEnsemble
module load cray-python/3.11.7

if [ -d "${SW_DIR}/adastra/gpu/venvs/warpx-adastra" ]
then
  source ${SW_DIR}/adastra/gpu/venvs/warpx-adastra/bin/activate
fi

# fix system defaults: do not escape $ with a \ on tab completion
shopt -s direxpand

# make output group-readable by default
umask 0027

# an alias to request an interactive batch node for one hour
# for paralle execution, start on the batch node: srun <command>
alias getNode="salloc --account=$proj --job-name=warpx --constraint=MI250 --nodes=1 --ntasks-per-node=8 --cpus-per-task=8 --gpus-per-node=8 --threads-per-core=1 --exclusive --time=01:00:00"
# note: to access a compute note it is required to get its name (look at the `NODELIST` column)
#    $ squeue -u $USER
# and then to ssh into the node:
#    $ ssh node_name

# GPU-aware MPI
export MPICH_GPU_SUPPORT_ENABLED=1

# optimize ROCm/HIP compilation for MI250X
export AMREX_AMD_ARCH=gfx90a

# compiler environment hints
export CC=$(which cc)
export CXX=$(which CC)
export FC=$(which ftn)
```

Edit the 2nd line of this script, which sets the `export proj=""` variable using a text editor
such as `nano`, `emacs`, or `vim` (all available by default on Adastra login nodes) and
uncomment the 3rd line (which sets `$proj` as the active project).

#### IMPORTANT
Now, and as the first step on future logins to Adastra, activate these environment settings:

```bash
source $HOME/adastra_warpx.profile
```

Finally, since Adastra does not yet provide software modules for some of our dependencies, install them once:

```bash
bash $WORKDIR/src/warpx/Tools/machines/adastra-cines/install_dependencies.sh
source $WORKDIR/sw/adastra/gpu/venvs/warpx-adastra/bin/activate
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


# Check: ######################################################################
#
#   Was perlmutter_gpu_warpx.profile sourced and configured correctly?
if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your adastra_warpx.profile file! Please edit its line 2 to continue!"; exit 1; fi


# Remove old dependencies #####################################################
#
SW_DIR="${WORKDIR}/sw/adastra/gpu"
SRC_DIR="${WORKDIR}/src"
rm -rf ${SW_DIR}
mkdir -p ${SW_DIR}

# remove common user mistakes in python, located in .local instead of a venv
python3 -m pip uninstall -qq -y pywarpx
python3 -m pip uninstall -qq -y warpx
python3 -m pip uninstall -qqq -y mpi4py 2>/dev/null || true


# General extra dependencies ##################################################
#

# define how many threads are used for compilation
PARALLEL=16

# BLAS++ (for PSATD+RZ)
if [ -d ${SRC_DIR}/blaspp ]
then
  cd ${SRC_DIR}/blaspp
  git fetch --prune
  git checkout v2024.05.31
  cd -
else
  git clone -b v2024.05.31 https://github.com/icl-utk-edu/blaspp.git ${SRC_DIR}/blaspp
fi
rm -rf ${SRC_DIR}/blaspp-adastra-gpu-build
cmake -S ${SRC_DIR}/blaspp -B ${SRC_DIR}/blaspp-adastra-gpu-build -Duse_openmp=OFF -Dgpu_backend=hip -DGPU_TARGETS=gfx90a  -DCMAKE_CXX_STANDARD=20 -DCMAKE_INSTALL_PREFIX=${SW_DIR}/blaspp-2024.05.31
cmake --build ${SRC_DIR}/blaspp-adastra-gpu-build --target install --parallel ${PARALLEL}
rm -rf ${SRC_DIR}/blaspp-adastra-gpu-build

# LAPACK++ (for PSATD+RZ)
if [ -d ${SRC_DIR}/lapackpp ]
then
  cd ${SRC_DIR}/lapackpp
  git fetch --prune
  git checkout v2024.05.31
  cd -
else
  git clone -b v2024.05.31 https://github.com/icl-utk-edu/lapackpp.git ${SRC_DIR}/lapackpp
fi
rm -rf ${SRC_DIR}/lapackpp-adastra-gpu-build
cmake -S ${SRC_DIR}/lapackpp -B ${SRC_DIR}/lapackpp-adastra-gpu-build -DGPU_TARGETS=gfx90a  -DCMAKE_CXX_STANDARD=20 -Dbuild_tests=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX=${SW_DIR}/lapackpp-2024.05.31
cmake --build ${SRC_DIR}/lapackpp-adastra-gpu-build --target install --parallel ${PARALLEL}
rm -rf ${SRC_DIR}/lapackpp-adastra-gpu-build

# Boost (QED tables)
rm -rf ${SRC_DIR}/boost-temp
mkdir -p ${SRC_DIR}/boost-temp
curl -Lo ${SRC_DIR}/boost-temp/boost.tar.gz https://archives.boost.io/release/1.88.0/source/boost_1_88_0.tar.gz
tar -xzf ${SRC_DIR}/boost-temp/boost.tar.gz -C ${SRC_DIR}/boost-temp
cd ${SRC_DIR}/boost-temp/boost_1_88_0
./bootstrap.sh -with-libraries=math  --prefix=${SW_DIR}/boost-1.88.0
./b2 --with-math cxxflags="-std=c++17" install -j ${PARALLEL}
cd -
rm -rf ${SRC_DIR}/boost-temp

# c-blosc2 (I/O compression, for OpenPMD)
if [ -d ${SRC_DIR}/c-blosc2 ]
then
  # git repository is already there
  :
else
  git clone -b v2.23.0 https://github.com/Blosc/c-blosc2.git ${SRC_DIR}/c-blosc2
fi
rm -rf ${SRC_DIR}/c-blosc2-ad-build
cmake -S ${SRC_DIR}/c-blosc2 -B ${SRC_DIR}/c-blosc2-ad-build -DBUILD_TESTS=OFF -DBUILD_EXAMPLES=OFF  -DBUILD_FUZZERS=OFF -DBUILD_BENCHMARKS=OFF -DDEACTIVATE_AVX2=OFF -DCMAKE_INSTALL_PREFIX=${SW_DIR}/c-blosc-2.23.0
cmake --build ${SRC_DIR}/c-blosc2-ad-build --target install --parallel ${PARALLEL}
rm -rf ${SRC_DIR}/c-blosc2-ad-build

# ADIOS2 v. 2.10.2 (for OpenPMD)
if [ -d ${SRC_DIR}/adios2 ]
then
  cd ${SRC_DIR}/adios2
  git fetch --prune
  git checkout v2.10.2
  cd -
else
  git clone -b v2.10.2 https://github.com/ornladios/ADIOS2.git ${SRC_DIR}/adios2
fi
rm -rf ${SRC_DIR}/adios2-ad-build
cmake -S ${SRC_DIR}/adios2 -B ${SRC_DIR}/adios2-ad-build -DADIOS2_USE_Blosc2=ON -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DADIOS2_USE_ZeroMQ=OFF -DCMAKE_INSTALL_PREFIX=${SW_DIR}/adios2-2.10.2
cmake --build ${SRC_DIR}/adios2-ad-build --target install -j ${PARALLEL}
rm -rf ${SRC_DIR}/adios2-ad-build


# Python ######################################################################
#
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade virtualenv
python3 -m pip cache purge
rm -rf ${SW_DIR}/venvs/warpx-adastra
python3 -m venv ${SW_DIR}/venvs/warpx-adastra
source ${SW_DIR}/venvs/warpx-adastra/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade build
python3 -m pip install --upgrade packaging
python3 -m pip install --upgrade wheel
python3 -m pip install --upgrade setuptools[core]
python3 -m pip install --upgrade cython
python3 -m pip install --upgrade numpy
python3 -m pip install --upgrade pandas
python3 -m pip install --upgrade scipy
python3 -m pip install --upgrade jupyter
MPICC="cc -shared" python3 -m pip install --upgrade mpi4py --no-cache-dir --no-build-isolation --no-binary mpi4py
python3 -m pip install --upgrade openpmd-api
python3 -m pip install --upgrade matplotlib
python3 -m pip install --upgrade yt
python3 -m pip install --upgrade openpmd-viewer
python3 -m pip install --upgrade adios2
# install or update WarpX dependencies such as picmistandard
python3 -m pip install --upgrade -r ${SRC_DIR}/warpx/requirements.txt
# optional: for optimas (based on libEnsemble & ax->botorch->gpytorch->pytorch)
python3 -m pip install 'optimas[all]'
# optional: for lasy
python3 -m pip install --upgrade lasy
```

<a id="building-adastra-compilation"></a>

## Compilation

Use the following [cmake commands](../cmake.md#install-build-cmake) to compile the application executable:

```bash
cd $WORKDIR/src/warpx
rm -rf build_adastra

cmake -S . -B build_adastra -DWarpX_COMPUTE=HIP -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_DIMS="1;2;RZ;3" -DWarpX_QED_TABLES_GEN_OMP=OFF
cmake --build build_adastra -j 16
```

The WarpX application executables are now in `$WORKDIR/src/warpx/build_adastra/bin/`.
Additionally, the following commands will install WarpX as a Python module:

```bash
rm -rf build_adastra_py

cmake -S . -B build_adastra_py -DWarpX_COMPUTE=HIP -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_APP=OFF -DWarpX_PYTHON=ON -DWarpX_DIMS="1;2;RZ;3" -DWarpX_QED_TABLES_GEN_OMP=OFF
cmake --build build_adastra_py -j 16 --target pip_install
```

#### NOTE
Enabling openMP support for QED lookup tables generation in WarpX while compiling WarpX for GPUs on Adastra does not work. It is recommended to generate QED lookup tables locally (e.g., with the standalone tool) and then transfer them to Adastra. The usage of the standalone tool for QED lookup tables generation is documented in the usage/workflows section of the documentation.

Now, you can [submit Adstra compute jobs](#running-cpp-adastra) for WarpX [Python (PICMI) scripts](../../usage/python.md#usage-picmi) ([example scripts](../../usage/examples.md#usage-examples)).
Or, you can use the WarpX executables to submit Adastra jobs ([example inputs](../../usage/examples.md#usage-examples)).
For executables, you can reference their location in your [job script](#running-cpp-adastra) .

<a id="building-adastra-update"></a>

## Update WarpX & Dependencies

If you already installed WarpX in the past and want to update it, start by getting the latest source code:

```bash
cd $WORKDIR/src/warpx

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

- [update the adastra_warpx.profile file](#building-adastra-preparation),
- log out and into the system, activate the now updated environment profile as usual,
- [execute the dependency install scripts](#building-adastra-preparation).

As a last step, clean the build directory `rm -rf $WORKDIR/src/warpx/build_adastra` and rebuild WarpX.

<a id="running-cpp-adastra"></a>

## Running

<a id="running-cpp-adastra-mi250x-gpus"></a>

### MI250X GPUs (2x64 GB)

In non-interactive runs:

```bash
#!/bin/bash
#SBATCH --account=<account_to_charge>
#SBATCH --job-name=warpx
#SBATCH --constraint=MI250
#SBATCH --nodes=2
#SBATCH --exclusive
#SBATCH --output=%x-%j.out
#SBATCH --time=00:10:00

module purge

# A CrayPE environment version
module load cpe/25.09
# An architecture
module load craype-accel-amd-gfx90a craype-x86-trento
# A compiler to target the architecture
module load rocm/6.4.3
module load PrgEnv-amd
# The MPI library
module load cray-mpich/9.0.1

date
module list

export MPICH_GPU_SUPPORT_ENABLED=1

# note
# this environment setting is needed to avoid that rocFFT writes a cache in
# the home directory, which does not scale.
export ROCFFT_RTC_CACHE_PATH=/dev/null

export OMP_NUM_THREADS=1
export WARPX_NMPI_PER_NODE=8
export TOTAL_NMPI=$(( ${SLURM_JOB_NUM_NODES} * ${WARPX_NMPI_PER_NODE} ))
srun -N${SLURM_JOB_NUM_NODES} -n${TOTAL_NMPI} --ntasks-per-node=${WARPX_NMPI_PER_NODE} \
     --cpus-per-task=8 --threads-per-core=1 --gpu-bind=closest \
    ./warpx inputs > output.txt
```

<a id="post-processing-adastra"></a>

## Post-Processing

#### NOTE
TODO: Document any Jupyter or data services.

<a id="known-adastra-issues"></a>

## Known System Issues

#### WARNING
May 16th, 2022:
There is a caching bug in Libfabric that causes WarpX simulations to occasionally hang on on more than 1 node.

As a work-around, please export the following environment variable in your job scripts until the issue is fixed:

```bash
#export FI_MR_CACHE_MAX_COUNT=0  # libfabric disable caching
# or, less invasive:
export FI_MR_CACHE_MONITOR=memhooks  # alternative cache monitor
```

#### WARNING
Sep 2nd, 2022:
rocFFT in ROCm 5.1-5.3 tries to [write to a cache](https://rocfft.readthedocs.io/en/latest/#runtime-compilation) in the home area by default.
This does not scale, disable it via:

```bash
export ROCFFT_RTC_CACHE_PATH=/dev/null
```

#### WARNING
January, 2023:
We discovered a regression in AMD ROCm, leading to 2x slower current deposition (and other slowdowns) in ROCm 5.3 and 5.4.
Reported to AMD and fixed for the next release of ROCm.

Stay with the ROCm 5.2 module to avoid.

#### WARNING
April 30th, 2025:
We observed several issues that can cause WarpX simulations to hang or crash on releases `25.02` and `25.03`.

Releases `<=25.01` and `>=25.04` are currently working.

#### WARNING
August 2025:
We observed a heavy node memory increase over time when using module `cray-mpich` versions `8.1.28` and `8.1.30`, which
causes simulations to slow down and eventually crash.

While no `cray-mpich` version `>8.1.30` is available on Adastra, stay with version `8.1.26` to avoid this issue.
