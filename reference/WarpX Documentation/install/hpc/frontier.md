<a id="building-frontier"></a>

# Frontier (OLCF)

The [Frontier cluster](https://www.olcf.ornl.gov/frontier/) is located at OLCF.

On Frontier, each compute node provides four AMD MI250X GPUs, each with two Graphics Compute Dies (GCDs) for a total of 8 GCDs per node.
You can think of the 8 GCDs as 8 separate GPUs, each having 64 GB of high-bandwidth memory (HBM2E).

## Introduction

If you are new to this system, **please see the following resources**:

* [Frontier user guide](https://docs.olcf.ornl.gov/systems/frontier_user_guide.html)
* Batch system: [Slurm](https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#running-jobs)
* [Filesystems](https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#data-and-storage):
  * `$HOME`: per-user directory, use only for inputs, source and scripts; backed up; mounted as read-only on compute nodes, that means you cannot run in it (50 GB quota)
  * `$PROJWORK/$proj/`: shared with all members of a project, purged every 90 days, Lustre (recommended)
  * `$MEMBERWORK/$proj/`: single user, purged every 90 days, Lustre (usually smaller quota, 50TB default quota)
  * `$WORLDWORK/$proj/`: shared with all users, purged every 90 days, Lustre (50TB default quota)

Note: the Orion Lustre filesystem on Frontier and the older Alpine GPFS filesystem on Summit are not mounted on each others machines.
Use [Globus](https://www.globus.org) to transfer data between them if needed.

<a id="building-frontier-preparation"></a>

## Preparation

Use the following commands to download the WarpX source code:

```bash
git clone https://github.com/BLAST-WarpX/warpx.git $HOME/src/warpx
```

We use system software modules, add environment hints and further dependencies via the file `$HOME/frontier_warpx.profile`.
Create it now:

```bash
cp $HOME/src/warpx/Tools/machines/frontier-olcf/frontier_warpx.profile.example $HOME/frontier_warpx.profile
```

### Script Details

```bash
# please set your project account
export proj=""  # change me!

# remembers the location of this script
export MY_PROFILE=$(cd $(dirname $BASH_SOURCE) && pwd)"/"$(basename $BASH_SOURCE)
if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your $MY_PROFILE file! Please edit its line 2 to continue!"; return; fi

# required dependencies
module load cmake/3.30.5
module load craype-accel-amd-gfx90a
module load rocm/6.2.4
module load cray-mpich/8.1.31
module load cce/18.0.1  # must be loaded after rocm
# https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#compatible-compiler-rocm-toolchain-versions

# Fix for OpenMP Runtime (OLCFHELP-21543)
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ROCM_PATH}/llvm/lib

# optional: faster builds
module load ccache
module load ninja

# optional: just an additional text editor
module load nano

# optional: for PSATD in RZ geometry support
export CMAKE_PREFIX_PATH=${HOME}/sw/frontier/gpu/blaspp-2024.05.31:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${HOME}/sw/frontier/gpu/lapackpp-2024.05.31:$CMAKE_PREFIX_PATH
export LD_LIBRARY_PATH=${HOME}/sw/frontier/gpu/blaspp-2024.05.31/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${HOME}/sw/frontier/gpu/lapackpp-2024.05.31/lib64:$LD_LIBRARY_PATH

# optional: for QED lookup table generation support
module load boost/1.85.0

# optional: for openPMD support
module load adios2/2.10.2-mpi
module load hdf5/1.14.3-mpi

# optional: for Python bindings or libEnsemble
module load cray-python/3.11.5

if [ -d "${HOME}/sw/frontier/gpu/venvs/warpx-frontier" ]
then
  source ${HOME}/sw/frontier/gpu/venvs/warpx-frontier/bin/activate
fi

# fix system defaults: do not escape $ with a \ on tab completion
shopt -s direxpand

# make output group-readable by default
umask 0027

# an alias to request an interactive batch node for one hour
#   for paralle execution, start on the batch node: srun <command>
alias getNode="salloc -A $proj -J warpx -t 01:00:00 -p batch -N 1"
# an alias to run a command on a batch node for up to 30min
#   usage: runNode <command>
alias runNode="srun -A $proj -J warpx -t 00:30:00 -p batch -N 1"

# GPU-aware MPI
export MPICH_GPU_SUPPORT_ENABLED=1

# optimize ROCm/HIP compilation for MI250X
export AMREX_AMD_ARCH=gfx90a

# compiler environment hints
export CC=$(which hipcc)
export CXX=$(which hipcc)
export FC=$(which ftn)
export CFLAGS="-I${ROCM_PATH}/include"
export CXXFLAGS="-I${ROCM_PATH}/include -Wno-pass-failed"
export LDFLAGS="-L${ROCM_PATH}/lib -lamdhip64 ${PE_MPICH_GTL_DIR_amd_gfx90a} -lmpi_gtl_hsa"
```

Edit the 2nd line of this script, which sets the `export proj=""` variable.
For example, if you are member of the project `aph114`, then run `vi $HOME/frontier_warpx.profile`.
Enter the edit mode by typing `i` and edit line 2 to read:

```bash
export proj="aph114"
```

Exit the `vi` editor with `Esc` and then type `:wq` (write & quit).

#### IMPORTANT
Now, and as the first step on future logins to Frontier, activate these environment settings:

```bash
source $HOME/frontier_warpx.profile
```

Finally, since Frontier does not yet provide software modules for some of our dependencies, install them once:

```bash
bash $HOME/src/warpx/Tools/machines/frontier-olcf/install_dependencies.sh
source $HOME/sw/frontier/gpu/venvs/warpx-frontier/bin/activate
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
#   Was frontier_warpx.profile sourced and configured correctly?
if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your frontier_warpx.profile file! Please edit its line 2 to continue!"; exit 1; fi


# Check $proj variable is correct and has a corresponding CFS directory #######
#
if [ ! -d "${PROJWORK}/${proj}/" ]
then
    echo "WARNING: The directory $PROJWORK/$proj/ does not exist!"
    echo "Is the \$proj environment variable of value \"$proj\" correctly set? "
    echo "Please edit line 2 of your frontier_warpx.profile file to continue!"
    exit
fi


# Remove old dependencies #####################################################
#
SW_DIR="${HOME}/sw/frontier/gpu"
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
rm -rf $HOME/src/blaspp-frontier-gpu-build
CXX=$(which CC) cmake -S $HOME/src/blaspp -B $HOME/src/blaspp-frontier-gpu-build -Duse_openmp=OFF -Dgpu_backend=hip -DCMAKE_CXX_STANDARD=20 -DCMAKE_INSTALL_PREFIX=${SW_DIR}/blaspp-2024.05.31
cmake --build $HOME/src/blaspp-frontier-gpu-build --target install --parallel 16
rm -rf $HOME/src/blaspp-frontier-gpu-build

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
rm -rf $HOME/src/lapackpp-frontier-gpu-build
CXX=$(which CC) CXXFLAGS="-DLAPACK_FORTRAN_ADD_" cmake -S $HOME/src/lapackpp -B $HOME/src/lapackpp-frontier-gpu-build -DCMAKE_CXX_STANDARD=20 -Dbuild_tests=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX=${SW_DIR}/lapackpp-2024.05.31
cmake --build $HOME/src/lapackpp-frontier-gpu-build --target install --parallel 16
rm -rf $HOME/src/lapackpp-frontier-gpu-build

# Python ######################################################################
#
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade virtualenv
python3 -m pip cache purge
rm -rf ${SW_DIR}/venvs/warpx-frontier
python3 -m venv ${SW_DIR}/venvs/warpx-frontier
source ${SW_DIR}/venvs/warpx-frontier/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade build
python3 -m pip install --upgrade packaging
python3 -m pip install --upgrade wheel
python3 -m pip install --upgrade setuptools[core]
python3 -m pip install --upgrade "cython>=3.0"
# cupy for ROCm
#   https://docs.cupy.dev/en/stable/install.html#building-cupy-for-rocm-from-source
#   https://github.com/cupy/cupy/issues/7830
CC=cc CXX=CC \
CUPY_INSTALL_USE_HIP=1  \
ROCM_HOME=${ROCM_PATH}  \
HCC_AMDGPU_TARGET=${AMREX_AMD_ARCH}  \
  python3 -m pip install -v git+https://github.com/cupy/cupy.git@e669b994f976565bf2da4b1f82de51e10b58fbe1
python3 -m pip install --upgrade numpy
python3 -m pip install --upgrade h5py
python3 -m pip install --upgrade pandas
python3 -m pip install --upgrade scipy
MPICC="cc -shared" python3 -m pip install --upgrade mpi4py --no-cache-dir --no-build-isolation --no-binary mpi4py
python3 -m pip install --upgrade openpmd-api
python3 -m pip install --upgrade matplotlib
python3 -m pip install --upgrade yt
# install or update WarpX dependencies such as picmistandard
python3 -m pip install --upgrade -r $HOME/src/warpx/requirements.txt
# optional: for optimas (based on libEnsemble & ax->botorch->gpytorch->pytorch)
#python3 -m pip install --upgrade torch --index-url https://download.pytorch.org/whl/rocm5.4.2
#python3 -m pip install -r $HOME/src/warpx/Tools/optimas/requirements.txt
```

<a id="building-frontier-compilation"></a>

## Compilation

Use the following [cmake commands](../cmake.md#install-build-cmake) to compile the application executable:

```bash
cd $HOME/src/warpx
rm -rf build_frontier

cmake -S . -B build_frontier -DWarpX_COMPUTE=HIP -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_frontier -j 16
```

The WarpX application executables are now in `$HOME/src/warpx/build_frontier/bin/`.
Additionally, the following commands will install WarpX as a Python module:

```bash
rm -rf build_frontier_py

cmake -S . -B build_frontier_py -DWarpX_COMPUTE=HIP -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_APP=OFF -DWarpX_PYTHON=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_frontier_py -j 16 --target pip_install
```

Now, you can [submit Frontier compute jobs](#running-cpp-frontier) for WarpX [Python (PICMI) scripts](../../usage/python.md#usage-picmi) ([example scripts](../../usage/examples.md#usage-examples)).
Or, you can use the WarpX executables to submit Frontier jobs ([example inputs](../../usage/examples.md#usage-examples)).
For executables, you can reference their location in your [job script](#running-cpp-frontier) or copy them to a location in `$PROJWORK/$proj/`.

<a id="building-frontier-update"></a>

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

- [update the frontier_warpx.profile file](#building-frontier-preparation),
- log out and into the system, activate the now updated environment profile as usual,
- [execute the dependency install scripts](#building-frontier-preparation).

As a last step, clean the build directory `rm -rf $HOME/src/warpx/build_frontier` and rebuild WarpX.

<a id="running-cpp-frontier"></a>

## Running

<a id="running-cpp-frontier-mi250x-gpus"></a>

### MI250X GPUs (2x64 GB)

After requesting an interactive node with the `getNode` alias above, run a simulation like this, here using 8 MPI ranks and a single node:

```bash
runNode ./warpx inputs
```

Or in non-interactive runs:

```bash
#!/usr/bin/env bash

#SBATCH -A <project id>
#SBATCH -J warpx
#SBATCH -o %x-%j.out
#SBATCH -t 00:10:00
#SBATCH -p batch
#SBATCH --ntasks-per-node=8
# Due to Frontier's Low-Noise Mode Layout only 7 instead of 8 cores are available per process
# https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#low-noise-mode-layout
#SBATCH --cpus-per-task=7
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=closest
#SBATCH -N 20

# load cray libs and ROCm libs
#export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}

# From the documentation:
# Each Frontier compute node consists of [1x] 64-core AMD EPYC 7A53
# "Optimized 3rd Gen EPYC" CPU (with 2 hardware threads per physical core) with
# access to 512 GB of DDR4 memory.
# Each node also contains [4x] AMD MI250X, each with 2 Graphics Compute Dies
# (GCDs) for a total of 8 GCDs per node. The programmer can think of the 8 GCDs
# as 8 separate GPUs, each having 64 GB of high-bandwidth memory (HBM2E).

# note (5-16-22 and 7-12-22)
# this environment setting is currently needed on Frontier to work-around a
# known issue with Libfabric (both in the May and June PE)
#export FI_MR_CACHE_MAX_COUNT=0  # libfabric disable caching
# or, less invasive:
export FI_MR_CACHE_MONITOR=memhooks  # alternative cache monitor

# Seen since August 2023
# OLCFDEV-1597: OFI Poll Failed UNDELIVERABLE Errors
# https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#olcfdev-1597-ofi-poll-failed-undeliverable-errors
export MPICH_SMP_SINGLE_COPY_MODE=NONE
export FI_CXI_RX_MATCH_MODE=software

# note (9-2-22, OLCFDEV-1079)
# this environment setting is needed to avoid that rocFFT writes a cache in
# the home directory, which does not scale.
export ROCFFT_RTC_CACHE_PATH=/dev/null

export OMP_NUM_THREADS=1
export WARPX_NMPI_PER_NODE=8
export TOTAL_NMPI=$(( ${SLURM_JOB_NUM_NODES} * ${WARPX_NMPI_PER_NODE} ))
srun -N${SLURM_JOB_NUM_NODES} -n${TOTAL_NMPI} --ntasks-per-node=${WARPX_NMPI_PER_NODE} \
    ./warpx inputs > output.txt
```

<a id="post-processing-frontier"></a>

## Post-Processing

For post-processing, most users use Python via OLCFs’s [Jupyter service](https://jupyter.olcf.ornl.gov) ([Docs](https://docs.olcf.ornl.gov/services_and_applications/jupyter/index.html)).

We usually just install our software on-the-fly on Frontier.
When starting up a post-processing session, run this in your first cells:

#### NOTE
The following software packages are installed only into a temporary directory.

```bash
# work-around for OLCFHELP-4242
!jupyter serverextension enable --py --sys-prefix dask_labextension

# next Jupyter cell: the software you want
!mamba install --quiet -c conda-forge -y openpmd-api openpmd-viewer ipympl ipywidgets fast-histogram yt

# restart notebook
```

<a id="known-frontier-issues"></a>

## Known System Issues

#### WARNING
May 16th, 2022 (OLCFHELP-6888):
There is a caching bug in Libfabric that causes WarpX simulations to occasionally hang on Frontier on more than 1 node.

As a work-around, please export the following environment variable in your job scripts until the issue is fixed:

```bash
#export FI_MR_CACHE_MAX_COUNT=0  # libfabric disable caching
# or, less invasive:
export FI_MR_CACHE_MONITOR=memhooks  # alternative cache monitor
```

#### WARNING
Sep 2nd, 2022 (OLCFDEV-1079):
rocFFT in ROCm 5.1-5.3 tries to [write to a cache](https://rocfft.readthedocs.io/en/latest/#runtime-compilation) in the home area by default.
This does not scale, disable it via:

```bash
export ROCFFT_RTC_CACHE_PATH=/dev/null
```

#### WARNING
January, 2023 (OLCFDEV-1284, AMD Ticket: ORNLA-130):
We discovered a regression in AMD ROCm, leading to 2x slower current deposition (and other slowdowns) in ROCm 5.3 and 5.4.

June, 2023:
Although a fix was planned for ROCm 5.5, we still see the same issue in this release and continue to exchange with AMD and HPE on the issue.

Stay with the ROCm 5.2 module to avoid a 2x slowdown.

#### WARNING
August, 2023 (OLCFDEV-1597, OLCFHELP-12850, OLCFHELP-14253):
With runs above 500 nodes, we observed issues in `MPI_Waitall` calls of the kind `OFI Poll Failed UNDELIVERABLE`.
According to the system known issues entry [OLCFDEV-1597](https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#olcfdev-1597-ofi-poll-failed-undeliverable-errors), we work around this by setting this environment variable in job scripts:

```bash
export MPICH_SMP_SINGLE_COPY_MODE=NONE
export FI_CXI_RX_MATCH_MODE=software
```

#### WARNING
Checkpoints and AMReX plotfile I/O at scale is very slow with the default Lustre filesystem configuration.
Using openPMD with ADIOS2 is the best performing output for regular diagnostics.

For checkpoint-restart, you have no choice and need to use AMReX plotfiles.
Please test checkpointing and I/O with short `#SBATCH -q debug` runs before running the full simulation.
Set [the following options for performance of plotfiles/checkpoints](https://github.com/AMReX-Codes/amrex/pull/4426):

```ini
warpx.field_io_nfiles = <1-per-node>
warpx.particle_io_nfiles = <1-per-node>

# These parameters are for a workaround needed for simulations on FRONTIER
# and enable checkpointing
vismf.noflushafterwrite = true
vismf.barrierafterlevel = true
```

Execute `lfs getstripe -d <dir>` to show the default progressive file layout.
For further tuning, consider using `lfs setstripe` to change the [striping](https://wiki.lustre.org/Configuring_Lustre_File_Striping) for new files **before** you submit the run.

```bash
mkdir /lustre/orion/proj-shared/<your-project>/<path/to/new/sim/dir>
cd <new/sim/dir/above>
# create your diagnostics directory first
mkdir diags
# change striping for new files before you submit the simulation
#   this is an example, striping 10 MB blocks onto 32 nodes
lfs setstripe -S 16M -c 1 $SLURM_SUBMIT_DIR  # or diags only or checkpoints only
```
