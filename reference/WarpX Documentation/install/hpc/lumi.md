<a id="building-lumi"></a>

# LUMI (CSC)

The [LUMI cluster](https://www.lumi-supercomputer.eu) is located at CSC (Finland).
Each node contains 4 AMD MI250X GPUs, each with 2 Graphics Compute Dies (GCDs) for a total of 8 GCDs per node.
You can think of the 8 GCDs as 8 separate GPUs, each having 64 GB of high-bandwidth memory (HBM2E).

## Introduction

If you are new to this system, **please see the following resources**:

* [Lumi user guide](https://docs.lumi-supercomputer.eu)
  * [Project Maintainance](https://my.lumi-supercomputer.eu) and [SSH Key management](https://mms.myaccessid.org)
  * [Quotas and projects](https://docs.lumi-supercomputer.eu/runjobs/lumi_env/dailymanagement/)
* Batch system: [Slurm](https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/slurm-quickstart/)
* [Data analytics and visualization](https://docs.lumi-supercomputer.eu/hardware/lumid/)
* [Production directories](https://docs.lumi-supercomputer.eu/storage/):
  * `$HOME`: single user, intended to store user configuration files and personal data (20GB default quota)
  * `/project/$proj`: shared with all members of a project, purged at the end of a project (50 GB default quota)
  * `/scratch/$proj`: temporary storage, main storage to be used for disk I/O needs when running simulations on LUMI, purged every 90 days (50TB default quota)

<a id="building-lumi-preparation"></a>

## Preparation

Use the following commands to download the WarpX source code:

```bash
git clone https://github.com/BLAST-WarpX/warpx.git $HOME/src/warpx
```

We use system software modules, add environment hints and further dependencies via the file `$HOME/lumi_warpx.profile`.
Create it now:

```bash
cp $HOME/src/warpx/Tools/machines/lumi-csc/lumi_warpx.profile.example $HOME/lumi_warpx.profile
```

### Script Details

```bash
# please set your project account
#export proj="project_..."

# required dependencies
module load LUMI/25.03  partition/G
module load rocm/6.3.4
module load buildtools/25.03

# optional: just an additional text editor
module load nano

# optional: for PSATD in RZ geometry support
SW_DIR="${HOME}/sw/lumi/gpu"
export CMAKE_PREFIX_PATH=${SW_DIR}/blaspp-2024.05.31:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/lapackpp-2024.05.31:$CMAKE_PREFIX_PATH
export LD_LIBRARY_PATH=${SW_DIR}/blaspp-2024.05.31/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/lapackpp-2024.05.31/lib64:$LD_LIBRARY_PATH

# optional: for QED lookup table generation support
module load Boost/1.88.0-cpeCray-25.03

# optional: for openPMD support
export CMAKE_PREFIX_PATH=${SW_DIR}/c-blosc-1.21.1:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/hdf5-1.14.1.2:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/adios2-2.10.2:$CMAKE_PREFIX_PATH
export PATH=${SW_DIR}/hdf5-1.14.1.2/bin:${PATH}
export PATH=${SW_DIR}/adios2-2.10.2/bin:${PATH}

# optional: for Python bindings or libEnsemble
module load cray-python/3.11.7

if [ -d "${SW_DIR}/venvs/warpx-lumi" ]
then
  source ${SW_DIR}/venvs/warpx-lumi/bin/activate
fi

# an alias to request an interactive batch node for one hour
#   for paralle execution, start on the batch node: srun <command>
alias getNode="salloc -A $proj -J warpx -t 01:00:00 -p dev-g -N 1 --ntasks-per-node=8 --gpus-per-task=1 --gpu-bind=closest"
# an alias to run a command on a batch node for up to 30min
#   usage: runNode <command>
alias runNode="srun -A $proj -J warpx -t 00:30:00 -p dev-g -N 1 --ntasks-per-node=8 --gpus-per-task=1 --gpu-bind=closest"

# GPU-aware MPI
export MPICH_GPU_SUPPORT_ENABLED=1

# optimize ROCm/HIP compilation for MI250X
export AMREX_AMD_ARCH=gfx90a

# compiler environment hints
# Warning: using the compiler wrappers cc and CC
#          instead of amdclang and amdclang++
#          currently results in a significant
#          loss of performances
export CC=$(which amdclang)
export CXX=$(which amdclang++)
export FC=$(which amdflang)
export CFLAGS="-I${ROCM_PATH}/include"
export CXXFLAGS="-I${ROCM_PATH}/include -Wno-pass-failed"
export LDFLAGS="-L${ROCM_PATH}/lib -lamdhip64 ${PE_MPICH_GTL_DIR_amd_gfx90a} -lmpi_gtl_hsa"
```

Edit the 2nd line of this script, which sets the `export proj="project_..."` variable using a text editor
such as `nano`, `emacs`, or `vim` (all available by default on LUMI login nodes).
You can find out your project name by running `lumi-ldap-userinfo` on LUMI.
For example, if you are member of the project `project_465000559`, then run `nano $HOME/lumi_impactx.profile` and edit line 2 to read:

```bash
export proj="project_465000559"
```

Exit the `nano` editor with `Ctrl` + `O` (save) and then `Ctrl` + `X` (exit).

#### IMPORTANT
Now, and as the first step on future logins to LUMI, activate these environment settings:

```bash
source $HOME/lumi_warpx.profile
```

Finally, since LUMI does not yet provide software modules for some of our dependencies, install them once:

```bash
bash $HOME/src/warpx/Tools/machines/lumi-csc/install_dependencies.sh
source $HOME/sw/lumi/gpu/venvs/warpx-lumi/bin/activate
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
#   Was lumi_warpx.profile sourced and configured correctly?
if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your lumi_warpx.profile file! Please edit its line 2 to continue!"; exit 1; fi


# Remove old dependencies #####################################################
#
SRC_DIR="${HOME}/src"
SW_DIR="${HOME}/sw/lumi/gpu"
rm -rf ${SW_DIR}
mkdir -p ${SW_DIR}
mkdir -p ${SRC_DIR}

# remove common user mistakes in python, located in .local instead of a venv
python3 -m pip uninstall -qq -y pywarpx
python3 -m pip uninstall -qq -y warpx
python3 -m pip uninstall -qqq -y mpi4py 2>/dev/null || true


# General extra dependencies ##################################################
#

# tmpfs build directory: avoids issues often seen with $HOME and is faster
build_dir=$(mktemp -d)

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
rm -rf ${build_dir}/blaspp-lumi-gpu-build
CXX=$(which CC)                              \
cmake -S ${SRC_DIR}/blaspp                   \
      -B ${build_dir}/blaspp-lumi-gpu-build  \
      -Duse_openmp=OFF                       \
      -Dgpu_backend=hip                      \
      -DCMAKE_CXX_STANDARD=20                \
      -DCMAKE_INSTALL_PREFIX=${SW_DIR}/blaspp-2024.05.31
cmake --build ${build_dir}/blaspp-lumi-gpu-build --target install --parallel 16
rm -rf ${build_dir}/blaspp-lumi-gpu-build

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
rm -rf ${build_dir}/lapackpp-lumi-gpu-build
CXX=$(which CC) CXXFLAGS="-DLAPACK_FORTRAN_ADD_" \
cmake -S ${SRC_DIR}/lapackpp                     \
      -B ${build_dir}/lapackpp-lumi-gpu-build    \
      -DCMAKE_CXX_STANDARD=20                    \
      -Dbuild_tests=OFF                          \
      -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON     \
      -DCMAKE_INSTALL_PREFIX=${SW_DIR}/lapackpp-2024.05.31
cmake --build ${build_dir}/lapackpp-lumi-gpu-build --target install --parallel 16
rm -rf ${build_dir}/lapackpp-lumi-gpu-build

# c-blosc (I/O compression, for openPMD)
if [ -d ${SRC_DIR}/c-blosc ]
then
  cd ${SRC_DIR}/c-blosc
  git fetch --prune
  git checkout v1.21.1
  cd -
else
  git clone -b v1.21.1 https://github.com/Blosc/c-blosc.git ${SRC_DIR}/c-blosc
fi
rm -rf ${build_dir}/c-blosc-lu-build
cmake -S ${SRC_DIR}/c-blosc             \
      -B ${build_dir}/c-blosc-lu-build  \
      -DBUILD_TESTS=OFF                 \
      -DBUILD_BENCHMARKS=OFF            \
      -DDEACTIVATE_AVX2=OFF             \
      -DCMAKE_INSTALL_PREFIX=${HOME}/sw/lumi/gpu/c-blosc-1.21.1
cmake --build ${build_dir}/c-blosc-lu-build --target install --parallel 16
rm -rf ${build_dir}/c-blosc-lu-build

# HDF5 (for openPMD)
if [ -d ${SRC_DIR}/hdf5 ]
then
  cd ${SRC_DIR}/hdf5
  git fetch --prune
  git checkout hdf5-1_14_1-2
  cd -
else
  git clone -b hdf5-1_14_1-2 https://github.com/HDFGroup/hdf5.git ${SRC_DIR}/hdf5
fi
rm -rf ${build_dir}/hdf5-lu-build
cmake -S ${SRC_DIR}/hdf5          \
      -B ${build_dir}/hdf5-lu-build  \
      -DBUILD_TESTING=OFF         \
      -DHDF5_ENABLE_PARALLEL=ON   \
      -DCMAKE_INSTALL_PREFIX=${SW_DIR}/hdf5-1.14.1.2
cmake --build ${build_dir}/hdf5-lu-build --target install --parallel 10
rm -rf ${build_dir}/hdf5-lu-build

# ADIOS2 (for openPMD)
if [ -d ${SRC_DIR}/adios2 ]
then
  cd ${SRC_DIR}/adios2
  git fetch --prune
  git checkout v2.10.2
  cd -
else
  git clone -b v2.10.2 https://github.com/ornladios/ADIOS2.git ${SRC_DIR}/adios2
fi
rm -rf ${build_dir}/adios2-lu-build
cmake -S ${SRC_DIR}/adios2             \
      -B ${build_dir}/adios2-lu-build  \
      -DADIOS2_USE_Blosc=ON            \
      -DADIOS2_USE_Fortran=OFF         \
      -DADIOS2_USE_HDF5=OFF            \
      -DADIOS2_USE_Python=OFF          \
      -DADIOS2_USE_ZeroMQ=OFF          \
      -DCMAKE_INSTALL_PREFIX=${HOME}/sw/lumi/gpu/adios2-2.10.2
cmake --build ${build_dir}/adios2-lu-build --target install -j 16
rm -rf ${build_dir}/adios2-lu-build


# Python ######################################################################
#
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade virtualenv
python3 -m pip cache purge
rm -rf ${SW_DIR}/venvs/warpx-lumi
python3 -m venv ${SW_DIR}/venvs/warpx-lumi
source ${SW_DIR}/venvs/warpx-lumi/bin/activate
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
# install or update WarpX dependencies such as picmistandard
python3 -m pip install --upgrade -r ${SRC_DIR}/warpx/requirements.txt
# optional: for optimas (based on libEnsemble & ax->botorch->gpytorch->pytorch)
python3 -m pip install -r ${SRC_DIR}/warpx/Tools/optimas/requirements.txt
#python3 -m pip install --upgrade torch --index-url https://download.pytorch.org/whl/rocm5.4.2
#python3 -m pip install -r ${SRC_DIR}/warpx/Tools/optimas/requirements.txt
```

<a id="building-lumi-compilation"></a>

## Compilation

Use the following [cmake commands](../cmake.md#install-build-cmake) to compile the application executable:

```bash
cd $HOME/src/warpx
rm -rf build_lumi

cmake -S . -B build_lumi -DWarpX_COMPUTE=HIP -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_QED_TABLES_GEN_OMP=OFF -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_lumi -j 16
```

The WarpX application executables are now in `$HOME/src/warpx/build_lumi/bin/`.
Additionally, the following commands will install WarpX as a Python module:

```bash
rm -rf build_lumi_py

cmake -S . -B build_lumi_py -DWarpX_COMPUTE=HIP -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_QED_TABLES_GEN_OMP=OFF -DWarpX_APP=OFF -DWarpX_PYTHON=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_lumi_py -j 16 --target pip_install
```

<a id="building-lumi-update"></a>

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

- [update the lumi_warpx.profile file](#building-lumi-preparation),
- log out and into the system, activate the now updated environment profile as usual,
- [execute the dependency install scripts](#building-lumi-preparation).

As a last step, clean the build directory `rm -rf $HOME/src/warpx/build_lumi` and rebuild WarpX.

<a id="running-cpp-lumi"></a>

## Running

<a id="running-cpp-lumi-mi250x-gpus"></a>

### MI250X GPUs (2x64 GB)

The GPU partition on the supercomputer LUMI at CSC has up to [2978 nodes](https://docs.lumi-supercomputer.eu/hardware/lumig/), each with 8 Graphics Compute Dies (GCDs).
WarpX runs one MPI rank per Graphics Compute Die.

For interactive runs, simply use the aliases `getNode` or `runNode ...`.

The batch script below can be used to run a WarpX simulation on multiple nodes (change `-N` accordingly).
Replace descriptions between chevrons `<>` by relevant values, for instance `<project id>` or the concete inputs file.
Copy the executable or point to it via `EXE` and adjust the path for the `INPUTS` variable accordingly.

```bash
#!/bin/bash -l

#SBATCH -A <project id>
#SBATCH --job-name=warpx
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --partition=standard-g
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --gpus-per-node=8
#SBATCH --time=00:10:00

date

# note (12-12-22)
# this environment setting is currently needed on LUMI to work-around a
# known issue with Libfabric
#export FI_MR_CACHE_MAX_COUNT=0  # libfabric disable caching
# or, less invasive:
export FI_MR_CACHE_MONITOR=memhooks  # alternative cache monitor

# Seen since August 2023 seen on OLCF (not yet seen on LUMI?)
# OLCFDEV-1597: OFI Poll Failed UNDELIVERABLE Errors
# https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#olcfdev-1597-ofi-poll-failed-undeliverable-errors
#export MPICH_SMP_SINGLE_COPY_MODE=NONE
#export FI_CXI_RX_MATCH_MODE=software

# note (9-2-22, OLCFDEV-1079)
# this environment setting is needed to avoid that rocFFT writes a cache in
# the home directory, which does not scale.
export ROCFFT_RTC_CACHE_PATH=/dev/null

# Seen since August 2023
# OLCFDEV-1597: OFI Poll Failed UNDELIVERABLE Errors
# https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#olcfdev-1597-ofi-poll-failed-undeliverable-errors
export MPICH_SMP_SINGLE_COPY_MODE=NONE
export FI_CXI_RX_MATCH_MODE=software

# LUMI documentation suggests using the following wrapper script
# to set the ROCR_VISIBLE_DEVICES to the value of SLURM_LOCALID
# see https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/lumig-job/
cat << EOF > select_gpu
#!/bin/bash

export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID
exec \$*
EOF

chmod +x ./select_gpu

sleep 1

# LUMI documentation suggests using the following CPU bind
# in order to have 6 threads per GPU (blosc compression in adios2 uses threads)
# see https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/lumig-job/
#
# WARNING: the following CPU_BIND options don't work on the dev-g partition.
#          If you want to run your simulation on dev-g, please comment them
#          out and replace them with CPU_BIND="map_cpu:49,57,17,25,1,9,33,41"
#
CPU_BIND="mask_cpu:7e000000000000,7e00000000000000"
CPU_BIND="${CPU_BIND},7e0000,7e000000"
CPU_BIND="${CPU_BIND},7e,7e00"
CPU_BIND="${CPU_BIND},7e00000000,7e0000000000"

export OMP_NUM_THREADS=6

export MPICH_GPU_SUPPORT_ENABLED=1

srun --cpu-bind=${CPU_BIND} ./select_gpu ./warpx inputs | tee outputs.txt
rm -rf ./select_gpu
```

To run a simulation, copy the lines above to a file `lumi.sbatch` and run

```bash
sbatch lumi.sbatch
```

to submit the job.

<a id="post-processing-lumi"></a>

## Post-Processing

#### NOTE
TODO: Document any Jupyter or data services.

<a id="known-lumi-issues"></a>

## Known System Issues

#### WARNING
December 12th, 2022:
There is a caching bug in libFabric that causes WarpX simulations to occasionally hang on LUMI on more than 1 node.

As a work-around, please export the following environment variable in your job scripts until the issue is fixed:

```bash
#export FI_MR_CACHE_MAX_COUNT=0  # libfabric disable caching
# or, less invasive:
export FI_MR_CACHE_MONITOR=memhooks  # alternative cache monitor
```

#### WARNING
January, 2023:
We discovered a regression in AMD ROCm, leading to 2x slower current deposition (and other slowdowns) in ROCm 5.3 and 5.4.

June, 2023:
Although a fix was planned for ROCm 5.5, we still see the same issue in this release and continue to exchange with AMD and HPE on the issue.

Stay with the ROCm 5.2 module to avoid a 2x slowdown.

#### WARNING
May 2023:
rocFFT in ROCm 5.1-5.3 tries to [write to a cache](https://rocfft.readthedocs.io/en/latest/#runtime-compilation) in the home area by default.
This does not scale, disable it via:

```bash
export ROCFFT_RTC_CACHE_PATH=/dev/null
```
