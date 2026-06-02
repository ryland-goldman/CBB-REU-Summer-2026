<a id="building-dane"></a>

# Dane (LLNL)

The [Dane Intel CPU cluster](https://hpc.llnl.gov/hardware/compute-platforms/dane) is located at LLNL.

## Introduction

If you are new to this system, **please see the following resources**:

* [LLNL user account](https://lc.llnl.gov) (login required)
* [Jupyter service](https://lc.llnl.gov/jupyter) ([documentation](https://lc.llnl.gov/confluence/display/LC/JupyterHub+and+Jupyter+Notebook), login required)
* [Production directories](https://hpc.llnl.gov/hardware/file-systems):
  * `/p/lustre1/$(whoami)` and `/p/lustre2/$(whoami)`: personal directory on the parallel filesystem
  * Note that the `$HOME` directory and the `/usr/workspace/$(whoami)` space are NFS mounted and *not* suitable for production quality data generation.

<a id="building-dane-preparation"></a>

## Preparation

Use the following commands to download the WarpX source code.
Note that these commands and the shell scripts all assume the bash shell.
This downloads WarpX into the workspace directory, which is recommended.
WarpX can be downloaded elsewhere if that doesn’t work with your directory structure, but note that the commands shown below refer to WarpX in the workspace directory.

```bash
git clone https://github.com/BLAST-WarpX/warpx.git /usr/workspace/${USER}/dane/src/warpx
```

The system software modules, environment hints, and further dependencies are setup via the file `$HOME/dane_warpx.profile` which is copied from the WarpX source.
Set it up now:

```bash
cp /usr/workspace/${USER}/dane/src/warpx/Tools/machines/dane-llnl/dane_warpx.profile.example $HOME/dane_warpx.profile
```

### Script Details

```bash
# please set your project account
#export proj="<yourProjectNameHere>"  # edit this and comment in

# required dependencies
module load cmake/3.26.3
module load clang/14.0.6-magic
module load mvapich2/2.3.7

# optional: for PSATD support
module load fftw/3.3.10

# optional: for QED lookup table generation support
module load boost/1.80.0

# optional: for openPMD support
module load hdf5-parallel/1.14.0

if [ -z "${WARPX_SW_DIR+x}" ]; then
    WARPX_SW_DIR="/usr/workspace/${USER}/dane"
fi

export CMAKE_PREFIX_PATH=${WARPX_SW_DIR}/install/c-blosc-1.21.6:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${WARPX_SW_DIR}/install/adios2-2.10.2:$CMAKE_PREFIX_PATH
export PATH=${WARPX_SW_DIR}/install/adios2-2.10.2/bin:${PATH}

# optional: for PSATD in RZ geometry support
export CMAKE_PREFIX_PATH=${WARPX_SW_DIR}/install/blaspp-2024.10.26:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${WARPX_SW_DIR}/install/lapackpp-2024.10.26:$CMAKE_PREFIX_PATH
export LD_LIBRARY_PATH=${WARPX_SW_DIR}/install/blaspp-2024.10.26/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${WARPX_SW_DIR}/install/lapackpp-2024.10.26/lib64:$LD_LIBRARY_PATH

# optional: for Python bindings
#module load python/3.12.2 # Note this version uses a too new GLIBCXX and breaks the adios build
module load python/3.11.5

if [ -d "${WARPX_SW_DIR}/venvs/warpx-dane" ]
then
    source ${WARPX_SW_DIR}/venvs/warpx-dane/bin/activate
fi

# optional: an alias to request an interactive node for two hours
alias getNode="srun --time=0:30:00 --nodes=1 --ntasks-per-node=2 --cpus-per-task=56 -p pdebug --pty bash"
# an alias to run a command on a batch node for up to 30min
#   usage: runNode <command>
alias runNode="srun --time=0:30:00 --nodes=1 --ntasks-per-node=2 --cpus-per-task=56 -p pdebug"

# fix system defaults: do not escape $ with a \ on tab completion
shopt -s direxpand

# optimize CPU microarchitecture for Intel Sapphire Rapids
# note: the cc/CC/ftn wrappers below add those
export CXXFLAGS="-march=sapphirerapids"
export CFLAGS="-march=sapphirerapids"

# compiler environment hints
export CC=$(which clang)
export CXX=$(which clang++)
export FC=$(which gfortran)
```

Edit the 2nd line of this script, which sets the `export proj=""` variable.
For example, if you are member of the project `tps`, then run `vi $HOME/dane_warpx.profile`.
Enter the edit mode by typing `i` and edit line 2 to read:

```bash
export proj="tps"
```

Exit the `vi` editor with `Esc` and then type `:wq` (write & quit).

#### IMPORTANT
Now, and as the first step on future logins to Dane, activate these environment settings by executing the file:

```bash
source $HOME/dane_warpx.profile
```

Finally, since Dane does not yet provide software modules for some of our dependencies, WarpX provides a script to install them.
This is done executed now.
They are by default installed in the workspace directory (which is recommended), but can be installed elsewhere by setting the environment variable `WARPX_SW_DIR`.
The second command activates the Python virtual environment.
This would normally be done by the `dane_warpx.profile` script, but the environment is created by the install script and so wasn’t created yet when the profile was run above.
So the activation needs to be done this way only this one time.

```bash
bash /usr/workspace/${USER}/dane/src/warpx/Tools/machines/dane-llnl/install_dependencies.sh
source /usr/workspace/${USER}/dane/venvs/warpx-dane/bin/activate
```

### Script Details

```bash
#!/bin/bash
#
# Copyright 2024 The WarpX Community
#
# This file is part of WarpX.
#
# Author: Axel Huebl, David Grote
# License: BSD-3-Clause-LBNL

# Exit on first error encountered #############################################
#
set -eu -o pipefail


# Check: ######################################################################
#
#   Was dane_warpx.profile sourced and configured correctly?
if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your dane_warpx.profile file! Please edit its line 2 to continue!"; exit 1; fi

# Make sure that a virtual environment is not already activated
if declare -F deactivate &>/dev/null; then
  deactivate
fi

# remove common user mistakes in python, located in .local instead of a venv
python3 -m pip uninstall -qq -y pywarpx
python3 -m pip uninstall -qq -y warpx
python3 -m pip uninstall -qqq -y mpi4py 2>/dev/null || true

# Setup the directories where the packages will be installed
if [ -z "${WARPX_SW_DIR+x}" ]; then
    WARPX_SW_DIR="/usr/workspace/${USER}/dane"
fi
rm -rf ${WARPX_SW_DIR}/install
mkdir -p ${WARPX_SW_DIR}/install

# General extra dependencies ##################################################
#

# tmpfs build directory: avoids issues often seen with ${HOME} and is faster
build_dir=$(mktemp -d)

# c-blosc (I/O compression)
if [ -d ${WARPX_SW_DIR}/src/c-blosc ]
then
  cd ${WARPX_SW_DIR}/src/c-blosc
  git fetch --prune
  git checkout v1.21.6
  cd -
else
  git clone -b v1.21.6 https://github.com/Blosc/c-blosc.git ${WARPX_SW_DIR}/src/c-blosc
fi
cmake -S ${WARPX_SW_DIR}/src/c-blosc -B ${build_dir}/c-blosc-dane-build -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF -DDEACTIVATE_AVX2=OFF -DCMAKE_INSTALL_PREFIX=${WARPX_SW_DIR}/install/c-blosc-1.21.6
cmake --build ${build_dir}/c-blosc-dane-build --target install --parallel 6

# ADIOS2
if [ -d ${WARPX_SW_DIR}/src/adios2 ]
then
  cd ${WARPX_SW_DIR}/src/adios2
  git fetch --prune
  git checkout v2.10.2
  cd -
else
  git clone -b v2.10.2 https://github.com/ornladios/ADIOS2.git ${WARPX_SW_DIR}/src/adios2
fi
cmake -S ${WARPX_SW_DIR}/src/adios2 -B ${build_dir}/adios2-dane-build -DBUILD_TESTING=OFF -DADIOS2_BUILD_EXAMPLES=OFF -DADIOS2_USE_Blosc=ON -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DADIOS2_USE_SST=OFF -DADIOS2_USE_ZeroMQ=OFF -DCMAKE_INSTALL_PREFIX=${WARPX_SW_DIR}/install/adios2-2.10.2
cmake --build ${build_dir}/adios2-dane-build --target install -j 6

# BLAS++ (for PSATD+RZ)
if [ -d ${WARPX_SW_DIR}/src/blaspp ]
then
  cd ${WARPX_SW_DIR}/src/blaspp
  git fetch --prune
  git checkout v2024.10.26
  cd -
else
  git clone -b v2024.10.26 https://github.com/icl-utk-edu/blaspp.git ${WARPX_SW_DIR}/src/blaspp
fi
cmake -S ${WARPX_SW_DIR}/src/blaspp -B ${build_dir}/blaspp-dane-build -Duse_openmp=ON -Duse_cmake_find_blas=ON -DCMAKE_CXX_STANDARD=20 -DCMAKE_INSTALL_PREFIX=${WARPX_SW_DIR}/install/blaspp-2024.10.26
cmake --build ${build_dir}/blaspp-dane-build --target install --parallel 6

# LAPACK++ (for PSATD+RZ)
if [ -d ${WARPX_SW_DIR}/src/lapackpp ]
then
  cd ${WARPX_SW_DIR}/src/lapackpp
  git fetch --prune
  git checkout v2024.10.26
  cd -
else
  git clone -b v2024.10.26 https://github.com/icl-utk-edu/lapackpp.git ${WARPX_SW_DIR}/src/lapackpp
fi
CXXFLAGS="-DLAPACK_FORTRAN_ADD_" cmake -S ${WARPX_SW_DIR}/src/lapackpp -B ${build_dir}/lapackpp-dane-build -Duse_cmake_find_lapack=ON -DCMAKE_CXX_STANDARD=20 -Dbuild_tests=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX=${WARPX_SW_DIR}/install/lapackpp-2024.10.26
cmake --build ${build_dir}/lapackpp-dane-build --target install --parallel 6


# Python ######################################################################
#
# Create a virtual environment and install the Python packages there.
rm -rf ${WARPX_SW_DIR}/venvs/warpx-dane
python3 -m venv ${WARPX_SW_DIR}/venvs/warpx-dane
source ${WARPX_SW_DIR}/venvs/warpx-dane/bin/activate
python3 -m pip install --upgrade pip
#python3 -m pip cache purge
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
SCRIPT_PATH="$(realpath ${BASH_SOURCE[0]})"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
python3 -m pip install --upgrade -r ${SCRIPT_DIR}/../../../requirements.txt

# ML dependencies
python3 -m pip install --upgrade torch


# remove build temporary directory ############################################
#
rm -rf ${build_dir}
```

<a id="building-dane-compilation"></a>

## Compilation

Use the following [cmake commands](../cmake.md#install-build-cmake) to compile the application executable.
The options should be modified to suit your needs, for example only building for the dimensions needed.

```bash
cd /usr/workspace/${USER}/dane/src/warpx
rm -rf build_dane

cmake -S . -B build_dane -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_dane -j 6
```

The WarpX application executables are now in `/usr/workspace/${USER}/dane/src/warpx/build_dane/bin/`.
Additionally, the following commands will install WarpX as a Python module:

```bash
rm -rf build_dane_py

cmake -S . -B build_dane_py -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_APP=OFF -DWarpX_PYTHON=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_dane_py -j 6 --target pip_install
```

Now, you can [submit Dane compute jobs](#running-cpp-dane) for WarpX [Python (PICMI) scripts](../../usage/python.md#usage-picmi) ([example scripts](../../usage/examples.md#usage-examples)).
Or, you can use the WarpX executables to submit Dane jobs ([example inputs](../../usage/examples.md#usage-examples)).
For executables, you can reference their location in your [job script](#running-cpp-dane) or copy them to a location in `$PROJWORK/$proj/`.

<a id="building-dane-update"></a>

## Update WarpX & Dependencies

If you already installed WarpX in the past and want to update it, start by getting the latest source code:

```bash
cd /usr/workspace/${USER}/dane/src/warpx

# read the output of this command - does it look ok?
git status

# get the latest WarpX source code
git pull

# read the output of these commands - do they look ok?
git status
git log     # press q to exit
```

And, if needed,

- [update the dane_warpx.profile file](#building-dane-preparation),
- log out and into the system, activate the now updated environment profile as usual,
- [execute the dependency install scripts](#building-dane-preparation).

As a last step, clean the build directory `rm -rf /usr/workspace/${USER}/dane/src/warpx/build_dane` and rebuild WarpX.

<a id="running-cpp-dane"></a>

## Running

<a id="running-cpp-dane-cpus"></a>

### Intel Sapphire Rapids CPUs

The batch script below can be used to run a WarpX simulation on 2 nodes on the supercomputer Dane at LLNL.
Replace descriptions between chevrons `<>` by relevant values, for instance `<input file>` could be `plasma_mirror_inputs`.

```bash
#!/bin/bash -l

# Just increase this number of you need more nodes.
#SBATCH -N 2
#SBATCH -t 24:00:00
#SBATCH -A <allocation ID>

#SBATCH -J WarpX
#SBATCH -q pbatch
#SBATCH --qos=normal
#SBATCH --license=lustre1,lustre2
#SBATCH --export=ALL
#SBATCH -e error.txt
#SBATCH -o output.txt
# one MPI rank per half-socket (see below)
#SBATCH --tasks-per-node=2
# request all logical (virtual) cores per half-socket
#SBATCH --cpus-per-task=112


# each Dane node has 2 sockets of Intel Sapphire Rapids with 56 cores each
export WARPX_NMPI_PER_NODE=2

# each MPI rank per half-socket has 56 physical cores
#   or 112 logical (virtual) cores
# over-subscribing each physical core with 2x
#   hyperthreading led to a slight (3.5%) speedup on Cori's Intel Xeon E5-2698 v3,
#   so we do the same here
# the settings below make sure threads are close to the
#   controlling MPI rank (process) per half socket and
#   distribute equally over close-by physical cores and,
#   for N>9, also equally over close-by logical cores
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_NUM_THREADS=112

EXE="<path/to/executable>"  # e.g. ./warpx

srun --cpu_bind=cores -n $(( ${SLURM_JOB_NUM_NODES} * ${WARPX_NMPI_PER_NODE} )) ${EXE} <input file>
```

To run a simulation, copy the lines above to a file `dane.sbatch` and run

```bash
sbatch dane.sbatch
```

to submit the job.
