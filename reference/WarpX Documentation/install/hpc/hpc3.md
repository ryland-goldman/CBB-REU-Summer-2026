<a id="building-hpc3"></a>

# HPC3 (UCI)

The [HPC3 supercomputer](https://rcic.uci.edu/hpc3/index.html) is located at University of California, Irvine.

## Introduction

If you are new to this system, **please see the following resources**:

* [HPC3 user guide](https://rcic.uci.edu/hpc3/index.html)
* Batch system: [Slurm](https://rcic.uci.edu/hpc3/slurm.html) ([notes](https://rcic.uci.edu/hpc3/examples.html#_submit_different_job_types))
* [Jupyter service](https://rcic.uci.edu/hpc3/examples.html#jupyterhub-portal)
* [Filesystems](https://rcic.uci.edu/storage/beegfs-howtos.html):
  * `$HOME`: per-user directory, use only for inputs, source and scripts; backed up (40GB)
  * `/pub/$USER`: per-user production directory; fast and larger storage for parallel jobs (1TB default quota)
  * `/dfsX/<lab-path>` lab group quota (based on PI’s purchase allocation). The storage owner (PI) can specify what users have read/write capability on the specific filesystem.

<a id="building-hpc3-preparation"></a>

## Preparation

Use the following commands to download the WarpX source code:

```bash
git clone https://github.com/BLAST-WarpX/warpx.git $HOME/src/warpx
```

On HPC3, you recommend to run on the [fast GPU nodes with V100 GPUs](https://rcic.uci.edu/hpc3/slurm.html#memmap).

We use system software modules, add environment hints and further dependencies via the file `$HOME/hpc3_gpu_warpx.profile`.
Create it now:

```bash
cp $HOME/src/warpx/Tools/machines/hpc3-uci/hpc3_gpu_warpx.profile.example $HOME/hpc3_gpu_warpx.profile
```

### Script Details

```bash
# please set your project account
export proj=""  # change me! GPU projects must end in "..._g"

# remembers the location of this script
export MY_PROFILE=$(cd $(dirname $BASH_SOURCE) && pwd)"/"$(basename $BASH_SOURCE)
if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your $MY_PROFILE file! Please edit its line 2 to continue!"; return; fi

# required dependencies
module load cmake/3.30.2
module load gcc/11.2.0
module load cuda/11.7.1
module load openmpi/4.1.2/gcc.11.2.0

# optional: for QED support with detailed tables
module load boost/1.78.0/gcc.11.2.0

# optional: for openPMD and PSATD+RZ support
module load OpenBLAS/0.3.21
module load hdf5/1.13.1/gcc.11.2.0-openmpi.4.1.2
export CMAKE_PREFIX_PATH=${HOME}/sw/hpc3/gpu/c-blosc-1.21.1:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${HOME}/sw/hpc3/gpu/adios2-2.10.2:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${HOME}/sw/hpc3/gpu/blaspp-2024.05.31:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${HOME}/sw/hpc3/gpu/lapackpp-2024.05.31:$CMAKE_PREFIX_PATH

export LD_LIBRARY_PATH=${HOME}/sw/hpc3/gpu/c-blosc-1.21.1/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${HOME}/sw/hpc3/gpu/adios2-2.10.2/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${HOME}/sw/hpc3/gpu/blaspp-2024.05.31/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${HOME}/sw/hpc3/gpu/lapackpp-2024.05.31/lib64:$LD_LIBRARY_PATH

export PATH=${HOME}/sw/hpc3/gpu/adios2-2.10.2/bin:${PATH}

# optional: CCache
#module load ccache  # missing

# optional: for Python bindings
# Any supported Python >=3.11 module should work
module load python/3.14.3

if [ -d "${HOME}/sw/hpc3/gpu/venvs/warpx-gpu" ]
then
  source ${HOME}/sw/hpc3/gpu/venvs/warpx-gpu/bin/activate
fi

# an alias to request an interactive batch node for one hour
#   for parallel execution, start on the batch node: srun <command>
alias getNode="salloc -N 1 -t 0:30:00 --gres=gpu:V100:1 -p free-gpu"
# an alias to run a command on a batch node for up to 30min
#   usage: runNode <command>
alias runNode="srun -N 1 -t 0:30:00 --gres=gpu:V100:1 -p free-gpu"

# optimize CUDA compilation for V100
export AMREX_CUDA_ARCH=7.0

# compiler environment hints
export CXX=$(which g++)
export CC=$(which gcc)
export FC=$(which gfortran)
export CUDACXX=$(which nvcc)
export CUDAHOSTCXX=${CXX}
```

Edit the 2nd line of this script, which sets the `export proj=""` variable.
For example, if you are member of the project `plasma`, then run `vi $HOME/hpc3_gpu_warpx.profile`.
Enter the edit mode by typing `i` and edit line 2 to read:

```bash
export proj="plasma"
```

Exit the `vi` editor with `Esc` and then type `:wq` (write & quit).

#### IMPORTANT
Now, and as the first step on future logins to HPC3, activate these environment settings:

```bash
source $HOME/hpc3_gpu_warpx.profile
```

Finally, since HPC3 does not yet provide software modules for some of our dependencies, install them once:

```bash
bash $HOME/src/warpx/Tools/machines/hpc3-uci/install_gpu_dependencies.sh
source $HOME/sw/hpc3/gpu/venvs/warpx-gpu/bin/activate
```

### Script Details

```bash
#!/bin/bash
#
# Copyright 2023 The WarpX Community
#
# This file is part of WarpX.
#
# Authors: Axel Huebl, Victor Flores
# License: BSD-3-Clause-LBNL

# Exit on first error encountered #############################################
#
set -eu -o pipefail


# Check: ######################################################################
#
#   Was hpc3_gpu_warpx.profile sourced and configured correctly?
if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your hpc3_gpu_warpx.profile file! Please edit its line 2 to continue!"; exit 1; fi


# Check $proj variable is correct and has a corresponding project directory ####
#
#if [ ! -d "/dfsX/${proj}/" ]
#then
#    echo "WARNING: The directory /dfsX/${proj}/ does not exist!"
#    echo "Is the \$proj environment variable of value \"$proj\" correctly set? "
#    echo "Please edit line 2 of your hpc3_gpu_warpx.profile file to continue!"
##    exit
fi


# Remove old dependencies #####################################################
#
SW_DIR="${HOME}/sw/hpc3/gpu"
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
cmake --build $HOME/src/c-blosc-pm-gpu-build --target install --parallel 8
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
cmake -S $HOME/src/adios2 -B $HOME/src/adios2-pm-gpu-build -DBUILD_TESTING=OFF -DADIOS2_BUILD_EXAMPLES=OFF -DADIOS2_USE_Blosc=ON -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_HDF5=OFF -DADIOS2_USE_Python=OFF -DADIOS2_USE_ZeroMQ=OFF -DCMAKE_INSTALL_PREFIX=${SW_DIR}/adios2-2.10.2
cmake --build $HOME/src/adios2-pm-gpu-build --target install --parallel 8
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
cmake -S $HOME/src/blaspp -B $HOME/src/blaspp-pm-gpu-build -Duse_openmp=OFF -Dgpu_backend=cuda -DCMAKE_CXX_STANDARD=20 -DCMAKE_INSTALL_PREFIX=${SW_DIR}/blaspp-2024.05.31
cmake --build $HOME/src/blaspp-pm-gpu-build --target install --parallel 8
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
CXXFLAGS="-DLAPACK_FORTRAN_ADD_" cmake -S $HOME/src/lapackpp -B $HOME/src/lapackpp-pm-gpu-build -DCMAKE_CXX_STANDARD=20 -Dbuild_tests=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX=${SW_DIR}/lapackpp-2024.05.31
cmake --build $HOME/src/lapackpp-pm-gpu-build --target install --parallel 8
rm -rf $HOME/src/lapackpp-pm-gpu-build


# Python ######################################################################
#
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade virtualenv
python3 -m pip cache purge
rm -rf ${SW_DIR}/venvs/warpx-gpu
python3 -m venv ${SW_DIR}/venvs/warpx-gpu
source ${SW_DIR}/venvs/warpx-gpu/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade build
python3 -m pip install --upgrade packaging
python3 -m pip install --upgrade wheel
python3 -m pip install --upgrade setuptools[core]
python3 -m pip install --upgrade pipx
python3 -m pip install --upgrade cython
python3 -m pip install --upgrade numpy
python3 -m pip install --upgrade pandas
python3 -m pip install --upgrade scipy
python3 -m pip install --upgrade mpi4py --no-cache-dir --no-build-isolation --no-binary mpi4py
python3 -m pip install --upgrade openpmd-api
python3 -m pip install --upgrade matplotlib
python3 -m pip install --upgrade yt
# install or update WarpX dependencies such as picmistandard
python3 -m pip install --upgrade -r $HOME/src/warpx/requirements.txt
# optional: for libEnsemble
python3 -m pip install -r $HOME/src/warpx/Tools/LibEnsemble/requirements.txt
# optional: for optimas (based on libEnsemble & ax->botorch->gpytorch->pytorch)
python3 -m pip install --upgrade torch  # CUDA 11.7 compatible wheel
python3 -m pip install -r $HOME/src/warpx/Tools/optimas/requirements.txt
```

<a id="building-hpc3-compilation"></a>

## Compilation

Use the following [cmake commands](../cmake.md#install-build-cmake) to compile the application executable:

```bash
cd $HOME/src/warpx
rm -rf build

cmake -S . -B build -DWarpX_COMPUTE=CUDA -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build -j 8
```

The WarpX application executables are now in `$HOME/src/warpx/build/bin/`.
Additionally, the following commands will install WarpX as a Python module:

```bash
rm -rf build_py

cmake -S . -B build_py -DWarpX_COMPUTE=CUDA -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_APP=OFF -DWarpX_PYTHON=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_py -j 8 --target pip_install
```

Now, you can [submit HPC3 compute jobs](#running-cpp-hpc3) for WarpX [Python (PICMI) scripts](../../usage/python.md#usage-picmi) ([example scripts](../../usage/examples.md#usage-examples)).
Or, you can use the WarpX executables to submit HPC3 jobs ([example inputs](../../usage/examples.md#usage-examples)).
For executables, you can reference their location in your [job script](#running-cpp-hpc3) or copy them to a location in `$PSCRATCH`.

<a id="building-hpc3-update"></a>

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

- [update the hpc3_gpu_warpx.profile file](#building-hpc3-preparation),
- log out and into the system, activate the now updated environment profile as usual,
- [execute the dependency install scripts](#building-hpc3-preparation).

As a last step, clean the build directory `rm -rf $HOME/src/warpx/build` and rebuild WarpX.

<a id="running-cpp-hpc3"></a>

## Running

The batch script below can be used to run a WarpX simulation on multiple nodes (change `-N` accordingly) on the supercomputer HPC3 at UCI.
This partition as up to [32 nodes](https://rcic.uci.edu/hpc3/slurm.html#memmap) with four V100 GPUs (16 GB each) per node.

Replace descriptions between chevrons `<>` by relevant values, for instance `<proj>` could be `plasma`.
Note that we run one MPI rank per GPU.

```bash
#!/bin/bash -l

# Copyright 2023 The WarpX Community
#
# This file is part of WarpX.
#
# Authors: Axel Huebl, Victor Flores
# License: BSD-3-Clause-LBNL

#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH -J WarpX
#S BATCH -A <proj>
# V100 GPU options: gpu, free-gpu, debug-gpu
#SBATCH -p free-gpu
# use all four GPUs per node
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:V100:4
#SBATCH --cpus-per-task=10
#S BATCH --mail-type=begin,end
#S BATCH --mail-user=<your-email>@uci.edu
#SBATCH -o WarpX.o%j
#SBATCH -e WarpX.e%j

# executable & inputs file or python interpreter & PICMI script here
EXE=./warpx.rz
INPUTS=inputs_rz

# OpenMP threads
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# run
mpirun -np ${SLURM_NTASKS} bash -c "
    export CUDA_VISIBLE_DEVICES=\${SLURM_LOCALID};
    ${EXE} ${INPUTS}" \
  > output.txt
```

To run a simulation, copy the lines above to a file `hpc3_gpu.sbatch` and run

```bash
sbatch hpc3_gpu.sbatch
```

to submit the job.

<a id="post-processing-hpc3"></a>

## Post-Processing

UCI provides a pre-configured [Jupyter service](https://rcic.uci.edu/hpc3/examples.html#jupyterhub-portal) that can be used for data-analysis.

We recommend to install at least the following `pip` packages for running Python3 Jupyter notebooks on WarpX data analysis:
`h5py ipympl ipywidgets matplotlib numpy openpmd-viewer openpmd-api pandas scipy yt`
