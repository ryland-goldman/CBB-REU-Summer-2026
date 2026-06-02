<a id="building-greatlakes"></a>

# Great Lakes (UMich)

The [Great Lakes cluster](https://arc.umich.edu/greatlakes/) is located at University of Michigan.
The cluster has various partitions, including [GPU nodes and CPU nodes](https://arc.umich.edu/greatlakes/configuration/).

## Introduction

If you are new to this system, **please see the following resources**:

* [Great Lakes user guide](https://arc.umich.edu/greatlakes/)
* Batch system: [Slurm](https://arc.umich.edu/greatlakes/slurm-user-guide/)
* [Jupyter service](https://greatlakes.arc-ts.umich.edu) ([documentation](https://arc.umich.edu/greatlakes/user-guide/#document-2))
* [Filesystems](https://arc.umich.edu/greatlakes/user-guide/#document-1):
  * `$HOME`: per-user directory, use only for inputs, source and scripts; backed up (80GB)
  * `/scratch`: per-project [production directory](https://arc.umich.edu/greatlakes/user-guide/#scratchpolicies); very fast for parallel jobs; purged every 60 days (10TB default)

<a id="building-greatlakes-preparation"></a>

## Preparation

Use the following commands to download the WarpX source code:

```bash
git clone https://github.com/BLAST-WarpX/warpx.git $HOME/src/warpx
```

On Great Lakes, you can run either on GPU nodes with [fast V100 GPUs (recommended), the even faster A100 GPUs (only a few available) or CPU nodes](https://arc.umich.edu/greatlakes/configuration/).

### V100 GPUs

We use system software modules, add environment hints and further dependencies via the file `$HOME/greatlakes_v100_warpx.profile`.
Create it now:

```bash
cp $HOME/src/warpx/Tools/machines/greatlakes-umich/greatlakes_v100_warpx.profile.example $HOME/greatlakes_v100_warpx.profile
```

### Script Details

```bash
# please set your project account
export proj=""  # change me!

# remembers the location of this script
export MY_PROFILE=$(cd $(dirname $BASH_SOURCE) && pwd)"/"$(basename $BASH_SOURCE)
if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your $MY_PROFILE file! Please edit its line 2 to continue!"; return; fi

# required dependencies
module purge
module load gcc/10.3.0
module load cuda/12.1.1
module load cmake/3.26.3
module load openblas/0.3.23
module load openmpi/4.1.6-cuda

# optional: for QED support
module load boost/1.78.0

# optional: for openPMD and PSATD+RZ support
module load phdf5/1.12.1

SW_DIR="${HOME}/sw/greatlakes/v100"
export CMAKE_PREFIX_PATH=${SW_DIR}/c-blosc2-2.14.4:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/adios2-2.10.2:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/blaspp-2024.05.31:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/lapackpp-2024.05.31:$CMAKE_PREFIX_PATH

export LD_LIBRARY_PATH=${SW_DIR}/c-blosc2-2.14.4/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/adios2-2.10.2/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/blaspp-2024.05.31/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/lapackpp-2024.05.31/lib64:$LD_LIBRARY_PATH

export PATH=${SW_DIR}/adios2-2.10.2/bin:${PATH}

# optional: for Python bindings or libEnsemble
module load python/3.12.1

if [ -d "${SW_DIR}/venvs/warpx-v100" ]
then
  source ${SW_DIR}/venvs/warpx-v100/bin/activate
fi

# an alias to request an interactive batch node for one hour
#   for parallel execution, start on the batch node: srun <command>
alias getNode="salloc -N 1 --partition=gpu --ntasks-per-node=2 --cpus-per-task=20 --gpus-per-task=v100:1 -t 1:00:00 -A $proj"
# an alias to run a command on a batch node for up to 30min
#   usage: runNode <command>
alias runNode="srun -N 1 --partition=gpu --ntasks-per-node=2 --cpus-per-task=20 --gpus-per-task=v100:1 -t 1:00:00 -A $proj"

# optimize CUDA compilation for V100
export AMREX_CUDA_ARCH=7.0

# optimize CPU microarchitecture for Intel Xeon Gold 6148
export CXXFLAGS="-march=skylake-avx512"
export CFLAGS="-march=skylake-avx512"

# compiler environment hints
export CC=$(which gcc)
export CXX=$(which g++)
export FC=$(which gfortran)
export CUDACXX=$(which nvcc)
export CUDAHOSTCXX=${CXX}
```

Edit the 2nd line of this script, which sets the `export proj=""` variable.
For example, if you are member of the project `iloveplasma`, then run `nano $HOME/greatlakes_v100_warpx.profile` and edit line 2 to read:

```bash
export proj="iloveplasma"
```

Exit the `nano` editor with `Ctrl` + `O` (save) and then `Ctrl` + `X` (exit).

#### IMPORTANT
Now, and as the first step on future logins to Great Lakes, activate these environment settings:

```bash
source $HOME/greatlakes_v100_warpx.profile
```

Finally, since Great Lakes does not yet provide software modules for some of our dependencies, install them once:

```bash
bash $HOME/src/warpx/Tools/machines/greatlakes-umich/install_v100_dependencies.sh
source ${HOME}/sw/greatlakes/v100/venvs/warpx-v100/bin/activate
```

### Script Details

```bash
#!/bin/bash
#
# Copyright 2024 The WarpX Community
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
#   Was greatlakes_v100_warpx.profile sourced and configured correctly?
if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your greatlakes_v100_warpx.profile file! Please edit its line 2 to continue!"; exit 1; fi


# Remove old dependencies #####################################################
#
echo "Cleaning up prior installation directory... This may take several minutes."
SW_DIR="${HOME}/sw/greatlakes/v100"
rm -rf ${SW_DIR}
mkdir -p ${SW_DIR}

# remove common user mistakes in python, located in .local instead of a venv
python3 -m pip uninstall -qq -y pywarpx
python3 -m pip uninstall -qq -y warpx
python3 -m pip uninstall -qqq -y mpi4py 2>/dev/null || true


# General extra dependencies ##################################################
#

# tmpfs build directory: avoids issues often seen with $HOME and is faster
build_dir=$(mktemp -d)

# c-blosc (I/O compression)
if [ -d $HOME/src/c-blosc2 ]
then
  cd $HOME/src/c-blosc2
  git fetch --prune
  git checkout v2.14.4
  cd -
else
  git clone -b v2.14.4 https://github.com/Blosc/c-blosc2.git $HOME/src/c-blosc2
fi
rm -rf $HOME/src/c-blosc2-v100-build
cmake -S $HOME/src/c-blosc2 -B ${build_dir}/c-blosc2-v100-build -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF -DBUILD_EXAMPLES=OFF -DBUILD_FUZZERS=OFF -DDEACTIVATE_AVX2=OFF -DCMAKE_INSTALL_PREFIX=${SW_DIR}/c-blosc2-2.14.4
cmake --build ${build_dir}/c-blosc2-v100-build --target install --parallel 8
rm -rf ${build_dir}/c-blosc2-v100-build

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
rm -rf $HOME/src/adios2-v100-build
cmake                                \
  -S $HOME/src/adios2                \
  -B ${build_dir}/adios2-v100-build  \
  -DADIOS2_USE_Blosc2=ON             \
  -DADIOS2_USE_Campaign=OFF          \
  -DADIOS2_USE_Fortran=OFF           \
  -DADIOS2_USE_Python=OFF            \
  -DADIOS2_USE_ZeroMQ=OFF            \
  -DCMAKE_INSTALL_PREFIX=${SW_DIR}/adios2-2.10.2
cmake --build ${build_dir}/adios2-v100-build --target install -j 8
rm -rf ${build_dir}/adios2-v100-build

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
rm -rf $HOME/src/blaspp-v100-build
cmake -S $HOME/src/blaspp -B ${build_dir}/blaspp-v100-build -Duse_openmp=OFF -Dgpu_backend=cuda -DCMAKE_CXX_STANDARD=20 -DCMAKE_INSTALL_PREFIX=${SW_DIR}/blaspp-2024.05.31
cmake --build ${build_dir}/blaspp-v100-build --target install --parallel 8
rm -rf ${build_dir}/blaspp-v100-build

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
rm -rf $HOME/src/lapackpp-v100-build
CXXFLAGS="-DLAPACK_FORTRAN_ADD_" cmake -S $HOME/src/lapackpp -B ${build_dir}/lapackpp-v100-build -DCMAKE_CXX_STANDARD=20 -Dbuild_tests=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX=${SW_DIR}/lapackpp-2024.05.31
cmake --build ${build_dir}/lapackpp-v100-build --target install --parallel 8
rm -rf ${build_dir}/lapackpp-v100-build


# Python ######################################################################
#
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade virtualenv
python3 -m pip cache purge
rm -rf ${SW_DIR}/venvs/warpx-v100
python3 -m venv ${SW_DIR}/venvs/warpx-v100
source ${SW_DIR}/venvs/warpx-v100/bin/activate
python3 -m pip install --upgrade pip
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
# install or update WarpX dependencies
python3 -m pip install --upgrade -r $HOME/src/warpx/requirements.txt
python3 -m pip install --upgrade cupy-cuda12x  # CUDA 12 compatible wheel
# optimas (based on libEnsemble & ax->botorch->gpytorch->pytorch)
python3 -m pip install --upgrade torch  # CUDA 12 compatible wheel
python3 -m pip install --upgrade optimas[all]


# remove build temporary directory
rm -rf ${build_dir}
```

### A100 Nodes

#### NOTE
This section is TODO.

### CPU Nodes

#### NOTE
This section is TODO.

<a id="building-greatlakes-compilation"></a>

## Compilation

Use the following [cmake commands](../cmake.md#install-build-cmake) to compile the application executable:

### V100 GPUs

```bash
cd $HOME/src/warpx
rm -rf build_v100

cmake -S . -B build_v100 -DWarpX_COMPUTE=CUDA -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_v100 -j 8
```

The WarpX application executables are now in `$HOME/src/warpx/build_v100/bin/`.
Additionally, the following commands will install WarpX as a Python module:

```bash
cd $HOME/src/warpx
rm -rf build_v100_py

cmake -S . -B build_v100_py -DWarpX_COMPUTE=CUDA -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_APP=OFF -DWarpX_PYTHON=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_v100_py -j 8 --target pip_install
```

### A100 Nodes

#### NOTE
This section is TODO.

### CPU Nodes

#### NOTE
This section is TODO.

Now, you can [submit Great Lakes compute jobs](#running-cpp-greatlakes) for WarpX [Python (PICMI) scripts](../../usage/python.md#usage-picmi) ([example scripts](../../usage/examples.md#usage-examples)).
Or, you can use the WarpX executables to submit greatlakes jobs ([example inputs](../../usage/examples.md#usage-examples)).
For executables, you can reference their location in your [job script](#running-cpp-greatlakes) or copy them to a location in `/scratch`.

<a id="building-greatlakes-update"></a>

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

- [update the greatlakes_v100_warpx.profile file](#building-greatlakes-preparation),
- log out and into the system, activate the now updated environment profile as usual,
- [execute the dependency install scripts](#building-greatlakes-preparation).

As a last step, clean the build directory `rm -rf $HOME/src/warpx/build_*` and rebuild WarpX.

<a id="running-cpp-greatlakes"></a>

## Running

### V100 (16GB) GPUs

The batch script below can be used to run a WarpX simulation on multiple nodes (change `-N` accordingly) on the supercomputer Great Lakes at University of Michigan.
This partition has [20 nodes, each with two V100 GPUs](https://arc.umich.edu/greatlakes/configuration/).

Replace descriptions between chevrons `<>` by relevant values, for instance `<input file>` could be `plasma_mirror_inputs`.
Note that we run one MPI rank per GPU.

```bash
#!/bin/bash -l

# Copyright 2024 The WarpX Community
#
# Author: Axel Huebl
# License: BSD-3-Clause-LBNL

#SBATCH -t 00:10:00
#SBATCH -N 1
#SBATCH -J WarpX
#SBATCH -A <proj>
#SBATCH --partition=gpu
#SBATCH --exclusive
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=20
#SBATCH --gpus-per-task=v100:1
#SBATCH --gpu-bind=single:1
#SBATCH -o WarpX.o%j
#SBATCH -e WarpX.e%j

# executable & inputs file or python interpreter & PICMI script here
EXE=./warpx
INPUTS=inputs

# threads for OpenMP and threaded compressors per MPI rank
#   per node are 2x 2.4 GHz Intel Xeon Gold 6148
#   note: the system seems to only expose cores (20 per socket),
#         not hyperthreads (40 per socket)
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# GPU-aware MPI optimizations
GPU_AWARE_MPI="amrex.use_gpu_aware_mpi=1"

# run WarpX
srun --cpu-bind=cores \
  ${EXE} ${INPUTS} ${GPU_AWARE_MPI} \
  > output.txt
```

To run a simulation, copy the lines above to a file `greatlakes_v100.sbatch` and run

```bash
sbatch greatlakes_v100.sbatch
```

to submit the job.

### A100 (80GB) GPUs

This partition has [2 nodes, each with four A100 GPUs](https://arc.umich.edu/greatlakes/configuration/) that provide 80 GB HBM per A100 GPU.
To the user, each node will appear as if it has 8 A100 GPUs with 40 GB memory each.

#### NOTE
This section is TODO.

### CPU Nodes

The Great Lakes CPU partition as up to [455 nodes](https://arc.umich.edu/greatlakes/configuration/), each with 2x Intel Xeon Gold 6154 CPUs and 180 GB RAM.

#### NOTE
This section is TODO.

<a id="post-processing-greatlakes"></a>

## Post-Processing

For post-processing, many users prefer to use the online [Jupyter service](https://greatlakes.arc-ts.umich.edu) ([documentation](https://arc.umich.edu/greatlakes/user-guide/#document-2)) that is directly connected to the cluster’s fast filesystem.

#### NOTE
This section is a stub and contributions are welcome.
We can document further details, e.g., which recommended post-processing Python software to install or how to customize Jupyter kernels here.
