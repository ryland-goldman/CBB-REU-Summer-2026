<a id="building-karolina"></a>

# Karolina (IT4I)

The [Karolina cluster](https://docs.it4i.cz/karolina/introduction/) is located at [IT4I, Technical University of Ostrava](https://www.it4i.cz/en).

## Introduction

If you are new to this system, **please see the following resources**:

* [IT4I user guide](https://docs.it4i.cz)
* Batch system: [SLURM](https://docs.it4i.cz/general/job-submission-and-execution/)
* Jupyter service: not provided/documented (yet)
* [Filesystems](https://docs.it4i.cz/karolina/storage/):
  * `$HOME`: per-user directory, use only for inputs, source and scripts; backed up (25GB default quota)
  * `/scratch/`: [production directory](https://docs.it4i.cz/karolina/storage/#scratch-file-system); very fast for parallel jobs (10TB default)
  * `/mnt/proj<N>/<proj>`: per-project work directory, used for long term data storage (20TB default)

<a id="building-karolina-preparation"></a>

## Installation

We show how to install from scratch all the dependencies using [Spack](https://spack.io).

For size reasons it is not advisable to install WarpX in the `$HOME` directory, it should be installed in the “work directory”. For this purpose we set an environment variable `$WORK` with the path to the “work directory”.

On Karolina, you can run either on GPU nodes with fast A100 GPUs (recommended) or CPU nodes.

### Profile file

One can use the pre-prepared `karolina_warpx.profile` script below,
which you can copy to `${HOME}/karolina_warpx.profile`, edit as required and then `source`.

### Script Details

```bash
# please set your project account, ie DD-N-N
export proj="<proj_id>"  # change me!

# Name and Path of this Script ################### (DO NOT change!)
export MY_PROFILE=$(cd $(dirname $BASH_SOURCE) && pwd)"/"$(basename $BASH_SOURCE)

if [ -z ${proj-} ]; then
    echo "WARNING: The 'proj' variable is not yet set in your $MY_PROFILE file!"
    echo "Please edit its line 2 to continue!"
    return
fi

# set env variable storing the path to the work directory
# please check if your project ID belongs to proj1, proj2, proj3 etc
export WORK="/mnt/proj<N>/${proj,,}/${USER}"  # change me!
mkdir -p WORK

# clone warpx
# you can also clone your own fork here, eg git@github.com:<user>/WarpX.git
if [ ! -d "$WORK/src/warpx" ]
then
    git clone https://github.com/BLAST-WarpX/warpx.git $WORK/src/warpx
fi

# load required modules
module purge
module load OpenMPI/4.1.4-GCC-11.3.0-CUDA-11.7.0

source $WORK/spack/share/spack/setup-env.sh && spack env activate warpx-karolina-cuda && {
    echo "Spack environment 'warpx-karolina-cuda' activated successfully."
} || {
    echo "Failed to activate Spack environment 'warpx-karolina-cuda'. Please run install_dependencies.sh."
}

# Text Editor for Tools ########################## (edit this line)
# examples: "nano", "vim", "emacs -nw" or without terminal: "gedit"
#export EDITOR="nano"  # change me!

# allocate an interactive shell for one hour
# usage: getNode 2  # allocates two interactive nodes (default: 1)
function getNode() {
    if [ -z "$1" ] ; then
        numNodes=1
    else
        numNodes=$1
    fi
    export OMP_NUM_THREADS=16
    srun --time=1:00:00 --nodes=$numNodes --ntasks=$((8 * $numNodes)) --ntasks-per-node=8 --cpus-per-task=16 --exclusive --gpus-per-node=8 -p qgpu -A $proj --pty bash
}

# Environment #####################################################
# optimize CUDA compilation for A100
export AMREX_CUDA_ARCH="8.0"
export SCRATCH="/scratch/project/${proj,,}/${USER}"

# optimize CPU microarchitecture for AMD EPYC 7763 (zen3)
export CFLAGS="-march=znver3"
export CXXFLAGS="-march=znver3"

# compiler environment hints
export CC=$(which gcc)
export CXX=$(which g++)
export FC=$(which gfortran)
export CUDACXX=$(which nvcc)
export CUDAHOSTCXX=${CXX}
```

To have the environment activated on every login, add the following line to `${HOME}/.bashrc`:

```bash
source $HOME/karolina_warpx.profile
```

To install the `spack` environment and Python packages:

```bash
bash $WORK/src/warpx/Tools/machines/karolina-it4i/install_dependencies.sh
```

### Script Details

```bash
#!/bin/bash
#
# Copyright 2023 The WarpX Community
#
# This file is part of WarpX.
#
# Author: Axel Huebl, Andrei Berceanu
# License: BSD-3-Clause-LBNL

# Exit on first error encountered #################################
#
set -eu -o pipefail

# Check: ##########################################################
#
# Was karolina_warpx.profile sourced and configured correctly?
if [ -z ${proj-} ]; then
    echo "WARNING: The 'proj' variable is not yet set in your karolina_warpx.profile file!"
    echo "Please edit its line 2 to continue!"
    return
fi

# download and activate spack
# this might take about ~ 1 hour
if [ ! -d "$WORK/spack" ]
then
    git clone -c feature.manyFiles=true -b v0.21.0 https://github.com/spack/spack.git $WORK/spack
    source $WORK/spack/share/spack/setup-env.sh
else
    # If the directory exists, checkout v0.21.0 branch
    cd $WORK/spack
    git checkout v0.21.0
    git pull origin v0.21.0
    source $WORK/spack/share/spack/setup-env.sh

    # Delete spack env if present
    spack env deactivate || true
    spack env rm -y warpx-karolina-cuda || true

    cd -
fi

# create and activate the spack environment
spack env create warpx-karolina-cuda $WORK/src/warpx/Tools/machines/karolina-it4i/spack-karolina-cuda.yaml
spack env activate warpx-karolina-cuda
spack install

# Python ##########################################################
#
python -m pip install --user --upgrade pandas
python -m pip install --user --upgrade matplotlib
# optional
#python -m pip install --user --upgrade yt

# install or update WarpX dependencies
python -m pip install --user --upgrade picmistandard==0.34.0
python -m pip install --user --upgrade lasy

# optional: for optimas (based on libEnsemble & ax->botorch->gpytorch->pytorch)
# python -m pip install --user --upgrade -r $WORK/src/warpx/Tools/optimas/requirements.txt
```

<a id="building-karolina-compilation"></a>

## Compilation

Use the following [cmake commands](../cmake.md#install-build-cmake) to compile the application executable:

```bash
cd $WORK/src/warpx
rm -rf build_gpu

cmake -S . -B build_gpu -DWarpX_COMPUTE=CUDA -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_gpu -j 48
```

The WarpX application executables are now in `$WORK/src/warpx/build_gpu/bin/`.
Additionally, the following commands will install WarpX as a Python module:

```bash
cd $WORK/src/warpx
rm -rf build_gpu_py

cmake -S . -B build_gpu_py -DWarpX_COMPUTE=CUDA -DWarpX_FFT=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_APP=OFF -DWarpX_PYTHON=ON -DWarpX_DIMS="1;2;RZ;3"
cmake --build build_gpu_py -j 48 --target pip_install
```

Now, you can [submit Karolina compute jobs](#running-cpp-karolina) for WarpX [Python (PICMI) scripts](../../usage/python.md#usage-picmi) ([example scripts](../../usage/examples.md#usage-examples)).
Or, you can use the WarpX executables to submit Karolina jobs ([example inputs](../../usage/examples.md#usage-examples)).
For executables, you can reference their location in your [job script](#running-cpp-karolina) or copy them to a location in `/scratch/`.

<a id="running-cpp-karolina"></a>

## Running

The batch script below can be used to run a WarpX simulation on multiple GPU nodes (change `#SBATCH --nodes=` accordingly) on the supercomputer Karolina at IT4I.
This partition has up to [72 nodes](https://docs.it4i.cz/karolina/hardware-overview/).
Every node has 8x A100 (40GB) GPUs and 2x AMD EPYC 7763, 64-core, 2.45 GHz processors.

Replace descriptions between chevrons `<>` by relevant values, for instance `<proj>` could be `DD-23-83`.
Note that we run one MPI rank per GPU.

```bash
#!/bin/bash -l

# Copyright 2023 The WarpX Community
#
# This file is part of WarpX.
#
# Authors: Axel Huebl, Andrei Berceanu
# License: BSD-3-Clause-LBNL

#SBATCH --account=<proj>
#SBATCH --partition=qgpu
#SBATCH --time=00:10:00
#SBATCH --job-name=WarpX
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=8
#SBATCH --gpu-bind=single:1

#SBATCH --mail-type=ALL
# change me!
#SBATCH --mail-user=someone@example.com
#SBATCH --chdir=/scratch/project/<proj>/it4i-<user>/runs/warpx

#SBATCH -o stdout_%j
#SBATCH -e stderr_%j

# set user rights to u=rwx;g=r-x;o=---
umask 0027

# OpenMP threads per MPI rank
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# executable & inputs file or python interpreter & PICMI script here
EXE=./warpx.rz
INPUTS=./inputs_rz

# run
srun -K1 ${EXE} ${INPUTS}
```

To run a simulation, copy the lines above to a file `karolina_gpu.sbatch` and run

```bash
sbatch karolina_gpu.sbatch
```

to submit the job.

<a id="post-processing-karolina"></a>

## Post-Processing

#### NOTE
This section was not yet written.
Usually, we document here how to use a Jupyter service.
