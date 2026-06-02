<a id="building-taurus"></a>

# Taurus (ZIH)

The [Taurus cluster](https://doc.zih.tu-dresden.de/jobs_and_resources/overview) is located at [ZIH (TU Dresden)](https://doc.zih.tu-dresden.de).

The cluster has multiple partitions, this section describes how to use the [AMD Rome CPUs + NVIDIA A100¶](https://doc.zih.tu-dresden.de/jobs_and_resources/hardware_overview/#amd-rome-cpus-nvidia-a100).

## Introduction

If you are new to this system, **please see the following resources**:

* [ZIH user guide](https://docs.nersc.gov/)
* Batch system: [Slurm](https://docs.nersc.gov/systems/perlmutter/#running-jobs)
* Jupyter service: Missing?
* [Production directories](https://docs.nersc.gov/filesystems/perlmutter-scratch/):
  * `$PSCRATCH`: per-user production directory, purged every 30 days (<TBD>TB)
  * `/global/cscratch1/sd/m3239`: shared production directory for users in the project `m3239`, purged every 30 days (50TB)
  * `/global/cfs/cdirs/m3239/`: community file system for users in the project `m3239` (100TB)

## Installation

Use the following commands to download the WarpX source code and switch to the correct branch:

```bash
git clone https://github.com/BLAST-WarpX/warpx.git $HOME/src/warpx
```

We use the following modules and environments on the system (`$HOME/taurus_warpx.profile`).

```bash
# please set your project account
#export proj="<yourProject>"  # change me

# required dependencies
module load modenv/hiera
module load foss/2021b
module load CUDA/11.8.0
module load CMake/3.27.6

# optional: for QED support with detailed tables
#module load Boost  # TODO

# optional: for openPMD and PSATD+RZ support
module load HDF5/1.13.1

# optional: for Python bindings or libEnsemble
#module load python  # TODO
#
#if [ -d "$HOME/sw/taurus/venvs/warpx" ]
#then
#  source $HOME/sw/taurus/venvs/warpx/bin/activate
#fi

# an alias to request an interactive batch node for one hour
#   for parallel execution, start on the batch node: srun <command>
alias getNode="salloc --time=2:00:00 -N1 -n1 --cpus-per-task=6 --mem-per-cpu=2048 --gres=gpu:1 --gpu-bind=single:1 -p alpha-interactive --pty bash"
# an alias to run a command on a batch node for up to 30min
#   usage: runNode <command>
alias runNode="srun --time=2:00:00 -N1 -n1 --cpus-per-task=6 --mem-per-cpu=2048 --gres=gpu:1 --gpu-bind=single:1 -p alpha-interactive --pty bash"

# optimize CUDA compilation for A100
export AMREX_CUDA_ARCH=8.0

# compiler environment hints
#export CC=$(which gcc)
#export CXX=$(which g++)
#export FC=$(which gfortran)
#export CUDACXX=$(which nvcc)
#export CUDAHOSTCXX=${CXX}
```

We recommend to store the above lines in a file, such as `$HOME/taurus_warpx.profile`, and load it into your shell after a login:

```bash
source $HOME/taurus_warpx.profile
```

Then, `cd` into the directory `$HOME/src/warpx` and use the following commands to compile:

```bash
cd $HOME/src/warpx
rm -rf build

cmake -S . -B build -DWarpX_DIMS="1;2;3" -DWarpX_COMPUTE=CUDA
cmake --build build -j 16
```

The general [cmake compile-time options](../cmake.md#install-build-cmake) apply as usual.

<a id="running-cpp-taurus"></a>

## Running

<a id="running-cpp-taurus-a100-gpus"></a>

### A100 GPUs (40 GB)

The alpha partition has 34 nodes with 8 x NVIDIA A100-SXM4 Tensor Core-GPUs and 2 x AMD EPYC CPU 7352 (24 cores) @ 2.3 GHz (multithreading disabled) per node.

The batch script below can be used to run a WarpX simulation on multiple nodes (change `-N` accordingly).
Replace descriptions between chevrons `<>` by relevant values, for instance `<input file>` could be `plasma_mirror_inputs`.
Note that we run one MPI rank per GPU.

```bash
#!/bin/bash -l

# Copyright 2023 Axel Huebl, Thomas Miethlinger
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

#SBATCH -t 00:10:00
#SBATCH -N 1
#SBATCH -J WarpX
#SBATCH -p alpha
#SBATCH --exclusive
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=2048
#SBATCH --gres=gpu:1
#SBATCH --gpu-bind=single:1
#SBATCH -o WarpX.o%j
#SBATCH -e WarpX.e%j

# executable & inputs file or python interpreter & PICMI script here
EXE=./warpx
INPUTS=inputs_small

# run
srun ${EXE} ${INPUTS} \
  > output.txt
```

To run a simulation, copy the lines above to a file `taurus.sbatch` and run

```bash
sbatch taurus.sbatch
```

to submit the job.
