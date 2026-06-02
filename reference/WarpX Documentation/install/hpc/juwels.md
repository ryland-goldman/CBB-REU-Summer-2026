<a id="building-juwels"></a>

# Juwels (JSC)

#### NOTE
For the moment, WarpX doesn’t run on Juwels with MPI_THREAD_MULTIPLE.
Please compile with this compilation flag: `MPI_THREAD_MULTIPLE=FALSE`.

The [Juwels supercomputer](https://www.fz-juelich.de/ias/jsc/EN/Expertise/Supercomputers/JUWELS/JUWELS_node.html) is located at JSC.

## Introduction

If you are new to this system, **please see the following resources**:

See [this page](https://apps.fz-juelich.de/jsc/hps/juwels/quickintro.html) for a quick introduction.
(Full [user guide](http://www.fz-juelich.de/ias/jsc/EN/Expertise/Supercomputers/JUWELS/UserInfo/UserInfo_node.html)).

* Batch system: [Slurm](https://apps.fz-juelich.de/jsc/hps/juwels/quickintro.html#batch-system-on-system-name)
* [Production directories](https://apps.fz-juelich.de/jsc/hps/juwels/environment.html?highlight=scratch#available-filesystems):
  * `$SCRATCH`: Scratch filesystem for [temporary data](http://www.fz-juelich.de/ias/jsc/EN/Expertise/Supercomputers/JUWELS/FAQ/juwels_FAQ_node.html#faq1495160) (90 day purge)
  * `$FASTDATA/`: Storage location for large data (backed up)
  * Note that the `$HOME` directory is not designed for simulation runs and producing output there will impact performance.

## Installation

Use the following commands to download the WarpX source code and switch to the correct branch:

```bash
git clone https://github.com/BLAST-WarpX/warpx.git $HOME/src/warpx
```

We use the following modules and environments on the system.

```bash
# please set your project account
#export proj=<yourProject>

# required dependencies
module load ccache
module load CMake
module load GCC
module load CUDA/11.3
module load OpenMPI
module load FFTW
module load HDF5
module load Python

# JUWELS' job scheduler may not map ranks to GPUs,
# so we give a hint to AMReX about the node layout.
# This is usually done in Make.<supercomputing center> files in AMReX
# but there is no such file for JSC yet.
export GPUS_PER_SOCKET=2
export GPUS_PER_NODE=4

# optimize CUDA compilation for V100 (7.0) or for A100 (8.0)
export AMREX_CUDA_ARCH=8.0
```

Note that for now WarpX must rely on OpenMPI instead of the recommended MPI implementation on this platform MVAPICH2.

We recommend to store the above lines in a file, such as `$HOME/juwels_warpx.profile`, and load it into your shell after a login:

```bash
source $HOME/juwels_warpx.profile
```

Then, `cd` into the directory `$HOME/src/warpx` and use the following commands to compile:

```bash
cd $HOME/src/warpx
rm -rf build

cmake -S . -B build -DWarpX_DIMS="1;2;3" -DWarpX_COMPUTE=CUDA -DWarpX_FFT=ON -DWarpX_MPI_THREAD_MULTIPLE=OFF
cmake --build build -j 16
```

The other [general compile-time options](../cmake.md#install-build-cmake) apply as usual.

**That’s it!**
A 3D WarpX executable is now in `build/bin/` and [can be run](#running-cpp-juwels) with a [3D example inputs file](../../usage/examples.md#usage-examples).
Most people execute the binary directly or copy it out to a location in `$SCRATCH`.

#### NOTE
Currently, if you want to use HDF5 output with openPMD, you need to add

```bash
export OMPI_MCA_io=romio321
```

in your job scripts, before running the `srun` command.

<a id="running-cpp-juwels"></a>

## Running

### Queue: gpus (4 x Nvidia V100 GPUs)

The [Juwels GPUs](https://apps.fz-juelich.de/jsc/hps/juwels/configuration.html) are V100 (16GB) and A100 (40GB).

An example submission script reads

```bash
#!/bin/bash -l

#SBATCH -A $proj
#SBATCH --partition=booster
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --time=00:05:00
#SBATCH --job-name=warpx
#SBATCH --output=warpx-%j-%N.txt
#SBATCH --error=warpx-%j-%N.err

export OMP_NUM_THREADS=1
export OMPI_MCA_io=romio321  # for HDF5 support in openPMD

# you can comment this out if you sourced the warpx.profile
# files before running sbatch:
module load GCC
module load OpenMPI
module load CUDA/11.3
module load HDF5
module load Python

srun -n 8 --cpu_bind=sockets $HOME/src/warpx/build/bin/warpx.3d.MPI.CUDA.DP.OPMD.QED inputs
```

### Queue: batch (2 x Intel Xeon Platinum 8168 CPUs, 24 Cores + 24 Hyperthreads/CPU)

*todo*

See the [data analysis section](../../dataanalysis/formats.md#dataanalysis-formats) for more information on how to visualize the simulation results.
