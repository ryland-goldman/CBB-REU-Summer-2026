# LUME-Impact Overview

LUME-Impact wraps the IMPACT-T and IMPACT-Z particle tracking codes in a Python interface. Install via conda-forge:

```bash
# Single-core
conda create -n impact -c conda-forge lume-impact

# OpenMPI
conda create -n impact -c conda-forge lume-impact impact-t=*=mpi_openmpi* impact-z=*=mpi_openmpi*

# MPICH
conda create -n impact -c conda-forge lume-impact impact-t=*=mpi_mpich* impact-z=*=mpi_mpich*

conda activate impact
```

After installation, `ImpactTexe` / `ImpactTexe-mpi` are available on PATH.

## IMPACT-T Quick Start

```python
from impact import Impact

I = Impact("/path/to/ImpactT.in", verbose=True)

# Modify header parameters
I.header["Np"] = 10000
I.header["Nx"] = 32
I.header["Ny"] = 32
I.header["Nz"] = 32

I.run()
I.plot()

# Archive all input/output to HDF5
I.archive("test.h5")

# Phase-space plot using openPMD-beamphysics
I.particles["final_particles"].plot("z", "pz")
```

## IMPACT-Z Quick Start

```python
from impact.z import ImpactZ

I = ImpactZ("/path/to/ImpactZ.in", verbose=True)
I.run()
I.plot()
```
