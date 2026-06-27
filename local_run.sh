#!/bin/bash
# Local launcher for the Xopt/CNSGA optimizer on this Mac (NOT the CLASSE cluster -- that path is
# launch_opt.sge + cluster_env.sh). config/xopt.yaml must have `executor: process`. The chain runs
# single-core per stage (OMP=1); max_workers in xopt.yaml caps concurrent full-chain evals.
#   bash local_run.sh   (long-running; CNSGA checkpoints to logs/opt/data.csv -- stop anytime, no resume)
source ~/miniforge3/bin/activate CBB
export OMP_NUM_THREADS=1                 # each stage is single-core (MLMG is memory-bandwidth bound)
export HDF5_USE_FILE_LOCKING=FALSE
export PYTHONNOUSERSITE=1                # ~/.local lume/openpmd would shadow the CBB env otherwise
export LINACSIM_RUNS_DIR=/tmp            # per-eval sandboxes under /tmp/linac_runs/<hash> (auto-cleaned on success)
export PATH="/Applications/G4beamline-3.08.app/Contents/MacOS:$PATH"   # converter g4bl (self-sets Geant4 data)
export PYTHONPATH="/Users/rylandgoldman/Documents/Coding/Cornell2:${PYTHONPATH:-}"
cd /Users/rylandgoldman/Documents/Coding/Cornell2 || exit 1
echo "host=$(hostname) python=$(command -v python) g4bl=$(command -v g4bl)"
exec python sim/optimize.py
