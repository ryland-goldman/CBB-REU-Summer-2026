# Site environment for running the linac sim chain on the CLASSE SGE cluster.
# ONE source of truth: `launch_opt.sge` and the dask-sge worker jobs (sim/optimize.make_executor)
# both `source` this, so the controller and every worker share an identical environment.
# Edit the paths below for your checkout. NOT used for local runs (it's CLASSE-specific).

# ---- paths (edit for your scratch layout) ----
export LINACSIM_BASE=/nfs/acc/temp/rjg343/Cornell2     # this repo (REPO_ROOT) on shared /nfs scratch
CBB_ENV=/nfs/acc/temp/rjg343/cbb                       # scratch conda env (full stack + xopt + dask-jobqueue)
G4BL_DIR=/nfs/acc/temp/rjg343/G4beamline-3.08          # G4beamline 3.08 install
G4DATA=/nfs/acc/temp/rjg343/Geant4Data                 # Geant4 11.0 datasets (versions g4bl 3.08 needs)

# ---- conda env ----
source ~/miniforge3/etc/profile.d/conda.sh
conda activate "$CBB_ENV"

# ---- WarpX/Impact-T/HDF5 (mirrors injphase/submit.sge) ----
export OMP_NUM_THREADS=1                # each eval is single-core (MLMG is memory-bandwidth bound)
export HDF5_USE_FILE_LOCKING=FALSE
# NB: do NOT set OPENPMD_DEFER_ITERATION_PARSING=1 -- buildfields.py reads chk.iterations[0].meshes
# eagerly, which deferred parsing leaves empty (IndexError: Key 'B' does not exist).
export PYTHONNOUSERSITE=1               # ~/.local lume/openpmd shadow the env otherwise
export LINACSIM_RUNS_DIR=/tmp               # per-eval sandboxes under /tmp/linac_runs/<hash> on NODE-LOCAL
                                           # /tmp, never /nfs: each eval writes hundreds of openPMD diag
                                           # dumps; on shared NFS that saturates the server and the whole
                                           # job stalls in disk-sleep
export PYTHONPATH="$LINACSIM_BASE:${PYTHONPATH:-}"   # dask-sge worker jobs must import sim.optimize.evaluate
                                           # (the controller adds it at runtime, but a bare worker does not)

# ---- G4beamline + Geant4 11.0 data (converter stage) ----
# g4bldata's GUI downloader can't run headless, so the dataset paths are set directly. Geant4 pins
# exact versions: g4bl 3.08 == Geant4 11.0 == G4EMLOW8.0 / G4ENSDFSTATE2.3 / G4PARTICLEXS4.0 / ...
export PATH="$G4BL_DIR/bin:$PATH"
# ONE-TIME DEPLOY STEP (not an env var): the g4bl binary spawns $G4BL_DIR/bin/g4bldata (a Qt GUI data
# downloader) on every launch and BLOCKS waiting for it -- headless, that GUI hangs and g4bl deadlocks
# (confirmed via top). Replace g4bldata with a no-op so g4bl proceeds straight to the sim; the data
# comes from the G4*DATA vars below, not g4bldata:
#   cd $G4BL_DIR/bin && mv g4bldata g4bldata.real && printf '#!/bin/bash\nexit 0\n' > g4bldata && chmod +x g4bldata
export G4LEDATA="$G4DATA/G4EMLOW8.0"
export G4ENSDFSTATEDATA="$G4DATA/G4ENSDFSTATE2.3"
export G4PARTICLEXSDATA="$G4DATA/G4PARTICLEXS4.0"
export G4LEVELGAMMADATA="$G4DATA/PhotonEvaporation5.7"
export G4RADIOACTIVEDATA="$G4DATA/RadioactiveDecay5.6"
export G4SAIDXSDATA="$G4DATA/G4SAIDDATA2.0"
export G4ABLADATA="$G4DATA/G4ABLA3.1"
export G4INCLDATA="$G4DATA/G4INCL1.0"
export G4PIIDATA="$G4DATA/G4PII1.3"
export G4REALSURFACEDATA="$G4DATA/RealSurface2.2"
export G4NEUTRONHPDATA="$G4DATA/G4NDL4.6"
