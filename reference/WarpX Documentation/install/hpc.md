<a id="install-hpc"></a>

# HPC Systems

On selected high-performance computing (HPC) systems, WarpX has documented or even pre-build installation routines.
Follow the guide here instead of the generic installation routines for optimal stability and best performance.

<a id="install-hpc-profile"></a>

## warpx.profile

Use a `warpx.profile` file to set up your software environment without colliding with other software.
Ideally, store that file directly in your `$HOME/` and source it after connecting to the machine:

```bash
source $HOME/warpx.profile
```

We list example `warpx.profile` files below, which can be used to set up WarpX on various HPC systems.

<a id="install-hpc-machines"></a>

## HPC Machines

This section documents quick-start guides for a selection of supercomputers that WarpX users are active on.

* [Adastra (CINES)](hpc/adastra.md)
* [Aurora (ALCF)](hpc/aurora.md)
* [Crusher (OLCF)](hpc/crusher.md)
* [Frontier (OLCF)](hpc/frontier.md)
* [Fugaku (Riken)](hpc/fugaku.md)
* [Great Lakes (UMich)](hpc/greatlakes.md)
* [HPC3 (UCI)](hpc/hpc3.md)
* [Juwels (JSC)](hpc/juwels.md)
* [Karolina (IT4I)](hpc/karolina.md)
* [Lawrencium (LBNL)](hpc/lawrencium.md)
* [Leonardo (CINECA)](hpc/leonardo.md)
* [Lonestar6 (TACC)](hpc/lonestar6.md)
* [LUMI (CSC)](hpc/lumi.md)
* [LXPLUS (CERN)](hpc/lxplus.md)
* [Ookami (Stony Brook)](hpc/ookami.md)
* [Perlmutter (NERSC)](hpc/perlmutter.md)
* [Pitzer (OSC)](hpc/pitzer.md)
* [Polaris (ALCF)](hpc/polaris.md)
* [Dane (LLNL)](hpc/dane.md)
* [Taurus (ZIH)](hpc/taurus.md)
* [Tuolumne (LLNL)](hpc/tuolumne.md)

<a id="install-hpc-batch"></a>

## Batch Systems

HPC systems use a scheduling (“batch”) system for time sharing of computing resources.
The batch system is used to request, queue, schedule and execute compute jobs asynchronously.
The individual HPC machines above document job submission example scripts, as templates for your modifications.

In this section, we document a quick reference guide (or cheat sheet) to interact in more detail with the various batch systems that you might encounter on different systems.

### Slurm

Slurm is a modern and very popular batch system.
Slurm is used at NERSC, OLCF Frontier, among others.

#### Job Submission

* `sbatch your_job_script.sbatch`

#### Job Control

* interactive job:
  * `salloc --time=1:00:00 --nodes=1 --ntasks-per-node=4 --cpus-per-task=8`
    * e.g. `srun "hostname"`
  * GPU allocation on most machines require additional flags, e.g. `--gpus-per-task=1` or `--gres=...`
* details for my jobs:
  * `scontrol -d show job 12345` all details for job with <job id> `12345`
  * `squeue -u $(whoami) -l` all jobs under my user name
* details for queues:
  * `squeue -p queueName -l` list full queue
  * `squeue -p queueName --start` (show start times for pending jobs)
  * `squeue -p queueName -l -t R` (only show running jobs in queue)
  * `sinfo -p queueName` (show online/offline nodes in queue)
  * `sview` (alternative on taurus: `module load llview` and `llview`)
  * `scontrol show partition queueName`
* communicate with job:
  * `scancel <job id>` abort job
  * `scancel -s <signal number> <job id>` send signal or signal name to job
  * `scontrol update timelimit=4:00:00 jobid=12345` change the walltime of a job
  * `scontrol update jobid=12345 dependency=afterany:54321` only start job `12345` after job with id `54321` has finished
  * `scontrol hold <job id>` prevent the job from starting
  * `scontrol release <job id>` release the job to be eligible for run (after it was set on hold)

#### References

* [https://slurm.schedmd.com/documentation.html](https://slurm.schedmd.com/documentation.html)

### Flux

Flux is a modern batch system and resource manager framework.
Flux is used at LLNL LC, among others.

#### Job Submission

* `flux batch your_job_script.flux`

#### Job Control

* [interactive job](https://flux-framework.readthedocs.io/projects/flux-core/en/latest/man1/flux-submit.html):
  * `flux submit --time-limit=1:00:00 --nodes=1 --tasks-per-node=4 --cores-per-task=8`
    * e.g. `flux submit "hostname"`
  * GPU allocation requires additional flags, e.g. `--gpus-per-task=1`
* details for my jobs:
  * `flux jobs` all jobs under my user name
  * `flux job info abc123 jobspec` all details for job with <job id> `abc123`
  * `flux job info 12345 eventlog` history of events for job with <job id> `12345`
* details for queues:
  * `flux queue list` list all queues
  * `flux queue status` show status of queues
  * *unclear/TODO* show start times for pending jobs
  * `sinfo -p queueName` show online/offline nodes in queue
* communicate with job:
  * `flux cancel <job id>` abort job
  * `flux job kill --signal=<signal number> <job id>` send signal or signal name to job
  * *unclear/TODO* change the walltime of a job
  * *unclear/TODO* only start job `12345` after job with id `54321` has finished
  * `flux job urgency <job id> hold` prevent the job from starting
  * `flux job urgency <job id> default` release the job to be eligible for run (after it was set on hold)

#### References

* [Flux commands](https://flux-framework.readthedocs.io/projects/flux-core/en/latest/man1/index.html)

### LSF

LSF (for *Load Sharing Facility*) is an IBM batch system.
It is used at OLCF Summit, LLNL Lassen, and other IBM systems.

#### Job Submission

* `bsub your_job_script.bsub`

#### Job Control

* interactive job:
  * `bsub -P $proj -W 2:00 -nnodes 1 -Is /bin/bash`
* [details for my jobs](https://docs.olcf.ornl.gov/systems/summit_user_guide.html#monitoring-jobs):
  * `bjobs 12345` all details for job with <job id> `12345`
  * `bjobs [-l]` all jobs under my user name
  * `jobstat -u $(whoami)` job eligibility
  * `bjdepinfo 12345` job dependencies on other jobs
* details for queues:
  * `bqueues` list queues
* communicate with job:
  * `bkill <job id>` abort job
  * `bpeek [-f] <job id>` peek into `stdout`/`stderr` of a job
  * `bkill -s <signal number> <job id>` send signal or signal name to job
  * `bchkpnt` and `brestart` checkpoint and restart job (untested/unimplemented)
  * `bmod -W 1:30 12345` change the walltime of a job (currently not allowed)
  * `bstop <job id>` prevent the job from starting
  * `bresume <job id>` release the job to be eligible for run (after it was set on hold)

#### References

* [https://www.ibm.com/docs/en/spectrum-lsf](https://www.ibm.com/docs/en/spectrum-lsf)

### PBS

PBS (for *Portable Batch System*) is a popular HPC batch system.
The OpenPBS project is related to [PBS, PBS Pro and TORQUE](https://en.wikipedia.org/wiki/Portable_Batch_System).

#### Job Submission

* `qsub your_job_script.qsub`

#### Job Control

* interactive job:
  * `qsub -I`
* details for my jobs:
  * `qstat -f 12345` all details for job with <job id> `12345`
  * `qstat -u $(whoami)` all jobs under my user name
* details for queues:
  * `qstat -a queueName` show all jobs in a queue
  * `pbs_free -l` compact view on free and busy nodes
  * `pbsnodes` list all nodes and their detailed state (free, busy/job-exclusive, offline)
* communicate with job:
  * `qdel <job id>` abort job
  * `qsig -s <signal number> <job id>` send signal or signal name to job
  * `qalter -lwalltime=12:00:00 <job id>` change the walltime of a job
  * `qalter -Wdepend=afterany:54321 12345` only start job `12345` after job with id `54321` has finished
  * `qhold <job id>` prevent the job from starting
  * `qrls <job id>` release the job to be eligible for run (after it was set on hold)

#### References

* [https://www.openpbs.org](https://www.openpbs.org)

### PJM

PJM (probably for *Parallel Job Manager*?) is a Fujitsu batch system
It is used at RIKEN Fugaku and on other Fujitsu systems.

#### NOTE
This section is a stub and improvements to complete the `(TODO)` sections are welcome.

#### Job Submission

* `pjsub your_job_script.pjsub`

#### Job Control

* interactive job:
  * `pjsub --interact`
* details for my jobs:
  * `pjstat` status of all jobs
  * (TODO) all details for job with <job id> `12345`
  * (TODO) all jobs under my user name
* details for queues:
  * (TODO) show all jobs in a queue
  * (TODO) compact view on free and busy nodes
  * (TODO) list all nodes and their detailed state (free, busy/job-exclusive, offline)
* communicate with job:
  * `pjdel <job id>` abort job
  * (TODO) send signal or signal name to job
  * (TODO) change the walltime of a job
  * (TODO) only start job `12345` after job with id `54321` has finished
  * `pjhold <job id>` prevent the job from starting
  * `pjrls <job id>` release the job to be eligible for run (after it was set on hold)

#### References

* [https://www.bsc.es/user-support/arm.php#ToC-runningjobs](https://www.bsc.es/user-support/arm.php#ToC-runningjobs)
* [https://www.cc.kyushu-u.ac.jp/scp/eng/system/ITO/02-2_batch.html](https://www.cc.kyushu-u.ac.jp/scp/eng/system/ITO/02-2_batch.html)
* [https://www.r-ccs.riken.jp/en/fugaku/user-guide/](https://www.r-ccs.riken.jp/en/fugaku/user-guide/)
