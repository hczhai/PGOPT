
### start:deriv.master
#MSUB -A @PROJ
#MSUB -l nodes=@NNODES:ppn=@NPROCS
#MSUB -l walltime=@TIME
#MSUB -N @NAME
#MSUB -o LOG.%j
#MSUB -j oe
#MSUB -V
### end:deriv.master

### start:prelude.master
export JOBID=${SLURM_JOBID}
source /etc/profile.d/modules.bash
module purge
module load pnnl_env
source ~/.bashrc
python --version
### end:prelude.master

### ::deriv.worker.para = deriv.master
### ::prelude.worker.para = prelude.master
### ::deriv.worker.seq = deriv.master
### ::prelude.worker.seq = prelude.master
### ::deriv.worker.x = deriv.master
### ::prelude.worker.x = prelude.master
### ::run.worker.x.para-like

### start:vasp.prelude.main
module load intel/16.1.150
module load impi/5.1.2.150
### end:vasp.prelude.main

### start:run.main
$ST --srun < "$FNAME.in.$ORDER_ID" > "$FNAME.out.$ORDER_ID"
### end:run.main