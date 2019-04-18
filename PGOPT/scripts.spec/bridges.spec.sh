
### start:deriv.master
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t @TIME
#SBATCH -J @NAME
#SBATCH --ntasks-per-node 1
#SBATCH -o LOG.%j
### end:deriv.master

### start:deriv.worker.para
#SBATCH -N @NNODES
#SBATCH -p RM-shared
#SBATCH -t @TIME
#SBATCH -J @NAME
#SBATCH --ntasks-per-node @NPROCS
#SBATCH -o LOG.%j
### end:deriv.worker.para

### start:prelude.master
export JOBID=`echo ${SLURM_JOB_ID} | cut -f1 -d'.'`
source ~/.bashrc
python --version
export WORKDIR=$SCRATCH
### end:prelude.master

### start:run.master
python -u $ACNNHOME/main.py para.json >> master.out.$OUTORDER 2>&1
### end:run.master

### start:prelude.worker.para
export JOBID=`echo ${SLURM_JOB_ID} | cut -f1 -d'.'`
if [ "$1" != "" ]; then
    export JOBID=${JOBID}.${1}
    export XBLOCK=${2}
    export XCORNER=${3}
fi
source ~/.bashrc
python --version
export WORKDIR=$SCRATCH
### end:prelude.worker.para

### ::deriv.worker.seq = deriv.worker.para
### ::prelude.worker.seq = prelude.master

### start:block.main
sed -i "1 i\\% mirablocks=${XBLOCK},${XCORNER}," "$FNAME.in.$ORDER_ID"
### end:block.main

### start:run.worker.seq
$PGOPTHOME/scripts/torun-parallel.sh $RUNTYPE $XX $XY
### end:run.worker.seq

### start:vasp.prelude.main
### end:vasp.prelude.main

### start:run.main
if [ "$XBLOCK" != "" ]; then
    $ST --mpirunst < "$FNAME.in.$ORDER_ID" > "$FNAME.out.$ORDER_ID"
else
    $ST < "$FNAME.in.$ORDER_ID" > "$FNAME.out.$ORDER_ID"
fi
### end:run.main