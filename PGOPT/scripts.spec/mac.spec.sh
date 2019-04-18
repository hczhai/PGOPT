
### start:deriv.master
#$ --time=@TIME --name=@NAME --procs=1
### end:deriv.master

### start:prelude.master
source ~/.bashrc
python --version
### end:prelude.master

### start:deriv.worker.para
#$ --time=@TIME --name=@NAME --procs=@NPROCS
### end:deriv.worker.para

### start:prelude.worker.para
export JOBID=12345
source ~/.bashrc
python --version
### end:prelude.worker.para

### start:deriv.worker.seq
#$ -t @LOWER-@UPPER:1 --time=@TIME --name=@NAME --procs=@NPROCS
### end:deriv.worker.seq

### ::prelude.worker.seq = prelude.worker.para

### start:assign.worker.seq
X=`expr $SGE_TASK_ID '-' 1`
RUNTYPE=@RUNT
### end:assign.worker.seq

### start:run.worker.seq
$PGOPTHOME/scripts/torun-single.sh $RUNTYPE $X
### end:run.worker.seq

### ::deriv.worker.x = deriv.worker.seq

### start:prelude.worker.x
export JOBID=${JOB_ID}.${SGE_TASK_ID}
source ~/.bashrc
python --version
### end:prelude.worker.x

### start:assign.worker.x
X=`expr $SGE_TASK_ID '-' 1`
RUNTYPE=@RUNT
NUML=3
PWX=`pwd`
PROCX="PROC`addzero $X $NUML`"
cd ../../procs/$RUNTYPE/$PROCX
### end:assign.worker.x

### ::run.worker.x.seq-like

### start:scratch.main
SCRATCH_DIR=$WORKDIR
echo "SCRATCH: $SCRATCH_DIR"
### end:scratch.main

### start:hosts.main
NPX=2
ulimit -s hard
### end:hosts.main

### start:vasp.prelude.main
### end:vasp.prelude.main
