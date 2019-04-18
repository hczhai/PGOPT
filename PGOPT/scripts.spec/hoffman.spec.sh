
### start:deriv.master
#$ -l h_data=4G,h_rt=@TIME
#$ -pe shared 1
#$ -cwd
#$ -o LOG.$JOB_ID
#$ -e ERR.$JOB_ID
#$ -N @NAME
### end:deriv.master

### start:prelude.master
source ~/.bashrc
python --version
### end:prelude.master

### start:deriv.worker.para
#$ -l h_data=4G,h_rt=@TIME
#$ -pe dc* @NPROCS
#$ -cwd
#$ -o LOG.$JOB_ID
#$ -e ERR.$JOB_ID
### end:deriv.worker.para

### start:prelude.worker.para
export JOBID=${JOB_ID}
source ~/.bashrc
python --version
### end:prelude.worker.para

### start:deriv.worker.seq
#$ -l h_data=4G,h_rt=@TIME
#$ -t @LOWER-@UPPER:1
#$ -pe shared 1
#$ -cwd
#$ -o LOG.$JOB_ID
#$ -e ERR.$JOB_ID
### end:deriv.worker.seq

### ::prelude.worker.seq = prelude.worker.para

### start:assign.worker.seq
X=`expr $SGE_TASK_ID '-' 1`
RUNTYPE=@RUNT
### end:assign.worker.seq

### start:run.worker.seq
$PGOPTHOME/scripts/torun-single.sh $RUNTYPE $X
### end:run.worker.seq

### start:deriv.worker.x
#$ -l h_data=4G,h_rt=@TIME
#$ -t @LOWER-@UPPER:1
#$ -pe dc* @NPROCS
#$ -cwd
#$ -o LOG.$JOB_ID
#$ -e ERR.$JOB_ID
### end:deriv.worker.x

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
INII=${MYNAME:0:1}
SCRATCH_DIR=/u/scratch/$INII/$MYNAME
echo "SCRATCH: $SCRATCH_DIR"
### end:scratch.main

### start:hosts.main
cat $PE_HOSTFILE | awk "{print \$1\":\"int(\$2/$NTHREADS)}" | sort -t : -k 2nr | grep -v ":0$" > hosts.txt.$ORDER_ID
NP=`cat hosts.txt.$ORDER_ID | awk -F : '{s+=$2} END {printf "%.0f\n", s}'`
NPX=`expr $NP '*' $NTHREADS`
sed -i "s/:1\$//g" hosts.txt.$ORDER_ID
NPOPT=`cat hosts.txt.$ORDER_ID | tr '\n' ',' | sed 's/,$//g'`
### end:hosts.main

### start:vasp.prelude.main
source /u/local/Modules/default/init/modules.sh
# module load intel/16.0.2
module load intel/17.0.1
module load intelmpi/5.0.0
### end:vasp.prelude.main
