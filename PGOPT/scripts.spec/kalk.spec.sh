
### start:deriv.master
#PBS -l nodes=1:ppn=1:xeon
#PBS -l walltime=@TIME
#PBS -l mem=2500mb
#PBS -N @NAME
#PBS -j oe
### end:deriv.master

### start:prelude.master
export JOBID=`echo ${PBS_JOBID} | cut -f1 -d'.'`
cd $PBS_O_WORKDIR
source /software/intel/compilers_and_libraries_2017/linux/bin/compilervars.sh intel64
source ~/.bashrc
python --version
export WORKDIR=$SCRATCH
### end:prelude.master

### start:run.master
python -u $ACNNHOME/main.py para.json >> master.out.$OUTORDER 2>&1
### end:run.master

### start:deriv.gpu
#PBS -l nodes=@NNODES:ppn=8:xeon
#PBS -l walltime=@TIME
#PBS -l mem=7520mb
#PBS -j oe
### end:deriv.gpu

### start:prelude.gpu
export JOBID=`echo ${PBS_JOBID} | cut -f1 -d'.'`
cd $PBS_O_WORKDIR
source /software/intel/compilers_and_libraries_2017/linux/bin/compilervars.sh intel64
source ~/.bashrc
python --version
### end:prelude.gpu

### start:run.gpu
python -u $ACNNHOME/main.py net.json > net.out.$OUTORDER 2>&1
### end:run.gpu

### start:deriv.worker.para
#PBS -l nodes=@NNODES:ppn=8:xeon
#PBS -l walltime=@TIME
#PBS -l mem=7520mb
#PBS -N @NAME
#PBS -j oe
### end:deriv.worker.para

### ::prelude.worker.para = prelude.master
### ::deriv.worker.seq = deriv.master
### ::prelude.worker.seq = prelude.master

### start:run.worker.seq
$PGOPTHOME/scripts/torun-parallel.sh $RUNTYPE $XX $XY
### end:run.worker.seq

### start:vasp.prelude.main
### end:vasp.prelude.main

### start:run.main
$ST < "$FNAME.in.$ORDER_ID" > "$FNAME.out.$ORDER_ID"
### end:run.main