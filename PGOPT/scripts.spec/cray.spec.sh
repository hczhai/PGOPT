
### start:deriv.master
#PBS -A @PROJ
#PBS -l select=@NNODES:ncpus=@NPROCS:mpiprocs=@NPROCS
#PBS -q @QUEUE
#PBS -l walltime=@TIME
#PBS -l ccm=1
#PBS -N @NAME
#PBS -j oe
### end:deriv.master

### start:prelude.master
export JOBID=`echo ${PBS_JOBID} | cut -f1 -d'.'`
cd $PBS_O_WORKDIR
. $MODULESHOME/init/bash
module load ccm
if [ "$BC_HOST" = "lightning" ]; then
    module swap PrgEnv-cray PrgEnv-intel/5.2.40
elif [ "$BC_HOST" = "garnet" ]; then
    module swap PrgEnv-pgi PrgEnv-intel/5.2.82
elif [ "$BC_HOST" = "conrad" ] || [ "$BC_HOST" = "gordon" ] || \
    [ "$BC_HOST" = "shepard" ] || [ "$BC_HOST" = "armstrong" ]; then
    module swap PrgEnv-cray PrgEnv-intel/5.2.40
    module swap intel intel/16.0.2.181
elif [ "$BC_HOST" = "excalibur" ]; then
    . $MODULESHOME/init/bash
    module unload cray-mpich
    module unload PrgEnv-intel
    module unload intel
    module unload cray-libsci
    module load cray-mpich/7.2.4
    module load intel/16.0.0.109
    module load cray-libsci/13.1.0
    module load PrgEnv-intel/5.2.82
elif [ "$BC_HOST" = "onyx" ]; then
    module unload PrgEnv-cray
    module load intel/16.0.2.181
    module load PrgEnv-intel/6.0.4
fi
source ~/.bashrc
python --version
### end:prelude.master

### start:run.master
ccmrun python -u $ACNNHOME/main.py para.json >> master.out.$OUTORDER 2>&1
### end:run.master

### start:deriv.gpu
#PBS -A @PROJ
#PBS -l select=1:ncpus=@NPROCS:accelerator_model=Tesla_K40s
#PBS -q gpu
#PBS -l walltime=@TIME
#PBS -N GBATCH
#PBS -j oe
### end:deriv.gpu

### start:prelude.gpu
export JOBID=`echo ${PBS_JOBID} | cut -f1 -d'.'`
cd $PBS_O_WORKDIR
. $MODULESHOME/init/bash
module load cudatoolkit
source ~/.bashrc
python --version
### end:prelude.gpu

### start:run.gpu
aprun -n 1 python -u $ACNNHOME/main.py net.json > net.out.$OUTORDER 2>&1
### end:run.gpu

### ::deriv.worker.para = deriv.master
### ::prelude.worker.para = prelude.master
### ::deriv.worker.seq = deriv.master
### ::prelude.worker.seq = prelude.master

### start:run.worker.seq
ccmrun $PGOPTHOME/scripts/torun-parallel.sh $RUNTYPE $XX $XY
### end:run.worker.seq

### start:vasp.prelude.main
if [ "$BC_HOST" = "conrad" ] || [ "$BC_HOST" = "gordon" ] || [ "$BC_HOST" = "shepard" ]; then
    LDBAK=$LD_LIBRARY_PATH
    LD_LIBRARY_PATH=
    module swap PrgEnv-cray PrgEnv-intel/5.2.40
    module swap intel intel/16.0.2.181
    LD_LIBRARY_PATH=$LDBAK
elif [ "$BC_HOST" = "armstrong" ]; then
    module swap PrgEnv-cray PrgEnv-intel/5.2.40
    module swap intel intel/16.0.2.181
elif [ "$BC_HOST" = "lightning" ]; then
    module swap PrgEnv-cray PrgEnv-intel/5.2.40
    module swap intel intel/13.1.2.183
elif [ "$BC_HOST" = "garnet" ]; then
    module swap PrgEnv-pgi PrgEnv-intel/5.2.82
    module swap intel intel/14.0.2.144
elif [ "$BC_HOST" = "excalibur" ]; then
    . $MODULESHOME/init/bash
    module unload cray-mpich
    module unload PrgEnv-intel
    module unload intel
    module unload cray-libsci
    module load cray-mpich/7.2.4
    module load intel/16.0.0.109
    module load cray-libsci/13.1.0
    module load PrgEnv-intel/5.2.82
elif [ "$BC_HOST" = "onyx" ]; then
    module unload PrgEnv-cray
    module load intel/16.0.2.181
    module load PrgEnv-intel/6.0.4
fi
### end:vasp.prelude.main

### start:run.main
if [ "$APROG" = "VASP" ]; then
    $ST --aprun < "$FNAME.in.$ORDER_ID" > "$FNAME.out.$ORDER_ID"
else
    if [ "$FNAME" != "freq" ]; then
        ccmrun $ST < "$FNAME.in.$ORDER_ID" > "$FNAME.out.$ORDER_ID"
    else
        $ST < "$FNAME.in.$ORDER_ID" > "$FNAME.out.$ORDER_ID"
    fi
fi
### end:run.main