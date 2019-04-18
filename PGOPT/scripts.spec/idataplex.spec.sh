
### start:deriv.master
#PBS -A @PROJ
#PBS -l select=@NNODES:ncpus=@NPROCS:mpiprocs=@NPROCS
#PBS -q @QUEUE
#PBS -l walltime=@TIME
#PBS -N @NAME
#PBS -j oe
### end:deriv.master

### start:prelude.master
export JOBID=`echo ${PBS_JOBID} | cut -f1 -d'.'`
cd $PBS_O_WORKDIR
source ~/.bashrc
export MPI_REMSH=$I_MPI_MPD_RSH
python --version
### end:prelude.master

### start:deriv.gpu
#PBS -A @PROJ
#PBS -l select=1:ncpus=28:mpiprocs=28:ngpus=1
#PBS -q standard
#PBS -l walltime=@TIME
#PBS -N GBATCH
#PBS -j oe
### end:deriv.gpu

### start:prelude.gpu
export JOBID=`echo ${PBS_JOBID} | cut -f1 -d'.'`
cd $PBS_O_WORKDIR
. $MODULESHOME/init/bash 2> /dev/null
module load cuda/6.5
module swap mpi/intelmpi mpi/intelmpi/16.0.0
module load compiler/intel/16.0.0
source ~/.bashrc
python --version
### end:prelude.gpu

### ::deriv.worker.para = deriv.master
### ::prelude.worker.para = prelude.master
### ::deriv.worker.seq = deriv.master
### ::prelude.worker.seq = prelude.master

### start:hostname.main
HOSTXN=`hostname -f | grep '\.' || hostname -a`
echo "$HOSTXN" | grep haise && HOSTNN=haise
echo "$HOSTXN" | grep "ice-x.erdc" && HOSTNN=topaz
echo "$HOSTXN" | grep "icex.afrl" && HOSTNN=spirit
echo "$HOSTXN" | grep "ice-x.afrl" && HOSTNN=thunder
echo "$HOSTXN" | grep "arlu.arl" && HOSTNN=centennial
### end:hostname.main

### start:turbomole.prelude.main
if [ "$HOSTNN" = "thunder" ]; then
    source /app/.bashrc
    . $MODULESHOME/init/bash
    module unload mpt
    module load intel-mpi-15/5.0.3.048
    module load intel-compilers/15.3.187
fi
### end:turbomole.prelude.main

### start:vasp.prelude.main
if [ "$HOSTNN" = "topaz" ]; then
    module unload compiler/intel
    module unload mpi/sgimpt
    module load compiler/intel/16.0.0
    module load mpi/intelmpi/16.0.0
elif [ "$HOSTNN" = "haise" ]; then
    . $MODULESHOME/init/bash # 2> /dev/null
    module swap mpi/intel/impi mpi/intel/impi/5.0.2
    module load compiler/intel/15.1
    module load mkl/15.1
    export MKLROOT=${MKLPATH%:*}/../..
    export LD_LIBRARY_PATH=$MKLROOT/lib/intel64:$LD_LIBRARY_PATH
elif [ "$HOSTNN" = "thunder" ]; then
    source /app/.bashrc
    . $MODULESHOME/init/bash
    module swap mpt intel-mpi-16/5.1.3.210
    module unload intel-compilers
    module load intel-compilers/16.0.4
elif [ "$HOSTNN" = "centennial" ]; then
    . $MODULESHOME/init/bash 2> /dev/null
    module unload compiler/intel
    module unload mpi/sgimpt
    module load compiler/intel/16.0.4.258
    module load mpi/intelmpi/2017.1.132
fi
### end:vasp.prelude.main
