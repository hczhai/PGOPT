#!/bin/bash

### ::deriv.worker.seq

### ::prelude.worker.seq

# addzero 1 3 ==> 003
function addzero {
  ix="$1"
  while [ ${#ix} -lt $2 ]; do
    ix="0$ix"
  done
  echo $ix
}

### start:assign.worker.seq
X=@X
NN=`expr @NPROCS '*' @NNODES`
RUNTYPE=@RUNT

XX=`expr $X '*' $NN`
XY=`expr $XX '+' $NN '-' 1`
### end:assign.worker.seq

### start:run.worker.seq
$PGOPTHOME/scripts/torun-parallel.sh $RUNTYPE $XX $XY
### end:run.worker.seq