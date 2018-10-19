#!/bin/bash

### ::deriv.worker.x

### ::prelude.worker.x

# addzero 1 3 ==> 003
function addzero {
    ix="$1"
    while [ ${#ix} -lt $2 ]; do
        ix="0$ix"
    done
    echo $ix
}

### start:assign.worker.x
RUNTYPE=@RUNT
XX=`expr @LOWER - 1`
XY=`expr @UPPER - 1`
### end:assign.worker.x

### start:run.worker.x.seq-like
$PGOPTHOME/scripts/torun-parallel.sh $RUNTYPE $XX $XY
### end:run.worker.x.seq-like

### start:run.worker.x.para-like
$PGOPTHOME/scripts/timestamp.sh &

II=0
while true; do
    while [ ! -f REQUEST ]; do
        sleep 5
    done
    if [ "`cat REQUEST`" = "FINISH" ]; then
        mv "REQUEST" "REQUEST.final"
        echo "FINISHED ..."
        break
    fi
    echo "RUN JOB $II ..."
    ./run-all.sh
    echo "END JOB $II ..."
    II=`expr $II + 1`
    if [ "@ONCE" = "TRUE" ]; then
        break
    fi
done
### end:run.worker.x.para-like