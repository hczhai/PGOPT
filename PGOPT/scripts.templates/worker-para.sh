#!/bin/bash

### ::deriv.worker.para

### ::prelude.worker.para

# addzero 1 3 ==> 003
function addzero {
    ix="$1"
    while [ ${#ix} -lt $2 ]; do
        ix="0$ix"
    done
    echo $ix
}

X=@X
RUNTYPE=@RUNT
NUML=3

PWX=`pwd`
PROCX="PROC`addzero $X $NUML`"
cd ../../procs/$RUNTYPE/$PROCX

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
