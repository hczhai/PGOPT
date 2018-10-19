#!/bin/bash

# addzero 1 3 ==> 003
function addzero {
  ix="$1"
  while [ ${#ix} -lt $2 ]; do
    ix="0$ix"
  done
  echo $ix
}

X=$2
RUNTYPE=$1
NUML=3
export XBLOCK=$3
export XCORNER=$4
export XSHAPE=$5

PWX=`pwd`
PROCX="PROC`addzero $X $NUML`"
cd ../../procs/$RUNTYPE/$PROCX
export JOBID=$JOBID.$X

if [ "$XBLOCK" = "" ]; then
  $PGOPTHOME/scripts/timestamp.sh &
else
  sleep `expr $RANDOM % 3`.`expr $RANDOM % 10`
  echo 'BLOCK' > BLOCK_RUNNING
fi

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
  echo "#$X:RUN JOB $II ..."
  ./run-all.sh
  echo "#$X:END JOB $II ..."
  II=`expr $II + 1`
done
