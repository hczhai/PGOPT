#!/bin/bash

### ::deriv.master

### ::prelude.master

PROCTYPE=@PROT
RUNTYPE=@RUNT
RUNSTAGE=@RUNS

export PG_BATCH_MODE=ON

if [ -f ../../DIRECTORIES ]; then
  $PGOPTHOME/scripts/autosync.sh &
  $PGOPTHOME/scripts/autocheck.sh $PROCTYPE &
  STED=${PWD##*/}
  cd `cat ../../DIRECTORIES | tail -n 1`/tomaster/$STED
fi

cd ..
for i in @MULTS; do
  cd $RUNTYPE-$RUNSTAGE.$i
  ./run-master.sh
  [ $? -ne 0 ] && exit 1
  cd ..
done

cd ../procs/$PROCTYPE
for i in *; do
  cd $i
  echo 'FINISH' > REQUEST
  cd ..
done
cd ../..
[ -f ./DIRECTORIES ] && ${PGOPTHOME:+$PGOPTHOME/}pgopt sync