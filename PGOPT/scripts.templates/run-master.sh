#!/bin/bash

### ::deriv.master

RUNTYPE=@RUNT

if [ "$PG_BATCH_MODE" != "ON" ]; then
    ### ::prelude.master
    if [ -f ../../DIRECTORIES ]; then
        $PGOPTHOME/scripts/autosync.sh &
        $PGOPTHOME/scripts/autocheck.sh $RUNTYPE &
    fi
fi

OUTORDER=`echo master.out.* | tr ' ' '\n' | grep '\*$' -v | wc -l`
echo "## `hostname -f | grep '\\.' || hostname -a`" > master.out.$OUTORDER
### start:run.master
python -u $ACNNHOME/main.py para.json >> master.out.$OUTORDER 2>&1
### end:run.master
[ $? -ne 0 ] && exit 1

if [ "$PG_BATCH_MODE" != "ON" ]; then
  cd ../../procs/$RUNTYPE
  for i in *; do
    cd $i
    echo 'FINISH' > REQUEST
    cd ..
  done
  cd ../..
  [ -f ./DIRECTORIES ] && ${PGOPTHOME:+$PGOPTHOME/}pgopt sync
fi