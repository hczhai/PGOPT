#!/bin/bash

RUNTYPE=$1
XX=$2
XY=$3
LISTX=($(seq $XX $XY))

for l in ${LISTX[@]}; do
  $PGOPTHOME/scripts/torun-single.sh $RUNTYPE $l &
done

wait
