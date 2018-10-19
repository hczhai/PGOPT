#!/bin/bash

ccmrun ./worker-main.pre.sh @NPX
[ "$?" != 0 ] && exit 4

PWDX=`pwd`
TJ=`tail -n 1 tmpdir.txt`
ORDER_ID=`tail -n 1 orderid.txt`
FNAME=`tail -n 1 fname.txt`

cd $TJ
ccmrun $STMOLE_HOME/SVASP --aprun --pre < "$FNAME.in.$ORDER_ID"
cd `cat ./RUNDIR`
aprun -n @NPX `cat ../PROGRAM` > "../$FNAME.out.$ORDER_ID"
if [ "$?" != 0 ]; then
    echo "Error termination" >> "../$FNAME.out.$ORDER_ID"
else
    echo "Normal termination" >> "../$FNAME.out.$ORDER_ID"
fi
cd ..
ccmrun $STMOLE_HOME/SVASP --aprun --post < "$FNAME.in.$ORDER_ID" >> "$FNAME.out.$ORDER_ID"

cd $PWDX
ccmrun ./worker-main.post.sh
