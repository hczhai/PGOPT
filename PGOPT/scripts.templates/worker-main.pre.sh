#!/bin/bash

if [ ! -f REQUEST ]; then
    exit 4
elif [ "`cat REQUEST | cut -d ' ' -f 1`" = "SLEEP" ]; then
    mv REQUEST REQUEST.slp
    touch SLEEPING
    exit 10
elif [ ! -f REQUEST-INDEX ]; then
    mv REQUEST REQUEST.nid
    exit 4
elif [ "`cat REQUEST | cut -d ' ' -f 2`" = "INIT" ]; then
    RUNFIRST=1
else
    RUNFIRST=0
fi

# -- analyze request
RUNTYPE=`cat REQUEST | cut -d ' ' -f 1`
OPTNUM=`cat REQUEST | cut -d ' ' -f 3`
AMETHOD=`cat REQUEST | cut -d ' ' -f 4`
AMULTI=`cat REQUEST | cut -d ' ' -f 5`
ACHARGE=`cat REQUEST | cut -d ' ' -f 6`
APROG=`cat REQUEST | cut -d ' ' -f 7`

NTHREADS=1
NPX=$1
NPOPT="unknown:$NPX"

case $RUNTYPE in
    RELAX)
        FNAME=relax;;
    FREQ)
        FNAME=freq;;
    NEB)
        FNAME=neb;;
esac

if [ ! -f orderid.txt ]; then
    ORDER_ID=1
else
    LO=`cat orderid.txt | wc -l`
    ORDER_ID=`expr $LO '+' 1`
fi

if [ -f tmpdir.txt ]; then
    OTJ=`cat tmpdir.txt | tail -1`
    OOID=`cat orderid.txt | tail -1`
else
    OTJ=""
    OOID=0
fi

HN=`hostname`
MYNAME=`whoami`
echo "MYNAME: $MYNAME"

### start:scratch.main
SCRATCH_DIR=$WORKDIR
echo "SCRATCH: $SCRATCH_DIR"
### end:scratch.main

### ::hosts.main
export OMP_NUM_THREADS=$NTHREADS

TMPORDER=`echo $SCRATCH_DIR/$JOBID.* | tr ' ' '\n' | grep '\*$' -v | wc -l`
export TMPORDER=`expr $TMPORDER '+' 1`
TJ="$SCRATCH_DIR/$JOBID.$TMPORDER"
[ ! -x "$TJ" ] && mkdir "$TJ"

echo "$TJ" >> tmpdir.txt
echo $ORDER_ID >> orderid.txt
echo $FNAME > fname.txt

### ::hostname.main

# -- program preludes
if [ "$APROG" = "TURBOMOLE" ]; then
    ST=$STMOLE_HOME/STMOLE
    if [ "$TM_PRE" != "OK" ]; then
        ### ::turbomole.prelude.main
        export TM_PRE=OK
    fi
elif [ "$APROG" = "VASP" ]; then
    ST=$STMOLE_HOME/SVASP
    if [ "$VASP_PRE" != "OK" ]; then
        ### ::vasp.prelude.main
        export VASP_PRE=OK
    fi
fi

echo "START AT: `date`" >> info.$ORDER_ID.log

# -- prepare input file
cp "template.in" "$FNAME.in.$ORDER_ID"
GPP=`grep "^#" "$FNAME.in.$ORDER_ID"`
GPP="$GPP coords=coord.xyz"
[ "$RUNFIRST" != "1" ] && GPP="$GPP continue"
[ "$FNAME" = "relax" ] && GPP="$GPP opt(iter=$OPTNUM)"
[ "$FNAME" = "freq" ] && GPP="$GPP freq"
[ "$FNAME" = "neb" ] && GPP="$GPP neb(iter=$OPTNUM)"
sed -i "1 i\\% workers=$NPOPT" "$FNAME.in.$ORDER_ID"
### ::block.main
sed -i "/^#/c\\$GPP" "$FNAME.in.$ORDER_ID"
sed -i "s|@LLTT|${AMETHOD//;/ }|g" "$FNAME.in.$ORDER_ID"
sed -i "s|@MULT|${AMULTI}|g" "$FNAME.in.$ORDER_ID"
sed -i "s|@CHAR|${ACHARGE}|g" "$FNAME.in.$ORDER_ID"

# -- copy files to temparary directory
mv REQUEST-INDEX REQUEST-INDEX-AC
PWDX=`pwd`
cp -r *.$ORDER_ID "$TJ"
[ -f coord.xyz ] && cp coord.xyz "$TJ"

WTST=`date +%s`
echo "$WTST" > STARTTIME

cd "$TJ"
if [ $RUNFIRST -eq 0 ]; then
    file "$PWDX/restart.tar.gz" | grep 'gzip'
    if [ "$?" = "0" ]; then
        tar xzvf "$PWDX/restart.tar.gz"
    else
        tar xvf "$PWDX/restart.tar.gz"
    fi
fi
