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
NPX=@NPX
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
cd "$TJ"
if [ $RUNFIRST -eq 0 ]; then
    file "$PWDX/restart.tar.gz" | grep 'gzip'
    if [ "$?" = "0" ]; then
        tar xzvf "$PWDX/restart.tar.gz"
    else
        tar xvf "$PWDX/restart.tar.gz"
    fi
fi

# -- run program
echo "!! $FNAME"
WTST=`date +%s`
### start:run.main
$ST < "$FNAME.in.$ORDER_ID" > "$FNAME.out.$ORDER_ID"
### end:run.main
wait
RESTIME=`expr \`date +%s\` '-' $WTST`
echo "TOTAL TIME IN SECONDS: $RESTIME" >> "$FNAME.out.$ORDER_ID"

FLIST="energy gradient vib_normal_modes vibspectrum final.xyz"
FLIST="$FLIST trajectory.xyz cell.txt bader.json trajectory.zip forces.xyz"

# --clean previous files
for i in $FLIST; do
    [ -f $PWDX/$i ] && rm $PWDX/$i
done

# -- copy files back
cp *.out.$ORDER_ID "$PWDX"
echo *.out.$ORDER_ID >> "$PWDX/FILELIST"
echo *.in.$ORDER_ID >> "$PWDX/FILELIST"
[ -L "$PWDX/restart.tar.gz" ] && rm "$PWDX/restart.tar.gz"
if [ -f *.chk/WAVECAR ]; then
    tar cvf "$PWDX/restart.tar.gz" *.chk
else
    tar czvf "$PWDX/restart.tar.gz" *.chk --exclude='*.chk/slave*.output'
fi
for i in $FLIST; do
    if [ -f *.chk/$i ]; then
        cp *.chk/$i "$PWDX/$i.$ORDER_ID"
        echo "$i.$ORDER_ID" >> "$PWDX/FILELIST"
    fi
done
rm -r *.chk
cd "$PWDX"
rm -r "$TJ"

# -- analyze the output
CNVED=0
if [ "$RUNTYPE" = "RELAX" ]; then
    VRF=`grep "Normal termination" "relax.out.$ORDER_ID" | wc -l`
    VRT=`grep "timeout error" "relax.out.$ORDER_ID" | wc -l`
    VRI=`grep "ired error" "relax.out.$ORDER_ID" | wc -l`
    VRN=`grep "OPTIMIZATION DID NOT CONVERGE WITHIN" "relax.out.$ORDER_ID" | wc -l`
    VRC=`grep "GEO_OPT_CONVERGED" "relax.out.$ORDER_ID" | wc -l`
    VRFAT=`grep "ERROR: Module dscf failed" "relax.out.$ORDER_ID" | wc -l`
    VRNC=`grep "ERROR: your energy calculation did not converge" "relax.out.$ORDER_ID" | wc -l`
    VRND=`grep "ATTENTION: dscf did not converge" "relax.out.$ORDER_ID" | wc -l`
    VRZ=`grep "Call to ZHEGV failed" "relax.out.$ORDER_ID" | wc -l`
    if [ $VRT != 0 ]; then
        echo 'FAILED REPEAT' > RESPONSE
    elif [ $VRFAT != 0 ]; then
        echo 'FAILED FATAL' > RESPONSE
    elif [ $VRNC != 0 ] || [ $VRND != 0 ]; then
        echo 'FAILED NOT_CONVERGED' > RESPONSE
    elif [ $VRI != 0 ]; then
        echo 'FAILED BAD_STRUCTURE' > RESPONSE
    elif [ $VRF != 0 ] && [ $VRN != 0 ] && [ $OPTNUM != 0 ]; then
        cp trajectory.xyz.$ORDER_ID trajectory.xyz
        echo 'SUCCESS STEP' > RESPONSE
    elif [ $VRF != 0 ] && [ $VRC != 0 ] && [ $OPTNUM != 0 ]; then
        CNVED=1
        cp trajectory.xyz.$ORDER_ID trajectory.xyz
        [ -f cell.txt.$ORDER_ID ] && cp cell.txt.$ORDER_ID cell.txt
        echo 'SUCCESS CONVERGED' > RESPONSE
    elif [ $VRF != 0 ]; then
        cp trajectory.xyz.$ORDER_ID trajectory.xyz
        [ -f forces.xyz.$ORDER_ID ] && cp forces.xyz.$ORDER_ID forces.xyz
        [ -f bader.json.$ORDER_ID ] && cp bader.json.$ORDER_ID bader.json
        echo 'SUCCESS ENERGY' > RESPONSE
    elif [ $VRZ != 0 ]; then
        echo 'FAILED BAD_STRUCTURE' > RESPONSE
    else
        echo 'FAILED UNKNOWN' > RESPONSE
    fi
elif [ "$RUNTYPE" = "NEB" ]; then
    VRF=`grep "Normal termination" "neb.out.$ORDER_ID" | wc -l`
    VRT=`grep "timeout error" "neb.out.$ORDER_ID" | wc -l`
    VRI=`grep "ired error" "neb.out.$ORDER_ID" | wc -l`
    VRN=`grep "NEB DID NOT CONVERGE WITHIN" "neb.out.$ORDER_ID" | wc -l`
    VRC=`grep "NEB_CONVERGED" "neb.out.$ORDER_ID" | wc -l`
    VRNC=`grep "ERROR: your energy calculation did not converge" "neb.out.$ORDER_ID" | wc -l`
    VRZ=`grep "Call to ZHEGV failed" "neb.out.$ORDER_ID" | wc -l`
    if [ $VRT != 0 ]; then
        echo 'FAILED REPEAT_NEB' > RESPONSE
    elif [ $VRNC != 0 ]; then
        echo 'FAILED NOT_CONVERGED_NEB' > RESPONSE
    elif [ $VRI != 0 ]; then
        echo 'FAILED BAD_STRUCTURE_NEB' > RESPONSE
    elif [ $VRF != 0 ] && [ $VRN != 0 ]; then
        cp final.xyz.$ORDER_ID final.xyz
        cp trajectory.zip.$ORDER_ID trajectory.zip
        echo 'SUCCESS STEP_NEB' > RESPONSE
    elif [ $VRF != 0 ] && [ $VRC != 0 ]; then
        CNVED=1
        cp final.xyz.$ORDER_ID final.xyz
        cp trajectory.zip.$ORDER_ID trajectory.zip
        echo 'SUCCESS CONVERGED_NEB' > RESPONSE
    elif [ $VRZ != 0 ]; then
        echo 'FAILED BAD_STRUCTURE_NEB' > RESPONSE
    else
        echo 'FAILED UNKNOWN_NEB' > RESPONSE
    fi
elif [ "$RUNTYPE" = "FREQ" ]; then
    VRQ=`grep "Normal termination" "freq.out.$ORDER_ID" | wc -l`
    VRT=`grep "timeout error" "freq.out.$ORDER_ID" | wc -l`
    VRFAT=`grep "ERROR: Module dscf failed" "freq.out.$ORDER_ID" | wc -l`
    VRI=`grep "ired error" "freq.out.$ORDER_ID" | wc -l`
    VRF=`grep "force ended abnormally" "freq.out.$ORDER_ID" | wc -l`
    if [ $VRT != 0 ]; then
        echo 'FAILED REPEAT_FREQ' > RESPONSE
    elif [ $VRFAT != 0 ]; then
        echo 'FAILED FATAL_FREQ' > RESPONSE
    elif [ $VRI != 0 ] || [ $VRF != 0 ]; then
        echo 'FAILED BAD_FREQ' > RESPONSE
    elif [ $VRQ != 0 ]; then
        CNVED=1
        cp vibspectrum.$ORDER_ID vibspectrum
        cp vib_normal_modes.$ORDER_ID vib_normal_modes
        echo 'SUCCESS FREQ' > RESPONSE
    else
        echo 'FAILED UNKNOWN_FREQ' > RESPONSE
    fi
fi
echo "$RESTIME" > RESPONSE-TIME
[ -f coord.xyz ] && mv coord.xyz coord.xyz.old.$ORDER_ID
mv REQUEST REQUEST.old.$ORDER_ID

echo "END AT: `date`" >> info.$ORDER_ID.log

[ -f "coord.xyz.old.$ORDER_ID" ] && echo "coord.xyz.old.$ORDER_ID" >> "$PWDX/FILELIST"
echo "REQUEST.old.$ORDER_ID" >> "$PWDX/FILELIST"
echo "info.$ORDER_ID.log" >> "$PWDX/FILELIST"

if [ "`expr $ORDER_ID '%' 10`" = "0" ] || [ $CNVED -eq 1 ]; then
    tar czvf archive.tar.gz.$ORDER_ID --files-from FILELIST --remove-files
    : > FILELIST
fi
