#!/bin/bash
### ::deriv.block

RUNTYPE=@RUNT

### ::prelude.block

OUTORDER=`echo master.out.* | tr ' ' '\n' | grep '\*$' -v | wc -l`
echo "## `hostname -f | grep '\\.' || hostname -a`" > master.out.$OUTORDER

XMPIPR=@MPIPR
XBLKPR=@BLKPR
if [ "$XMPIPR" = "16" ]; then
    XSHAPE="1x1x1x1x1"
elif [ "$XMPIPR" = "32" ]; then
    XSHAPE="2x1x1x1x1"
elif [ "$XMPIPR" = "64" ]; then
    XSHAPE="2x2x1x1x1"
elif [ "$XMPIPR" = "128" ]; then
    XSHAPE="2x2x2x1x1"
elif [ "$XMPIPR" = "256" ]; then
    XSHAPE="2x2x2x2x1"
elif [ "$XMPIPR" = "512" ]; then
    XSHAPE="2x2x2x2x2"
elif [ "$XMPIPR" = "1024" ]; then
    XSHAPE="4x2x2x2x2"
elif [ "$XMPIPR" = "2048" ]; then
    XSHAPE="4x4x2x2x2"
elif [ "$XMPIPR" = "4096" ]; then
    XSHAPE="4x4x4x2x2"
elif [ "$XMPIPR" = "8192" ]; then
    XSHAPE="4x4x4x4x2"
else
    echo "SHAPE cannot be calculated for $XMPIPR !"
    exit 3
fi

BLOCKS=`get-bootable-blocks --size 512 $COBALT_PARTNAME`
for BLOCK in $BLOCKS; do
    boot-block --block $BLOCK &
done; wait
IBX=0
for BLOCK in $BLOCKS; do
    CORS=`/soft/cobalt/bgq_hardware_mapper/get-corners.py $BLOCK $XSHAPE`
    for COR in $CORS; do
        $PGOPTHOME/scripts/torun-single.sh $RUNTYPE $IBX $BLOCK $COR $XSHAPE &
        IBX=`expr $IBX + 1`
    done
    if [ $IBX -ge $XBLKPR ]; then
        break
    fi
done

### start:run.block
python -u $ACNNHOME/main.py para.json >> master.out.$OUTORDER 2>&1
### end:run.block

cd ../../procs/$RUNTYPE
for i in *; do
    cd $i
    echo 'FINISH' > REQUEST
    cd ..
done
cd ../..

[ -f ./DIRECTORIES ] && ${PGOPTHOME:+$PGOPTHOME/}pgopt sync

for BLOCK in $BLOCKS; do
    boot-block --block $BLOCK --free &
done; wait
