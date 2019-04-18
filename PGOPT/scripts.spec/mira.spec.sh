
### start:deriv.block
#COBALT -t @TIMEMIN -n @NPROCSTOT --disable_preboot -O @NAME
### end:deriv.block

### start:prelude.block
export JOBID=${COBALT_JOBID}
python --version
### end:prelude.block

### start:block.main
sed -i "1 i\\% mirablocks=${XBLOCK},${XCORNER},${XSHAPE}" "$FNAME.in.$ORDER_ID"
### end:block.main

### start:run.main
$ST --mira < "$FNAME.in.$ORDER_ID" > "$FNAME.out.$ORDER_ID"
### end:run.main
