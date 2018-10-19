#!/bin/bash

### ::deriv.gpu

if [ "$PG_BATCH_MODE" != "ON" ]; then
    ### ::prelude.gpu
fi

OUTORDER=`echo net.out.* | tr ' ' '\n' | grep '\*$' -v | wc -l`
echo "## `hostname -f | grep '\\.' || hostname -a`" > net.out.$OUTORDER
### start:run.gpu
python -u $ACNNHOME/main.py net.json > net.out.$OUTORDER 2>&1
### end:run.gpu
[ $? -ne 0 ] && exit 1
