#!/bin/bash

cat << EOF > rv.in
name @GNAME
program      ADF
charge       0
composition   @GNAMEX
reoptimize    0
matings       0
3D @GNUMBER
energy_resolution -1
repulsion_resolution -1
energy_window -1
EOF

$CLUSTER_KANTERS_HOME/cluster.py rv.in > rv.out 2> /dev/null
mv Pt6_all.xyz rv_structs.xyz