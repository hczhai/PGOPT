#!/bin/bash

while true; do
    ${PGOPTHOME:+$PGOPTHOME/}pgopt check --runt=$1 --maxt=600 --delete
    sleep 600
done

