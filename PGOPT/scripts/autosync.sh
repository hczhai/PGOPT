#!/bin/bash

while true; do
    ${PGOPTHOME:+$PGOPTHOME/}pgopt sync
    sleep 120
done
