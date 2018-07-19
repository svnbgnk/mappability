#!/bin/bash

# 0,36 1,24 2,36 2,50 3,75 4,100
#for i in 0,36,18 0,36,24 1,24,4 1,24,2 2,36,2 2,50,4 2,50,10; do
#for i in 0,30 0,40 0,50 0,60 0,70 0,80 0,90 0,95 1,5 1,10 1,15 1,20 1,25 1,30 1,35 1,40 1,45 1,50 2,0 2,5 2,10 2,15; do
#for i in 3,3 3,5 3,10 3,15 4,3 4,5 4,10 4,15; do
E=4
for K in 24 36 50 75 100 150 200 ; do
    for OPERC in 0.05 0.1 0.15 0.2 0.25 0.3 0.5 0.75 ; do
        #IFS=',' read O K <<< "${i}"
        Ofloat=$(echo $K*$OPERC | bc)
        O=${Ofloat%.*}
        echo "### K = ${K}, O = ${OPERC} (${O}) ###"
        /usr/bin/time -f "%E, %e, %P" ./mappability -I /srv/public/cpockrandt/chr13-index/index -O /dev/shm/chr13_wo_k -E ${E} -K ${K} -o ${O} -m -t 20
        echo "# KNUT #"
        /usr/bin/time -f "%E, %e, %P" ./mappability -I /srv/public/cpockrandt/chr13-index/index -O /dev/shm/chr13_wi_k -E ${E} -K ${K} -o ${O} -m -t 20 --knut
    done
done

