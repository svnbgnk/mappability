#!/bin/bash

E=$1

for K in 24 36 50 75 100 150 200 ; do
    for OPERC in 0.05 0.1 0.15 0.2 0.25 0.3 0.5 0.75 ; do
        Ofloat=$(echo $K*$OPERC | bc)
        O=${Ofloat%.*}
        echo "### K = ${K}, O = ${OPERC} (${O}) ###"
        /usr/bin/time -f "%E, %e, %P" ~/Dev/mappability-build/release/mappability -I /srv/public/cpockrandt/chr13-index/index -O /srv/public/cpockrandt/gmapp_results/chr13 -E ${E} -K ${K} -o ${O} -m -t 20
        echo "# KNUT #"
        /usr/bin/time -f "%E, %e, %P" ~/Dev/mappability-build/release/mappability -I /srv/public/cpockrandt/chr13-index/index -O /srv/public/cpockrandt/gmapp_results/chr13 -E ${E} -K ${K} -o ${O} -m -t 20 --knut
        echo "# KNUT2 #"
        /usr/bin/time -f "%E, %e, %P" ~/Dev/mappability-build/release/mappability -I /srv/public/cpockrandt/chr13-index/index -O /srv/public/cpockrandt/gmapp_results/chr13 -E ${E} -K ${K} -o ${O} -m -t 20 --knut2
    done
done

