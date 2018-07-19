#!/bin/bash

FILE=$1

DATA=$(cat $FILE | sed '/^\[/ d' | sed '/^For/ d' | sed '/^Ski/ d' | sed '/^# KNUT/ d' | sed -e 's/### K = //g' -e 's/###//g' -e 's/, O = /,/g' -e 's/ (/,/g' -e 's/)//g' -e 's/\(.*\), \(.*\), \([0-9]*\)%/\2,\3/g' -e 's/,/;/g' -e 's/\./,/g')

#echo -e "${DATA}"

FIRST_K=1
K=0

PERC_GLOBAL=""

TIME_WITHOUT_GLOBAL=""
#TIME_WITH_GLOBAL=""
CPU_WITHOUT_GLOBAL=""
#CPU_WITH_GLOBAL=""

while read -r line; do
    IFS=";" read CURRENT_K PERC O <<< "${line}"
    if [ "$K" -eq "$CURRENT_K" ]
    then
        if [ "$FIRST_K" -eq 1 ]
        then
            PERC_GLOBAL+=";${PERC}"
        fi
    else
        if [ "$K" -ne 0 ]
        then
            TIME_WITHOUT_GLOBAL+="${TIME_WITHOUT}${TIME_WITH}\n"
            #TIME_WITH_GLOBAL+="${TIME_WITH}\n"
            CPU_WITHOUT_GLOBAL+="${CPU_WITHOUT}${CPU_WITH}\n"
            #CPU_WITH_GLOBAL+="${CPU_WITH}\n"
            FIRST_K=0
        else
            PERC_GLOBAL+=";${PERC}"
        fi
        K=$CURRENT_K
        TIME_WITHOUT="${K}"
        TIME_WITH=";"
        CPU_WITHOUT="${K}"
        CPU_WITH=";"
    fi
    read -r line
    IFS=";" read _TIME _CPU <<< "${line}"
    TIME_WITHOUT+=";${_TIME}"
    CPU_WITHOUT+=";${_CPU}"

    read -r line
    IFS=";" read _TIME _CPU <<< "${line}"
    TIME_WITH+=";${_TIME}"
    CPU_WITH+=";${_CPU}"

done <<< "$DATA"

echo -e "${PERC_GLOBAL};${PERC_GLOBAL}\n${TIME_WITHOUT_GLOBAL}${TIME_WITHOUT}${TIME_WITH}"
echo ""
echo -e "${PERC_GLOBAL};${PERC_GLOBAL}\n${CPU_WITHOUT_GLOBAL}${CPU_WITHOUT}${CPU_WITH}"

