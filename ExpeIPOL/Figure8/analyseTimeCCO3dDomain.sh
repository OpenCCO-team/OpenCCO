#!/bin/bash


function usage(){
    echo $0 commandLine
    exit 1
}

STATFILE="statTimeCCO3DimplicitDomain.dat"
echo "" > ${STATFILE}
for ((i=5; i<1000; i=i+5))
do
    echo "process $i.."
    a=$(time -p (generateTree3D -n $i -a 20000 -d ../../.ipol/data/maskLiver05.vol  >/dev/null 2>&1)  2>&1)
    a=${a##*user}
    a=${a%%s*}
    a=$(echo $a | tr -d '\n')
    echo "Res=${a}=" 
    echo "$i $a" >> ${STATFILE}
done
