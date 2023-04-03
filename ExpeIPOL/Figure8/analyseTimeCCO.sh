#!/bin/bash


function usage(){
    echo $0 commandLine
    exit 1
}

STATFILE="statTime.dat"
echo "" > ${STATFILE}
for ((i=0; i<1000; i=i+5))
do
    echo "process $i.."
    a=$(time -p (generateTree3D -n $i -a 200000   -x graphExport.xml >/dev/null 2>&1)  2>&1)
    a=${a##*user}
    a=${a%%s*}
    a=$(echo $a | tr -d '\n')
    echo "Res=${a}=" 
    echo "$i $a" >> ${STATFILE}
done
