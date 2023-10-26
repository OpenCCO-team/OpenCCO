#!/bin/bash


PARAMFILE="paramStatTime.txt"


function generateParamFile(){
N=$1
echo "SUPPLY_MAP: paramS.txt" > $PARAMFILE
echo "OXYGENATION_MAP: paramOx.txt" >> $PARAMFILE
echo "RANDOM_SEED: 1" >> $PARAMFILE
echo "PERF_POINT: 0 50 50" >> $PARAMFILE
echo "PERF_PRESSURE: 200000" >> $PARAMFILE
echo "TERM_PRESSURE: 83000" >> $PARAMFILE
echo "PERF_FLOW: 8.33" >> $PARAMFILE
echo "RHO: 0.036" >> $PARAMFILE
echo "GAMMA: 3" >> $PARAMFILE
echo "LAMBDA: 2" >> $PARAMFILE
echo "MU: 1" >> $PARAMFILE
echo "MIN_DISTANCE: 1" >> $PARAMFILE
echo "NUM_NODES: $N" >> $PARAMFILE
echo "VOXEL_WIDTH: 0.04" >> $PARAMFILE
echo "CLOSEST_NEIGHBOURS: 5" >> $PARAMFILE
}

function usage(){
    echo $0 commandLine
    exit 1
}

STATFILE="statTimeVascu.dat"


for ((i=5; i<1000; i=i+5))
do
    echo "process $i.."
    generateParamFile $i
    a=$(time -p (VascuSynth paramMainTime.txt imageNameListParamTime.txt 0.04 >/dev/null 2>&1)  2>&1)
    a=${a##*user}
    a=${a%%s*}
    a=$(echo $a | tr -d '\n')
    echo "Res=${a}=" 
    echo "$i $a" >> ${STATFILE}
done
