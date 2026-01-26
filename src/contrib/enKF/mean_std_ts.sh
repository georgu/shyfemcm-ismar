#!/bin/bash
#
# Copyright (C) 2021, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
#
# Make mean and std of an ensemble of timeseries: time,value.
# files must have the same length and time-records. No check on times.
#

#----------------------------------------------------------

# This finds the path of the current script
SCRIPT=$(realpath $0)
SCRIPTPATH=$(dirname $SCRIPT)

FEMDIR=$SCRIPTPATH/..   # fem directory
SIMDIR=$(pwd)           # current dir

Usage()
{
  echo "Usage: ./mean_std_ts.sh [basename]"
  exit 0
}

###################################################################
[[ "$#" -ne "1" ]] && Usage

bfile=$1

k=0
for tfile in $(ls ${bfile}*); do
    if [ ! -s "$tfile" ]; then
	    echo "Bad $tfile"
	    exit 1
    fi
    echo "Ensemble member: $tfile"
    k=$((k+1))
done
echo
echo "Number of ensemble member: $k"
[[ "$k" -lt "3" ]] && echo "Stopping, too few ens members." && exit 1

awk -v nr=$k '{T[FNR]=$1; v1[FNR]+=$2; v2[FNR]+=$2^2} jmax<FNR {jmax=FNR} END {for (j=1; j<=jmax; j++) 
{df=v2[j]/nr-(v1[j]/nr)^2; print T[j],v1[j]/nr,sqrt(df+0.000000000001)}}' ${bfile}* > ${bfile}mean_std.dat

