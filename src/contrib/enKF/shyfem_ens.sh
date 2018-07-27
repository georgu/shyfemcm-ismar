#!/bin/bash
#
# Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
#
# Run many shyfem simulations using different str files with a common
# basename. This program uses "parallel".
#

Usage(){
 echo
 echo "Run many shyfem simulations using different str files with a common"
 echo "basename. This program uses parallel."
 echo
 echo "Usage: shyfem_ens.sh [n. of threads] [basename]"
 echo "with basename the base name of the str files"
 exit 0
}

#-------------------------------------------------------

Check_file()
{
  if [ ! -s $1 ]; then
     echo "File $1 does not exist or has zero size."
     exit 1
  fi
}


#-------------------------------------------------------

Check_exec(){
  command -v parallel > /dev/null 2>&1 || \
   { echo "parallel it's not installed.  Aborting." >&2; exit 1; }
}

#-------------------------------------------------------

Make_sim(){
  basen=$(basename $1 .str)
  $2/shyfem $1 > $basen.log
}

#-------------------------------------------------------
#-------------------------------------------------------

# This finds the path of the current script
SCRIPT=$(realpath $0)
SCRIPTPATH=$(dirname $SCRIPT)

FEMDIR=$SCRIPTPATH/..   # fem directory

if [ $2 ]; then
   nth=$1
   fbasename=$2
else
   Usage
fi

Check_exec

for strf in $(ls ${fbasename}*.str); do
  Check_file $strf
done
strfiles=$(ls ${fbasename}*.str)

# with nthreads=0 uses the maximum number
export -f Make_sim
parallel --no-notice -P $nth Make_sim ::: $strfiles ::: $FEMDIR/fem3d

exit 0
