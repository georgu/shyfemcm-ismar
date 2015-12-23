#!/bin/bash
#
# Runs the Ensemble Kalman Filter analysis. It needs a configuration file with 
# these inputs:
# - number of ensemble members (nrens);
# - basename of nrens str-files from which to start the simulations;
# - number of analysis steps (nranl);
# - basename of nranl obs-files;
# - nranl obs-times.
#
# You must prepare also nrens restart files to be read in the initial str-files
# and nranl observation files. Each observation file must include all the observations
# to be used in that analysis step. Each row of the file must contain:
# time x y z obs_id value error value_1 ... value_nrens
#
#----------------------------------------------------------

FEMDIR=$HOME/shyfem

# min/max number of ensemble members
ensmin=10
ensmax=50
# min/max number of analysis steps
anmin=1
anmax=100

#----------------------------------------------------------

Usage()
{
  echo "Usage: enKF.sh [opt-file]"
  exit 0
}

#----------------------------------------------------------

Check_file()
{
  if [ ! -s $1 ]; then
     echo "File $1 does not exist or has zero size."
     exit 1
  fi
}

#----------------------------------------------------------

Check_num()
{
  nint='^[0-9]+$'
  nreal='^-?[0-9]+([.][0-9]+)?$'

  if [ $3 = 'int' ]; then
     if ! [[  $4 =~ $nint ]]; then
        echo "$4 is not an integer number"
        exit 1
     fi
  elif [ $3 = 'real' ]; then
     if ! [[ $4 =~ $nreal ]]; then
        echo "$4 is not a real number"
        exit 1
     fi
  fi
  if [[ $(echo "$4 < $1" | bc) = 1 ]] || [[ $(echo "$4 > $2" | bc) = 1 ]]; then
     echo "Number $4 out of range"
  fi 
}

#----------------------------------------------------------

Read_conf()
{
  echo "Reading configuration file: $1"
  Check_file $1

  nrows=0
  while read line
  do
     # nr of ens members
     if [ $nrows = 0 ]; then
        Check_num $ensmin $ensmax 'int' $line
        nrens=$line
     # basename for the str-files
     elif [ $nrows = 1 ]; then
        basestr=$line
     # nr of analysis steps
     elif [ $nrows = 2 ]; then
        Check_num $anmin $anmax 'int' $line
        nranl=$line
     # basename for the obs-files
     elif [ $nrows = 3 ]; then
        baseobs=$line
     # times of observations
     elif [ $nrows -gt 3 ] & [ $nrows -le $((nranl+3)) ]; then
        Check_num -100000000 100000000 'real' $line
        timeo[$((nrows - 4))]=$line
     else
        echo "Too many rows"
        exit 1
     fi
     nrows=$((nrows + 1))
  done < $1
}

#----------------------------------------------------------

if [ $1 ]; then
   Read_conf $1
else
   Usage
fi

# Assimilation sims
for (( na = 1; na <= $nranl; na++ )); do
   # run nrens sims before the obs
   # make the analysis 
done
# run last sims from the last analysis and save restarts

