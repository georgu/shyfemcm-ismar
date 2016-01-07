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
# Example of str filename: basename_an000_en001.str
# Initial rst names: an000_en001.rst
#
#----------------------------------------------------------

FEMDIR=$HOME/git_shyfem/shyfem

# min/max number of ensemble members
ensmin=2
ensmax=50
# min/max number of analysis steps
anmin=1
anmax=100

#----------------------------------------------------------

Usage()
{
  echo "Usage: enKF.sh [conf-file]"
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
     exit 1
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

Make_str()
{
  nens=$1; nan=$2; atime=$3; basestr=$4

  nensl=$(printf "%03d" $nens)
  nanl=$(printf "%03d" $nan)

  nanold=$((nan - 1))
  nanoldl=$(printf "%03d" $nanold)

  strold=${basestr}_an${nanoldl}_en${nensl}.str
  Check_file $strold
  strnew=${basestr}_an${nanl}_en${nensl}.str

  namesim="an${nanl}_en${nensl}"
  namesimold="an${nanoldl}_en${nensl}"
  itanf=$($FEMDIR/fembin/strparse.pl -value=itend $strold)
  itend=$atime
  itrst=$itanf
  idtrst=$(echo "$itend - $itanf" | bc)
  rstname="${namesimold}.rst"
  Check_file $rstname

  $FEMDIR/fembin/strparse.pl -value=itanf -replace=$itanf $strold
  $FEMDIR/fembin/strparse.pl -value=itend -replace=$itend replace.str
  $FEMDIR/fembin/strparse.pl -value=itrst -replace=$itrst replace.str
  $FEMDIR/fembin/strparse.pl -value=idtrst -replace=$idtrst replace.str
  # TMP
  # namesim
  # rstname
  cat replace.str | \
      sed -e "s/$namesimold/$namesim/" | \
      sed -e "s/restrt.\+/restrt = \'$rstname\'/" > replace1.str
  mv -f replace1.str $strnew; rm -f replace.str
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
   for (( ne = 1; ne <= $nrens; ne++ )); do
      Make_str $ne $na ${timeo[$na]} $basestr
      #Make_sim $ne $na #TODO
   done
   # make the analysis
   #Make_analysis $na ${timeo[$na]} #TODO
done

exit 0
