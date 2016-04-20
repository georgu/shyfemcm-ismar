#!/bin/bash
#
# Runs the Ensemble Kalman Filter analysis. It needs a configuration file with 
# these inputs:
# - number of ensemble members (nrens);
# - basename of nrens str-files from which to start the simulations;
# - number of analysis steps (nran);
# - basename of nran obs-files;
# - nran obs-times.
#
# You must prepare also nrens restart files to be read in the initial str-files
# and nran observation files.
#
# Example of str filename: basename_an000_en001.str
#
# Initial rst names: an000_en001.rst
#
# Observations (see read_obs):
# The file of observations must be written in ascii format, separated by spaces.
# One file for each analysis step must be used, which includes all the observations 
# for that specific time.
# observation filenames: myobs_an000.obs
# format:
# time n_of_records
# ...
# obs_type(5 chars) x y z value stand_dev
# obs_type(5 chars) x y z value stand_dev
# ...
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
        nran=$line
     # basename for the obs-files
     elif [ $nrows = 3 ]; then
        baseobs=$line
     # times of observations
     elif [ $nrows -gt 3 ] & [ $nrows -le $((nran+3)) ]; then
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

  if [ $nanold = 0} ]; then
     rstname="${namesimold}.rst"
  else
     rstname="an_${namesimold}.rst"
  fi
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

Make_sim()
{
  basen=$(basename $1 .str)
  $FEMDIR/fem3d/shyfem $1 > $basen.log
}

#----------------------------------------------------------

Make_analysis()
{
  echo; echo "*** Assimilation step $na of $nran ***"

  #basname=$($FEMDIR/fembin/strparse.pl -value=basin  $4) TODO
  basname=$(sed -n '8p' < $4) #TMP
  Check_file $basname.bas

  sdate=$($FEMDIR/fembin/strparse.pl -value=date $4)
  stime=$($FEMDIR/fembin/strparse.pl -value=time $4)
  [[ $stime = '(unknown)' ]] && stime=0
    
  rm -f analysis.info
  echo $basname > analysis.info	# name of the basin
  echo $date >> analysis.info	# sim date0
  echo $time >> analysis.info	# sim time0
  echo $3 >> analysis.info		# obs time
  echo $1 >> analysis.info		# nr of ens members
  echo $na >> analysis.info		# n of analysis step
  echo $baseobs >> analysis.info	# basename of obs files
  for (( ne = 1; ne <= $1; ne++ )); do
      nanl=$(printf "%03d" $2); nensl=$(printf "%03d" $ne)
      rstname="an${nanl}_en${nensl}.rst"
      Check_file $rstname
      $rstname >> analysis.info	# names of the rst files
  done
  $FEMDIR/enKF/enKF_analysis
}

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------

if [ $1 ]; then
   Read_conf $1
else
   Usage
fi

# Assimilation sims
for (( na = 1; na <= $nran; na++ )); do

   # run nrens sims before the obs
   echo; echo "*** Running $nrens ensemble runs. Step $na of $nran ***"
   for (( ne = 1; ne <= $nrens; ne++ )); do
      Make_str $ne $na ${timeo[$na]} $basestr
      Make_sim $strnew 
   done

   # make the analysis
   Make_analysis $nrens $na ${timeo[$na]} $basestr

done

exit 0
