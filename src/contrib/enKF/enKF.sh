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
        nrens=$line
        Check_num $ensmin $ensmax 'int' $nrens
     # str-file of the base simulation. First time must be the same of rst-files
     elif [ $nrows = 1 ]; then
        strfile=$line
	Check_file $strfile
     # Basename of the rst-files
     elif [ $nrows = 2 ]; then
        baserst=$line
     # name of the obs-file_list, containing rows like: time obs_file
     elif [ $nrows = 3 ]; then
        obslist=$line
	Check_file $obslist
     else
        echo "Too many rows"
        exit 1
     fi
     nrows=$((nrows + 1))
  done < $1
}

#----------------------------------------------------------

Read_obs_list()
{
  echo "Reading list of observation files"
  Check_file $obslist
  
  nrows=0
  while read line
  do
    timeo[$nrows]=$(echo $line | cut -d " " -f 1)
    Check_num -100000000 100000000 'real' $timeo[$nrows]
    obsfile[$nrows]=$(echo $line | cut -d " " -f 2)
    Check_file $obsfile[$nrows]
    nrows=$((nrows + 1))
  done
  nran=$nrows
}

#----------------------------------------------------------

Make_str()
{
  nens=$1; nan=$2; obstime=$3; strfile=$4
  basestr=$(basename $strfile .str)

  nensl=$(printf "%03d" $nens)
  nanl=$(printf "%03d" $nan)
  strnew=${basestr}_st${nanl}_en${nensl}.str

  namesim="st${nanl}_en${nensl}"
  itanf=$timeo[$((nan - 1))]
  itend=$timeo[$nan]
  itrst=$itanf
  idtrst=$(echo "$itend - $itanf" | bc)
  #rstname="an_${namesimold}.rst"
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

Write_info_file()
{
  nanl=$(printf "%03d" $na)

  #basname=$($FEMDIR/fembin/strparse.pl -value=basin  $4) TODO
  basname=$(sed -n '8p' < $4) #TMP
  Check_file $basname.bas

  sdate=$($FEMDIR/fembin/strparse.pl -value=date $4)
  stime=$($FEMDIR/fembin/strparse.pl -value=time $4)
  [[ $stime = '(unknown)' ]] && stime=0
    
  # Write analysis.info, file for the fortran program
  rm -f analysis.info
  echo $date > analysis.info		# sim date0
  echo $time >> analysis.info		# sim time0
  echo $2 >> analysis.info		# obs time
  echo $1 >> analysis.info		# nr of ens members
  echo $basname >> analysis.info	# name of the basin
  echo "${rstfile}_st${nanl}" >> analysis.info	# basname of the restart files
  echo $obsfile[$na] >> analysis.info	# obs file
}

#----------------------------------------------------------

Run_ensemble_analysis()
{
# Run fortran program
make enKF_analysis
./enKF_analysis

# Check restart files
for (( ne = 1; ne <= $nrens; ne++ )); do
	nensl=$(printf "%03d" $nens)
	filename="${rstfile}_st${nanl}_en${nensl}_an.rst"
	Check_file $filename
done

filename="${rstfile}_st${nanl}_avr_bk.rst"	#Average
Check_file $filename
filename="${rstfile}_st${nanl}_avr_an.rst"	#Average
Check_file $filename
}

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------

if [ $1 ]; then
   Read_conf $1
else
   Usage
fi

# Reading obs list file
Read_obs_list

# Assimilation cycle for every analysis time step
for (( na = 0; na < $nran; na++ )); do

   if [ $na != 0 ]; then # Suppose to start from t=first obs
    # run nrens sims before the obs
    echo; echo "*** Running $nrens ensemble runs ***"
    for (( ne = 1; ne <= $nrens; ne++ )); do
      Make_str $ne $na ${timeo[$na]} $basestr
      Make_sim $strnew 
    done
   fi

   # make the analysis
   echo; echo "*** Analysis step $na of $nran ***"
   Write_info_file $nrens ${timeo[$na]} $basestr
   Run_ensemble_analysis

done

exit 0
