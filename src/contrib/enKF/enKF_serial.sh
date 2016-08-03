#!/bin/bash
#
# Runs the Ensemble Kalman Filter analysis. It needs a configuration file with 
# these inputs:
# - name of the bas file with the extension;
# - skel-file for the str (it is possible to give a list of nens skel-files, using a standard name "skel_list.txt");
# - name of the list of the restart files. This determines the n of ens members;
# - name of the list of the observation files, containing rows like: time obs_file.
#
# You must prepare the restart files to be read at the initial time
# and the observation files, one at each time.
#
# Observations (see read_obs):
# The file of observations must be written in ascii format, separated by spaces.
# format:
# n_of_rows
# time obs_type(5 chars) x y z value stand_dev
# time obs_type(5 chars) x y z value stand_dev
# ...
#
#----------------------------------------------------------

# This finds the path of the current script
SCRIPT=$(realpath $0)
SCRIPTPATH=$(dirname $SCRIPT)

FEMDIR=$SCRIPTPATH/..	# fem directory
SIMDIR=$(pwd)		# current dir

#----------------------------------------------------------

Usage()
{
  echo "Usage: enKF_serial.sh [conf-file]"
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
     # Name of the bas file
     if [ $nrows = 0 ]; then
        bas_file=$line
        Check_file $bas_file
     # Name of the skel file to make the str files
     elif [ $nrows = 1 ]; then
        sk_file=$line
	Check_file $sk_file
     # List of the initial restart files
     elif [ $nrows = 2 ]; then
        rst_file_list=$line
	Check_file $rst_file_list
     # List of the observation files
     elif [ $nrows = 3 ]; then
        obs_file_list=$line
        Check_file $obs_file_list
     else
        echo "Too many rows"
        exit 1
     fi
     nrows=$((nrows + 1))
  done < $1

  if [ $sk_file = 'skel_list.txt' ]; then
     many_skels=1
     nrows=0
     while read line
     do
        Check_file $line
        skel_file[$nrows]=$line
        nrows=$((nrows + 1))
     done < $sk_file
     echo "Using different ensemble skel-files for the ens simulations"
  else
     many_skels=0
  fi
}

#----------------------------------------------------------

Read_rst_list(){
# Reads the list of the initial restart files and links them in order
# to have a standard name (an000_en***b.rst) and determines
# the size of the ensemble problem (nrens)
  echo "Reading list of restart files for the initial ensemble"
  
  rm -f an001_en*b.rst

  nrow=0
  while read line
  do
    Check_file $line
    nel=$(printf "%03d" $nrow)
    ln -s $line an001_en${nel}b.rst
    nrow=$((nrow + 1))
  done < $rst_file_list
  nrens=$nrow
  echo "Number of ensemble members: $nrens"
}


#----------------------------------------------------------

Read_obs_list(){
# Reads the list of the observation files and determines
# the number of analysis steps (nran)
  echo "Reading list of observation files"
  Check_file $obs_file_list
  
  # timeo nad obsfile starts from 1 as the analysis steps
  nrow=1
  while read line
  do
    timeo[$nrow]=$(echo $line | cut -d " " -f 1)
    Check_num -100000000 100000000 'real' ${timeo[$nrow]}
    obsfile[$nrow]=$(echo $line | cut -d " " -f 2)
    Check_file ${obsfile[$nrow]}
    nrow=$((nrow + 1))
  done < $obs_file_list
  nran=$((nrow - 1))
  echo "Number of analysis steps: $nran"
}

#----------------------------------------------------------

SkelStr(){
# Makes a str file from a skel file
inamesim=$1; iitrst=$2; iitanf=$3; iitend=$4; iidtout=$5; idragco=$6
inomp=$7; irestrt=$8; ibasin=$9; iskelname=${10}; istrname=${11}        #needs brackets!

[ $iitrst = 'none' ] && iitrst="''"
[ $irestrt = 'none' ] && irestrt="" 

if [ ! -s $iskelname ]; then
        echo "File $iskelname does not exist"
        exit 1
fi

irestrt="'$irestrt'"

cat $iskelname | sed -e "s/NAMESIM/$inamesim/g" |  sed -e "s/ITRST/$iitrst/g" | \
        sed -e "s/ITANF/$iitanf/g" | sed -e "s/ITEND/$iitend/g" | \
        sed -e "s/IDTOUT/$iidtout/g" |  sed -e "s/DRAGCO/$idragco/g" | \
        sed -e "s/NOMP/$inomp/g" |  sed -e "s/RESTRT/$irestrt/g" | \
        sed -e "s/BASIN/$ibasin/g" >  $istrname
}

#----------------------------------------------------------

Make_sim()
{
  basen=$(basename $1 .str)
  $2/shyfem $1 > $basen.log 
}

#----------------------------------------------------------

Write_info_file(){
# Write a file with informations for the fortran analysis program

  na=$1

  rm -f analysis.info
  echo $nrens > analysis.info		# nr of ens members
  echo $na >> analysis.info		# analysis step
  echo $bas_file >> analysis.info	# name of the basin
  echo ${obsfile[$na]} >> analysis.info	# obs file
}

#----------------------------------------------------------

Run_ensemble_analysis()
{
nanl=$(printf "%03d" $1)

# Run fortran program
cd $FEMDIR/enKF
make cleanall > make.log
make enkf_analysis > make.log
rm -f make.log

cd $SIMDIR
$FEMDIR/enKF/enkf_analysis

# Check restart files
for (( ne = 0; ne < $nrens; ne++ )); do
	nensl=$(printf "%03d" $nens)
	filename="an${nanl}_en${nensl}a.rst"
	Check_file $filename
done

filename="an${nanl}_enavrb.rst"	#Average
Check_file $filename
filename="an${nanl}_enavra.rst"	#Average
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

# Reading rst file list
Read_rst_list

# Reading obs file list
Read_obs_list

# Assimilation cycle for every analysis time step
for (( na = 1; na <= $nran; na++ )); do

   # make the analysis
   echo; echo "			ANALYSIS STEP $na OF $nran"; echo
   Write_info_file $na
   Run_ensemble_analysis $na

   if [ $na != $nran ]; then # not the last one

    # run nrens sims before the obs
    echo; echo "       Running $nrens ensemble runs..."
    for (( ne = 0; ne < $nrens; ne++ )); do

      # Use different skel files if a list is provided
      if [ $many_skels = 1 ]; then
         ens_skel_file=${skel_file[$ne]}
      else
         ens_skel_file=$sk_file
      fi
      Check_file $ens_skel_file

      nel=$(printf "%03d" $ne); nal=$(printf "%03d" $na)
      naa=$((na + 1)); naal=$(printf "%03d" $naa)
      name_sim="an${naal}_en${nel}b"; itrst=${timeo[$na]}; itanf=${timeo[$na]}
      itend=${timeo[$naa]}; idtout=300; wdrag=0.0025; nomp=1
      rstfile="an${nal}_en${nel}a.rst"; strnew="${name_sim}.str"
      SkelStr $name_sim $itrst $itanf $itend $idtout $wdrag $nomp $rstfile $bas_file $ens_skel_file $strnew
      Make_sim $strnew $FEMDIR/fem3d

    done
    echo "       Done"; echo

   fi

done

exit 0
