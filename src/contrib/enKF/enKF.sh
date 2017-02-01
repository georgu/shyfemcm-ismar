#!/bin/bash
#
# Ensemble Kalman Filter for SHYFEM. 
# 
# Marco Bajo, ISMAR-CNR Venice - 2016
#
# See the README file to set-up
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
  echo "Usage: enKF.sh [n. of threads] [conf-file]"
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
  echo "Read the configuration file: $1"
  Check_file $1

  nrows=0
  while read line
  do

     # Name of the bas file
     if [ $nrows = 0 ]; then
        bas_file=$line
        Check_file $bas_file

     # List of the skel files
     elif [ $nrows = 1 ]; then
        skel_file_list=$line
	Check_file $skel_file_list

     # List of the initial restart files
     elif [ $nrows = 2 ]; then
        rst_file_list=$line
	Check_file $rst_file_list

     # List of the observation files
     elif [ $nrows = 3 ]; then
        obs_file_list=$line
        Check_file $obs_file_list

     # List of times with observations
     elif [ $nrows = 4 ]; then
        obs_time_list=$line
        Check_file $obs_time_list

     # Dimension of the state vector: nkn nel nlv
     elif [ $nrows = 5 ]; then
        sdim=$line
        
     # If 1 makes a new ens of initial states from 1
     elif [ $nrows = 6 ]; then
        is_new_ens=$line
        Check_num 0 1 'int' $is_new_ens

     # If 1 uses an augmented state with the model errors
     elif [ $nrows = 7 ]; then
        is_mod_err=$line
        Check_num 0 1 'int' $is_mod_err

     else

        echo "Too many rows"
        exit 1

     fi
     nrows=$((nrows + 1))
  done < $1
}

#----------------------------------------------------------

Compile_enkf(){

  nkn=$1
  nel=$2
  nlv=$3

  cd $FEMDIR/enKF

  cat mod_dimensions.skel | sed -e "s/NKN/$nkn/" | sed -e "s/NEL/$nel/" | \
	sed -e "s/NLV/$nlv/" > mod_dimensions.F90

  make cleanall > $SIMDIR/make.log
  make enkf_analysis >> $SIMDIR/make.log
  cd $SIMDIR
}


#----------------------------------------------------------

Check_exec(){
  echo "Check the exec programs"
  command -v parallel > /dev/null 2>&1 || { echo "parallel it's not installed.  Aborting." >&2; exit 1; }
  [ ! -s $FEMDIR/fem3d/shyfem ] && echo "shyfem exec does not exist. Compile the model first." && exit 1
  # Make here the mod_dimensions and compile enkf_analysis
  
  [ ! -s $FEMDIR/enKF/enkf_analysis ] && echo "enkf_analysis exec does not exist. Compile the enKF first." && exit 1
}

#----------------------------------------------------------

Read_skel_list(){
# Reads the list of skel files used to make the str files in the
# ensemble data assimilation window
  echo "Read the list of skel files"

  nrow=0
  while read line
  do
     Check_file $line
     skel_file[$nrow]=$line
     nrow=$((nrow + 1))
  done < $skel_file_list
  nrens=$nrow
  echo ""; echo "Number of ensemble members: $nrens"; echo ""
}

#----------------------------------------------------------

Read_rst_list(){
# Reads the list of the initial restart files and links them in order
# to have a standard name (an000_en***b.rst) and determines
# the size of the ensemble problem (nrens)
  echo "Read the list of restart files"
  
  rm -f an001_en*b.rst

  nrow=0
  while read line
  do
    Check_file $line
    nel=$(printf "%03d" $nrow)
    ln -s $line an001_en${nel}b.rst
    nrow=$((nrow + 1))
  done < $rst_file_list
  if [ $nrow -ne $nrens ]; then
     echo "The number of restart files differs from the number of skel files"
     echo "You must specify the same number, which is the dimension of the ensemble"
     exit 1
  fi
}

#----------------------------------------------------------

Read_obs_time_list(){
# Reads the list of the observation times and determines
# the number of analysis steps (nran)
  echo "Read the list of observation times"
  
  # timeo starts from 1 as the analysis steps
  nrow=1
  while read line
  do
    timeo[$nrow]=$line
    Check_num -100000000 100000000 'real' ${timeo[$nrow]}
    nrow=$((nrow + 1))
  done < $obs_time_list
  nran=$((nrow - 1))
  echo ""; echo "Number of analysis steps: $nran"
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
        sed -e "s/BASIN/$ibasin/g" | sed -e "s/IDTRST/-1/" >  $istrname
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

  echo $nrens > analysis.info		# nr of ens members
  echo $na >> analysis.info		# analysis step
  echo $bas_file >> analysis.info	# name of the basin
  echo ${timeo[$na]} >> analysis.info	# current time
  echo $obs_file_list >> analysis.info	# obs file list
  echo $is_new_ens >> analysis.info	# if to make a new ens of states
  echo $is_mod_err >> analysis.info	# if to use an augmented state with mod err
}

#----------------------------------------------------------

Run_ensemble_analysis()
{
nanl=$(printf "%03d" $1)

# Run fortran program
#cd $FEMDIR/enKF
#make cleanall > make.log
#make enkf_analysis > make.log
#rm -f make.log

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

Make_ens_str(){
strfiles=""
for (( ne = 0; ne < $nrens; ne++ )); do

   ens_skel_file=${skel_file[$ne]}
   Check_file $ens_skel_file

   nel=$(printf "%03d" $ne); nal=$(printf "%03d" $na)
   naa=$((na + 1)); naal=$(printf "%03d" $naa)
   name_sim="an${naal}_en${nel}b"; itrst=${timeo[$na]}; itanf=${timeo[$na]}
   itend=${timeo[$naa]}; idtout=300; wdrag=0.0025; nomp=1
   rstfile="an${nal}_en${nel}a.rst"; strnew="${name_sim}.str"
   SkelStr $name_sim $itrst $itanf $itend $idtout $wdrag $nomp $rstfile $bas_file $ens_skel_file $strnew
   strfiles="$strfiles $strnew"

done
}


#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------

if [ $2 ]; then
   nthreads=$1
   Read_conf $2
else
   Usage
fi

# Compiles the enKF code with the right total dimensions
Compile_enkf $sdim

# Checking the executable programs
Check_exec

# Reading skel file list
Read_skel_list

# Reading rst file list
Read_rst_list

# Reading obs file list
Read_obs_time_list

# Assimilation cycle for every analysis time step
for (( na = 1; na <= $nran; na++ )); do

   # make the analysis
   echo; echo "			ANALYSIS STEP $na OF $nran"; echo
   Write_info_file $na
   Run_ensemble_analysis $na

   if [ $na != $nran ]; then # not the last one

      # Makes nrens str files for the simulations
      Make_ens_str

      # run nrens sims before the obs
      echo; echo "       Running $nrens ensemble simulations..."

      # with nthreads=0 uses the maximum number
      export -f Make_sim
      #parallel --no-notice -j -k $nthreads Make_sim ::: $strfiles ::: $FEMDIR/fem3d
      parallel --no-notice -P $nthreads Make_sim ::: $strfiles ::: $FEMDIR/fem3d

      echo "       Done"; echo

   fi

done

exit 0
