#!/bin/bash
#
# Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
#
# Ensemble Kalman Filter for SHYFEM. 
# 
# 2016 first version
# 2018-06 important updates
#
# See the README file
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
  echo "Usage: enKF.sh [n] [out]"
  echo
  echo "n = n. of threads"
  echo "out = 0 only mean and std, 1 all members (for enKS)"
  echo
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
  cfile='assconf.dat'
  echo "Read the configuration file: $cfile"
  Check_file $cfile

  nrows=0
  while read line
  do

     # Name of the bas file
     if [ "$nrows" -eq "0" ]; then
        bas_file=$line
        Check_file $bas_file

     # Dimension of the state vector: nkn nel nlv
     elif [ "$nrows" -eq "1" ]; then
        sdim=$line
        
     # idtrst
     elif [ "$nrows" -eq "2" ]; then
        idtrst=$line

     # If 1 makes a new ens of initial states from 1
     elif [ "$nrows" -eq "3" ]; then
        is_new_ens=$line
        Check_num 0 1 'int' $is_new_ens

     else

        echo "Too many rows"
        exit 1

     fi
     nrows=$((nrows + 1))
  done < $cfile

  ens_file_list='ens_list.txt'
  Check_file $ens_file_list  
  obs_file_list='obs_list.txt'
  Check_file $obs_file_list
  an_time_list='antime_list.txt'
  Check_file $an_time_list
}

#----------------------------------------------------------

Compile_enkf(){

  nkn=$1
  nel=$2
  nlv=$3

  cd $FEMDIR/enKF

  cat mod_dimensions.skel | sed -e "s/NKN/$nkn/" | sed -e "s/NEL/$nel/" | \
	sed -e "s/NLV/$nlv/" > mod_dimensions.F90

  make clean > $SIMDIR/make.log
  make main >> $SIMDIR/make.log

  cd $SIMDIR
}


#----------------------------------------------------------

Check_exec(){
  echo "Check the exec programs"
  command -v parallel > /dev/null 2>&1 || { echo "parallel it's not installed.  Aborting." >&2; exit 1; }
  [ ! -s $FEMDIR/fem3d/shyfem ] && echo "shyfem exec does not exist. Compile the model first." && exit 1
  # Make here the mod_dimensions and compile main
  
  [ ! -s $FEMDIR/enKF/main ] && echo "main exec does not exist. Compile the enKF first." && exit 1
}

#----------------------------------------------------------

Read_ens_list(){
# Reads the list of skel and restart files of the ensemble
  echo "Reading the ensemble list"

  rm -f an00001_en*b.rst

  nrow=0
  while read line
  do
     skelf=$(echo $line | cut -d " " -f 1)
     rstf=$(echo $line | cut -d " " -f 2)
     Check_file $skelf
     Check_file $rstf

     skel_file[$nrow]=$skelf

     if [ "$is_new_ens" -eq "0" ] || [ "$nrow" -eq "0" ]; then
        nel=$(printf "%05d" $nrow)
        ln -fs $rstf an00001_en${nel}b.rst
     fi

     nrow=$((nrow + 1))

  done < $ens_file_list

  nrens=$nrow
  echo ""; echo "Number of ensemble members: $nrens"; echo ""
}

#----------------------------------------------------------

Read_an_time_list(){
# Reads the list of the analysis times and determines
# the number of analysis steps (nran)
  echo "Read the list of analysis steps"
  
  # timeo starts from 1 as the analysis steps
  nrow=1
  while read line
  do
    timeo[$nrow]=$line
    nrow=$((nrow + 1))
  done < $an_time_list
  nran=$((nrow - 1))
  echo ""; echo "Number of analysis steps: $nran"
}


#----------------------------------------------------------

SkelStr(){
# Makes a str file from a skel file
inamesim=$1; iitanf=$2; iitend=$3; irestrt=$4; iidtrst=$5; iskelname=$6; istrname=$7

if [ ! -s $iskelname ]; then
        echo "File $iskelname does not exist"
        exit 1
fi

cat $iskelname | sed -e "s/NAMESIM/$inamesim/g" |  sed -e "s/ITANF/$iitanf/g" \
               | sed -e "s/ITEND/$iitend/g" | sed -e "s/RESTRT/$irestrt/" \
               | sed -e "s/IDTRST/$iidtrst/" \
               >  $istrname
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
}

#----------------------------------------------------------

Run_ensemble_analysis()
{
nanl=$(printf "%05d" $1)

cd $SIMDIR
$FEMDIR/enKF/main
if [ "$?" -ne "0" ]; then
          echo "Errors while running main."
          exit 1
fi

# Check restart files
for (( ne = 0; ne < $nrens; ne++ )); do
	nensl=$(printf "%05d" $nens)
	filename="an${nanl}_en${nensl}a.rst"
	Check_file $filename
done

#filename="an${nanl}_mean_state_b.rst"	#Average
#Check_file $filename
#filename="an${nanl}_mean_state_a.rst"	#Average
#Check_file $filename
}

#----------------------------------------------------------

Make_ens_str(){
strfiles=""
for (( ne = 0; ne < $nrens; ne++ )); do

   if [ "$is_new_ens" -eq "0" ]; then
      ens_skel_file=${skel_file[$ne]}
   else
      ens_skel_file=${skel_file[0]}
   fi
   Check_file $ens_skel_file

   nel=$(printf "%05d" $ne); nal=$(printf "%05d" $na)
   naa=$((na + 1)); naal=$(printf "%05d" $naa)

   itanf=${timeo[$na]}

   if [ "$na" -ne "$nran" ]; then
        name_sim="an${naal}_en${nel}b"
	itend=${timeo[$naa]}
        rstfile="an${nal}_en${nel}a.rst" 
        strnew="${name_sim}.str"

        SkelStr $name_sim $itanf $itend $rstfile $idtrst $ens_skel_file $strnew
        strfiles="$strfiles $strnew"
   fi

done
}


#----------------------------------------------------------
#----------------------------------------------------------
#	MAIN
#----------------------------------------------------------
#----------------------------------------------------------

if [ $2 ]; then
   nthreads=$1
   Read_conf
   out_verb=$2
else
   Usage
fi

# Export num of threads for main. Warning! Use the max num of cpu, not threads.
export OMP_NUM_THREADS=$nthreads

# Compiles the enKF code with the right total dimensions
Compile_enkf $sdim

# Checking the executable programs
Check_exec

# Reading skel file list
Read_ens_list

# Reading obs file list
Read_an_time_list

# Assimilation cycle for every analysis time step
rm -f X5*.uf backKF_*.rst analKF_*.rst 	# old files
for (( na = 1; na <= $nran; na++ )); do

   # make the analysis
   echo; echo "			ANALYSIS STEP $na OF $nran"; echo
   Write_info_file $na
   Run_ensemble_analysis $na

   # Makes nrens str files for the simulations
   Make_ens_str

   if [ "$na" -ne "$nran" ]; then # not the last one

      # run nrens sims before the obs
      echo; echo "       running $nrens ensemble simulations..."

      # with nthreads=0 uses the maximum number
      export -f Make_sim
      parallel --no-notice -P $nthreads Make_sim ::: $strfiles ::: $FEMDIR/fem3d

   fi

   # merge the rst files
   nanl=$(printf "%05d" $na)
   for (( ne = 0; ne < $nrens; ne++ )); do
        nel=$(printf "%05d" $ne)
        filename1="an${nanl}_en${nel}b.rst"
        filename2="an${nanl}_en${nel}a.rst"
        Check_file $filename1
        Check_file $filename2
	# If not verbose, saves only the last rst for each member
	if [ "$out_verb" -eq "1" ]; then
		[[ "$na" -gt "1" ]] && cat $filename1 >> backKF_en$nel.rst
		cat $filename2 >> analKF_en$nel.rst
	else
	        [[ "$na" -eq "$nran" ]] && mv -f $filename2 analKF_en$nel.rst
	fi
	rm -f $filename1 $filename2
   done
   filename1="an${nanl}_mean_a.rst"
   filename2="an${nanl}_std_a.rst"
   Check_file $filename1
   Check_file $filename2
   cat $filename1 >> analKF_mean.rst
   cat $filename2 >> analKF_std.rst
   rm -f $filename1 $filename2
   rm -f an*_en*b.inf an*_en*.log an*_en*b.str 

done
rm -f X5col.dat X5row.dat X5.uf make.log

exit 0
