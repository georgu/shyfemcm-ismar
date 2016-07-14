#/bin/bash
# Merges all the ens output files extracted from ext or shy files in the analysis period.
# Note: shyplot is still not able to plot shy files with some records with the same time,
# as the merged one.
#----------------------------------------------------------

# This finds the path of the current script
SCRIPT=$(realpath $0)
SCRIPTPATH=$(dirname $SCRIPT)

FEMDIR=$SCRIPTPATH/..   # fem directory
SIMDIR=$(pwd)           # current dir

#----------------------------------------------------------

#----------------------------------------------------------

Usage()
{
  echo "Usage: merge_ens_output.sh [conf-file] [ file_type ]"
  echo "With file_type: ext shy"
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
        skel_file=$line
        Check_file $skel_file
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
}

#----------------------------------------------------------

Read_rst_list(){
# Reads the list of the initial restart files and links them in order
# to have a standard name (an000_en***b.rst) and determines
# the size of the ensemble problem (nrens)
  echo "Reading list of restart files for the initial ensemble"

  nrow=0
  while read line
  do
    Check_file $line
    nel=$(printf "%03d" $nrow)
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
#----------------------------------------------------------
#----------------------------------------------------------

if [ $2 ]; then
   Read_conf $1
else
   Usage
fi

file_type=$2

# Reading rst file list
Read_rst_list

# Reading obs file list
Read_obs_list

# Assimilation cycle for every analysis time step
rm -f *_z*.dat
for (( ne = 0; ne < $nrens; ne++ )); do

  nensl=$(printf "%03d" $ne)
  filout_ext="en${nensl}_z1_tot.dat"
  filout_shy="en${nensl}_tot.hydro.shy"
  if [ $file_type = 'ext' ]; then
        rm -f $filout_ext
        echo "File output $filout_ext"
  elif [ $file_type = 'shy' ]; then
        rm -f $filout_shy
        echo "File output $filout_shy"
  else
        echo "File type unknown"
        exit 1
  fi

  for (( na = 2; na <= $nran; na++ )); do

    nanl=$(printf "%03d" $na)

    if [ $file_type = 'ext' ]; then
       filename="an${nanl}_en${nensl}b.ext"
       Check_file $filename
       memory -s $filename > pp
       $FEMDIR/fem3d/splitext > pp
       cat z.1 >> $filout_ext
       rm -f pp z.1 u.1 v.1 m.1

    elif [ $file_type = 'shy' ]; then
       filename="an${nanl}_en${nensl}b.hydro.shy"
       Check_file $filename
       cat $filename >> $filout_shy

    fi

  done
done
