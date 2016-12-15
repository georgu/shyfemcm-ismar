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
  echo "Usage: merge_ens_output.sh [file_type] [ens_member]"
  echo "file_type: ext ets shy"
  echo "ens_member: number of ens member"
  exit 0
}

#----------------------------------------------------------

Merge_files()
{
   type=$1
   files=$(ls an*_en${ens_member}b.${type})
   rm -f z.* u.* v.* m.* total_${ens_member}_*
   for fil in $files; do
      echo "Processing file: $fil"
      memory -s $fil > log
      $FEMDIR/fembin/split${type} >> log
      for flev in $(ls z.*); do
          cat $flev >> total_${ens_member}_$flev
      done
      for flev in $(ls u.*); do
          cat $flev >> total_${ens_member}_$flev
      done
      for flev in $(ls v.*); do
          cat $flev >> total_${ens_member}_$flev
      done
   done
}
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------

if [ $2 ]; then
   file_type=$1
   ens_member=$2
else
   Usage
fi

if [ $file_type = 'ets' ]; then
   Merge_files 'ets'
elif [ $file_type = 'ext' ]; then
   Merge_files 'ext'
elif [ $file_type = 'shy' ]; then
   files=$(ls an*_en${ens_member}b.hydro.shy)
   for fil in $files; do
      echo "Processing file: $fil"
      cat $fil >> total_$fil
   done
else
   echo "Unknown file type"
   exit 1
fi

