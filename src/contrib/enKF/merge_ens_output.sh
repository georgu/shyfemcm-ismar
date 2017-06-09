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
  echo "Usage: merge_ens_output.sh [file_type]"
  echo "file_type: ext ets"
  exit 0
}

#----------------------------------------------------------

Merge_files()
{
   type=$1
   files=$(ls an*_en${ens_member}b.${type})
   rm -f z.* u.* v.* m.* tot${ens_member}_*
   for fil in $files; do
      echo "Processing file: $fil"
      memory -s $fil > log
      $FEMDIR/fembin/split${type} >> log
      for flev in $(ls z.*); do
          cat $flev >> tot${ens_member}_$flev
      done
      for flev in $(ls u.*); do
          cat $flev >> tot${ens_member}_$flev
      done
      for flev in $(ls v.*); do
          cat $flev >> tot${ens_member}_$flev
      done
   done
}
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------

if [ $1 ]; then
   file_type=$1
else
   Usage
fi

k=0
for efile in $(ls an002_en*.${file_type}); do
    ens_member=${efile:8:3}
    Merge_files ${file_type}
done
