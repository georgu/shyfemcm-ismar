#/bin/bash
# Merges all the ens files extracted from ext or shy files in the 
# analysis period.
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
  echo "Usage: merge_ens.sh [file_type]"
  echo "file_type: ext shy"
  exit 0
}

#----------------------------------------------------------

Merge_timeseries()
{
   nen=$1
   ftype=$2
   files=$(ls an*_en${nen}b.${ftype})
   rm -f zeta.2d.* velx.2d.* vely.2d.* speed.2d.* merged_en${nen}_*
   for fil in $files; do
      echo "Processing file: $fil"
      $FEMDIR/fembin/shyelab -split ${fil} > log
      for flev in $(ls zeta.2d.*); do
          cat $flev >> merged_en${nen}_$flev
      done
      for flev in $(ls velx.2d.*); do
          cat $flev >> merged_en${nen}_$flev
      done
      for flev in $(ls vely.2d.*); do
          cat $flev >> merged_en${nen}_$flev
      done
   done
}

#----------------------------------------------------------

Merge_shy()
{
   nen=$1
   files=$(ls an*_en${nen}b.hydro.shy)
   $FEMDIR/fembin/shyelab -out -catmode +1 $files
   mv -f out.shy en${nen}.hydro.shy
}


#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------

if [ $1 ]; then
   file_type=$1
else
   Usage
fi

# loop on ens members
#
for efile in $(ls an002_en*.${file_type}); do

    ens_member=${efile:8:3}

    if [ ${file_type} = 'ext' ]; then

       Merge_timeseries ${ens_member} ${file_type}

    elif [ ${file_type} = 'shy' ]; then

       Merge_shy ${ens_member}

    else

       echo "not a valid file"
       exit 1

    fi

done


exit 0
