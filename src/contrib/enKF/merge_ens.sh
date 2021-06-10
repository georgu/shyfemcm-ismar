#/bin/bash
#
# Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
#
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

   # not working!
   #$FEMDIR/fembin/shyelab -catmode +1 -out an*_en${nen}b.${ftype}
   #mv -f out.ext analysis_en${nen}.${ftype}

   rm -f *_st*_en${nen}.ts *.2d.* *.3d.*
   for fil in $files; do
      [[ ! -s "$fil" ]] && echo "Error in file: $fil" && exit 1
      echo "Processing file: $fil"
      $FEMDIR/fembin/shyelab -split ${fil} &> /dev/null

      var2d=$(ls *.2d.1 2>/dev/null |cut -d '.' -f 1,2)
      var3d=$(ls *.3d.1 2>/dev/null |cut -d '.' -f 1,2)
      for vv in $(echo "$var2d $var3d"); do
         for flev in $(ls $vv.*); do
	    idst=$(echo $flev | cut -d '.' -f 3)
            cat $flev |head -n -1|tail -n +2 >> ${vv}_st${idst}_en${nen}.ts
	    rm $flev
         done
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
for efile in $(ls an00002_en*.${file_type}); do

    ens_member=${efile:10:5}

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
