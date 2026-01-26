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
  echo
  echo "Usage: merge_ens.sh [file_type] [output]"
  echo
  echo "file_type: ext shy"
  echo "output: small (z), medium (z,t,s), full (all). Only for the ext files."
  exit 0
}

#----------------------------------------------------------

Merge_timeseries()
{
   nen=$1
   ftype=$2
   outt=$3

   files=$(ls an*_en${nen}b.${ftype})

   # not working!
   #$FEMDIR/fembin/shyelab -catmode +1 -out an*_en${nen}b.${ftype}
   #mv -f out.ext analysis_en${nen}.${ftype}

   allvars='all.2d dir.2d salt.2d speed.2d temp.2d velx.2d vely.2d zeta.2d dir.3d salt.3d speed.3d temp.3d velx.3d vely.3d'
   if [ "$outt" = "small" ]; then
	   vars='zeta.2d'
   elif [ "$outt" = "medium" ]; then
	   vars='zeta.2d temp.2d salt.2d temp.3d salt.3d'
   elif [ "$outt" = "full" ]; then
	   vars=$allvars
   else
	   echo "Bad input"
	   exit 1
   fi

   rm -f *_st*_en${nen}.ts *.2d.* *.3d.*
   for fil in $files; do
      [[ ! -s "$fil" ]] && echo "Error in file: $fil" && exit 1
      echo "Processing file: $fil"
      $FEMDIR/fembin/shyelab -split ${fil} &> /dev/null

      for vv in $vars; do
         for flev in $(ls $vv.*); do
	    idst=$(echo $flev | cut -d '.' -f 3)
            cat $flev |head -n -1|tail -n +2 >> ${vv}_st${idst}_en${nen}.ts
	    rm $flev
         done
      done
   done
   rm -f *.2d.* *.3d.*
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

if [ $2 ]; then
   file_type=$1
   outt=$2
else
   Usage
fi

# loop on ens members
#
for efile in $(ls an00002_en*.${file_type}); do

    ens_member=${efile:10:5}

    if [ ${file_type} = 'ext' ]; then

       Merge_timeseries ${ens_member} ${file_type} $outt

    elif [ ${file_type} = 'shy' ]; then

       Merge_shy ${ens_member}

    else

       echo "not a valid file"
       exit 1

    fi

done


exit 0
