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
   vars='zeta.2d velx.2d vely.2d velx.3d vely.3d speed.2d speed.3d dir.2d dir.3d all.2d temp.2d temp.3d salt.2d salt.3d'
   rm -f *_st*_en${nen}.ts
   for fil in $files; do
      echo "Processing file: $fil"
      $FEMDIR/fembin/shyelab -split ${fil} > log
       
      for vv in $vars; do
	if [ -e "$vv.1" ]; then
          for flev in $(ls $vv.*); do
	    idst=$(echo $flev | cut -d '.' -f 3)
            cat $flev |head -n -1|tail -n +2 >> ${vv}_st${idst}_en${nen}.ts
	    rm $flev
          done
	fi
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
