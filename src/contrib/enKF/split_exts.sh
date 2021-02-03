#/bin/bash
#
# Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
#
# Just split ext files 
#----------------------------------------------------------

# This finds the path of the current script
SCRIPT=$(realpath $0)
SCRIPTPATH=$(dirname $SCRIPT)

FEMDIR=$SCRIPTPATH/..   # fem directory

Usage()
{
  echo "Usage: split_exts.sh [basename]"
  exit 0
}

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------

if [ $1 ]; then
   fbase=$1
else
   Usage
fi

for efile in $(ls ${fbase}*.ext); do
    echo "file: $efile"
    basefile=$(basename $efile .ext)
    vars='zeta.2d velx.2d vely.2d velx.3d vely.3d speed.2d speed.3d dir.2d dir.3d all.2d temp.2d temp.3d salt.2d salt.3d'

    $FEMDIR/fembin/shyelab -split $efile > ext.log
    rm -f ${basefile}_*.ts
    for vv in $vars; do
	if [ -e "$vv.1" ]; then
           for fl in $(ls $vv.*); do
               mv -f $fl ${basefile}_${fl}.ts
           done
	fi
    done
done
rm -f ext.log
