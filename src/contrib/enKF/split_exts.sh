#/bin/bash
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

    rm -f zeta.2d.* velx.2d.* vely.2d.* temp.2d.* salt.2d.*
    $FEMDIR/fembin/shyelab -split $efile > ext.log

    for fl in $(ls zeta.2d.*); do
        mv -f $fl ${basefile}_${fl}
    done
    for fl in $(ls velx.2d.*); do
        mv -f $fl ${basefile}_${fl}
    done
    for fl in $(ls vely.2d.*); do
        mv -f $fl ${basefile}_${fl}
    done
    if [ -e temp.2d.1 ]; then
     for fl in $(ls temp.2d.*); do
        mv -f $fl ${basefile}_${fl}
     done
    fi
    if [ -e salt.2d.1 ]; then
     for fl in $(ls salt.2d.*); do
        mv -f $fl ${basefile}_${fl}
     done
    fi
done
rm -f ext.log
rm -f zeta.2d.* velx.2d.* vely.2d.* temp.2d.* salt.2d.*
