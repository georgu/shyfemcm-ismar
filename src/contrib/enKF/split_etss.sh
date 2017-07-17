#/bin/bash
# Just split ets files 
#----------------------------------------------------------

# This finds the path of the current script
SCRIPT=$(realpath $0)
SCRIPTPATH=$(dirname $SCRIPT)

FEMDIR=$SCRIPTPATH/..   # fem directory

Usage()
{
  echo "Usage: split_etss.sh [basename]"
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

for efile in $(ls ${fbase}*.ets); do
    echo "file: $efile"
    basefile=$(basename $efile .ets)

    rm -f z.* u.* v.* t.* s.*
    memory -s $efile > ets.log
    $FEMDIR/fembin/splitets >> ets.log

    for fl in $(ls z.*); do
        mv -f $fl ${basefile}_${fl}
    done
    for fl in $(ls u.*); do
        mv -f $fl ${basefile}_${fl}
    done
    for fl in $(ls v.*); do
        mv -f $fl ${basefile}_${fl}
    done
    for fl in $(ls t.*); do
        mv -f $fl ${basefile}_${fl}
    done
    for fl in $(ls s.*); do
        mv -f $fl ${basefile}_${fl}
    done
done
rm -f ets.log
rm -f z.* u.* v.* t.* s.*
