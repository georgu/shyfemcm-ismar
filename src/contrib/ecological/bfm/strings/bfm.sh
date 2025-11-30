#!/bin/sh
#
# creates bfm_strings.f90
#
#------------------------------------------------------

if [ "$1" = "-nocheck" ]; then
  options="-nocheck"
  shift
elif [ "$1" = "-check" ]; then
  options="-check"
  shift
else
  options="-check"
fi

if [ $# -ne 1 ]; then
  echo "Usage: ./bfm.sh [-check|-nocheck] what"
  echo "  what can be var_num_2_name cbms bfm afm bfm_d afm_d"
  exit 1
fi

#------------------------------------------------------

what=$1
file="${what}_strings.f90"
echo "creating file $file"

MODDIR=$HOME/shyfemcm/lib/mod
echo "moddir is $MODDIR"

echo "compiling parse_bfm.f90"
make compile
[ $? -ne 0 ] && echo "error in compilation" && exit 1

echo "parsing shyfem_vars.py"
a.out $options $what < shyfem_vars.py
status=$?
#echo "status = $status"
if [ $status -eq 77 ]; then
  echo "*** error stop: entries not unique"
  echo "(to avoind this error run with -nocheck)"
  exit 77
fi
[ $status -ne 99 ] && echo "error in parsing" && exit 3

echo "compiling $file"
gfortran -c -J$MODDIR $file
[ $? -ne 0 ] && echo "error in compilation" && exit 5

echo "successful compilation"
echo "file with strings is in $file"

#------------------------------------------------------

