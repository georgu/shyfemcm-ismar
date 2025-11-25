#!/bin/sh
#
# creates bfm_strings.f90
#
#------------------------------------------------------

what='bfm'
what='afm'

file=${what}"_strings.f90"
echo "treating $what with file $file"

MODDIR=$HOME/shyfemcm/lib/mod
echo "moddir is $MODDIR"

echo "compiling parse_bfm.f90"
gfortran parse_bfm.f90 convert.f90
[ $? -ne 0 ] && echo "error in compilation" && exit 1

echo "parsing shyfem_vars.py"
a.out < shyfem_vars.py

echo "compiling $file"
gfortran -c -J$MODDIR $file
[ $? -ne 0 ] && echo "error in compilation" && exit 1

echo "successful compilation"

#------------------------------------------------------

