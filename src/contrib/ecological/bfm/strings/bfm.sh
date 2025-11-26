#!/bin/sh
#
# creates bfm_strings.f90
#
#------------------------------------------------------

file="new_strings.f90"
echo "creating file $file"

MODDIR=$HOME/shyfemcm/lib/mod
echo "moddir is $MODDIR"

echo "compiling parse_bfm.f90"
gfortran parse_bfm.f90 convert.f90
[ $? -ne 0 ] && echo "error in compilation" && exit 1

echo "parsing shyfem_vars.py"
a.out < shyfem_vars.py
[ $? -ne 99 ] && echo "error in parsing" && exit 3

echo "compiling $file"
gfortran -c -J$MODDIR $file
[ $? -ne 0 ] && echo "error in compilation" && exit 5

echo "successful compilation"

#------------------------------------------------------

