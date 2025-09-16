#!/bin/sh
#
# creates bfm_strings.f90
#
#------------------------------------------------------

MODDIR=$HOME/shyfemcm/lib/mod
echo "moddir is $MODDIR"

echo "compiling parse_bfm.f90"
gfortran parse_bfm.f90 convert.f90
[ $? -ne 0 ] && echo "error in compilation" && exit 1

echo "parsing shyfem_vars.py"
a.out < shyfem_vars.py

echo "compiling bfm_strings.f90"
gfortran -c -J$MODDIR bfm_strings.f90 
[ $? -ne 0 ] && echo "error in compilation" && exit 1

echo "successful compilation"

#------------------------------------------------------

