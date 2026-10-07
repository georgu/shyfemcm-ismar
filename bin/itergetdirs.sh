#!/bin/sh
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# create list of all subdirectories
#
#if [ $# -lt 1 ]; then
#  echo "Usage: itergetdirs.sh [dir(s)]"
#  exit 1
#fi
#
#------------------------------------------------------------------------

home=`pwd`
if [ $# -eq 0 ]; then
  dirs=.
else
  dirs=$*
fi

#------------------------------------------------------------------------

founddirs=""

for dir in $dirs
do
  >&2 echo "starting directory: $dir"
  fdirs=$( find $dir -type d )
  founddirs="$founddirs $fdirs"
done

echo "$founddirs"

#------------------------------------------------------------------------

