#!/bin/bash
#
#--------------------------------------------------------------------------
#
#    Copyright (C) 1985-2026  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#--------------------------------------------------------------------------
#
# checks all files with code
#
#---------------------------------------------------------

bindir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
actdir=$( pwd )

exts="f90 f sh pl"
exclude="tmp netcdf arc"
exclude="tmp arc"

#--------------------------------------------------------

#subdirs=$( iterdir.sh -dir pwd . )
subdirs=$( find . -type d )

for dir in $subdirs
do
  last=$( basename $dir )
  #echo "last: $last"
  for excl in $exclude
  do
    [ "$last" = $excl ] && continue 2
  done
  cd $dir
  for ext in $exts
  do
    files=$( ls *.$ext 2> /dev/null )
    status=$?
    [ $status -ne 0 ] && continue
    echo "checking extension $ext in $dir"
    $bindir/git-revlog -nowrite $files
    status=$?
    if [ $status -ne 0 ]; then
      :
      #echo "errors found... aborting"
      #exit 1
    fi
  done
  cd $actdir
done

#---------------------------------------------------------

