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

exts="f90 f F90 sh pl"
exclude="tmp arc obsolete"

option="-nowrite"
quiet="NO"
quiet="YES"

#--------------------------------------------------------

SearchSubdirs()
{
  local subdirs=$*

  local actdir=$( pwd )
  local dir last files

  for dir in $subdirs
  do
    last=$( basename $dir )
    for excl in $exclude
    do
      [ "$last" = $excl ] && continue 2
      [[ "$dir" == *"/$excl/"* ]] && continue 2
    done
    local options="$option -dir $dir"
    cd $dir
    for ext in $exts
    do
      files=$( ls *.$ext 2> /dev/null )
      status=$?
      [ $status -ne 0 ] && continue
      [ $quiet = "NO" ] && echo "checking extension $ext in $dir"
      $bindir/git-revlog $options $files
      status=$?
    done
    cd $actdir
  done
}

#---------------------------------------------------------

subdirs=$( find . -type d )

SearchSubdirs $subdirs

#---------------------------------------------------------

