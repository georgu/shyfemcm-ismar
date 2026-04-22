#!/bin/bash
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2026  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# sets new path with shyfemcm
#
# "shypath"        just shows actual path
# "shypath new"    sets new path to shyfemcm
# "shypath all"    sets new path to shyfemcm and all subdirs
# "shypath [dirs]" sets new path to shyfemcm and subdirs dirs
# "shypath old"    resets path to old shyfem distribution
#
# in order to work this script must be sourced (called as ". shypath.sh"
#
#-----------------------------------------------------------

all="git fem3d grd shyfem_util gps"

has_dot="NO"
set_path="NO"
debug="YES"
debug="NO"

quiet="NO"
[ "$1" = "-quiet" ] && quiet="YES"
[ "$1" = "-silent" ] && quiet="YES"
[ $quiet = "YES" ] && shift

#-----------------------------------------------------------

AddPath()
{
  for p in $*
  do
    path="$path:$bindir/$p"
  done
}

AddAbsolutePath()
{
  for p in $*
  do
    path="$path:$p"
  done
}

CleanPath()
{
  newpath=""
  declare -A paths

  for p in $*
  do
    if [[ $p == *shyfem* ]]; then
      continue
    elif [[ $p == *games* ]]; then
      continue
    elif [[ $p == '.' ]]; then
      has_dot="YES"
      continue
    fi
    
    if [ ${paths[$p]+_} ]; then
      [ $debug = "YES" ] && echo "$p already inserted"
      :
    else 
      [ $debug = "YES" ] && echo "$p not yet inserted"
      newpath="$newpath $p"
      paths[$p]="YES"
    fi
  done
}

MakeBinDir()
{
  # in order for this to work we must be in the directory tree of the repo

  local execdir=$( dirname ${BASH_SOURCE[0]} | pwd -P )
  local dir=$execdir

  [ $quiet = "NO" ] && echo "execdir = $execdir"

  while [ -n "$dir" ]
  do
    #echo "looking for bin in $dir"
    if [ -d $dir/bin ]; then
      bindir=$dir/bin
      #echo "bindir found: $bindir"
      return
    else
      dir=$( dirname $dir )
    fi
  done

  bindir="/home/shyfem/bin"		# last resort
}

#-----------------------------------------------------------

actdir=$( pwd )
bindir=
execdir=$( dirname ${BASH_SOURCE[0]} | pwd -P )
#echo "execdir = $execdir"
MakeBinDir
bindir=$SHYFEMDIR/bin

shyfemdir=$( echo $bindir | sed -E 's/\/[^\/]*$//' )

path=$( echo $PATH | tr ':' ' ' )
CleanPath $path
path=$( echo $newpath | tr ' ' ':' )

[ $debug = "YES" ] && echo "bindir: $bindir"

if [ $quiet = "NO" ]; then
  #echo "actdir: $actdir"
  echo "bindir: $bindir"
  echo "shyfemdir: $shyfemdir"
  echo "SHYFEMDIR: $SHYFEMDIR"
fi

if [ $# -eq 0 ]; then
  echo "PATH: $PATH"
elif [ $1 = "all" ]; then
  set_path="YES"
  AddAbsolutePath $bindir
  AddPath $all
elif [ $1 = "new" ]; then
  set_path="YES"
  AddAbsolutePath $bindir
elif [ $1 = "old" ]; then
  set_path="YES"
  AddAbsolutePath $HOME/shyfem/bin
  shyfemdir=$HOME/shyfem
else
  set_path="YES"
  AddAbsolutePath $bindir
  AddPath $*
fi

if [ $has_dot = "YES" ]; then
  path="$path:."
fi

if [ $set_path = "YES" ]; then
  [ $quiet = "NO" ] && echo "new PATH: $path"
  export PATH=$path
  export SHYFEMDIR=$shyfemdir
fi

#cd $actdir

#-----------------------------------------------------------

