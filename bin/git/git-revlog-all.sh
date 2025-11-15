#!/bin/bash
#
# checks all files with code
#
#---------------------------------------------------------

bindir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
actdir=$( pwd )

exts="f90 f"
exclude="tmp netcdf"

#--------------------------------------------------------

subdirs=$( iterdir.sh -dir pwd . )

for dir in $subdirs
do
  echo "$dir"
done

for dir in $subdirs
do
  last=$( basename $dir )
  #echo "last: $last"
  [ "$last" = "tmp" ] && continue
  [ "$last" = "netcdf" ] && continue
  [ "$last" = "arc" ] && continue
  cd $dir
  for ext in $exts
  do
    files=$( ls *.$ext 2> /dev/null )
    status=$?
    [ $status -ne 0 ] && continue
    echo "checking extension $ext in $dir"
    $bindir/git-revlog $files
    status=$?
    if [ $status -ne 0 ]; then
      echo "errors found... aborting"
  #    exit 1
    fi
  done
  cd $actdir
done

#---------------------------------------------------------

