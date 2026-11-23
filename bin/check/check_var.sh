#!/bin/bash
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# checks various things
#
#--------------------------------------------------

except="Makefile INFO_INSTALL README CR copyright_notice.txt"

CheckExe()
{
  pushd $1 > /dev/null
  if [ $? -ne 0 ]; then
    echo "*** cannot cd to directory $1 ...aborting"
    return
  else
    echo "checking directory `pwd`"
  fi

  files=`ls`

  for file in $files
  do
    [ -x $file ] && continue
    [ -d $file ] && continue
    is_special="NO"
    for special in $except
    do
      [ $file = $special ] && is_special="YES"
    done
    if [ $is_special = "NO" ]; then
      echo "*** file is not executable: $file"
    fi
  done

  popd > /dev/null
}

#--------------------------------------------------

if [ ! -f VERSION ]; then
  echo "must be in base directory to run this script... aborting"
  exit 1
fi

CheckExe bin
CheckExe bin/check

#--------------------------------------------------

