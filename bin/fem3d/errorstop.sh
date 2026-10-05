#!/bin/bash
#
# find incomplete error stop specifications
#
#-----------------------------------------------

FindErrorStop()
{
  cd $1

  local files=$( ls *.f90 2>/dev/null )
  [ -z "$files" ] && cd $actdir && return

  grep "error stop" $files | sed -e 's/\.f90: */.f90 /' \
				| grep -v ": " \
				| grep -E -v 'f90\s*\!' \
				| grep -E -v 'f90\s*write'

  echo "---------------------------------------"

  grep "error stop:" $files | sed -e 's/\.f90: */.f90 /' \
				| grep -E -v 'f90\s*\!' \
				| grep -E -v 'f90\s*write'

  echo "---------------------------------------"

  grep -E "\s+stop" $files | grep -v "error *stop" \
				| grep -E -v 'f90\s*\!' \
				| grep -E -v 'f90\s*write'

  cd $actdir
}

#-----------------------------------------------

actdir=$( pwd )

subdirs=$( ~/shyfem/bin/itergetdirs.sh $* )

#-----------------------------------------------

for dir in $subdirs
do
  [[ $dir == */tmp ]] && continue			#skip tmp directories

  echo "======================================="
  echo "$dir"
  echo "======================================="

  FindErrorStop $dir
done

#-----------------------------------------------

