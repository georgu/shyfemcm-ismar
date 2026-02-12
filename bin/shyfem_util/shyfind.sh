#!/bin/sh
#
# finds file or pattern in source directories
#
#-----------------------------------------------

base="$HOME/shyfemcm"
base="$SHYFEMDIR"

alldir=$( find $base/src -type d )

actdir=$( pwd )

#-----------------------------------------------

Usage()
{
  echo "Usage: shyfind.sh [-h|-help] [-options] {file|pattern}"
}

FullUsage()
{
  Usage
  
  echo ""
  echo "looks for file names or pattterns in standard locations"
  echo ""
  echo "Available options:"
  echo "  -h|-help         this help"
  echo "  -verbose         be verbose"
  echo "  -quiet           be quiet"
  echo "  -i               ignore case for pattern"
  echo "  -file            look for pattern in file names"
  echo "  -pattern         look for pattern in files"
  echo "  -all             look for pattern in all files (not only f and f90)"
  echo "  -obsolete        look also in obsolete dirs"
  echo "  -full            show also full path with file name"
  echo "  -fullonly        show only full path with file name"
  echo "  -dirs            show directories where to look"
  echo ""
}

SkipObsolete()
{
  echo $dir | grep -q obsolete 

  if [ $? -eq 0 -a $show_obsolete = "NO" ]; then
    return 0
  else
    return 1
  fi
}

Do_File()
{
  echo "looking for file $*"
  for dir in $alldir
  do
    if SkipObsolete; then continue; fi
    cd $dir
    result=$( ls $* 2> /dev/null )
    if [ -n "$result" ]; then
      if [ $full_name = "YES" ]; then
        ls -1 $PWD/$* 2> /dev/null
      else
        echo "---- $dir ----" >&2
        ls -1 $* 2> /dev/null
      fi
    fi
    cd $actdir
  done
}

Do_Pattern()
{
  echo "looking for pattern $*"
  for dir in $alldir
  do
    if SkipObsolete; then continue; fi
    cd $dir
    result=$( grep $options "$*" *.f90 *.f 2> /dev/null )
    if [ -n "$result" ]; then
      if [ $full_name = "YES" ]; then
        if [ $full_only = "NO" ]; then
          echo "---- $dir ----" >&2
          grep $options "$*" *.f90 *.f 2> /dev/null
	fi
        grep -l $options "$*" $PWD/*.f90 $PWD/*.f 2> /dev/null
      else
        echo "---- $dir ----" >&2
        grep "$*" *.f90 *.f 2> /dev/null
      fi
    fi
    cd $actdir
  done
}

Do_Pattern_All()
{
  echo "looking for pattern $* in all files"
  for dir in $alldir
  do
    if SkipObsolete; then continue; fi
    cd $dir
    result=$( grep $options "$*" * 2> /dev/null )
    if [ -n "$result" ]; then
      if [ $full_name = "YES" ]; then
        if [ $full_only = "NO" ]; then
          echo "---- $dir ----" >&2
          grep $options "$*" * 2> /dev/null
	fi
        grep -l $options "$*" $PWD/* 2> /dev/null
      else
        echo "---- $dir ----" >&2
        grep "$*" * 2> /dev/null
      fi
    fi
    cd $actdir
  done
}

#-----------------------------------------------

do_file=NO
do_pattern=NO
do_pattern_all=NO
do_dirs=NO
quiet=NO
verbose=NO
show_obsolete="NO"
full_name="NO"
full_only="NO"
options=""

while [ -n "$1" ]
do
   case $1 in
        -quiet)         quiet="YES";;
        -verbose)       verbose="YES";;
        -file)          do_file="YES";;
        -i)             options="$options -i";;
        -pattern)       do_pattern="YES";;
        -all)           do_pattern_all="YES";;
        -obsolete)      show_obsolete="YES";;
        -full)          full_name="YES";;
        -fullonly)      full_only="YES";full_name="YES";;
        -dirs)          do_dirs="YES";;
        -h|-help)       FullUsage; exit 0;;
        -*)             echo "No such option: $1"; exit 1;;
        *)              break;;
   esac
   shift
done

if [ $do_dirs = "YES" ]; then
  echo $alldir
  exit 0
fi

[ $# -eq 0 ] && Usage && exit 0

echo "looking for $*"

#-----------------------------------------------

if [ $do_file = "YES" ]; then
  Do_File $*
elif [ $do_pattern_all = "YES" ]; then
  alldir=$( find $base -type d )
  Do_Pattern_All $*
elif [ $do_pattern = "YES" ]; then
  Do_Pattern $*
else
  Do_File $*
  Do_Pattern $*
fi

#-----------------------------------------------

