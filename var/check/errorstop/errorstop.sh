#!/bin/sh
#
# tests error stop in fortran
#
#------------------------------------------------------

error=0
. ../../../bin/colors.sh

Run()
{
  local expected=$1
  shift
  local what=$*

  echo "======================="
  echo "running $what..."
  echo "======================="
  ./a.out
  status=$?
  echo "exit code stop: $status"
  [ $status -ne 0 ] && status=1
  if [ $status -ne $expected ]; then
    echo "${red}error in status code: expected $expected got $status${normal}"
    error=$(( error + 1 ))
  else
    echo "${green}status code ok: expected $expected got $status${normal}"
  fi
}

#------------------------------------------------------

option="-fno-backtrace -g"
gfortran -v > /dev/null 2>&1
if [ $? -eq 0 ] ; then
  gfortran $option stop.f90
  Run 0 stop gfortran
  gfortran $option errorstop.f90
  Run 1 errorstop gfortran
fi

option="-g"
ifort -v > /dev/null 2>&1
if [ $? -eq 0 ] ; then
  ifort $option stop.f90
  Run 0 stop intel
  ifort $option errorstop.f90
  Run 1 errorstop intel
fi

echo "-----------------------------"
if [ $error -eq 0 ]; then
  echo "${green}successfull completetion${normal}"
else
  echo "${red}$error errors detected${normal}"
fi
echo "-----------------------------"

#------------------------------------------------------


