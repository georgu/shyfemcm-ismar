#!/bin/sh
#
# test random generator
#
#----------------------------------------------------

OBJS="random.f90 test_random.f90"

. ../../../bin/colors.sh

#----------------------------------------------------

#==========================================
# gfortran
#==========================================

gfortran -v > /dev/null 2>&1
if [ $? -eq 0 ] ; then
  echo "running with gfortran"
  gfortran $OBJS
  [ $? -ne 0 ] && echo "*** compilation error..." && exit 1
  ./a.out
  mv fort.66 gfortran.tmp
fi

#==========================================
# intel
#==========================================

ifort -v > /dev/null 2>&1
if [ $? -eq 0 ] ; then
  echo "running with ifort"
  ifort $OBJS
  [ $? -ne 0 ] && echo "*** compilation error..." && exit 1
  ./a.out
  mv fort.66 ifort.tmp
fi

#----------------------------------------------------

status=0
if [ -f ifort.tmp -a -f gfortran.tmp ]; then
  echo "comparing results..."
  cmp gfortran.tmp ifort.tmp
  status=$?
elif [ -f ifort.tmp ]; then
  echo "cannot compare results... gfortran compiler missing"
else
  echo "cannot compare results... intel compiler missing"
fi

echo "-----------------------------"
if [ $status -eq 0 ]; then
  echo "${green}successfull completetion${normal}"
else
  echo "${red}$error differences detected${normal}"
fi
echo "-----------------------------"

#----------------------------------------------------

