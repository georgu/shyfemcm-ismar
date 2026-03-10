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
# creates *.tex documentation files from *.f files

FEMDIR=..
SHYDIR=$FEMDIR/src/shyfem
BIODIR=$FEMDIR/src/contrib/ecological/weutro
UTIDIR=$FEMDIR/src/utils/shyutil
UTMDIR=$FEMDIR/src/utils/shyutilmpi
TOLDIR=$FEMDIR/src/tools/shyplot

Proc()
{
  if [ ! -f $1 ]; then
    echo "*** no such file: $1\n";
    exit 3
  fi
  ./femdoc.pl $1
  if [ $? -ne 0 ]; then
    echo "*** (femdoc.pl) Error in processing file $1"
    exit 1
  fi
}

#---------------------------------------------------------

Proc $SHYDIR/def_para.f90
Proc $TOLDIR/def_para_post.f90
Proc $UTMDIR/bnd_admin.f90
Proc $SHYDIR/meteo_forcing.f90
Proc $SHYDIR/chezy.f90
Proc $UTIDIR/version.f90
Proc $BIODIR/bio3d.f90
Proc $SHYDIR/sedim_admin.f90
Proc $SHYDIR/waves_admin.f90
Proc $SHYDIR/tidef.f90

#---------------------------------------------------------

