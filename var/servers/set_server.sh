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
# this script handles server specific settings
#
#--------------------------------------------------------------

debug="YES"
debug="NO"

#--------------------------------------------------------------

FindServer()
{
  host=$( hostname )

  if [ "$host" = galileo -o "$HPC_SYSTEM" = galileo ]; then
    server=galileo
  elif [ "$host" = fluxus -o "$SGE_CLUSTER_NAME" = fluxus ]; then
    server=fluxus
  else
    server=$host
  fi
}

#--------------------------------------------------------------

Usage()
{
  echo "Usage: set_server.sh [-h|-help] host"
}

FullUsage()
{
  Usage
}

ParseOptions()
{
  what="set"
  while [ -n "$1" ]
  do
    case $1 in
        -server)        what="server";;
        -show)          what="show";;
        -load)          what="load"; setting=$2; shift;;
        -h|-help)       FullUsage; exit 0;;
        -*)             echo "no such option: $1"; exit 1;;
        *)              break;;
    esac
    shift
  done

  compiler=$1
}

Sourced()
{
  [ $what = "show" ] && return

  if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "script is not sourced"
    echo "please run this script as \". ${0}\""
  fi
}

SetVars()
{
  [ $what = "show" ] && return

  if [ $server = "Caesium" ]; then
    export NETCDF_C_HOME=
    export NETCDF_FORTRAN_HOME=
    export NETCDFDIR=
    export METIS_HOME=$HOME/lib/metis
    export PARMETIS_HOME=
    export PETSC_HOME=$HOME/lib/petsc
    export BFM_HOME=$HOME/work/shyfem_repo/bfm-for-shyfem-private
  elif [ $server = "lagoon" ]; then
    export NETCDF_C_HOME=
    export NETCDF_FORTRAN_HOME=
    export NETCDFDIR=
    export METIS_HOME=$HOME/lib/metis
    export PARMETIS_HOME=
    export PETSC_HOME=$HOME/lib/petsc
    export BFM_HOME=$HOME/lib/bfm
  elif [ $server = "tide" ]; then
    export NETCDF_C_HOME=
    export NETCDF_FORTRAN_HOME=
    export NETCDFDIR=
    export METIS_HOME=$HOME/georg/lib/metis
    export PARMETIS_HOME=
    export PETSC_HOME=$HOME/georg/lib/petsc
    export BFM_HOME=$HOME/georg/lib/bfm
  elif [ $server = "storm" ]; then
    export NETCDF_C_HOME=
    export NETCDF_FORTRAN_HOME=
    export NETCDFDIR=
    export METIS_HOME=$HOME/georg/lib/metis
    export PARMETIS_HOME=$HOME/georg/lib/parmetis
    export PETSC_HOME=$HOME/georg/lib/petsc
    export BFM_HOME=$HOME/georg/lib/bfm
  elif [ $server = "stream" ]; then
    export NETCDF_C_HOME=
    export NETCDF_FORTRAN_HOME=
    export NETCDFDIR=
    export METIS_HOME=$HOME/georg/lib/metis
    export PARMETIS_HOME=$HOME/georg/lib/parmetis
    export PETSC_HOME=$HOME/georg/lib/petsc
    export BFM_HOME=$HOME/georg/lib/bfm
  elif [ $server = "nemunas.ku.lt" ]; then
    export NETCDF_C_HOME=
    export NETCDF_FORTRAN_HOME=
    export NETCDFDIR=/opt/sw/
    export METIS_HOME=$HOME/georg/lib/metis
    export PARMETIS_HOME=$HOME/georg/lib/parmetis
    export PETSC_HOME=$HOME/georg/lib/petsc
    export BFM_HOME=$HOME/georg/lib/bfm
  else
    echo "unknown server: $server"
  fi
}

#--------------------------------------------------------------

ParseOptions $*
FindServer

Sourced

SetVars

if [ $what = "show" ]; then
    echo " server = $server"
    echo " NETCDF_C_HOME = $NETCDF_C_HOME"
    echo " NETCDF_FORTRAN_HOME = $NETCDF_FORTRAN_HOME"
    echo " NETCDFDIR = $NETCDFDIR"
    echo " METIS_HOME = $METIS_HOME"
    echo " PARMETIS_HOME = $PARMETIS_HOME"
    echo " PETSC_HOME = $PETSC_HOME"
    echo " BFM_HOME = $BFM_HOME"
fi

#--------------------------------------------------------------

