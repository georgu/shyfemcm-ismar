#!/bin/bash
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

ctables=""

centries=$( grep _data colormap.dat )

for entry in $centries
do
  if [[ $entry == _* ]]; then
    #echo $entry
    name=$( echo $entry | sed -e 's/_data//' | sed -e 's/_//' ) 
    echo $name
  fi
done

#----------------------------------------------

