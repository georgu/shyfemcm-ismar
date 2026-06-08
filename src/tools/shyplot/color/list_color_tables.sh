#!/bin/bash
#
#----------------------------------------------

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

