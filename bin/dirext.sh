#!/bin/bash
#
# lists files ordered by extension
#
#------------------------------------------------------------

bindir=~/shyfemcm/bin
#. $bindir/bashutil.sh

files=$( ls )

declare -A exts
declare -A count

for file in $files
do
  [ -d $file ] && continue
  filename=$(basename -- "$file")
  extension="${filename##*.}"
  [[ "$filename" != *.* ]] && extension="none"
  filename="${filename%.*}"
  #echo $filename $extension
  exts[$extension]="${exts[$extension]} $file"
  count[$extension]=$(( ${count[$extension]} + 1 ))
done

nc=0
for key in "${!exts[@]}"
do
  #echo "$key: ${count[$key]}"
  echo "  ${exts[$key]}"
  nc=$(( nc + 1 ))
done

>&2 echo "$nc different extensions found"

for key in "${!exts[@]}"
do
  >&2 echo "$key: ${count[$key]}"
done

#------------------------------------------------------------

