#!/bin/bash


Usage()
{
  echo
  echo "Computes n numbers of a random uniform distribution between two extremes and its average value."
  echo "Then, substitutes these numbers to a specified label in a skel-file, creating n skel-files."
  echo
  echo "Usage: randomize_var.sh [val-min] [val-max] [n of values] [skel-file] [label]"
  exit 0
}

#------------------

if [ $5 ]; then
   minval=$1
   maxval=$2
   nval=$3
   skelfile=$4
   string=$5
else
   Usage
fi

[[ ! -s $skelfile ]] && echo "File $skelfile does not exist" && exit 1
skelbase=$(basename $skelfile .skel)

# Average val
av_val=$(echo "scale=6; ($maxval + $minval)/2" | bc)
cat $skelfile | sed -e "s/$string/$av_val/" > ${skelbase}_000.skel

# Ensemble val
for (( n=1; n<=$nval; n++)); do
  nl=$(printf "%03d" $n)

  rand_num=$(echo "scale=6; $RANDOM/32767" | bc)
  val=$(echo "scale=6; $rand_num * ($maxval - $minval) + $minval" | bc)

  cat $skelfile | sed -e "s/$string/$val/" > ${skelbase}_${nl}.skel
done
