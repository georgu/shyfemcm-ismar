#!/bin/sh
#
# extracts horizontal limits from fem file and writes it to grd file
#
#----------------------------------------------------

file=$1
[ $# -ne 1 ] && echo "Usage: femlimit.sh fem-file" && exit 1
[ ! -f $file ] && echo "*** no such file: $file" && exit 1

xy0=$( femelab $file | grep 'x0,y0:' | sed -E 's/^.*: +//' )
xy1=$( femelab $file | grep 'x1,y1:' | sed -E 's/^.*: +//' )

echo "xy0: $xy0"
echo "xy1: $xy1"

x0=$( echo $xy0 | sed -E 's/ +.*$//' )
y0=$( echo $xy0 | sed -E 's/^.* +//' )
x1=$( echo $xy1 | sed -E 's/ +.*$//' )
y1=$( echo $xy1 | sed -E 's/^.* +//' )

echo "x0=$x0  y0=$y0"
echo "x1=$x1  y1=$y1"

echo "0 limits of file $file"		 > limit.grd
echo "1 1 0 $x0 $y0"			>> limit.grd
echo "1 2 0 $x1 $y0"			>> limit.grd
echo "1 3 0 $x1 $y1"			>> limit.grd
echo "1 4 0 $x0 $y1"			>> limit.grd
echo "3 1 0 5 1 2 3 4 1"		>> limit.grd

echo "limits have been written to file limit.grd"

#----------------------------------------------------

