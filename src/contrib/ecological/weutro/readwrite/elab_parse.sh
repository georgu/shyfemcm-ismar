
cp ../weutro.f90 .
cp ../weutro.h .

./parse_weutro.pl weutro.f90
[ $? -ne 0 ] && echo "error parsing weutro.f90" && exit 1

gfortran read_nml.f90
[ $? -ne 0 ] && echo "error running read_nml.f90" && exit 1

a.out

