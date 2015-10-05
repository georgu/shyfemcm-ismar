#!/bin/ksh
#Extract module dependencies
#KAL -- * is also comment mark ....
#egrep -v '(^C|^c|^!)' MODEL.CPP *.F *.F90 |\
egrep -v '(^\*|^C|^c|^!)' MODEL.CPP *.F *.F90 |\
sed -e 's/.*:#/#/' | g++ -P |\
egrep -i '\:[ ]*use [      ]*' |\
sed  -e 's/\.F90:/.o/' -e 's/\.F:/.o/'  -e 's/use//g'  -e 's/USE//g' -e 's/mydim/modules/g' |\
awk '{printf("%s:\t%s.o\n",$1,$2);}' 

#Extract include dependencies
egrep -v '(^C|^c|^!)' MODEL.CPP *.H *.F *.F90 |\
#sed -e 's/.*:#/#/' | /usr/lib/cpp -P |\
sed -e 's/.*:#/#/' | g++ -P |\
egrep -i '\:[ ]*include [      ]*' |\
sed  -e 's/\.F90:/.o/' -e 's/\.F:/.o/' -e 's/\"//g' -e s/\'//g -e '/mpif\.h/d'  |\
awk '{printf("%s:\t%s\n",$1,$3);}'


#Add include dependencies   .h: .H for preprocessor
ls *.H  |\
awk -F. '{
    CC[NR]=$1
}
END {
    for (i = 1; i <= NR; i++) printf("%s.h:\t%s.H\n", CC[i], CC[i])
}'

grep '#ifdef' *.F *.F90 *.H | cut -f1 -d: | sort -u | sed -e 's/.F90/.o/' -e 's/.F/.o/' -e 's/.H/.h/' |\
awk '{printf("%s:\tMODEL.CPP\n",$1);}'

