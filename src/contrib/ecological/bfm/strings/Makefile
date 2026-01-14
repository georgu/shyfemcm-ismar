
OBJS = convert.f90 clo.f90

default:

bfm: run
run:
	./bfm.sh

compile:
	gfortran $(OBJS) parse_bfm.f90

clean:
	-rm -f *.o *.mod a.out
	-rm -f *.tmp

cleanall: clean
	-rm -f *_strings.f90
	-rm -f strings_*.f90

