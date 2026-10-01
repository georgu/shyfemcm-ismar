
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------


DIR = dstruc

TOTAL = dsall.f

#------------------------------------------------------

EXES = stack queue hash tree sort spline heap iheap
F77 = gfortran

FFLAGS = -fcheck=all
FFLAGS = -fbounds-check -g -fbacktrace

#------------------------------------------------------

default: all

all: $(EXES)

stack: stack.o
	$(F77) -o $@ $@.f90

queue: queue.o
	$(F77) -o $@ $@.f90

hash: hash.o
	$(F77) -o $@ $@.f90

tree: tree.o
	$(F77) -o $@ $@.f90

sort: sort.o
	$(F77) -o $@ $@.f90

spline: spline.o
	$(F77) -o $@ $@.f90

heap: heap.o
	$(F77) -o $@ $@.f90

iheap: iheap.o
	$(F77) -o $@ $@.f90

#------------------------------------------------------

test: all
	@for prog in $(EXES); do ./$$prog; done

compact:
	@rm -f $(TOTAL)
	@for prog in $(EXES); do ./compact.pl $$prog.f >> $(TOTAL) ; done

#------------------------------------------------------

clean:
	-rm -f *.o
	-rm -f *.mod
	-rm -f *.tmp
	-rm -f *.bak
	-rm -f a.out
	-rm -f $(EXES)
	-rm -f ggg hhh
	-rm -f out*.ps
	-rm -f fort.*

cleanzip: clean
	-rm -f $(DIR).zip

cleanall: clean cleanzip
	-rm -f $(TOTAL)

save: zip
zip: cleanzip
	zip $(DIR) *

copy:
	cp $(DIR).zip $(HOME)/C

#------------------------------------------------------

.f.o:
	$(F77) -c $(FFLAGS) $<

.f90.o:
	$(F77) -c $(FFLAGS) $<

%.o: %.f90
	$(F77) $(FFLAGS) -c $<

#------------------------------------------------------

