# This Makefile is not necessary. Compile from the upstram directory.
##---------------------
#F95 = gfortran
F95 = ifort
F77 = $(F95)
FFLAGS = -O -Wall -pedantic

LINKER = $(F95)

CC = gcc
CFLAGS =
#---------------------

m_mean_preserving_rotation.o: m_mean_preserving_rotation.F90 m_randrot.o

analysis.o: analysis.F90 m_multa.o m_mean_preserving_rotation.o mod_anafunc.o m_randrot.o

#---------------------
# Compilation rules
.SUFFIXES: .F90 .f90 .for
.f.o:
	$(F77) -c $(FFLAGS) $<
.for.obj:
	$(FLC) $< , $(<,B).OBJ $(FLFLAGS)
.f90.o:
	$(F95) -c $(FFLAGS) $<
.F90.o:
	$(F95) -c $(FFLAGS) $<
.c.o:
	$(CC) -c $(CFLAGS) $<
#---------------------

clean:
	rm -fr *.o *.mod

all:
	make analysis.o
