# Include shyfem rules and redefines moddir
include ../../Rules.make

MODDIR = .

# Switch to double precision for real
ifeq ($(FORTRAN_COMPILER),GNU_GFORTRAN)
  R8FL = -fdefault-real-8
endif

ifeq ($(FORTRAN_COMPILER),INTEL)
  R8FL = -real-size 64
endif
override FFLAGS += $(R8FL)
override LFLAGS += $(R8FL)


OBJS = analysis.o

# Rules
analysis.o: analysis.F90 m_multa.o m_mean_preserving_rotation.o mod_anafunc.o m_randrot.o
m_mean_preserving_rotation.o: m_mean_preserving_rotation.F90 m_randrot.o

default:
	all

clean:
	rm -fr *.o *.mod

all: 	$(OBJS)

#---------------------
# Compilation rules
.SUFFIXES: .F90
.F90.o:
	$(F95) -c $(FFLAGS) $<
.c.o:
	$(CC) -c $(CFLAGS) $<
#---------------------
