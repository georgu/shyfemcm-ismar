#--------------------
# Include shyfem rules and redefines moddir
include ../../Rules.make
MODDIR = .
FGNU_GENERAL = -J$(MODDIR)
FFLAGS        = $(LFLAGS) $(FGNU_NOOPT) $(FGNU_WARNING) $(FGNU_GENERAL)

OBJS = m_mean_preserving_rotation.o analysis.o

# Rules
analysis.o: analysis.F90 m_multa.o m_mean_preserving_rotation.o mod_anafunc.o m_randrot.o
m_mean_preserving_rotation.o: m_mean_preserving_rotation.F90 m_randrot.o

default:
	all

clean:
	rm -fr *.o *.mod

all:
	$(OBJS)

#---------------------
# Compilation rules
.SUFFIXES: .F90
.F90.o:
	$(F95) -c $(FFLAGS) $<
.c.o:
	$(CC) -c $(CFLAGS) $<
#---------------------
