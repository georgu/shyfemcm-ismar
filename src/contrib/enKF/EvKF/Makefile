MODDIR = .

OBJS = analysis.o analysis6c.o analysis2.o


default:	all

all: 	$(OBJS)


clean:
	rm -fr *.o *.mod


# Rules
analysis2.o: analysis2.F90 m_multa.o
analysis6c.o: analysis6c.F90 m_multa.o m_randrot.o
analysis.o: analysis.F90 m_multa.o m_mean_preserving_rotation.o mod_anafunc.o m_randrot.o
m_mean_preserving_rotation.o: m_mean_preserving_rotation.F90 m_randrot.o


#---------------------
# Compilation rules
.SUFFIXES: .F90
.F90.o:
	$(F95) -c $(FFLAGS) $<
.c.o:
	$(CC) -c $(CFLAGS) $<
#---------------------
