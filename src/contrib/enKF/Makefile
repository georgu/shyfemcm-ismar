
#--------------------
# Include shyfem rules
include ../Rules.make

# Reset some env vars
DIRFEM    = $(FEMDIR)/fem3d
ENKDIR    = EvKF

#--------------------
# ht objects
HTOBJS = $(DIRFEM)/subrst.o $(DIRFEM)/subbas.o $(DIRFEM)/subdts.o $(DIRFEM)/mod_bnd.o \
	$(DIRFEM)/subreg.o $(DIRFEM)/subpar3.o $(DIRFEM)/suboutput.o $(DIRFEM)/subfil.o \
	$(DIRFEM)/mod_hydro.o $(DIRFEM)/debug.o $(DIRFEM)/subsss.o $(DIRFEM)/subbit10.o \
 	$(DIRFEM)/mod_geom_dynamic.o $(DIRFEM)/mod_ts.o $(DIRFEM)/newlevels.o \
	$(DIRFEM)/mod_conz.o $(DIRFEM)/mod_hydro_vel.o $(DIRFEM)/ecological_dummy.o \
	$(DIRFEM)/subdef.o $(DIRFEM)/subnev.o $(DIRFEM)/mod_geom.o $(DIRFEM)/subscn.o

#--------------------
# Ensemble analysis modules
ENSOBJ = m_randrot.o m_multa.o mod_anafunc.o m_mean_preserving_rotation.o m_random2.o
#mod_anafunc.o: m_mean_preserving_rotation.o
#analysis.o: $(ENSOBJ)

#--------------------
#mod_enKF.o: mod_enKF.mod
#mod_enKF.mod: mod_enKF.o
#enKF_analysis.o: mod_enKF.mod m_randrot.mod m_multa.mod mod_anafunc.mod \
	m_mean_preserving_rotation.mod m_random2.mod

enKF_analysis: enKF_analysis.o mod_enKF.o analysis.o $(ENSOBJ) $(HTOBJS)
	$(LINKER) $(LFLAGS) $@.o mod_enKF.o analysis.o $(ENSOBJ) $(HTOBJS) \ 
	-o $@

#--------------------
#--------------------
EXES = enKF_analysis 

default:
	all

all:
	make $(EXES)

clean:
	rm -fr $(EXES) *.o *.mod tags

cleanall:
	rm -fr $(EXES) *.o *.mod tags *.dat *.uf fort.*; cd $(ENKDIR); rm -fr *.o *.mod 


#---------------------
# Compilation rules
.SUFFIXES: .F90 .f90

.f.o:
	$(F77) -c $(FFLAGS) $<

.f90.o:
	$(F77) -c $(FFLAGS) $< 

.F90.o:
	$(F77) -c $(FFLAGS) $<
