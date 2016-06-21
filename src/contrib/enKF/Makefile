
#--------------------
# Include shyfem rules
include ../Rules.make

DIRFEM    = $(FEMDIR)/fem3d
ENKDIR    = EvKF
MODDIR    = $(FEMDIR)/femlib/mod

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
ENSOBJ = analysis.o m_randrot.o m_multa.o mod_anafunc.o m_mean_preserving_rotation.o

#--------------------
make_EvKF:
	cd $(ENKDIR); make; mv -f *.o ..; mv -f *.mod ../$(MODDIR)

enKF_analysis: make_EvKF
	make enKF_analysis.o mod_enKF.o $(ENSOBJ) $(HTOBJS)
	$(LINKER) $(LFLAGS) -llapack -lblas $@.o mod_enKF.o $(ENSOBJ) $(HTOBJS) -o $@

#--------------------
#--------------------
EXES = enKF_analysis 

all:
	make $(EXES)

default:
	all

# I don't know why but this is not necessary:
depend:
	$(FEMDIR)/fembin/mkdp.pl -moddir=$(MODDIR) *.f90

clean:
	rm -fr $(EXES) *.o *.mod tags

cleanall:
	rm -fr $(EXES) *.o *.mod tags *.dat *.uf fort.*; cd $(ENKDIR); rm -fr *.o *.mod 


#---------------------
# Compilation rules
.SUFFIXES: .f90 .f

.f.o:
	$(F77) -c $(FFLAGS) $<

.f90.o:
	$(F77) -c $(FFLAGS) $< 


# DO NOT DELETE THIS LINE -- make depend depends on it.
