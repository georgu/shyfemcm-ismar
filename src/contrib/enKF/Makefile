
#--------------------
# Include shyfem rules
include ../Rules.make

# Switch to double precision for real
ifeq ($(FORTRAN_COMPILER),GNU_GFORTRAN)
  R8FL = -fdefault-real-8
endif

ifeq ($(FORTRAN_COMPILER),INTEL)
  R8FL = -real-size 64
endif
override FFLAGS += $(R8FL)
override LFLAGS += $(R8FL)


LIBDIR = ../femlib
DIRFEM = ../fem3d
MODDIR = $(LIBDIR)/mod
ENKDIR = EvKF

#--------------------
# ht objects
HTOBJS = $(DIRFEM)/subrst.o $(DIRFEM)/subbas.o $(DIRFEM)/subdts.o $(DIRFEM)/mod_bnd.o \
	$(DIRFEM)/subreg.o $(DIRFEM)/subpar3.o $(DIRFEM)/suboutput.o $(DIRFEM)/subfil.o \
	$(DIRFEM)/mod_hydro.o $(DIRFEM)/debug.o $(DIRFEM)/subsss.o $(DIRFEM)/subbit10.o \
 	$(DIRFEM)/mod_geom_dynamic.o $(DIRFEM)/mod_ts.o $(DIRFEM)/newlevels.o \
	$(DIRFEM)/mod_conz.o $(DIRFEM)/mod_hydro_vel.o $(DIRFEM)/ecological_dummy.o \
	$(DIRFEM)/subdef.o $(DIRFEM)/subnev.o $(DIRFEM)/mod_geom.o $(DIRFEM)/subscn.o

#--------------------
# Ensemble analysis objects
ENSOBJ1 = analysis.o m_randrot.o m_multa.o mod_anafunc.o m_mean_preserving_rotation.o
#ENSOBJ1 = analysis6c.o m_randrot.o m_multa.o 
#ENSOBJ1 = analysis2.o m_randrot.o m_multa.o 
ENSOBJ2 = subenkf.o mod_enkf.o mod_dimensions.o mod_states.o mod_observations.o
ENSOBJ3 = m_sample2D.o m_pseudo2D.o mod_fftw3.o m_newton2D.o \
	m_newtonfunc2D.o m_fixsample2D.o m_random.o

#--------------------
# Libraries to link
LIBS = -llapack -lblas -lfftw3 -lpthread

#--------------------
make_EvKF:
	cd $(ENKDIR); make clean; make; cd -; mv -f $(ENKDIR)/*.o .
	mv -f $(ENKDIR)/*.mod $(MODDIR)

make_sample:
	cd $(ENKDIR)/pseudo_random_fields; make clean; make; cd -
	mv -f $(ENKDIR)/pseudo_random_fields/TMP/*.o .
	mv -f $(ENKDIR)/pseudo_random_fields/TMP/*.mod $(MODDIR)

enkf_analysis: make_EvKF
	make make_sample
	make enkf_analysis.o $(ENSOBJ2) 
	$(LINKER) $(LFLAGS) $(LIBS) $(HTOBJS) $@.o \
        $(ENSOBJ1) $(ENSOBJ2) $(ENSOBJ3) -o $@

make_ens_meteo: make_EvKF
	make make_sample
	make make_ens_meteo.o 
	$(LINKER) $(LFLAGS) $@.o $(LIBS) $(LIBDIR)/libfem.a \
	$(ENSOBJ1) $(ENSOBJ3) -o $@

#--------------------
#--------------------
EXES = enkf_analysis make_ens_meteo

all: 	$(EXES)

default: 
	make all

# I don't know why but this is not necessary:
depend:
	$(FEMDIR)/fembin/mkdp.pl -moddir=$(MODDIR) *.F90

clean:
	rm -fr $(EXES) *.o *.mod tags

cleanall:
	rm -fr $(EXES) *.o *.mod tags *.dat *.uf fort.*; cd $(ENKDIR); rm -fr *.o *.mod; rm -f *.bak


#---------------------
# Compilation rules

.SUFFIXES: .f .F90

.f.o:
	$(F77) -c $(FFLAGS) $<

.F90.o:
	$(F77) -c $(FFLAGS) $< 




# DO NOT DELETE THIS LINE -- make depend depends on it.

../femlib/mod/mod_enkf.mod: mod_enkf.o
../femlib/mod/mod_observations.mod: mod_observations.o
../femlib/mod/mod_states.mod: mod_states.o ../femlib/mod/mod_dimensions.mod
../femlib/mod/mod_dimensions.mod: mod_dimensions.o
../femlib/mod/m_sample2D.mod: m_sample2D.o ../femlib/mod/m_pseudo2D.mod \
		../femlib/mod/m_randrot.mod ../femlib/mod/m_fixsample2D.mod
../femlib/mod/m_fixsample2D.mod: m_fixsample2D.o
../femlib/mod/m_randrot.mod: m_randrot.o

../femlib/mod/m_pseudo2D.mod: m_pseudo2D.o ../femlib/mod/mod_fftw3.mod \
		../femlib/mod/m_newton2D.mod
../femlib/mod/mod_fftw3.mod: mod_fftw3.o
../femlib/mod/m_newton2D.mod: m_newton2D.o ../femlib/mod/m_newtonfunc2D.mod
../femlib/mod/m_newtonfunc2D.mod: m_newtonfunc2D.o

enkf_analysis.o: ../femlib/mod/mod_enkf.mod \
		../femlib/mod/mod_states.mod 
mod_enkf.o: ../femlib/mod/basin.mod ../femlib/mod/mod_conz.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_observations.mod \
		../femlib/mod/mod_states.mod \
		../femlib/mod/mod_ts.mod 
mod_states.o: ../femlib/mod/mod_dimensions.mod
subenkf.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/m_sample2D.mod \
		../femlib/mod/mod_conz.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_restart.mod \
		../femlib/mod/mod_ts.mod ../femlib/mod/shyfile.mod 

make_ens_meteo.o: ../femlib/mod/m_sample2D.mod
