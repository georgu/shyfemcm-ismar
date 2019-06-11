#
# Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
#
#--------------------
# Include shyfem rules

include ../Rules.make

# Switch to double precision for real
ifeq ($(FORTRAN_COMPILER),GNU_GFORTRAN)
  R8FL = -fdefault-real-8 -fdefault-double-8
endif

ifeq ($(FORTRAN_COMPILER),INTEL)
  R8FL = -real-size 64 -double-size 64
endif
override FFLAGS += $(R8FL)
override LFLAGS += $(R8FL)


LIBDIR = ../femlib
DIRFEM = ../fem3d
MODDIR = $(LIBDIR)/mod

#--------------------
# ht objects
HTOBJS = $(DIRFEM)/subrst.o $(DIRFEM)/subbas.o $(DIRFEM)/subdts.o $(DIRFEM)/mod_bnd.o \
	$(DIRFEM)/subreg.o $(DIRFEM)/subfind.o $(DIRFEM)/subpar3.o $(DIRFEM)/suboutput.o \
	$(DIRFEM)/subfil.o $(DIRFEM)/mod_hydro.o $(DIRFEM)/debug.o $(DIRFEM)/subsss.o \
	$(DIRFEM)/subbit10.o $(DIRFEM)/mod_geom_dynamic.o $(DIRFEM)/mod_ts.o \
	$(DIRFEM)/newlevels.o $(DIRFEM)/mod_conz.o $(DIRFEM)/mod_hydro_vel.o \
	$(DIRFEM)/ecological_dummy.o $(DIRFEM)/subdef.o $(DIRFEM)/subnev.o \
	$(DIRFEM)/mod_geom.o $(DIRFEM)/subscn.o $(DIRFEM)/suboutputd.o \
	$(DIRFEM)/subiso8601.o

#--------------------
# Ensemble analysis objects
ENSOBJ1 = mod_para.o subenkf.o mod_ens_state.o mod_mod_err.o mod_enkf.o \
	mod_init_enkf.o mod_dimensions.o mod_mod_states.o mod_obs_states.o \
	mod_manage_obs.o sublocan.o
ENSOBJ2 = analysis.o m_randrot.o m_multa.o mod_anafunc.o m_mean_preserving_rotation.o
#ENSOBJ2 = analysis6c.o m_randrot.o m_multa.o 
#ENSOBJ2 = analysis2.o m_randrot.o m_multa.o 
ENSOBJ3 = m_sample2D.o m_pseudo2D.o mod_fftw3.o m_newton2D.o \
	m_newtonfunc2D.o m_fixsample2D.o m_random.o

#--------------------
# Libraries to link
LIBS = -llapack -lblas -lfftw3 -lpthread

#--------------------
enkf_analysis: enkf_analysis.o $(ENSOBJ1) $(ENSOBJ2) $(ENSOBJ3)
	$(LINKER) $(LFLAGS) $(LIBS) $(HTOBJS) $@.o \
        $(ENSOBJ1) $(ENSOBJ2) $(ENSOBJ3) $(LIBDIR)/libfem.a \
	-o $@

make_ens_meteo: make_ens_meteo.o $(ENSOBJ3) m_randrot.o
	$(LINKER) $(LFLAGS) $@.o $(ENSOBJ3) m_randrot.o $(LIBS) $(LIBDIR)/libfem.a -o $@

make_eof_ens_meteo: make_eof_ens_meteo.o m_random.o
	$(LINKER) $(LFLAGS) $@.o m_random.o $(LIBS) $(LIBDIR)/libfem.a -llapack -lblas -o $@

enKF2enKS: enKF2enKS.o 
	$(LINKER) $(LFLAGS) $(LIBS) $(HTOBJS) $@.o \
	$(LIBDIR)/libfem.a -o $@

#--------------------
#--------------------
EXES = enkf_analysis make_ens_meteo make_eof_ens_meteo enKF2enKS

all: $(EXES)

default: all

clean:
	rm -f $(EXES) *.o *.mod tags

cleanall:
	rm -f $(EXES) *.o *.mod tags *.dat *.uf fort.* *.bak


#---------------------
# Compilation rules

.SUFFIXES: .f .F90

.f.o:
	$(F77) -c $(FFLAGS) $<

.F90.o:
	$(F77) -c $(FFLAGS) $< 
