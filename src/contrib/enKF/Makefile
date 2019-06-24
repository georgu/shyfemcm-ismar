#
# Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
#

FEMDIR = ..
FEMBIN = $(FEMDIR)/fembin
FEM3D = $(FEMDIR)/fem3d

include ../Rules.make

# Switch to real64
ifeq ($(FORTRAN_COMPILER),GNU_GFORTRAN)
  R8FL = -fdefault-real-8 -fdefault-double-8
endif

ifeq ($(FORTRAN_COMPILER),INTEL)
  R8FL = -real-size 64 -double-size 64
endif
override FFLAGS += $(R8FL)
override LFLAGS += $(R8FL)

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#--------------------
# ht objects
HTOBJS = $(FEM3D)/subrst.o $(FEM3D)/subbas.o $(FEM3D)/subdts.o $(FEM3D)/mod_bnd.o \
	$(FEM3D)/subreg.o $(FEM3D)/subfind.o $(FEM3D)/subpar3.o $(FEM3D)/suboutput.o \
	$(FEM3D)/subfil.o $(FEM3D)/mod_hydro.o $(FEM3D)/debug.o $(FEM3D)/subsss.o \
	$(FEM3D)/subbit10.o $(FEM3D)/mod_geom_dynamic.o $(FEM3D)/mod_ts.o \
	$(FEM3D)/newlevels.o $(FEM3D)/mod_conz.o $(FEM3D)/mod_hydro_vel.o \
	$(FEM3D)/ecological_dummy.o $(FEM3D)/subdef.o $(FEM3D)/subnev.o \
	$(FEM3D)/mod_geom.o $(FEM3D)/subscn.o $(FEM3D)/suboutputd.o \
	$(FEM3D)/subiso8601.o

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
main: main.o $(ENSOBJ1) $(ENSOBJ2) $(ENSOBJ3)
	$(LINKER) $(LFLAGS) $(LIBS) $(HTOBJS) $@.o \
        $(ENSOBJ1) $(ENSOBJ2) $(ENSOBJ3) $(DIRLIB)/libfem.a \
	-o $@

make_ens_meteo: make_ens_meteo.o $(ENSOBJ3) m_randrot.o
	$(LINKER) $(LFLAGS) $@.o $(ENSOBJ3) m_randrot.o $(LIBS) $(DIRLIB)/libfem.a -o $@

make_eof_ens_meteo: make_eof_ens_meteo.o m_random.o
	$(LINKER) $(LFLAGS) $@.o m_random.o $(LIBS) $(DIRLIB)/libfem.a -llapack -lblas -o $@

enKF2enKS: enKF2enKS.o 
	$(LINKER) $(LFLAGS) $(LIBS) $(HTOBJS) $@.o \
	$(DIRLIB)/libfem.a -o $@

depend:
	@$(MAKEDEPEND) -moddir=$(MODDIR) *.F90


#--------------------
#--------------------
EXES = main make_ens_meteo make_eof_ens_meteo enKF2enKS

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



# DO NOT DELETE THIS LINE -- make depend depends on it.

../femlib/mod/m_fixsample2d.mod: m_fixsample2D.o
../femlib/mod/m_mean_preserving_rotation.mod: \
		m_mean_preserving_rotation.o 
../femlib/mod/m_multa.mod: m_multa.o
../femlib/mod/m_newton2d.mod: m_newton2D.o
../femlib/mod/m_newtonfunc2d.mod: m_newtonfunc2D.o
../femlib/mod/m_pseudo2d.mod: m_pseudo2D.o
../femlib/mod/m_random.mod: m_random.o
../femlib/mod/m_randrot.mod: m_randrot.o
../femlib/mod/m_sample2d.mod: m_sample2D.o
../femlib/mod/mod_anafunc.mod: mod_anafunc.o
../femlib/mod/mod_dimensions.mod: mod_dimensions.o
../femlib/mod/mod_enkf.mod: mod_enkf.o
../femlib/mod/mod_ens_state.mod: mod_ens_state.o
../femlib/mod/mod_fftw3.mod: mod_fftw3.o
../femlib/mod/mod_init_enkf.mod: mod_init_enkf.o
../femlib/mod/mod_manage_obs.mod: mod_manage_obs.o
../femlib/mod/mod_mod_err.mod: mod_mod_err.o
../femlib/mod/mod_mod_states.mod: mod_mod_states.o
../femlib/mod/mod_obs_states.mod: mod_obs_states.o
../femlib/mod/mod_para.mod: mod_para.o
analysis.o: ../femlib/mod/m_multa.mod \
		../femlib/mod/mod_anafunc.mod 
analysis2.o: ../femlib/mod/m_multa.mod
analysis6c.o: ../femlib/mod/m_multa.mod \
		../femlib/mod/m_randrot.mod 
enKF2enKS.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_conz.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_restart.mod \
		../femlib/mod/mod_ts.mod 
m_mean_preserving_rotation.o: ../femlib/mod/m_randrot.mod
m_newton2D.o: ../femlib/mod/m_newtonfunc2d.mod
m_pseudo2D.o: ../femlib/mod/m_newton2d.mod \
		../femlib/mod/mod_fftw3.mod 
m_sample2D.o: ../femlib/mod/m_fixsample2d.mod \
		../femlib/mod/m_pseudo2d.mod \
		../femlib/mod/m_randrot.mod 
main.o: ../femlib/mod/mod_enkf.mod ../femlib/mod/mod_ens_state.mod \
		../femlib/mod/mod_mod_err.mod \
		../femlib/mod/mod_mod_states.mod \
		../femlib/mod/mod_para.mod 
make_ens_meteo.o: ../femlib/mod/m_random.mod \
		../femlib/mod/m_sample2d.mod 
make_eof_ens_meteo.o: ../femlib/mod/m_random.mod
mod_anafunc.o: ../femlib/mod/m_mean_preserving_rotation.mod \
		../femlib/mod/m_randrot.mod 
mod_enkf.o: ../femlib/mod/mod_ens_state.mod \
		../femlib/mod/mod_manage_obs.mod \
		../femlib/mod/mod_para.mod 
mod_ens_state.o: ../femlib/mod/basin.mod \
		../femlib/mod/mod_conz.mod \
		../femlib/mod/mod_dimensions.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_init_enkf.mod \
		../femlib/mod/mod_mod_states.mod \
		../femlib/mod/mod_para.mod \
		../femlib/mod/mod_ts.mod 
mod_init_enkf.o: ../femlib/mod/basin.mod ../femlib/mod/iso8601.mod \
		../femlib/mod/mod_dimensions.mod 
mod_manage_obs.o: ../femlib/mod/iso8601.mod \
		../femlib/mod/levels.mod \
		../femlib/mod/mod_ens_state.mod \
		../femlib/mod/mod_init_enkf.mod \
		../femlib/mod/mod_obs_states.mod \
		../femlib/mod/mod_para.mod 
mod_mod_err.o: ../femlib/mod/mod_ens_state.mod \
		../femlib/mod/mod_init_enkf.mod \
		../femlib/mod/mod_mod_states.mod 
mod_mod_states.o: ../femlib/mod/mod_dimensions.mod
subenkf.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/m_random.mod \
		../femlib/mod/m_sample2d.mod \
		../femlib/mod/mod_conz.mod \
		../femlib/mod/mod_dimensions.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_restart.mod \
		../femlib/mod/mod_ts.mod ../femlib/mod/shyfile.mod 
sublocan.o: ../femlib/mod/basin.mod ../femlib/mod/mod_enkf.mod \
		../femlib/mod/mod_ens_state.mod \
		../femlib/mod/mod_manage_obs.mod \
		../femlib/mod/mod_mod_states.mod \
		../femlib/mod/mod_para.mod 

