#
# Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
#

#-------------------- default values from shyfem
FEMDIR = ..
FEMBIN = $(FEMDIR)/fembin
FEM3D = $(FEMDIR)/fem3d
include ../Rules.make

#-------------------- switch from real/double-precision to real64 (real*8)
ifeq ($(FORTRAN_COMPILER),GNU_GFORTRAN)
  R8FL = -fdefault-real-8 -fdefault-double-8
endif
ifeq ($(FORTRAN_COMPILER),INTEL)
  R8FL = -real-size 64 -double-size 64
endif
override FFLAGS += $(R8FL)
override LFLAGS += $(R8FL)

#-------------------- run to make depends
MAKEDEPEND = $(FEMBIN)/mkdp.pl

#-------------------- shyfem objects
SHYOBJ = $(FEM3D)/subrst.o $(FEM3D)/subbas.o $(FEM3D)/subdts.o $(FEM3D)/mod_bnd.o \
	$(FEM3D)/subreg.o $(FEM3D)/subfind.o $(FEM3D)/subpar3.o $(FEM3D)/suboutput.o \
	$(FEM3D)/subfil.o $(FEM3D)/mod_hydro.o $(FEM3D)/debug.o $(FEM3D)/subsss.o \
	$(FEM3D)/subbit10.o $(FEM3D)/mod_geom_dynamic.o $(FEM3D)/mod_ts.o \
	$(FEM3D)/newlevels.o $(FEM3D)/mod_conz.o $(FEM3D)/mod_hydro_vel.o \
	$(FEM3D)/ecological_dummy.o $(FEM3D)/subdef.o $(FEM3D)/subnev.o \
	$(FEM3D)/mod_geom.o $(FEM3D)/subscn.o $(FEM3D)/suboutputd.o \
	$(FEM3D)/subiso8601.o $(FEM3D)/mercury_dummy.o $(FEM3D)/mod_gotm_aux.o

#-------------------- ens objects
# Ensemble analysis objects
#
ENKOBJ = analysis.o m_ensmean.o m_ensvar.o m_mean_preserving_rotation.o \
	 m_multa.o m_random.o m_randrot.o mod_anafunc.o

SAMPOBJ = m_pseudo2D.o m_fixsample2D.o m_newton2D.o m_newtonfunc2D.o \
	  m_sample2D.o mod_fftw3.o
          
ENSOBJ = $(ENKOBJ) $(SAMPOBJ) subenkf.o mod_mod_err.o mod_enkf.o \
	 mod_ens_state.o mod_init_enkf.o mod_manage_obs.o mod_mod_states.o \
	 mod_obs_states.o mod_para.o sublocan.o mod_dimensions.o

#-------------------- external libraries
LIBS = -llapack -lblas -lfftw3 -lpthread

#-------------------- ens modules (just to clean)
ENSMOD = $(MODDIR)/m_random.mod $(MODDIR)/mod_mod_err.mod $(MODDIR)/mod_enkf.mod \
	 $(MODDIR)/mod_manage_obs.mod $(MODDIR)/mod_ens_state.mod $(MODDIR)/mod_para.mod \
	 $(MODDIR)/mod_init_enkf.mod $(MODDIR)/mod_anafunc.mod \
	 $(MODDIR)/m_mean_preserving_rotation.mod $(MODDIR)/m_multa.mod $(MODDIR)/m_sample2d.mod \
	 $(MODDIR)/m_randrot.mod $(MODDIR)/m_pseudo2d.mod $(MODDIR)/mod_fftw3.mod \
	 $(MODDIR)/m_newton2d.mod $(MODDIR)/m_newtonfunc2d.mod $(MODDIR)/m_fixsample2d.mod \
	 $(MODDIR)/mod_obs_states.mod $(MODDIR)/mod_mod_states.mod $(MODDIR)/mod_dimensions.mod

#-------------------- custom files (just to clean)
CUSTF = mod_dimensions.F90

#-------------------- all the exec files
EXES = main make_eof_ens_meteo enKF2enKS perturbeBC make_antime_list


#------------------- targets
default: $(EXES)

all: $(EXES)

clean:
	rm -f $(EXES) *.o $(ENSMOD)

cleanall:
	rm -f $(EXES) *.o $(ENSMOD) $(CUSTF) tags *.dat *.uf fort.* *.bak

main: main.o $(ENSOBJ) 
	$(LINKER) $(SHYOBJ) $@.o \
        $(ENSOBJ) $(DIRLIB)/libfem.a $(LFLAGS) $(LIBS) -o $@

perturbeBC: perturbeBC.o m_sample2D.o m_pseudo2D.o m_newton2D.o \
	m_randrot.o m_fixsample2D.o m_random.o m_newtonfunc2D.o
	$(LINKER) $@.o m_sample2D.o m_pseudo2D.o m_newton2D.o \
	m_randrot.o m_fixsample2D.o m_random.o m_newtonfunc2D.o \
	$(DIRLIB)/libfem.a $(LFLAGS) $(LIBS) -o $@

make_antime_list: make_antime_list.o
	$(LINKER) $(LFLAGS) $@.o $(DIRLIB)/libfem.a -o $@

make_eof_ens_meteo: make_eof_ens_meteo.o m_random.o
	$(LINKER) $@.o m_random.o \
	$(DIRLIB)/libfem.a $(LFLAGS) $(LIBS) -o $@

enKF2enKS: enKF2enKS.o 
	$(LINKER) $(SHYOBJ) $@.o \
	$(DIRLIB)/libfem.a $(LFLAGS) $(LIBS) -o $@

depend:
	@$(MAKEDEPEND) -moddir=$(MODDIR) *.F90


#------------------- compilation rules
.SUFFIXES: .f .F90
.f.o:
	$(F77) -c $(FFLAGS) $<
.F90.o:
	$(F77) -c $(FFLAGS) $< 


# DO NOT DELETE THIS LINE -- make depend depends on it.

../femlib/mod/m_ensmean.mod: m_ensmean.o
../femlib/mod/m_ensvar.mod: m_ensvar.o
../femlib/mod/m_fixsample1d.mod: m_fixsample1D.o
../femlib/mod/m_fixsample2d.mod: m_fixsample2D.o
../femlib/mod/m_mean_preserving_rotation.mod: \
		m_mean_preserving_rotation.o 
../femlib/mod/m_multa.mod: m_multa.o
../femlib/mod/m_newton1d.mod: m_newton1D.o
../femlib/mod/m_newton2d.mod: m_newton2D.o
../femlib/mod/m_newtonfunc1d.mod: m_newtonfunc1D.o
../femlib/mod/m_newtonfunc2d.mod: m_newtonfunc2D.o
../femlib/mod/m_pseudo1d.mod: m_pseudo1D.o
../femlib/mod/m_pseudo2d.mod: m_pseudo2D.o
../femlib/mod/m_random.mod: m_random.o
../femlib/mod/m_randrot.mod: m_randrot.o
../femlib/mod/m_sample1d.mod: m_sample1D.o
../femlib/mod/m_sample2d.mod: m_sample2D.o
../femlib/mod/m_set_random_seed2.mod: m_set_random_seed2.o
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
analysis.o: ../femlib/mod/m_ensmean.mod ../femlib/mod/m_ensvar.mod \
		../femlib/mod/m_multa.mod \
		../femlib/mod/mod_anafunc.mod 
analysis2_EnOI.o: ../femlib/mod/m_multa.mod
enKF2enKS.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_conz.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_restart.mod \
		../femlib/mod/mod_ts.mod 
m_mean_preserving_rotation.o: ../femlib/mod/m_randrot.mod
m_newton1D.o: ../femlib/mod/m_newtonfunc1d.mod
m_newton2D.o: ../femlib/mod/m_newtonfunc2d.mod
m_pseudo1D.o: ../femlib/mod/m_newton1d.mod \
		../femlib/mod/m_random.mod \
		../femlib/mod/mod_fftw3.mod 
m_pseudo2D.o: ../femlib/mod/m_newton2d.mod \
		../femlib/mod/mod_fftw3.mod 
m_sample1D.o: ../femlib/mod/m_fixsample1d.mod \
		../femlib/mod/m_pseudo1d.mod 
m_sample2D.o: ../femlib/mod/m_fixsample2d.mod \
		../femlib/mod/m_pseudo2d.mod \
		../femlib/mod/m_randrot.mod 
main.o: ../femlib/mod/mod_enkf.mod ../femlib/mod/mod_ens_state.mod \
		../femlib/mod/mod_mod_err.mod \
		../femlib/mod/mod_mod_states.mod \
		../femlib/mod/mod_para.mod \
		../femlib/mod/mod_restart.mod 
make_antime_list.o: ../femlib/mod/iso8601.mod
make_eof_ens_meteo.o: ../femlib/mod/m_random.mod
mod_anafunc.o: ../femlib/mod/m_ensmean.mod \
		../femlib/mod/m_mean_preserving_rotation.mod \
		../femlib/mod/m_multa.mod \
		../femlib/mod/m_random.mod \
		../femlib/mod/m_randrot.mod 
mod_enkf.o: ../femlib/mod/levels.mod \
		../femlib/mod/mod_ens_state.mod \
		../femlib/mod/mod_manage_obs.mod \
		../femlib/mod/mod_para.mod 
mod_ens_state.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/mod_conz.mod \
		../femlib/mod/mod_dimensions.mod \
		../femlib/mod/mod_geom_dynamic.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_hydro_vel.mod \
		../femlib/mod/mod_init_enkf.mod \
		../femlib/mod/mod_mod_states.mod \
		../femlib/mod/mod_para.mod \
		../femlib/mod/mod_restart.mod \
		../femlib/mod/mod_ts.mod \
		../femlib/mod/mod_gotm_aux.mod
mod_init_enkf.o: ../femlib/mod/iso8601.mod
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
perturbeBC.o: ../femlib/mod/iso8601.mod ../femlib/mod/m_random.mod \
		../femlib/mod/m_sample2d.mod 
subenkf.o: ../femlib/mod/basin.mod ../femlib/mod/levels.mod \
		../femlib/mod/m_random.mod \
		../femlib/mod/m_sample2d.mod \
		../femlib/mod/mod_dimensions.mod \
		../femlib/mod/mod_hydro.mod \
		../femlib/mod/mod_restart.mod \
		../femlib/mod/shyfile.mod 
sublocan.o: ../femlib/mod/basin.mod ../femlib/mod/mod_enkf.mod \
		../femlib/mod/mod_ens_state.mod \
		../femlib/mod/mod_manage_obs.mod \
		../femlib/mod/mod_para.mod \
		../femlib/mod/mod_restart.mod 

