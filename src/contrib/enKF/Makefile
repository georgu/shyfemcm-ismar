#
# Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
#

#-------------------- default values from shyfem

FEMDIR = ../../..
FEMBIN = $(FEMDIR)/bin
FEM3D = $(FEMDIR)/fem3d
FEMSRC = $(FEMDIR)/src
SHYGENERIC = $(FEMSRC)/utils/generic
SHYUTIL = $(FEMSRC)/utils/shyutil
SHYMPI = $(FEMSRC)/utils/shympi
SHYUTILMPI = $(FEMSRC)/utils/shyutilmpi
SHYMOD = $(FEMSRC)/utils/modules
SHYSHYFEM = $(FEMSRC)/shyfem

include $(FEMDIR)/Rules.make

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

SHYOBJ = \
	$(SHYGENERIC)/dts.o \
	$(SHYGENERIC)/file.o \
	$(SHYGENERIC)/femfile.o \
	$(SHYGENERIC)/iso8601.o \
	$(SHYGENERIC)/para.o \
	$(SHYGENERIC)/bit10.o \
	$(SHYGENERIC)/convert.o \
	$(SHYGENERIC)/time_util.o \
	$(SHYGENERIC)/round.o \
	$(SHYGENERIC)/sort_util.o \
	$(SHYUTIL)/basin.o \
	$(SHYUTIL)/levels.o \
	$(SHYUTIL)/find.o \
	$(SHYUTIL)/regular.o \
	$(SHYUTIL)/fem_util.o \
	$(SHYUTIL)/geom.o \
	$(SHYUTIL)/evgeom.o \
	$(SHYUTIL)/geom_dynamic.o \
	$(SHYUTIL)/defnames.o \
	$(SHYMOD)/mod_hydro.o \
	$(SHYMOD)/mod_hydro_vel.o \
	$(SHYMOD)/mod_ts.o \
	$(SHYSHYFEM)/outputd.o \
	$(SHYSHYFEM)/restart.o \
	$(SHYUTILMPI)/mod_restart.o \
	$(SHYMPI)/mpi_dummy.o \
	$(SHYMPI)/mpi_general.o \
	$(SHYSHYFEM)/conz_common.o \
	$(SHYSHYFEM)/ecological_dummy.o \
	$(SHYSHYFEM)/gotm_aux.o \
	$(FEMSRC)/contrib/ecological/mercury/mercury_dummy.o \
	$(FEMSRC)/contrib/ecological/bfm/bfm_dummy.o \

SHYOBJENS = \
	$(SHYGENERIC)/dts.o \
	$(SHYGENERIC)/file.o \
	$(SHYGENERIC)/femfile.o \
	$(SHYGENERIC)/iso8601.o \
	$(SHYGENERIC)/convert.o \

SHYOBJ_OLD = \
	$(FEM3D)/subrst.o $(FEM3D)/subbas.o \
	$(FEM3D)/subdts.o $(FEM3D)/mod_bnd.o \
	$(FEM3D)/subreg.o $(FEM3D)/subfind.o $(FEM3D)/subpar3.o \
	$(FEM3D)/subfil.o $(FEM3D)/mod_hydro.o \
	$(FEM3D)/debug.o $(FEM3D)/subsss.o \
	$(FEM3D)/subbit10.o $(FEM3D)/mod_geom_dynamic.o $(FEM3D)/mod_ts.o \
	$(FEM3D)/newlevels.o $(FEM3D)/mod_conz.o $(FEM3D)/mod_hydro_vel.o \
	$(FEM3D)/ecological_dummy.o $(FEM3D)/subdef.o $(FEM3D)/subnev.o \
	$(FEM3D)/mod_geom.o $(FEM3D)/subscn.o $(FEM3D)/suboutputd.o \
	$(FEM3D)/subiso8601.o $(FEM3D)/mercury_dummy.o $(FEM3D)/mod_gotm_aux.o

#-------------------- ens objects
# Ensemble analysis objects

ENKOBJ = m_analysis.o m_ensmean.o m_ensvar.o \
	 m_mean_preserving_rotation.o \
	 m_multa.o mod_anafunc.o

SAMPOBJ = m_pseudo2D.o m_fixsample2D.o m_newton2D.o m_newtonfunc2D.o \
	  m_sample2D.o mod_fftw3.o m_set_random_seed2.o m_randrot.o \
	  m_random.o
          
ENSOBJ = $(ENKOBJ) $(SAMPOBJ) subenkf.o mod_mod_err.o mod_enkf.o \
	 mod_ens_state.o mod_init_enkf.o mod_manage_obs.o mod_mod_states.o \
	 mod_obs_states.o mod_para.o sublocan.o

#-------------------- external libraries

LIBS = -llapack -lblas -lfftw3 -lpthread

#-------------------- ens modules (just to clean)

ENSMOD = \
	 $(MODDIR)/m_random.mod \
	 $(MODDIR)/mod_mod_err.mod \
	 $(MODDIR)/mod_enkf.mod \
	 $(MODDIR)/mod_manage_obs.mod \
	 $(MODDIR)/mod_ens_state.mod \
	 $(MODDIR)/mod_para.mod \
	 $(MODDIR)/mod_init_enkf.mod \
	 $(MODDIR)/mod_anafunc.mod \
	 $(MODDIR)/m_randrot.mod \
	 $(MODDIR)/m_mean_preserving_rotation.mod \
	 $(MODDIR)/m_multa.mod \
	 $(MODDIR)/m_sample2d.mod \
	 $(MODDIR)/m_pseudo2d.mod \
	 $(MODDIR)/mod_fftw3.mod \
	 $(MODDIR)/m_newton2d.mod \
	 $(MODDIR)/m_newtonfunc2d.mod \
	 $(MODDIR)/m_fixsample2d.mod \
	 $(MODDIR)/mod_obs_states.mod \
	 $(MODDIR)/mod_mod_states.mod \
	 $(MODDIR)/m_set_random_seed2.mod \
	 $(MODDIR)/m_analysis.mod \

#-------------------- all the exec files

EXES = main make_eof_ens_meteo enKF2enKS perturbeBC make_antime_list perturbe_state

#------------------- targets

default: $(EXES)

all: $(EXES)

main: $(ENSOBJ) main.o
	$(LINKER) $(SHYOBJ) $@.o \
	$(ENSOBJ) $(LFLAGS) $(LIBS) -o $@

perturbeBC: perturbeBC.o $(SAMPOBJ)
	$(LINKER) $@.o $(SAMPOBJ) \
	$(SHYOBJENS) \
	$(LFLAGS) $(LIBS) -o $@

make_antime_list: make_antime_list.o 
	$(LINKER) $(LFLAGS) $(SHYOBJENS) $@.o -o $@

make_eof_ens_meteo: make_eof_ens_meteo.o m_random.o
	$(LINKER) $@.o m_random.o $(SHYOBJENS) \
	$(LFLAGS) $(LIBS) -o $@

enKF2enKS: enKF2enKS.o 
	$(LINKER) $(SHYOBJ) $@.o \
	$(LFLAGS) $(LIBS) -o $@

perturbe_state: perturbe_state.o subenkf.o $(SAMPOBJ)
	$(LINKER) $(SHYOBJ) subenkf.o $(SAMPOBJ) $@.o \
	$(LFLAGS) $(LIBS) -o $@

#------------------- var

depend:
	@$(MAKEDEPEND) -moddir=$(MODDIR) *.f90

clean:
	rm -f $(EXES) *.o $(ENSMOD)

cleanall:
	rm -f $(EXES) *.o $(ENSMOD) tags *.dat *.uf fort.* *.bak

#------------------- compilation rules

.SUFFIXES: .f .F90 .f90

.f.o:
	$(F77) -c $(FFLAGS) $<
.F90.o:
	$(F77) -c $(FFLAGS) $< 
.f90.o:
	$(F77) -c $(FFLAGS) $< 


# DO NOT DELETE THIS LINE -- make depend depends on it.

../../../lib/mod/m_ensmean.mod: m_ensmean.o
../../../lib/mod/m_ensvar.mod: m_ensvar.o
../../../lib/mod/m_fixsample1d.mod: m_fixsample1D.o
../../../lib/mod/m_fixsample2d.mod: m_fixsample2D.o
../../../lib/mod/m_mean_preserving_rotation.mod: \
		m_mean_preserving_rotation.o 
../../../lib/mod/m_multa.mod: m_multa.o
../../../lib/mod/m_newton1d.mod: m_newton1D.o
../../../lib/mod/m_newton2d.mod: m_newton2D.o
../../../lib/mod/m_newtonfunc1d.mod: m_newtonfunc1D.o
../../../lib/mod/m_newtonfunc2d.mod: m_newtonfunc2D.o
../../../lib/mod/m_pseudo1d.mod: m_pseudo1D.o
../../../lib/mod/m_pseudo2d.mod: m_pseudo2D.o
../../../lib/mod/m_random.mod: m_random.o
../../../lib/mod/m_randrot.mod: m_randrot.o
../../../lib/mod/m_sample1d.mod: m_sample1D.o
../../../lib/mod/m_sample2d.mod: m_sample2D.o
../../../lib/mod/m_set_random_seed2.mod: m_set_random_seed2.o
../../../lib/mod/mod_anafunc.mod: mod_anafunc.o
../../../lib/mod/mod_enkf.mod: mod_enkf.o
../../../lib/mod/mod_ens_state.mod: mod_ens_state.o
../../../lib/mod/mod_fftw3.mod: mod_fftw3.o
../../../lib/mod/mod_init_enkf.mod: mod_init_enkf.o
../../../lib/mod/mod_manage_obs.mod: mod_manage_obs.o
../../../lib/mod/mod_mod_err.mod: mod_mod_err.o
../../../lib/mod/mod_mod_states.mod: mod_mod_states.o
../../../lib/mod/mod_obs_states.mod: mod_obs_states.o
../../../lib/mod/mod_para.mod: mod_para.o
../../../lib/mod/m_analysis.mod: m_analysis.o
m_analysis.o: ../../../lib/mod/m_ensmean.mod \
		../../../lib/mod/m_ensvar.mod \
		../../../lib/mod/m_multa.mod \
		../../../lib/mod/mod_anafunc.mod 
analysis2_EnOI.o: ../../../lib/mod/m_multa.mod
m_mean_preserving_rotation.o: ../../../lib/mod/m_randrot.mod
m_newton1D.o: ../../../lib/mod/m_newtonfunc1d.mod
m_newton2D.o: ../../../lib/mod/m_newtonfunc2d.mod
m_pseudo1D.o: ../../../lib/mod/m_newton1d.mod \
		../../../lib/mod/m_random.mod \
		../../../lib/mod/mod_fftw3.mod 
m_pseudo2D.o: ../../../lib/mod/m_newton2d.mod \
		../../../lib/mod/mod_fftw3.mod 
m_sample1D.o: ../../../lib/mod/m_fixsample1d.mod \
		../../../lib/mod/m_pseudo1d.mod 
m_sample2D.o: ../../../lib/mod/m_fixsample2d.mod \
		../../../lib/mod/m_pseudo2d.mod \
		../../../lib/mod/m_randrot.mod 
main.o: ../../../lib/mod/mod_enkf.mod \
		../../../lib/mod/mod_mod_err.mod \
		../../../lib/mod/mod_restart.mod \
		../../../lib/mod/m_set_random_seed2.mod \
		../../../lib/mod/m_analysis.mod
make_antime_list.o: ../../../lib/mod/iso8601.mod
make_eof_ens_meteo.o: ../../../lib/mod/m_random.mod
mod_anafunc.o: ../../../lib/mod/m_ensmean.mod \
		../../../lib/mod/m_mean_preserving_rotation.mod \
		../../../lib/mod/m_multa.mod \
		../../../lib/mod/m_random.mod \
		../../../lib/mod/m_randrot.mod 
mod_enkf.o: ../../../lib/mod/levels.mod \
		../../../lib/mod/mod_ens_state.mod \
		../../../lib/mod/mod_manage_obs.mod \
		../../../lib/mod/mod_para.mod 
mod_ens_state.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/levels.mod \
		../../../lib/mod/mod_conz.mod \
		../../../lib/mod/mod_geom_dynamic.mod \
		../../../lib/mod/mod_gotm_aux.mod \
		../../../lib/mod/mod_hydro.mod \
		../../../lib/mod/mod_hydro_vel.mod \
		../../../lib/mod/mod_init_enkf.mod \
		../../../lib/mod/mod_mod_states.mod \
		../../../lib/mod/mod_para.mod \
		../../../lib/mod/mod_restart.mod \
		../../../lib/mod/mod_ts.mod \
		../../../lib/mod/shympi.mod 
mod_init_enkf.o: ../../../lib/mod/iso8601.mod \
		../../../lib/mod/mod_para.mod 
mod_manage_obs.o: ../../../lib/mod/iso8601.mod \
		../../../lib/mod/levels.mod \
		../../../lib/mod/mod_ens_state.mod \
		../../../lib/mod/mod_init_enkf.mod \
		../../../lib/mod/mod_obs_states.mod \
		../../../lib/mod/mod_para.mod 
mod_mod_err.o: ../../../lib/mod/mod_ens_state.mod \
		../../../lib/mod/mod_init_enkf.mod \
		../../../lib/mod/mod_mod_states.mod 
perturbeBC.o: ../../../lib/mod/iso8601.mod \
		../../../lib/mod/m_random.mod \
		../../../lib/mod/m_sample2d.mod 
subenkf.o: ../../../lib/mod/basin.mod ../../../lib/mod/levels.mod \
		../../../lib/mod/m_random.mod \
		../../../lib/mod/m_sample2d.mod \
		../../../lib/mod/mod_hydro.mod \
		../../../lib/mod/mod_para.mod \
		../../../lib/mod/mod_restart.mod \
		../../../lib/mod/shyfile.mod \
		../../../lib/mod/shympi.mod 
sublocan.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/mod_enkf.mod \
		../../../lib/mod/mod_ens_state.mod \
		../../../lib/mod/mod_manage_obs.mod \
		../../../lib/mod/mod_para.mod \
		../../../lib/mod/mod_restart.mod \
		../../../lib/mod/m_analysis.mod
enKF2enKS.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/levels.mod \
		../../../lib/mod/mod_conz.mod \
		../../../lib/mod/mod_geom_dynamic.mod \
		../../../lib/mod/mod_hydro.mod \
		../../../lib/mod/mod_hydro_vel.mod \
		../../../lib/mod/mod_gotm_aux.mod \
		../../../lib/mod/mod_restart.mod \
		../../../lib/mod/mod_ts.mod \
		../../../lib/mod/shympi.mod 
perturbe_state.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/levels.mod \
		../../../lib/mod/mod_conz.mod \
		../../../lib/mod/mod_geom_dynamic.mod \
		../../../lib/mod/mod_hydro.mod \
		../../../lib/mod/mod_hydro_vel.mod \
		../../../lib/mod/mod_gotm_aux.mod \
		../../../lib/mod/mod_restart.mod \
		../../../lib/mod/mod_ts.mod \
		../../../lib/mod/shympi.mod 
