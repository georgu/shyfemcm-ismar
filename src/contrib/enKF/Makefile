#
# Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
#
################################################################################
# EnKF Data Assimilation Module - Makefile
################################################################################

# --- Directory Paths ---
FEMDIR       = ../../..
FEMBIN       = $(FEMDIR)/bin
FEM3D        = $(FEMDIR)/fem3d
FEMSRC       = $(FEMDIR)/src
SHYGENERIC   = $(FEMSRC)/utils/generic
SHYUTIL      = $(FEMSRC)/utils/shyutil
SHYMPI       = $(FEMSRC)/utils/shympi
SHYUTILMPI   = $(FEMSRC)/utils/shyutilmpi
SHYMOD       = $(FEMSRC)/utils/modules
SHYSHYFEM    = $(FEMSRC)/shyfem

# Include global configuration rules from the main SHYFEM framework
include $(FEMDIR)/Rules.make

# --- Hybrid Parallelization Logic ---

# 1. Determine the correct OpenMP flag based on the active compiler defined in Rules.make
# This ensures OpenMP is ALWAYS enabled for EnKF, even if PARALLEL_OMP = false globally.
ENKF_OMP_FLAG =
ifeq ($(FORTRAN_COMPILER),GNU_GFORTRAN)
  ENKF_OMP_FLAG = -fopenmp
endif
ifeq ($(FORTRAN_COMPILER),INTEL)
  ENKF_OMP_FLAG = -qopenmp
  ifeq ($(IMV_LE_14),true)
    ENKF_OMP_FLAG = -openmp
  endif
endif
ifeq ($(FORTRAN_COMPILER),$(filter $(FORTRAN_COMPILER),PGI PORTLAND))
  ENKF_OMP_FLAG = -mp
endif
ifeq ($(FORTRAN_COMPILER),IBM)
  ENKF_OMP_FLAG = -qsmp=omp -qnosave
endif

# 2. Setup the Linker wrapper and MPI objects if multi-node execution is active
ifneq ($(PARALLEL_MPI), NONE)
  ifdef MPIFC
    LINKER      = $(MPIFC)
  else
    LINKER      = mpif90
  endif
  # Inject real MPI buffers from SHYFEM upstream to solve link interdependencies
  SHY_MPI_OBJS  = $(SHYMPI)/mpi_node.o $(SHYMPI)/mpi_node_internal.o $(SHYMPI)/mpi_buffer.o $(SHYMPI)/mpi_time.o
else
  LINKER        = $(FC)
  SHY_MPI_OBJS  = $(SHYMPI)/mpi_dummy.o
endif

# 3. Safe fallback in case global FFLAGS lacks the OpenMP flag
# We append our enforced OpenMP flag to ensure the compiler processes parallel directives
override FFLAGS += $(ENKF_OMP_FLAG)

# --- External Library Configurations ---
NC_FLAGS := $(shell nf-config --fflags)
NC_LIBS  := $(shell nf-config --flibs)

MAKEDEPEND = $(FEMBIN)/mkdp.pl

# --- SHYFEM Base Objects ---
SHYOBJ = \
	$(SHYGENERIC)/dts.o \
	$(SHYGENERIC)/file.o \
	$(SHYGENERIC)/compatibility.o \
	$(SHYGENERIC)/femfile.o \
	$(SHYGENERIC)/iso8601.o \
	$(SHYGENERIC)/para.o \
	$(SHYGENERIC)/bit10.o \
	$(SHYGENERIC)/convert.o \
	$(SHYGENERIC)/time_util.o \
	$(SHYGENERIC)/round.o \
	$(SHYGENERIC)/sort_util.o \
	$(SHYGENERIC)/subsort.o \
	$(SHYUTIL)/basin.o \
	$(SHYUTIL)/levels.o \
	$(SHYUTIL)/find.o \
	$(SHYUTIL)/regular.o \
	$(SHYUTIL)/fem_util.o \
	$(SHYUTIL)/geom.o \
	$(SHYUTIL)/evgeom.o \
	$(SHYUTIL)/geom_dynamic.o \
	$(SHYUTIL)/defnames.o \
	$(SHYUTIL)/depth.o \
	$(SHYUTIL)/sigma.o \
	$(SHYUTIL)/zadaptutil.o \
	$(SHYMOD)/mod_hydro.o \
	$(SHYMOD)/mod_area.o \
	$(SHYMOD)/mod_hydro_vel.o \
	$(SHYMOD)/mod_ts.o \
	$(SHYMOD)/mod_shyfem.o \
	$(SHYMOD)/mod_layer_thickness.o \
	$(SHYSHYFEM)/outputd.o \
	$(SHYSHYFEM)/restart.o \
	$(SHYUTILMPI)/mod_restart.o \
	$(SHYUTILMPI)/info_output.o \
	$(SHYUTILMPI)/depth_util.o \
	$(SHY_MPI_OBJS) \
	$(SHYMPI)/mpi_general.o \
	$(SHYSHYFEM)/elems_dealing.o \
	$(SHYSHYFEM)/conz_common.o \
	$(SHYSHYFEM)/ecological_dummy.o \
	$(SHYSHYFEM)/gotm_aux.o \
	$(FEMSRC)/contrib/check_dummy.o \
	$(FEMSRC)/contrib/ecological/mercury/mercury_dummy.o \
	$(FEMSRC)/contrib/ecological/bfm/bfm_dummy.o

# --- SHYFEM Ensemble Base Objects ---
SHYOBJENS = \
	$(SHYGENERIC)/dts.o \
	$(SHYGENERIC)/file.o \
	$(SHYGENERIC)/femfile.o \
	$(SHYGENERIC)/iso8601.o \
	$(SHYGENERIC)/convert.o

# --- EnKF Objects ---
ENKOBJ = m_analysis.o m_ensmean.o m_ensvar.o \
	 m_mean_preserving_rotation.o \
	 m_multa.o mod_anafunc.o

SAMPOBJ = m_pseudo2D.o m_fixsample2D.o m_newton2D.o m_newtonfunc2D.o \
	  m_sample2D.o mod_fftw3.o m_set_random_seed2.o m_randrot.o \
	  m_random.o
          
ENSOBJ = $(ENKOBJ) $(SAMPOBJ) subenkf.o mod_mod_err.o mod_enkf.o \
	 mod_ens_state.o mod_init_enkf.o mod_manage_obs.o mod_mod_states.o \
	 mod_obs_states.o mod_para.o sublocan.o mod_pestimation.o

# Executables list
EXES = main \
       make_antime_list \
       perturbe_fem_wind_press \
       perturbe_fem_heat \
       perturbe_fem_scalar \
       enKF2enKS \
       perturbe_state \
       perturbe_ts \
       fem2nc

MATH_LIBS = -llapack -lblas -lfftw3 -lpthread

# --- Build Targets & Compilation Rules ---

default: $(EXES)
all: $(EXES)

.SUFFIXES: .f .F90 .f90

# Phase 1: Object compilation applying OpenMP and NetCDF flags dynamically
%.o: %.f90
	$(FC) $(FFLAGS) $(ENKF_OMP_FLAG) $(INCLUDE) $(NC_FLAGS) -c $< -o $@

%.o: %.F90
	$(FC) $(FFLAGS) $(ENKF_OMP_FLAG) $(INCLUDE) $(NC_FLAGS) -c $< -o $@

%.o: %.f
	$(FC) $(FFLAGS) $(ENKF_OMP_FLAG) $(INCLUDE) -c $< -o $@

# Phase 2: Final Linking commands using structural recipes

main: main.o $(ENSOBJ) $(SHYOBJ)
	$(LINKER) $(FFLAGS) $(ENKF_OMP_FLAG) main.o $(ENSOBJ) $(SHYOBJ) $(LFLAGS) $(MATH_LIBS) $(NC_LIBS) -o $@

make_antime_list: make_antime_list.o $(SHYOBJENS)
	$(LINKER) $(FFLAGS) $(ENKF_OMP_FLAG) make_antime_list.o $(SHYOBJENS) $(LFLAGS) -o $@

perturbe_fem_wind_press: perturbe_fem_wind_press.o $(SAMPOBJ) $(SHYOBJENS)
	$(LINKER) $(FFLAGS) $(ENKF_OMP_FLAG) perturbe_fem_wind_press.o $(SAMPOBJ) $(SHYOBJENS) $(LFLAGS) $(MATH_LIBS) $(NC_LIBS) -o $@

perturbe_fem_heat: perturbe_fem_heat.o $(SAMPOBJ) $(SHYOBJENS)
	$(LINKER) $(FFLAGS) $(ENKF_OMP_FLAG) perturbe_fem_heat.o $(SAMPOBJ) $(SHYOBJENS) $(LFLAGS) $(MATH_LIBS) $(NC_LIBS) -o $@

perturbe_fem_scalar: perturbe_fem_scalar.o $(SAMPOBJ) $(SHYOBJENS)
	$(LINKER) $(FFLAGS) $(ENKF_OMP_FLAG) perturbe_fem_scalar.o $(SAMPOBJ) $(SHYOBJENS) $(LFLAGS) $(MATH_LIBS) $(NC_LIBS) -o $@

enKF2enKS: enKF2enKS.o $(SHYOBJ)
	$(LINKER) $(FFLAGS) $(ENKF_OMP_FLAG) enKF2enKS.o $(SHYOBJ) $(LFLAGS) $(MATH_LIBS) $(NC_LIBS) -o $@

perturbe_state: perturbe_state.o subenkf.o $(SAMPOBJ) $(SHYOBJ)
	$(LINKER) $(FFLAGS) $(ENKF_OMP_FLAG) perturbe_state.o subenkf.o $(SAMPOBJ) $(SHYOBJ) $(LFLAGS) $(MATH_LIBS) $(NC_LIBS) -o $@

perturbe_ts: perturbe_ts.o $(SHYOBJ)
	$(LINKER) $(FFLAGS) $(ENKF_OMP_FLAG) perturbe_ts.o $(SHYOBJ) $(LFLAGS) -o $@

fem2nc: fem2nc.o $(SHYOBJENS)
	$(LINKER) $(FFLAGS) $(ENKF_OMP_FLAG) fem2nc.o $(SHYOBJENS) $(LFLAGS) $(MATH_LIBS) $(NC_LIBS) -o $@

# --- Module Interdependencies ---
../../../lib/mod/m_ensmean.mod: m_ensmean.o
../../../lib/mod/m_ensvar.mod: m_ensvar.o
../../../lib/mod/m_fixsample2d.mod: m_fixsample2D.o
../../../lib/mod/m_mean_preserving_rotation.mod: m_mean_preserving_rotation.o 
../../../lib/mod/m_multa.mod: m_multa.o
../../../lib/mod/m_newton2d.mod: m_newton2D.o
../../../lib/mod/m_newtonfunc2d.mod: m_newtonfunc2D.o
../../../lib/mod/m_pseudo2d.mod: m_pseudo2D.o
../../../lib/mod/m_random.mod: m_random.o
../../../lib/mod/m_randrot.mod: m_randrot.o
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
../../../lib/mod/mod_pestimation.mod: mod_pestimation.o

m_analysis.o: ../../../lib/mod/m_ensmean.mod ../../../lib/mod/m_ensvar.mod ../../../lib/mod/m_multa.mod ../../../lib/mod/mod_anafunc.mod 
analysis2_EnOI.o: ../../../lib/mod/m_multa.mod
m_mean_preserving_rotation.o: ../../../lib/mod/m_randrot.mod
m_newton2D.o: ../../../lib/mod/m_newtonfunc2d.mod
m_pseudo2D.o: ../../../lib/mod/m_newton2d.mod ../../../lib/mod/mod_fftw3.mod 
m_sample2D.o: ../../../lib/mod/m_fixsample2d.mod ../../../lib/mod/m_pseudo2d.mod ../../../lib/mod/m_randrot.mod 

main.o: ../../../lib/mod/mod_enkf.mod ../../../lib/mod/mod_mod_err.mod ../../../lib/mod/mod_restart.mod ../../../lib/mod/m_set_random_seed2.mod ../../../lib/mod/m_analysis.mod ../../../lib/mod/mod_pestimation.mod
make_antime_list.o: ../../../lib/mod/iso8601.mod
mod_anafunc.o: ../../../lib/mod/m_ensmean.mod ../../../lib/mod/m_mean_preserving_rotation.mod ../../../lib/mod/m_multa.mod ../../../lib/mod/m_random.mod ../../../lib/mod/m_randrot.mod 

mod_enkf.o: ../../../lib/mod/levels.mod ../../../lib/mod/mod_ens_state.mod ../../../lib/mod/mod_manage_obs.mod ../../../lib/mod/mod_para.mod 
mod_ens_state.o: ../../../lib/mod/basin.mod ../../../lib/mod/levels.mod \
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

mod_init_enkf.o: ../../../lib/mod/iso8601.mod ../../../lib/mod/mod_para.mod 
mod_manage_obs.o: ../../../lib/mod/iso8601.mod ../../../lib/mod/levels.mod ../../../lib/mod/mod_ens_state.mod ../../../lib/mod/mod_init_enkf.mod ../../../lib/mod/mod_obs_states.mod ../../../lib/mod/mod_para.mod 
mod_mod_err.o: ../../../lib/mod/mod_ens_state.mod ../../../lib/mod/mod_init_enkf.mod ../../../lib/mod/mod_mod_states.mod 

subenkf.o: ../../../lib/mod/basin.mod ../../../lib/mod/levels.mod ../../../lib/mod/m_random.mod ../../../lib/mod/m_sample2d.mod ../../../lib/mod/mod_hydro.mod ../../../lib/mod/mod_para.mod ../../../lib/mod/mod_restart.mod ../../../lib/mod/shyfile.mod ../../../lib/mod/shympi.mod 
sublocan.o: ../../../lib/mod/basin.mod ../../../lib/mod/mod_enkf.mod ../../../lib/mod/mod_ens_state.mod ../../../lib/mod/mod_manage_obs.mod ../../../lib/mod/mod_para.mod ../../../lib/mod/mod_restart.mod ../../../lib/mod/m_analysis.mod

enKF2enKS.o: ../../../lib/mod/basin.mod ../../../lib/mod/levels.mod ../../../lib/mod/mod_conz.mod ../../../lib/mod/mod_geom_dynamic.mod ../../../lib/mod/mod_hydro.mod ../../../lib/mod/mod_hydro_vel.mod ../../../lib/mod/mod_gotm_aux.mod ../../../lib/mod/mod_restart.mod ../../../lib/mod/mod_ts.mod ../../../lib/mod/shympi.mod 
perturbe_state.o: ../../../lib/mod/basin.mod ../../../lib/mod/levels.mod ../../../lib/mod/mod_conz.mod ../../../lib/mod/mod_geom_dynamic.mod ../../../lib/mod/mod_hydro.mod ../../../lib/mod/mod_hydro_vel.mod ../../../lib/mod/mod_gotm_aux.mod ../../../lib/mod/mod_restart.mod ../../../lib/mod/mod_ts.mod ../../../lib/mod/shympi.mod 

clean:
	rm -f *.o *.mod $(EXES)

