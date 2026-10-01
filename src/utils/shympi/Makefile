
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#--------------------------------------------- 

DIR	= shympi

#--------------------------------------------- just for default
FEMDIR = ../../..
FEMBIN = $(FEMDIR)/bin
F77    = g77
FFLAGS =  -O -Wall -g
LINKER = $(F77)
#--------------------------------------------- include Rules.make

include $(FEMDIR)/Rules.make

#--------------------------------------------- libraries

LIBNAME = $(DIR)
LIBDIR  = $(FEMDIR)/lib

#--------------------------------------------- macros

METIS_OBJ = shyparts_dummy.o
ifeq ($(PARTS),SDDA)
  METIS_OBJ = shyparts_basutil.o shyparts_custom.o shyparts_util.o
endif
ifeq ($(PARTS),$(filter $(PARTS),METIS PARMETIS))
  METIS_OBJ = shyparts_basutil.o shyparts_metis.o shyparts_util.o
  LIBG_METIS = -L$(METISDIR)/lib -lmetis
endif

ifeq ($(PARALLEL_MPI),NODE)
  MOD = mpi_node.o mpi_node_internal.o
  MPI_FILE = MPI_NODE
else
  MOD = mpi_dummy.o
  MPI_FILE = MPI_NONE
endif

#---- next is needed if PARALLEL_MPI has changed -> recreate library

MPI_TARGET = $(shell [ ! -f $(MPI_FILE) ] && echo reset )

#--------------------------------------------- objects

OBJ  = \
	mpi_basin.o \
	mpi_buffer.o \
	mpi_custom.o \
	mpi_debug.o \
	mpi_general.o \
	mpi_ghost.o \
	mpi_time.o \
	mpi_tripple.o \
	mpi_tvd.o \
	mpi_check.o \

ALL = $(MOD) $(OBJ) $(METIS_OBJ)

SPECIAL =       Makefile

EXES = check_shympi_debug check_shympi_check

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#--------------------------------------------- EXES

LIB = $(LIBDIR)/lib$(LIBNAME).a

default: fem
all:	fem
fem: shyfem
shyfem: lib $(EXES) links
#fem: lib
lib: $(MPI_TARGET) $(LIB) nothing

nothing:
	@true

reset: cleanall
	@echo "resetting library to $(MPI_FILE)"
	touch $(MPI_FILE)

$(LIB): $(ALL)
	ar rvs $(LIB) $(ALL)

LIBUTIL     = util
LIBSHYUTIL  = shyutil

LIBLUT  = -L$(LIBDIR) -l$(LIBUTIL)
LIBFUT  = $(LIBDIR)/lib$(LIBUTIL).a

LIBLSU  = -L$(LIBDIR) -l$(LIBSHYUTIL)
LIBFSU  = $(LIBDIR)/lib$(LIBSHYUTIL).a

LIBL = $(LIBLSU) $(LIBLUT) $(LIBG_METIS)
LIBF = $(LIBFSU) $(LIBFUT)

links:
	@$(FEMBIN)/symlink $(FEMBIN) $(EXES)

check_shympi_debug: $(LIBF) check_shympi_debug.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $@.o $(LIBL) -o $@

check_shympi_check: $(LIBF) check_shympi_check.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $@.o $(LIBL) -o $@

#--------------------------------------------- test compile

test: lib main_dummy

main_dummy:
	@echo " write(6,*) 'all dependencies ok'\n end" > main_dummy.f90
	$(F77) $(FFLAGS) -c $@.f90 
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $(ALL) $(LIBL) $@.o -o $@
	@echo "...running main_dummy"
	@./main_dummy
	@-rm -f main_dummy main_dummy.f90 main_dummy.o

#--------------------------------------------- var targets

MPI_TARGET = $(shell [ ! -f MPI_NODE ] && echo reset )

mpi_node_internal.o: mpi_node_internal.f90
	$(F77) -c $(FFLAG_SPECIAL) mpi_node_internal.f90

ifeq ($(COMPILER_PROFILE),CHECK)
shyparts_metis.o: shyparts_metis.f90
	$(F77) -c $(FFLAG_SPECIAL) shyparts_metis.f90
endif

#--------------------------------------------- clean

cleanall: clean cleanlib
	-rm -f MPI_*

clean: cleanobj cleanvar cleanlinks

cleanobj:
	-rm -f *.o
	-rm -f *.mod
	-rm -f *_genmod.f90

cleanvar:
	-rm -f $(EXES)
	-rm -f main_dummy main_dummy.f90
	-rm -f a.out core
	-rm -f *.bak *.tmp
	-rm -f ggg hhh
	-rm -f mon.out nohup.out
	-rm -f fort.*
	-rm -f gmon.out
	-rm -f out.ps plot.ps
	-rm -f MPI_NODE MPI_NONE

cleanlinks:
	@$(FEMBIN)/symlink -remove $(FEMBIN) $(EXES)

cleanlib:
	-rm -f $(LIB)

#--------------------------------------------- various

depend:
	@-$(MAKEDEPEND) -moddir=$(MODDIR) *.f90

check:
	#~/shyfem/fem3d/bin/compile.sh -single *.f90
	~/shyfem/fem3d/bin/usemod.pl *.f90

list:
	pwd; ls

save:   zip
zip:	cleanall
	-zip $(DIR) $(SAVEFILES)

#--------------------------------------------- rules

.SUFFIXES: .F .f90

.c.o:
	$(CC) -c $(CFLAGS) $*.c

.f.o:
	$(F77) -c $(FFLAGS) $<

.f90.o:
	$(F77) -c $(FFLAGS) $<

.F.o:
	$(F77) -c -cpp $(FFLAGS) $<

.for.obj:
	$(FLC) $< , $(<,B).OBJ $(FLFLAGS)

#--------------------------------------------- makedepend

# DO NOT DELETE THIS LINE -- make depend depends on it.

../../../lib/mod/mod_check_shympi_check.mod: check_shympi_check.o
../../../lib/mod/mod_save_index.mod: shyparts_basutil.o
../../../lib/mod/mod_shympi_check.mod: mpi_check.o
../../../lib/mod/mod_shympi_debug.mod: check_shympi_debug.o
../../../lib/mod/shympi_debug.mod: mpi_debug.o
../../../lib/mod/shympi_internal.mod: mpi_node_internal.o
../../../lib/mod/shympi_time.mod: mpi_time.o
../../../lib/mod/shympi_tripple.mod: mpi_tripple.o
../../../lib/mod/shympi_tvd.mod: mpi_tvd.o
check_shympi_check.o: ../../../lib/mod/clo.mod
check_shympi_debug.o: ../../../lib/mod/clo.mod
mpi_basin.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/shympi.mod 
mpi_check.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/levels.mod \
		../../../lib/mod/shympi.mod 
mpi_custom.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/shympi.mod 
mpi_debug.o: ../../../lib/mod/shympi.mod
mpi_dummy.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/levels.mod 
mpi_general.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/shympi.mod 
mpi_ghost.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/shympi.mod 
mpi_node.o: ../../../lib/mod/basin.mod ../../../lib/mod/levels.mod \
		../../../lib/mod/mod_compatibility.mod 
mpi_node_internal.o: ../../../lib/mod/shympi.mod
mpi_tripple.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/levels.mod \
		../../../lib/mod/mod_geom.mod \
		../../../lib/mod/mod_shyfem.mod \
		../../../lib/mod/shympi.mod 
mpi_tvd.o: ../../../lib/mod/basin.mod ../../../lib/mod/levels.mod \
		../../../lib/mod/mod_quad_tree.mod \
		../../../lib/mod/mod_tvd.mod \
		../../../lib/mod/shympi.mod 
shyparts_basutil.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/grd.mod \
		../../../lib/mod/mod_connect.mod \
		../../../lib/mod/mod_depth.mod \
		../../../lib/mod/mod_geom.mod \
		../../../lib/mod/shympi.mod 
shyparts_custom.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/mod_color.mod \
		../../../lib/mod/mod_connect.mod \
		../../../lib/mod/mod_flood.mod \
		../../../lib/mod/mod_sort.mod \
		../../../lib/mod/queue.mod 
shyparts_dummy.o: ../../../lib/mod/basin.mod
shyparts_metis.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/mod_geom.mod \
		../../../lib/mod/shympi.mod 
shyparts_util.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/grd.mod 

