
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#--------------------------------------------- 

DIR	= shyutilmpi

#--------------------------------------------- just for default
FEMDIR = ../../..
FEMBIN = $(FEMDIR)/bin
F77    = g77
FFLAGS =  -O -Wall -g
LINKER = $(F77)
#--------------------------------------------- include Rules.make

include $(FEMDIR)/Rules.make

#--------------------------------------------- names

LIBNAME = $(DIR)
LIBDIR  = $(FEMDIR)/lib

INCLUDE = -I$(FEMDIR)/include
FFLAGS += $(INCLUDE)

#--------------------------------------------- macros

MOD  = \

OBJ  = \
	topo_update.o \
	bnd_admin.o \
	depth_util.o \
	shy_util.o \
	shy_util1.o \
	intp_fem_file.o \
	coordinates.o \
	io_off.o \
	trace_point.o \
	info_output.o \
	mod_restart.o \

OBSOLETE = \
	evgeom.o \
	find.o \

ALL = $(MOD) $(OBJ)

SPECIAL =       Makefile

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#--------------------------------------------- EXES

LIB = $(LIBDIR)/lib$(LIBNAME).a

default: lib
all:	lib
fem: shyfem
shyfem: lib
lib:	$(LIB) nothing

nothing:
	@true

$(LIB): $(MOD) $(OBJ)
	ar rvs $(LIB) $(MOD) $(OBJ)

#--------------------------------------------- libraries

LIBUTIL     = util
LIBSHYUTIL  = shyutil
LIBSHYMPI   = shympi

LIBLUT  = -L$(LIBDIR) -l$(LIBUTIL)
LIBFUT  = $(LIBDIR)/lib$(LIBUTIL).a

LIBLSU  = -L$(LIBDIR) -l$(LIBSHYUTIL)
LIBFSU  = $(LIBDIR)/lib$(LIBSHYUTIL).a

LIBLMP  = -L$(LIBDIR) -l$(LIBSHYMPI)
LIBFMP  = $(LIBDIR)/lib$(LIBSHYMPI).a

LIBL = $(LIBLMP) $(LIBLSU) $(LIBLUT)
LIBF = $(LIBFMP) $(LIBFSU) $(LIBFUT)

#--------------------------------------------- test compile

test: lib main_dummy

main_dummy:
	@echo " write(6,*) 'all dependencies ok'\n end" > main_dummy.f90
	$(F77) $(FFLAGS) -c $@.f90 
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $(ALL) $(LIBL) $@.o -o $@
	@echo "...running main_dummy"
	@./main_dummy
	@-rm -f main_dummy main_dummy.f90 main_dummy.o

test_iff: test_iff.o lib
	$(F77) $(FFLAGS) -c $@.f90 
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $(ALL) $(LIBL) $@.o -o $@

#--------------------------------------------- clean

cleanall: clean cleanlib

clean: cleanobj cleanvar

cleanobj:
	-rm -f *.o
	-rm -f *.mod
	-rm -f *_genmod.f90

cleanvar:
	-rm -f main_dummy main_dummy.f90
	-rm -f a.out core
	-rm -f *.bak *.tmp
	-rm -f ggg hhh
	-rm -f mon.out nohup.out
	-rm -f fort.*
	-rm -f gmon.out
	-rm -f out.ps plot.ps

cleanlib:
	-rm -f $(LIB)

#--------------------------------------------- various

depend:
	@-$(MAKEDEPEND) -moddir=$(MODDIR) *.f90

check:
	~/shyfem/fem3d/bin/compile.sh -single *.f90
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

../../../lib/mod/coordinates.mod: coordinates.o
../../../lib/mod/intp_fem_file.mod: intp_fem_file.o
../../../lib/mod/mod_info_output.mod: info_output.o
../../../lib/mod/mod_offline.mod: io_off.o
../../../lib/mod/mod_restart.mod: mod_restart.o
../../../lib/mod/mod_trace_point.mod: trace_point.o
../../../lib/mod/shyutil.mod: shy_util1.o
bnd_admin.o: ../../../lib/mod/mod_bnd.mod \
		../../../lib/mod/mod_bound_geom.mod \
		../../../lib/mod/nls.mod ../../../lib/mod/para.mod \
		../../../lib/mod/shympi.mod 
coordinates.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/shympi.mod 
depth_util.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/mod_depth.mod \
		../../../lib/mod/shympi.mod 
info_output.o: ../../../lib/mod/shympi.mod
intp_fem_file.o: ../../../lib/mod/shympi.mod
io_off.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/mod_trace_point.mod \
		../../../lib/mod/shympi.mod 
mod_restart.o: ../../../lib/mod/shympi.mod
shy_util.o: ../../../lib/mod/basin.mod ../../../lib/mod/levels.mod \
		../../../lib/mod/mod_trace_point.mod \
		../../../lib/mod/shyfile.mod \
		../../../lib/mod/shympi.mod \
		../../../lib/mod/simul.mod 
shy_util1.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/evgeom.mod \
		../../../lib/mod/levels.mod \
		../../../lib/mod/mod_depth.mod \
		../../../lib/mod/mod_sort.mod 
topo_update.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/evgeom.mod \
		../../../lib/mod/mod_geom.mod \
		../../../lib/mod/mod_geom_dynamic.mod \
		../../../lib/mod/shympi.mod 
trace_point.o: ../../../lib/mod/shympi.mod

