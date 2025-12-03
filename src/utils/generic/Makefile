
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2023  The SHYFEM team
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#--------------------------------------------- 

DIR	= util

#--------------------------------------------- just for default
FEMDIR = ../../..
F77    = g77
FFLAGS =  -O -Wall -g
LINKER = $(F77)
#--------------------------------------------- include Rules.make

include $(FEMDIR)/Rules.make

#--------------------------------------------- libraries

LIBNAME = $(DIR)
FEMBIN  = $(FEMDIR)/bin
BINDIR  = $(FEMDIR)/bin
LIBDIR  = $(FEMDIR)/lib
INCDIR = $(FEMDIR)/include
INCLUDE = -I$(INCDIR)

#--------------------------------------------- macros

MOD  = \
	para.o \

OBJ  = \
	chk_NaN.o \
	clo.o \
	convert.o \
	dts.o \
	file.o \
	geometrical.o \
	intp.o \
	iso8601.o \
	mat_util.o \
	nls.o \
	projection.o \
	random.o \
	random_gen.o \
	round.o \
	grd.o \
	subsort.o \
	sort_util.o \
	util.o \
	utility.o \
	vec_util.o \
	femfile.o \
	color.o \
	to_linear.o \
	dates.o \
	unitutils.o \
	tsfile.o \
	heat_util.o \
	connections.o \
	time_util.o \
	vert_intp.o \
	int_to_cart.o \
	bit10.o \
	nodelist.o \
	color_utils.o \
	flood.o \
	queue.o \
	stack.o \
	reallocate.o \
	histogram.o \
	errorstop.o \
	compatibility.o \
	spline.o \

OBSOLETE = \
	expand.o \
	strings.o \

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

$(LIB): $(ALL)
	ar rvs $(LIB) $(ALL)

LIBL =
LIBF =

#--------------------------------------------- test compile

test: lib main_dummy
        
main_dummy:
	@echo " write(6,*) 'all dependencies ok'\n end" > main_dummy.f90
	$(F77) $(FFLAGS) -c $@.f90 
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $(ALL) $(LIBL) $@.o -o $@
	@echo "...running main_dummy"
	@./main_dummy
	@-rm -f main_dummy main_dummy.f90 main_dummy.o

clo:
	$(F77) $@.f90

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
	$(FEMBIN)/fem3d/compile.sh -single *.f90
	$(FEMBIN)/fem3d/compile.sh *.f90
	$(FEMBIN)/fem3d/usemod.pl *.f90

check_utils:
	@$(FEMBIN)/fem3d/check_utils.sh

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
	$(F77) -c $(INCLUDE) $(FFLAGS) $<

.F.o:
	$(F77) -c -cpp $(FFLAGS) $<

.for.obj:
	$(FLC) $< , $(<,B).OBJ $(FLFLAGS)

#--------------------------------------------- makedepend

# DO NOT DELETE THIS LINE -- make depend depends on it.

../../../lib/mod/clo.mod: clo.o
../../../lib/mod/custom_dates.mod: dates.o
../../../lib/mod/dts.mod: dts.o
../../../lib/mod/fem_file.mod: femfile.o
../../../lib/mod/femtime.mod: time_util.o
../../../lib/mod/grd.mod: grd.o
../../../lib/mod/iso8601.mod: iso8601.o
../../../lib/mod/list.mod: list.o
../../../lib/mod/mod_color.mod: color.o
../../../lib/mod/mod_compatibility.mod: compatibility.o
../../../lib/mod/mod_connect.mod: connections.o
../../../lib/mod/mod_debug.mod: chk_NaN.o
../../../lib/mod/mod_error_stop.mod: errorstop.o
../../../lib/mod/mod_flood.mod: flood.o
../../../lib/mod/mod_histo.mod: histogram.o
../../../lib/mod/mod_random.mod: random.o
../../../lib/mod/mod_reallocate.mod: reallocate.o
../../../lib/mod/mod_sort.mod: subsort.o
../../../lib/mod/nls.mod: nls.o
../../../lib/mod/para.mod: para.o
../../../lib/mod/projection.mod: projection.o
../../../lib/mod/queue.mod: queue.o
../../../lib/mod/scan_string.mod: convert.o
../../../lib/mod/shyelab_unit.mod: unitutils.o
../../../lib/mod/stack.mod: stack.o
dates.o: ../../../lib/mod/iso8601.mod
flood.o: ../../../lib/mod/mod_connect.mod \
		../../../lib/mod/queue.mod 
nls.o: ../../../lib/mod/para.mod
time_util.o: ../../../lib/mod/iso8601.mod
tsfile.o: ../../../lib/mod/iso8601.mod

