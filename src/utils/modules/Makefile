
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
	mod_area.o \
	mod_bclfix.o \
	mod_diff_aux.o \
	mod_diffusion.o \
	mod_fluidmud.o \
	mod_flux.o \
	mod_hydro.o \
	mod_hydro_baro.o \
	mod_hydro_print.o \
	mod_hydro_vel.o \
	mod_internal.o \
	mod_layer_thickness.o \
	mod_meteo.o \
	mod_nohyd.o \
	mod_nudging.o \
	mod_roughness.o \
	mod_sedim.o \
	mod_sinking.o \
	mod_system_matrix.o \
	mod_ts.o \
	mod_waves.o \
	mod_shyfem.o \

OBJ  = \

OBSOLETE = \

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

../../../lib/mod/mod_area.mod: mod_area.o
../../../lib/mod/mod_bclfix.mod: mod_bclfix.o
../../../lib/mod/mod_diff_aux.mod: mod_diff_aux.o
../../../lib/mod/mod_diff_visc_fric.mod: mod_diffusion.o
../../../lib/mod/mod_fluidmud.mod: mod_fluidmud.o
../../../lib/mod/mod_flux.mod: mod_flux.o
../../../lib/mod/mod_hydro.mod: mod_hydro.o
../../../lib/mod/mod_hydro_baro.mod: mod_hydro_baro.o
../../../lib/mod/mod_hydro_print.mod: mod_hydro_print.o
../../../lib/mod/mod_hydro_vel.mod: mod_hydro_vel.o
../../../lib/mod/mod_internal.mod: mod_internal.o
../../../lib/mod/mod_layer_thickness.mod: mod_layer_thickness.o
../../../lib/mod/mod_meteo.mod: mod_meteo.o
../../../lib/mod/mod_nohyd.mod: mod_nohyd.o
../../../lib/mod/mod_nudging.mod: mod_nudging.o
../../../lib/mod/mod_roughness.mod: mod_roughness.o
../../../lib/mod/mod_sediment.mod: mod_sedim.o
../../../lib/mod/mod_sediment_para.mod: mod_sedim.o
../../../lib/mod/mod_shyfem.mod: mod_shyfem.o
../../../lib/mod/mod_sinking.mod: mod_sinking.o
../../../lib/mod/mod_system.mod: mod_system_matrix.o
../../../lib/mod/mod_system_interface.mod: mod_system_matrix.o
../../../lib/mod/mod_ts.mod: mod_ts.o
../../../lib/mod/mod_waves.mod: mod_waves.o

