
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2023  The SHYFEM team
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

DIR     = shyelab

#--------------------------------------------- just for default
FEMDIR = ../../..
FEMBIN = $(FEMDIR)/bin
LIBDIR = $(FEMDIR)/lib
SRCDIR = $(FEMDIR)/src
F77    = g77
FFLAGS =  -O -Wall -g
LINKER = $(F77)
LFLAGS:=
INCLUDE = -I$(FEMDIR)/include/
#---------------------------------------------

include $(FEMDIR)/Rules.make

#---------------------------------------------

LIB_NETCDF =
ifeq ($(NETCDF),true)
  LIB_NETCDF = -L$(NETCDFDIR)/lib -lnetcdf -L$(NETCDFFDIR)/lib -lnetcdff
endif

LIBUTIL         = util
LIBSHYUTIL      = shyutil
LIBSHYMPI       = shympi
LIBSHYUTILMPI   = shyutilmpi
LIBNETCDF       = shynetcdf

LIBL     = -L$(LIBDIR)

LIBL    += -l$(LIBNETCDF)
LIBF    += $(LIBDIR)/lib$(LIBNETCDF).a

LIBL    += -l$(LIBSHYUTILMPI)
LIBF    += $(LIBDIR)/lib$(LIBSHYUTILMPI).a

LIBL    += -l$(LIBSHYMPI)
LIBF    += $(LIBDIR)/lib$(LIBSHYMPI).a

LIBL    += -l$(LIBSHYUTIL)
LIBF    += $(LIBDIR)/lib$(LIBSHYUTIL).a

LIBL    += -l$(LIBUTIL)
LIBF    += $(LIBDIR)/lib$(LIBUTIL).a

LIBL    += $(LIB_NETCDF)

#$(info FFLAGS = $(FFLAGS) )
#$(info LFLAGS = $(LFLAGS) )

#---------------------------------------------

EXES = shyelab extelab flxelab lgrelab tselab femelab shyproj shyfile \
	shystrings femcombine dbg_flx dbg_box \
	femadd femdate

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#----------------------------------------------------------------------------

LFLAGS	+= $(LIB_NETCDF) $(LIB_MPI)

#----------------------------------------------------------------------------

OBJ := \
	elabutil.o \
	elabtime.o \
	elab_off.o \
	elab_gis.o \
	elab_regular.o \
	elab_resample.o \
	elab_check.o \
	shyelab_output.o \
	shyelab_util.o \
	shyelab_area.o \
	shyelab_average.o \
	shyelab_nodes.o \
	shyelab_proj.o \
	shyelab_extract.o \
	shyelab1.o \
	extelab1.o \
	flxelab1.o \
	tselab1.o \
	femelab1.o \
	lgrelab1.o \
	shy_util2.o \

#----------------------------------------------------------------------------

default: fem
fem: $(EXES) links

links:
	@$(FEMBIN)/symlink $(FEMBIN) $(EXES)

shyelab: $(OBJ) $(LIBF) shyelab.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $(OBJ) $@.o $(LIBL) -o $@

extelab: $(OBJ) $(LIBF) extelab.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $(OBJ) $@.o $(LIBL) -o $@

flxelab: $(OBJ) $(LIBF) flxelab.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $(OBJ) $@.o $(LIBL) -o $@

lgrelab: $(OBJ) $(LIBF) lgrelab.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $(OBJ) $@.o $(LIBL) -o $@

femelab: $(OBJ) $(LIBF) femelab.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $(OBJ) $@.o $(LIBL) -o $@

femcombine: $(OBJ) $(LIBF) femcombine.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $(OBJ) $@.o $(LIBL) -o $@

tselab: $(OBJ) $(LIBF) tselab.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $(OBJ) $@.o $(LIBL) -o $@

shyproj: $(LIBF) shyproj.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $@.o $(LIBL) -o $@

shyfile: $(LIBF) shyfile.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $@.o $(LIBL) -o $@

shystrings: $(LIBF) shystrings.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $@.o $(LIBL) -o $@

dbg_flx: $(LIBF) dbg_flx.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $@.o $(LIBL) -o $@

dbg_box: $(LIBF) dbg_box.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $@.o $(LIBL) -o $@

femadd: $(OBJ) $(LIBF) femadd.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $(OBJ) $@.o $(LIBL) -o $@

femdate: $(OBJ) $(LIBF) femdate.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $(OBJ) $@.o $(LIBL) -o $@

#----------------------------------------------------------------------------

clean: cleanobj cleanlinks cleanexe cleanvar

cleanvar:
	-rm -f *.tmp *.bak

cleanobj:
	-rm -f *.o
	-rm -f *.mod
	-rm -f *_genmod.f90

cleanlinks:
	@$(FEMBIN)/symlink -remove $(FEMBIN) $(EXES)

cleanexe:
	-rm -f $(EXES)

cleanall: clean

#----------------------------------------------------------------------------

depend:
	@-$(MAKEDEPEND) -moddir=$(MODDIR) *.f90

#----------------------------------------------------------------------------

.SUFFIXES: .F .f90

.c.o:
	$(CC) -c $(CFLAGS) $<

.f.o:
	$(F77) -c $(FFLAGS) $<

.f90.o:
	${F77} -c $(FFLAGS) $<

.F.o:
	$(F77) -c -cpp $(FFLAGS) $<

#----------------------------------------------------------------------------

# DO NOT DELETE THIS LINE -- make depend depends on it.

../../../lib/mod/elabtime.mod: elabtime.o
../../../lib/mod/elabutil.mod: elabutil.o
../../../lib/mod/extelab_strings.mod: extelab1.o
../../../lib/mod/mod_dbg_box.mod: dbg_box.o
../../../lib/mod/mod_dbg_flx.mod: dbg_flx.o
../../../lib/mod/shy_extract.mod: shyelab_extract.o
../../../lib/mod/shyelab_out.mod: shyelab_output.o
dbg_box.o: ../../../lib/mod/clo.mod \
		../../../lib/mod/mod_error_stop.mod 
dbg_flx.o: ../../../lib/mod/clo.mod \
		../../../lib/mod/mod_error_stop.mod 
elab_check.o: ../../../lib/mod/iso8601.mod
elab_gis.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/shyfem_strings.mod 
elab_off.o: ../../../lib/mod/basin.mod ../../../lib/mod/evgeom.mod \
		../../../lib/mod/levels.mod 
elab_regular.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/levels.mod 
elab_resample.o: ../../../lib/mod/basin.mod
elabutil.o: ../../../lib/mod/clo.mod \
		../../../lib/mod/shyfem_strings.mod 
extelab1.o: ../../../lib/mod/basin.mod ../../../lib/mod/clo.mod \
		../../../lib/mod/elabtime.mod \
		../../../lib/mod/elabutil.mod \
		../../../lib/mod/evgeom.mod \
		../../../lib/mod/levels.mod \
		../../../lib/mod/mod_depth.mod \
		../../../lib/mod/shyfem_strings.mod 
femadd.o: ../../../lib/mod/clo.mod ../../../lib/mod/fem_util.mod
femcombine.o: ../../../lib/mod/clo.mod \
		../../../lib/mod/fem_util.mod 
femdate.o: ../../../lib/mod/clo.mod ../../../lib/mod/fem_util.mod
femelab1.o: ../../../lib/mod/clo.mod ../../../lib/mod/elabtime.mod \
		../../../lib/mod/elabutil.mod \
		../../../lib/mod/iso8601.mod \
		../../../lib/mod/shyfem_strings.mod 
femelab2.o: ../../../lib/mod/clo.mod ../../../lib/mod/elabtime.mod \
		../../../lib/mod/elabutil.mod \
		../../../lib/mod/fem_util.mod \
		../../../lib/mod/shyfem_strings.mod 
flxelab1.o: ../../../lib/mod/basin.mod ../../../lib/mod/clo.mod \
		../../../lib/mod/elabtime.mod \
		../../../lib/mod/elabutil.mod \
		../../../lib/mod/evgeom.mod \
		../../../lib/mod/levels.mod \
		../../../lib/mod/mod_depth.mod \
		../../../lib/mod/shyfem_strings.mod 
lgrelab1.o: ../../../lib/mod/basin.mod ../../../lib/mod/clo.mod \
		../../../lib/mod/elabtime.mod \
		../../../lib/mod/elabutil.mod \
		../../../lib/mod/evgeom.mod \
		../../../lib/mod/levels.mod \
		../../../lib/mod/mod_depth.mod \
		../../../lib/mod/shyelab_out.mod \
		../../../lib/mod/shyfem_strings.mod \
		../../../lib/mod/shyfile.mod \
		../../../lib/mod/shympi.mod \
		../../../lib/mod/shyutil.mod 
shyelab.o: ../../../lib/mod/clo.mod ../../../lib/mod/elabutil.mod \
		../../../lib/mod/shyfile.mod 
shyelab1.o: ../../../lib/mod/basin.mod ../../../lib/mod/clo.mod \
		../../../lib/mod/custom_dates.mod \
		../../../lib/mod/elabtime.mod \
		../../../lib/mod/elabutil.mod \
		../../../lib/mod/evgeom.mod \
		../../../lib/mod/levels.mod \
		../../../lib/mod/mod_depth.mod \
		../../../lib/mod/shy_extract.mod \
		../../../lib/mod/shyfile.mod \
		../../../lib/mod/shympi.mod \
		../../../lib/mod/shyutil.mod 
shyelab_area.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/elabutil.mod 
shyelab_nodes.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/elabtime.mod \
		../../../lib/mod/elabutil.mod \
		../../../lib/mod/levels.mod \
		../../../lib/mod/mod_depth.mod \
		../../../lib/mod/shyfem_strings.mod 
shyelab_output.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/elabtime.mod \
		../../../lib/mod/elabutil.mod \
		../../../lib/mod/levels.mod \
		../../../lib/mod/mod_offline.mod \
		../../../lib/mod/shyfem_strings.mod \
		../../../lib/mod/shyfile.mod 
shyelab_proj.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/clo.mod \
		../../../lib/mod/coordinates.mod \
		../../../lib/mod/projection.mod \
		../../../lib/mod/shympi.mod 
shyelab_util.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/evgeom.mod \
		../../../lib/mod/levels.mod \
		../../../lib/mod/mod_depth.mod 
shyfile.o: ../../../lib/mod/clo.mod ../../../lib/mod/elabutil.mod \
		../../../lib/mod/shyfile.mod 
shyproj.o: ../../../lib/mod/basin.mod ../../../lib/mod/clo.mod \
		../../../lib/mod/grd.mod \
		../../../lib/mod/projection.mod 
shystrings.o: ../../../lib/mod/clo.mod \
		../../../lib/mod/shyfem_strings.mod 
tselab1.o: ../../../lib/mod/clo.mod ../../../lib/mod/elabtime.mod \
		../../../lib/mod/elabutil.mod \
		../../../lib/mod/shyfem_strings.mod 

