
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#--------------------------------------------- 

DIR	= shynetcdf

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

INCLUDE = -I$(FEMDIR)/include
FFLAGS += $(INCLUDE)

FLAGS   = $(FFLAGS) $(INCLUDE) $(LFLAGS)

#--------------------------------------------- macros

EXES = nc2fem ts2nc

LIB_NETCDF = 
NCDF_TARGETS = lib
NETCDF_OBJ = elab_nc_dummy.o

ifeq ($(NETCDF),true)
  LIB_NETCDF = -L$(NETCDFDIR)/lib -lnetcdf -L$(NETCDFFDIR)/lib -lnetcdff
  NCDF_TARGETS = lib $(EXES)
  NETCDF_OBJ = \
        elab_nc.o \
        netcdf.o \
        netcdf_util.o
endif

#---------------------------------------------
# the next section checks if the variable NETCDF
# has changed between one call to make and another
# if it has changed, we must clean before compiling
# otherwise there will be a conflict of objects
# already defined in the library
#---------------------------------------------

CLEAN =
ifeq ("$(wildcard elab_nc_dummy.o)","")
  NETCDF_OLD = false
  ifeq ($(NETCDF),false)
    CLEAN = cleanall
  endif
else
  NETCDF_OLD = true
  ifeq ($(NETCDF),true)
    CLEAN = cleanall
  endif
endif

#--------------------------------------------- files

NC2FEM_OBJ = nc_util.o nc_tutil.o nc_domain.o nc_dim_coords.o nc2fem.o

ALL = $(NETCDF_OBJ)

SPECIAL =       Makefile

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#--------------------------------------------- EXES

LIB = $(LIBDIR)/lib$(LIBNAME).a

default: fem
all:	fem
fem:	$(CLEAN) $(NCDF_TARGETS) links
lib:	$(CLEAN) $(LIB)

$(LIB): $(ALL)
	ar rvs $(LIB) $(ALL)

LIBUTIL         = util
LIBSHYUTIL      = shyutil
LIBSHYMPI       = shympi
LIBSHYUTILMPI   = shyutilmpi
LIBNETCDF       = shynetcdf

LIBLUT  = -L$(LIBDIR) -l$(LIBUTIL)
LIBFUT  = $(LIBDIR)/lib$(LIBUTIL).a

LIBLSU  = -L$(LIBDIR) -l$(LIBSHYUTIL)
LIBFSU  = $(LIBDIR)/lib$(LIBSHYUTIL).a

LIBLMP  = -L$(LIBDIR) -l$(LIBSHYMPI)
LIBFMP  = $(LIBDIR)/lib$(LIBSHYMPI).a

LIBLUM  = -L$(LIBDIR) -l$(LIBSHYUTILMPI)
LIBFUM  = $(LIBDIR)/lib$(LIBSHYUTILMPI).a

LIBLNC  = -L$(LIBDIR) -l$(LIBNETCDF)
LIBFNC  = $(LIBDIR)/lib$(LIBNETCDF).a

LIBL = $(LIBLNC) $(LIBLUM) $(LIBLMP) $(LIBLSU) $(LIBLUT) $(LIB_NETCDF)
LIBF = $(LIBFNC) $(LIBFUM) $(LIBFMP) $(LIBFSU) $(LIBFUT)

nc2fem: $(NC2FEM_OBJ) $(LIBF)
	$(F77) $(FLAGS) $(NC2FEM_OBJ) $(LIBL) -o $@

ts2nc: ts2nc.o
	$(F77) $(FLAGS) ts2nc.o $(LIBL) -o $@

links:
	@$(FEMBIN)/symlink $(FEMBIN) $(EXES)

show:
	@echo "NETCDF_OLD = $(NETCDF_OLD)"
	@echo "CLEAN      = $(CLEAN)"

#--------------------------------------------- test compile

test: lib main_dummy

main_dummy:
	@echo " write(6,*) 'all dependencies ok'\n end" > main_dummy.f90
	$(F77) $(FFLAGS) -c $@.f90 
	$(F77) $(FLAGS) $(ALL) $(LIBL) $@.o -o $@
	@echo "...running main_dummy"
	@./main_dummy
	@-rm -f main_dummy main_dummy.f90 main_dummy.o

#--------------------------------------------- clean

cleanall: clean cleanlib

clean: cleanobj cleanvar cleanexe cleanlinks

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

cleanlinks:
	@$(FEMBIN)/symlink -remove $(FEMBIN) $(EXES)

cleanexe:
	-rm -f $(EXES)

cleanlib:
	-rm -f $(LIB)

#--------------------------------------------- various

check:
	~/shyfem/fem3d/bin/compile.sh -single *.f90
	~/shyfem/fem3d/bin/usemod.pl *.f90

list:
	pwd; ls

save:   zip
zip:	cleanall
	-zip $(DIR) $(SAVEFILES)

depend:
	@-$(MAKEDEPEND) -moddir=$(MODDIR) *.f90

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
	$(FLC) $< , $(<,B).OBJ $(FFLAGS)

#--------------------------------------------- makedepend

# DO NOT DELETE THIS LINE -- make depend depends on it.

../../../lib/mod/nc_domain.mod: nc_domain.o
../../../lib/mod/nc_interpol.mod: nc_domain.o
../../../lib/mod/nc_time.mod: nc_tutil.o
../../../lib/mod/ncnames.mod: nc_dim_coords.o
../../../lib/mod/netcdf_out.mod: elab_nc.o
../../../lib/mod/netcdf_params.mod: netcdf.o
elab_nc.o: ../../../lib/mod/basin.mod ../../../lib/mod/levels.mod \
		../../../lib/mod/mod_depth.mod 
elab_nc_dummy.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/levels.mod 
nc2fem.o: ../../../lib/mod/clo.mod ../../../lib/mod/iso8601.mod \
		../../../lib/mod/nc_time.mod \
		../../../lib/mod/ncnames.mod \
		../../../lib/mod/shyfem_strings.mod 
nc_tutil.o: ../../../lib/mod/iso8601.mod
netcdf.o: ../../../lib/mod/shyfem_strings.mod netcdf.inc
netcdf_util.o: ../../../lib/mod/basin.mod
ts2nc.o: ../../../lib/mod/basin.mod ../../../lib/mod/evgeom.mod \
		../../../lib/mod/levels.mod \
		../../../lib/mod/mod_depth.mod \
		../../../lib/mod/simul.mod 

