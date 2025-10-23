
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2023  The SHYFEM team
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

DIR     = shyplot

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
LIBPOST         = post

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

LIBLPS  = -L$(LIBDIR) -l$(LIBPOST)
LIBFPS  = $(LIBDIR)/lib$(LIBPOST).a

#LIBL = $(LIBLNC) $(LIBLUM) $(LIBLMP) $(LIBLSU) $(LIBLUT) $(LIB_NETCDF)
LIBL = $(LIBLNC) $(LIBLUM) $(LIBLMP) $(LIBLSU) $(LIBLUT) $(LIBLPS)
LIBF = $(LIBFNC) $(LIBFUM) $(LIBFMP) $(LIBFSU) $(LIBFUT)

#---------------------------------------------

EXES = shyplot
EXELINKS =      $(FEMBIN) $(EXES)

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#----------------------------------------------------------------------------

LFLAGS += $(LIB_NETCDF) $(LIB_MPI)

#----------------------------------------------------------------------------

OBJ := \
	col.o \
	supcust.o \
	mod_hydro_plot.o \
	$(FEMDIR)/src/tools/shyelab/elabtime.o \
	$(FEMDIR)/src/tools/shyelab/elabutil.o \
	$(FEMDIR)/src/tools/shyelab/shy_util2.o \
	plotutil.o \
	supout.o \
	supbas.o \
	supano.o \
	supdep.o \
	supint.o \
	supiso.o \
	suplag.o \
	suplin.o \
	suplin1.o \
	supmpi.o \
	suprgf.o \
	supsim.o \
	suptim.o \
	supver.o \
	genutil.o \
	def_para_post.o \
	subnsa.o \
	plotgeneric.o \
	color_table.o \

WHAT = \
	util/suputi.o \
	util/supcol.o \

#----------------------------------------------------------------------------

default: fem
fem: $(EXES) links

links:
	@$(FEMBIN)/symlink $(EXELINKS)

$(EXES): $(OBJ) $(LIBF) $(EXES).o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $^ $(LIBL) -o $@

#----------------------------------------------------------------------------

clean: cleanobj cleanlinks cleanexe cleanvar
	cd util; make clean

cleanvar:
	-rm -f *.tmp *.bak

cleanobj:
	-rm -f *.o
	-rm -f *.mod
	-rm -f *_genmod.f90

cleanlinks:
	@$(FEMBIN)/symlink -remove $(EXELINKS)

cleanexe:
	-rm -f $(EXES)

cleanall: clean
	cd util; make cleanall

#----------------------------------------------------------------------------

depend:
	@$(MAKEDEPEND) -moddir=$(MODDIR) *.f90

#----------------------------------------------------------------------------

%.o:%.f
	${F77} $(FFLAGS) $(INCLUDE) -c $<

%.o:%.f90
	${F77} $(FFLAGS) $(INCLUDE) -c $<

%.o:%.c
	${CC} $(CFLAGS) $(INCLUDE) -c $<

#----------------------------------------------------------------------------

# DO NOT DELETE THIS LINE -- make depend depends on it.

../../../lib/mod/color.mod: col.o
../../../lib/mod/legend.mod: supano.o
../../../lib/mod/mod_bash.mod: supbas.o
../../../lib/mod/mod_hydro_plot.mod: mod_hydro_plot.o
../../../lib/mod/mod_plot.mod: supint.o
../../../lib/mod/plot_fonts.mod: supcust.o
../../../lib/mod/plotutil.mod: plotutil.o
../../../lib/mod/supout.mod: supout.o
../../../lib/mod/timlim.mod: suptim.o
color_table.o: ../../../lib/mod/color.mod \
		../../../lib/mod/para.mod 
def_para_post.o: ../../../lib/mod/para.mod
line_nodes.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/mod_geom.mod 
plotutil.o: ../../../lib/mod/basin.mod ../../../lib/mod/clo.mod \
		../../../lib/mod/nls.mod \
		../../../lib/mod/shyfile.mod 
shyplot.o: ../../../lib/mod/basin.mod ../../../lib/mod/clo.mod \
		../../../lib/mod/elabtime.mod \
		../../../lib/mod/evgeom.mod \
		../../../lib/mod/levels.mod \
		../../../lib/mod/mod_depth.mod \
		../../../lib/mod/mod_geom.mod \
		../../../lib/mod/mod_hydro.mod \
		../../../lib/mod/mod_hydro_plot.mod \
		../../../lib/mod/mod_hydro_print.mod \
		../../../lib/mod/mod_hydro_vel.mod \
		../../../lib/mod/plotutil.mod \
		../../../lib/mod/shyfile.mod \
		../../../lib/mod/shympi.mod \
		../../../lib/mod/shyutil.mod 
subnsa.o: ../../../lib/mod/simul.mod
supano.o: ../../../lib/mod/basin.mod ../../../lib/mod/color.mod \
		../../../lib/mod/mod_bash.mod \
		../../../lib/mod/plot_fonts.mod \
		../../../lib/mod/plotutil.mod \
		../../../lib/mod/simul.mod 
supbas.o: ../../../lib/mod/basin.mod ../../../lib/mod/mod_geom.mod \
		../../../lib/mod/mod_hydro_plot.mod \
		../../../lib/mod/plot_fonts.mod 
supdep.o: ../../../lib/mod/basin.mod ../../../lib/mod/evgeom.mod \
		../../../lib/mod/levels.mod \
		../../../lib/mod/mod_depth.mod \
		../../../lib/mod/mod_hydro.mod \
		../../../lib/mod/mod_hydro_plot.mod 
supint.o: ../../../lib/mod/levels.mod ../../../lib/mod/para.mod
supiso.o: ../../../lib/mod/basin.mod ../../../lib/mod/color.mod \
		../../../lib/mod/mod_hydro_plot.mod 
suplag.o: ../../../lib/mod/basin.mod ../../../lib/mod/color.mod \
		../../../lib/mod/mod_hydro_plot.mod \
		../../../lib/mod/plotutil.mod 
suplin.o: ../../../lib/mod/basin.mod ../../../lib/mod/color.mod \
		../../../lib/mod/levels.mod \
		../../../lib/mod/mod_hydro_print.mod 
suplin1.o: ../../../lib/mod/basin.mod ../../../lib/mod/evgeom.mod \
		../../../lib/mod/levels.mod \
		../../../lib/mod/mod_hydro.mod \
		../../../lib/mod/mod_hydro_plot.mod \
		../../../lib/mod/mod_hydro_print.mod \
		../../../lib/mod/mod_hydro_vel.mod 
supout.o: ../../../lib/mod/basin.mod ../../../lib/mod/levels.mod \
		../../../lib/mod/mod_bash.mod \
		../../../lib/mod/mod_depth.mod \
		../../../lib/mod/mod_hydro.mod \
		../../../lib/mod/mod_hydro_plot.mod \
		../../../lib/mod/mod_hydro_print.mod \
		../../../lib/mod/mod_hydro_vel.mod \
		../../../lib/mod/simul.mod 
supsim.o: ../../../lib/mod/basin.mod ../../../lib/mod/levels.mod \
		../../../lib/mod/mod_color.mod \
		../../../lib/mod/mod_connect.mod \
		../../../lib/mod/mod_depth.mod \
		../../../lib/mod/mod_geom.mod \
		../../../lib/mod/mod_hydro.mod \
		../../../lib/mod/mod_hydro_plot.mod \
		../../../lib/mod/mod_hydro_vel.mod \
		../../../lib/mod/plotutil.mod 
supver.o: ../../../lib/mod/basin.mod ../../../lib/mod/evgeom.mod \
		../../../lib/mod/levels.mod \
		../../../lib/mod/mod_hydro.mod \
		../../../lib/mod/mod_hydro_plot.mod \
		../../../lib/mod/mod_hydro_vel.mod 

