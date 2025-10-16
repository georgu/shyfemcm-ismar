
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2023  The SHYFEM team
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

DIR	= weutro

#--------------------------------------------- just for default
FEMDIR = ../../../..
FEMBIN = $(FEMDIR)/bin
LIBDIR = $(FEMDIR)/lib
F77    = g77
FFLAGS =  -O -Wall -g
LINKER = $(F77)
#---------------------------------------------

include $(FEMDIR)/Rules.make

#--------------------------------------------- macros

EXES	=	

ifeq ($(ECOLOGICAL),EUTRO)
  OBJ	     =	bio3d.o \
		bio3d_util.o \
		weutro.o \
		weutro_light.o \
		weutro_sedim.o \
		weutro_seed.o \
		weutro_shell.o \
		weutro_ulva.o \
		weutro_wdia.o \
		weutro_param.o \

endif

WEUTROOBJ := $(foreach item,$(OBJ),$(PWD)/$(item))

ALL	=	$(OBJ)

SPECIAL =       Makefile TODO README
SAVEFILES =     *.f90 $(SPECIAL)

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#--------------------------------------------- EXES

default: fem
fem:	$(OBJ) nothing

nothing:
	@true

objs:
	@echo $(WEUTROOBJ)

#--------------------------------------------- clean

cleanall: clean cleanzip cleangrd

clean: cleanobj cleanexe cleanvar cleanlinks

cleanobj:
	-rm -f *.o
	-rm -f *.stb
	-rm -f *.mod
	-rm -f *_genmod.f90

cleanlinks:
	@$(FEMBIN)/symlink -remove $(FEMBIN) $(EXES)

cleanexe:
	-rm -f $(EXES)
	-rm -f *.exe

cleanvar:
	-rm -f a.out core
	-rm -f *.bak *.tmp
	-rm -f ggg hhh
	-rm -f mon.out nohup.out
	-rm -f fort.*
	-rm -f errout.dat
	-rm -f gmon.out
	-rm -f out.ps plot.ps

cleanzip:
	-rm -f $(DIR).zip
	-rm -f save.zip

cleangrd:
	-rm -f *.grd

#--------------------------------------------- various

depend:
	@-$(MAKEDEPEND) -moddir=$(MODDIR) *.f90

list:
	pwd; ls

save:   zip
zip:	cleanall
	-zip $(DIR) $(SAVEFILES)

strip:
	-strip $(EXES)

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

../../../../lib/mod/eutro.mod: bio3d.o
bio3d.o: ../../../../lib/mod/basin.mod \
		../../../../lib/mod/levels.mod \
		../../../../lib/mod/mkonst.mod \
		../../../../lib/mod/mod_diff_visc_fric.mod 
bio3d_util.o: ../../../../lib/mod/basin.mod \
		../../../../lib/mod/levels.mod donata.h 
weutro.o: donata.h weutro.h
weutro_sedim.o: donata.h weutro.h
weutro_seed.o: ../../../../lib/mod/basin.mod \
		../../../../lib/mod/levels.mod 
weutro_shell.o: donata.h
weutro_ulva.o: donata.h
weutro_wdia.o: ../../../../lib/mod/eutro.mod weutro.h

