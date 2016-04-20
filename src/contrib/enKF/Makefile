
#--------------------
# Include shyfem rules
include ../Rules.make

#F77=ifort -r8
#F95=ifort -r8
#LINKER=ifort -r8

# Reset some env vars
FEMDIR    = ..
DIRLIB    = $(FEMDIR)/femlib
DIRFEM    = $(FEMDIR)/fem3d
ENKDIR    = Evensen_s_code

#--------------------
# ht objects
HTOBJS = $(DIRFEM)/subrst.o $(DIRFEM)/subbas.o $(DIRFEM)/subdts.o $(DIRFEM)/mod_bnd.o \
	$(DIRFEM)/subreg.o $(DIRFEM)/subpar3.o $(DIRFEM)/suboutput.o $(DIRFEM)/subfil.o \
	$(DIRFEM)/mod_hydro.o $(DIRFEM)/debug.o $(DIRFEM)/subsss.o $(DIRFEM)/subbit10.o \
 	$(DIRFEM)/mod_geom_dynamic.o $(DIRFEM)/mod_ts.o $(DIRFEM)/newlevels.o \
	$(DIRFEM)/mod_conz.o $(DIRFEM)/mod_hydro_vel.o $(DIRFEM)/ecological_dummy.o \
	$(DIRFEM)/subdef.o $(DIRFEM)/subnev.o $(DIRFEM)/mod_geom.o $(DIRFEM)/subscn.o

#--------------------
# Ensemble analysis modules
ENSOBJ = m_randrot.o m_multa.o mod_anafunc.o m_mean_preserving_rotation.o
mod_anafunc.o: m_mean_preserving_rotation.o
analysis_ev.o: $(ENSOBJ)

#--------------------
enKF_analysis.o: mod_enKF.o
enKF_analysis: enKF_analysis.o mod_enKF.o analysis_ev.o $(ENSOBJ) $(HTOBJS)
	$(LINKER) $(LFLAGS) -llapack -lblas $@.o mod_enKF.o \
	analysis_ev.o $(ENSOBJ) $(HTOBJS) -o $@

#--------------------
#--------------------
EXES = enKF_analysis 

links:
	ln -s ../fem3d/basin.h
	ln -s ../fem3d/conz.h
	ln -s ../fem3d/geom_dynamic.h
	ln -s ../fem3d/hydro.h
	ln -s ../fem3d/hydro_vel.h
	ln -s ../fem3d/nlevel.h
	ln -s ../fem3d/param.h
	ln -s ../fem3d/ts.h

all:
	make $(EXES)

clean:
	rm -fr $(EXES) *.o *.mod tags

cleanall:
	rm -fr $(EXES) *.o *.mod tags *.dat *.uf fort.*; cd $(ENKDIR); rm -fr *.o *.mod 


#---------------------
# Compilation rules
.SUFFIXES: .F90 .f90 .for
.f.o:
	$(F77) -c $(FFLAGS) $<
.for.obj:
	$(FLC) $< , $(<,B).OBJ $(FLFLAGS)
.f90.o:
	$(F95) -c $(FFLAGS) $< 
.F90.o:
	$(F95) -c $(FFLAGS) $<
.c.o:
	$(CC) -c $(CFLAGS) $<

%.o : %.mod
