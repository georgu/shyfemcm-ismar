
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
PARMOD = $(DIRFEM)/subpar3.o
OUTOBJS = $(PARMOD) $(DIRFEM)/subnls.o $(DIRFEM)/subsss.o $(DIRFEM)/subscn.o \
          $(DIRFEM)/subfil.o $(DIRFEM)/subdef.o $(DIRFEM)/subnsu.o $(DIRFEM)/subssv.o \
          $(DIRFEM)/subbas.o $(DIRFEM)/subsys.o $(DIRFEM)/subdim.o \
          $(DIRFEM)/subapn.o $(DIRFEM)/subnsa.o $(DIRFEM)/subdum.o \
          $(DIRFEM)/subbnd.o $(DIRFEM)/subdts.o $(DIRFEM)/subdep.o \
          $(DIRFEM)/debug.o $(DIRFEM)/subreg.o $(DIRFEM)/subnev.o \
          $(DIRFEM)/subfemfile.o $(DIRFEM)/subtsfile.o \
          $(DIRFEM)/subvintp.o $(DIRFEM)/subver.o $(DIRFEM)/sigmautil.o $(DIRFEM)/suboutput.o \
          $(DIRFEM)/subclo.o
RSTOBJS = $(DIRFEM)/subrst.o $(DIRFEM)/ecological_dummy.o $(DIRFEM)/newlevels.o

#--------------------
# Ensemble analysis modules
ENSOBJ = m_randrot.o m_multa.o mod_anafunc.o m_mean_preserving_rotation.o
mod_anafunc.o: m_mean_preserving_rotation.o
analysis_ev.o: $(ENSOBJ)

#--------------------
enKF_analysis: enKF_analysis.o analysis_ev.o $(ENSOBJ) subio.o
	$(LINKER) $(LFLAGS) -llapack -lblas $@.o \
	analysis_ev.o $(ENSOBJ) subio.o $(OUTOBJS) $(RSTOBJS) -o $@

#--------------------
make_obs_file: make_obs_file.o subio.o
	$(LINKER) $(LFLAGS) $@.o subio.o -o $@


#--------------------
#--------------------
EXES = enKF_analysis make_obs_file

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
	rm -fr $(EXES) *.o *.mod

cleanall:
	rm -fr $(EXES) *.o *.mod; cd $(ENKDIR); rm -fr *.o *.mod


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
