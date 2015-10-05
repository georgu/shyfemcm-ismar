FEMDIR    = ../..
DIRLIB    = $(FEMDIR)/femlib
DIRFEM    = $(FEMDIR)/fem3d
ENKDIR    = ../Evensen_s_code

#---------------------
#F95 = gfortran
F95 = ifort
F77 = $(F95)
FFLAGS = -O -Wall -pedantic $(FGNU_GENERAL)
FFLAGS = -O

LINKER = $(F95)

CC = gcc
CFLAGS = 
#---------------------

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
# Ensemble analysis objects
ENSOBJ = $(ENKDIR)/analysis.o $(ENKDIR)/m_randrot.o $(ENKDIR)/mod_anafunc.o \
	 $(ENKDIR)/m_multa.o $(ENKDIR)/m_mean_preserving_rotation.o

#--------------------
enKF_analysis: enKF_analysis.o subio.o
	cd $(ENKDIR); make clean; make all; cd -
	$(LINKER) $(LFLAGS) -llapack -lblas $@.o \
	subio.o $(OUTOBJS) $(RSTOBJS) $(ENSOBJ) -o $@

#--------------------
make_obs_file: make_obs_file.o subio.o
	$(LINKER) $(LFLAGS) $@.o subio.o -o $@


#--------------------
#--------------------
EXES = enKF_analysis make_obs_file

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
#---------------------



