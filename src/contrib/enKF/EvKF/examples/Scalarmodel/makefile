VPATH = .:RCS:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f .H .h  


#Compiler Options -- serial hands-off version
#LD = xlf90_r -qspillsize=2152 
#CF90 = xlf90_r  -qsuffix=f=f90 -qspillsize=2152 -qfree=f90
#CF77 = xlf90_r -qfixed -qspillsize=2152
#




#Compiler -- F90
#CF90 = /usr/bin/xlf95_r 
CF90 = /usr/bin/gfortran 

#Compiler -- F77
#CF77 = /usr/bin/xlf90_r 
CF77 = /usr/bin/gfortran 

#Linker 
#LD = /usr/bin/xlf95_r 
LD = /usr/bin/gfortran 

#Paralellization opts
#PARO =
PARO = -qsmp=omp
PARO_OE = 

#Size defaults
SIZEO = -b64 -q64 -qrealsize=8

#Arch opts
ARCHO= -qarch=auto -qtune=auto -qcache=auto

#Optimalization opts
OPTO= -O3 -qstrict
OPTO= -O2 -qstrict

#Inlining opts
#INLO= -Q -qipa
INLO=

# Diverse options
DIVO= -qmaxmem=-1 -qnosave




# Flags for debugging -- slow, gives location of SIGTRAP
# -g        -- turns on debugging ...
# -qfloat   -- options for floating point operations
# -qflttrap -- Floating point operation traps
# -qextchk  -- Checks if subroutine args are consistent + common block checks
# -C        -- Checks array bounds 
# qinitauto=FF initializes all vars to NaN -- So that uninitialized vars can be detected
#
DEBUG_FLAGS = -C -qflttrap=underflow:overflow:zerodivide:invalid:enable  \
              -qextchk -qinitauto=FF -qfloat=nofold:nohssngl:nans  -g


# version without qfloat shoul allow for stmt_func
#         without qextchk which allows for inconsistent interfaces 
#         qextchk notes the different declarations of a in analysis(A,...)
#DEBUG_FLAGS = -C -qflttrap=overflow:zerodivide:invalid:enable -qinitauto=FF -g
#DEBUG_FLAGS = -qflttrap=underflow:overflow:zerodivide:invalid:enable 
DEBUG_FLAGS = 


F77FLG = -qextname -qfixed                   
F90FLG = -qextname -qsuffix=f=f90 -qfree=f90 


FFLAGS    = $(SIZEO) $(OPTO) $(ARCHO) $(PARO)    $(INLO) $(DIVO) $(DEBUG_FLAGS) 
FFLAGS_OE = $(SIZEO) $(OPTO) $(ARCHO) $(PARO_OE) $(INLO) $(DIVO) $(DEBUG_FLAGS) 
LINKFLAGS = $(SIZEO) $(OPTO) $(PARO) $(INLO) $(DIVO) $(DEBUG_FLAGS)  -bloadmap:load.out





CPPARCH = -DPWR4 -DIBM
CPPMODEL  =  -DTOPAZ1
CPPFLAGS =  -P $(CPPARCH) $(CPPMODEL)
LIBS = -lessl /home/parallab/nersc/knutali/lib/liblapack64.a


#CPP = /usr/lib/cpp
CPP = /usr/bin/g++


# Rules for running cpp and updating files in TMP directory
.H.h:
	rm -f ./TMP/$*.h
	cat MODEL.CPP $*.H | $(CPP) $(CPPFLAGS) | st.sh > ./TMP/$*.h


.F90.o:
	rm -f ./TMP/$*.f90
	cat MODEL.CPP $*.F90 | $(CPP) $(CPPFLAGS)  | st.sh > ./TMP/$*.f90
	cd ./TMP ; $(CF90) -c $(FFLAGS) $(F90FLG) -o $*.o $*.f90  

.F.o:
	rm -f ./TMP/$*.f
	cat MODEL.CPP $*.F | $(CPP) $(CPPFLAGS) | st.sh > ./TMP/$*.f
	cd ./TMP ; $(CF77) -c $(FFLAGS) $(F77FLG) -o $*.o $*.f  



TARGET =  scalar

include source.files

INC2 =$(INC1:.H=.h)
FILES =$(F90FILES) $(F77FILES) $(MODULES)
FFILES =$(F90FILES:.F90=.f90) $(F77FILES:.F=.f) $(MODULES:.F90=.f90)
OBJECTS = $(F90FILES:.F90=.o) $(F77FILES:.F=.o) 
OMOD = $(MODULES:.F90=.o) $(MODULES77:.F=.o)


all: $(TARGET) 


$(TARGET): $(INC2) $(OMOD) $(OBJECTS)  $(OMP_EXCEPTION_OBJ77) $(OMP_EXCEPTION_OBJ90)
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET) $(OMOD) $(OBJECTS) $(LIBS)





clean:
	cd ./TMP ; rm *.f  *.o *.f90 *.h *.mod


new: source depend

source:
	mksource.sh > source.files

depend:
	mkdepend.pl | sort -u > depends.file

include depends.file

