VPATH = .:RCS:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f .H .h  


#LD = pgf90
#CF90 = pgf90#CF77 = pgf90

LD = ifort
CF90 = ifort
CF77 = ifort

#pgg FFLAGS =  -c -r8 -byteswapio
OPT_FFLAGS = -implicitnone -nozero -warn unused -c -O3 -real_size 64 -assume byterecl -assume nobuffered_io -convert big_endian

PAR = 
DEBUG_FFLAGS = -fpe0 -traceback -warn unused -CB -u -g -c -real_size 64 -assume byterecl -assume nobuffered_io -convert big_endian
FFLAGS =  $(DEBUG_FFLAGS)
FFLAGS =  $(OPT_FFLAGS)

# run ulimit -c unlimited

F77FLG =  

F90FLG =  

LINKFLAGS =


CPPARCH = 
CPPMODEL  =  


CPPFLAGS = -traditional -P $(CPPARCH) $(CPPMODEL)

#LIBDIR = /local/pgi/linux86/6.0/lib
#LIBDIR = /local/intel/compiler90/ia32/fc/lib   -L /h/a684814/lib/linux
#LIBDIR = /h/a684814/lib/linux  -L /usr/lib/gcc-lib/i386-redhat-linux/3.2.3 
LIBDIR = -L /h/a684814/lib/linux  -L /usr/lib/gcc/x86_64-redhat-linux/3.4.3/32 -L/h/a152128/lib/linux 

#LIBS =  -lg2c  -lblas -llapack -lfftw3 -lecl -lhash -lpthread
LIBS =  -lblas -llapack -lfftw3 -lpthread


BINDIR = $(HOME)/bin

CPP = /usr/bin/cpp


# Rules for running cpp and updating files in TMP directory
.H.h:
	rm -f ./TMP/$*.h
	cat MODEL.CPP $*.H | $(CPP) $(CPPFLAGS) > ./TMP/$*.h


.F90.o:
	@rm -f ./TMP/$*.f90
	@cat MODEL.CPP $*.F90 | $(CPP) $(CPPFLAGS) > ./TMP/$*.f90
	cd ./TMP ; $(CF90) $(FFLAGS) $(F90FLG) -o $*.o $*.f90  

.F.o:
	rm -f ./TMP/$*.f
	cat MODEL.CPP $*.F | $(CPP) $(CPPFLAGS) > ./TMP/$*.f
	cd ./TMP ; $(CF77) $(FFLAGS) $(F77FLG) -o $*.o $*.f  



include source.files
include target.mk

INC2 =$(INC1:.H=.h)
FILES =$(F90FILES) $(F77FILES) $(MODULES)
FFILES =$(F90FILES:.F90=.f90) $(F77FILES:.F=.f) $(MODULES:.F90=.f90)
OBJECTS = $(F90FILES:.F90=.o) $(F77FILES:.F=.o) 
OMOD = $(MODULES:.F90=.o) $(MODULES77:.F=.o)



all: $(TARGET) install


$(TARGET): $(INC2) $(OMOD) $(OBJECTS) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o $(TARGET) $(OMOD) $(OBJECTS) $(LIBDIR) $(LIBS) 
	mv TMP/$(TARGET) .

install:
	cp $(TARGET) $(BINDIR)


clean:
	rm -f $(TARGET)
	cd ./TMP ; rm -rf *.o *.mod *.f90 rii_files

new: source depend

source:
	mksource.sh > source.files

depend:
	mkdepend_linux.pl | sort -u > depends.file

tags: $(FILES)
	f90tags.sh

include depends.file
