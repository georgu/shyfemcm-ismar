Author: Geir.Evensen@nersc.no

The files in this directory illustrates a Fortran90 programming model which provides
full control of dependencies, interface consistency and automatic makefile generation.

It assumes that every subroutine or function is placed in a seperate file with a standardized
naming convention.    It allows for four types of source files,
   F77 code,                     (subroutine sub4, stored in sub4.F,     fixed form)
   F77 code stored in module,    (subroutine sub5, stored in m_sub5.F,   fixed form)
   F90 code,                     (subroutine sub3, stored in sub3.F90,   free form)
   F90 code stored in module,    (subroutine sub1, stored in m_sub1.F90, free form)
It is encourages that all subroutines and functions are store in modules, e.g.,

module m_sub1
contains
subroutine sub1()
   use m_sub2
   implicit none
   print *,'Sub1'
   call sub2
end subroutine sub1
end module m_sub1

This requires the subroutines to be compiled in the correct order (sub2 must be compiled
before sub1 because of the use m_sub2 statement).   The "use m_sub2" statement also
provides "sub1" with the interface of "sub2" and any mismatch in the call and header
information will be reported at compile time.

To simplify the organization of the code we also allow for pure module definitions to
be contained in files mod_xxxx.F90,  see e.g. the mod_dimensions.F90 and mod_states.F90
files.

The file MODEL.CPP contains define statements which are used during the compilation.
See the example in m_sub1.F90.

A generic makefile is used, which is totally model independent.
The make new command will generate two files, source.files and depends.file
which contains all the source files in the current directory and the dependencies
between the files.   Thus if a new subroutine is added, just type "make new"
to update the makefile.   
 - Note that only one main program should reside in the current directory.
 - The present makefile is configured for IBM Regatta systems.
 - The makefile stores and compiles all code in the TMP directory.
