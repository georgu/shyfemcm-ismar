
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017-2019  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

! revision log :
!
! 19.05.2020	ccf	started from scratch
! 24.05.2020	ggu	debug option added
! 28.05.2020	ggu	some more informational messages
! 15.07.2020	ggu	bug fix for counting elements
! 22.04.2021	ggu	resolve bound check error (not yet finished)
! 11.04.2022	ggu	prepare for online partitioning
! 12.04.2022	ggu	bug fix: iarnv and iarv were not saved
! 15.10.2024	ggu	set array tpwgts (which is real, not dble) (IFX bug)
! 21.11.2024	ggu	call check_partition() with bdebug flag
! 23.11.2024	ggu	code cleaned
! 29.11.2024	ggu	option to force contigous areas
!
!****************************************************************

        subroutine do_partition(nkn,nel,nen3v,nparts,npart,epart)

! partitions grd file using the METIS library
! requires metis-5.1.0.tar.gz

	use shympi

	implicit none

	integer nkn,nel
	integer nen3v(3,nel)
	integer nparts
	integer npart(nkn)
	integer epart(nel)

	integer		      :: ie,ii,l
	integer, allocatable  :: eptr(:) 	!index for eind
	integer, allocatable  :: eind(:) 	!nodelist in elements
	integer, allocatable  :: vwgt(:) 	!weights of vertices
        integer, allocatable  :: vsize(:)	!size of the nodes
        real, allocatable     :: tpwgts(:)	!desired partition weight 

        integer               :: objval		!edge-cut or total comm vol
	integer               :: options(40)	!metis options

	logical bcontig
	integer iu,nfill,i

	bcontig = .false.	!want contigous areas?

!-----------------------------------------------------------------
! initialiaze arrays
!-----------------------------------------------------------------

	npart = 0
	epart = 0

        allocate(eptr(nel+1))
        allocate(eind(nel*3))

        allocate(vwgt(nkn))
        allocate(vsize(nkn))
	vwgt = 1
	vsize = 1

        allocate(tpwgts(nparts))
	tpwgts = 1./dble(nparts)

!-----------------------------------------------------------------
! set up METIS eptr and eind arrays structures
!-----------------------------------------------------------------

        eptr(1)=1
        do ie = 1,nel
          do ii = 1,3
            l = eptr(ie) + ii - 1
            eind(l) = nen3v(ii,ie)
          end do
          eptr(ie+1) = eptr(ie) + 3
        enddo 

!-----------------------------------------------------------------
! Set METIS options
! For the full list and order of options see metis-5.1.0/include/metis.h
! and the documentation in manual/manual.pdf
!-----------------------------------------------------------------

        call METIS_SetDefaultOptions(options)

        !options(1) = 1		!PTYPE (0=rb,1=kway)
        !options(2) = 1		!OBJTYPE (0=cut,1=vol)
        !options(12) = 1	!CONTIG (0=defoult,1=force contiguous)
        options(18) = 1		!NUMBERING (0 C-style, 1 Fortran-style)

	if( bcontig ) options(12) = 1

!-----------------------------------------------------------------
! Call METIS for patitioning on nodes
!-----------------------------------------------------------------

	write(6,*) 'partitioning with METIS...'
        call METIS_PartMeshNodal(nel, nkn, eptr, eind, vwgt, vsize,  &
     &       nparts, tpwgts, options, objval, epart, npart)

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

        end

!*******************************************************************

        subroutine check_partition(npart,epart,bdebug,ierr1,ierr2)

! check partition

        use basin
        use mod_geom

        implicit none

        integer npart(nkn)
        integer epart(nel)
	logical bdebug
        integer ierr1,ierr2

        integer :: nsave(nkn)
        integer :: esave(nel)

        write(6,*) 'checking connectivity and connections...'

	ierr1 = 0
	ierr2 = 0

        call set_geom

        call link_set_stop(.false.)     !not stopping after error
        call link_set_write(.false.)    !do not write error

	nsave = iarnv
	esave = iarv

        iarnv = npart
        iarv = epart
        call check_connectivity(bdebug,ierr1)
        call check_connections(bdebug,ierr2)
        npart = iarnv
        epart = iarv

	iarnv = nsave
	iarv  = esave

        call link_set_stop(.true.)

	call mod_geom_release

        end

!****************************************************************

