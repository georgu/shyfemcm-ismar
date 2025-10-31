
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010,2019  Georg Umgiesser
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
! 23.03.2010	ggu	changed v6.1.1
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 21.10.2025	ggu	new routine read_timeseries()
! 30.10.2025	ggu	invert dimensions in values

!***********************************************************************

	subroutine read_ts(file,ndim,nvar,n,itime,array)

! reads time series file

	implicit none

	character*(*) file		!file name
	integer ndim			!dimension of arrays
	integer nvar			!number of expected variables in file
	integer n			!number of records read (return)
	integer itime(ndim)		!time column
	real array(ndim,nvar)		!value column(s)

	integer naux
	parameter(naux=100)

	integer iunit,i,it
	real aux(naux)
	integer ifileo

	if( nvar .gt. naux ) stop 'error stop read_ts: naux'

	iunit = 55
        iunit = ifileo(iunit,file,'form','old')
        if( iunit .le. 0 ) stop 'error stop read_ts: no such file'

	n = 0
    1	continue
	  read(iunit,*,end=2) it,(aux(i),i=1,nvar)
	  n = n + 1
	  if( n .gt. ndim ) stop 'error stop read_ts: ndim'
	  itime(n) = it
	  do i=1,nvar
	    array(n,i) = aux(i)
	  end do
	goto 1
    2	continue

	close(iunit)

	end

!***********************************************************************

        subroutine read_timeseries(iu,nmax,nval,times,values)

        implicit none

        integer iu
        integer nmax,nval
        double precision times(nmax)
        real values(nmax,nval)

	logical bstring
        integer i,n,nv
        integer ios,ioff,ierr
	double precision atime
	character*20 aline
	character*256 line
	double precision d(nval)

	integer iscand,istos,istod,iston,istot

	i = 0
	n = 0
        do
	  i = i + 1
          read(iu,'(a)',iostat=ios) line
          if( ios < 0 ) exit
          if( ios > 0 ) goto 99
          n = iscand(line,d,0)
          if( n /= 0 ) exit
        end do
	rewind(iu)
	if( n == 0 ) then	!nothing in file
	  nmax = 0
	  nval = 0
	  return
	else if ( n > 0 ) then	!relative time (number)
	  bstring = .false.
	else			!absolute time (string)
	  bstring = .true.
	end if

        i = 0
	nv = 0
        do
	  ioff = 1
          read(iu,'(a)',iostat=ios) line
          if( ios < 0 ) exit
          if( ios > 0 ) goto 99
          if( bstring ) then
	    n = istot(line,aline,ioff)
	    if( n == 0 ) cycle
	    if( n < 0 ) goto 95
	    call dts_string2time(aline,atime,ierr)
	    if( ierr /= 0 ) goto 96
          else
            n = istod(line,atime,ioff)
	    if( n == 0 ) cycle
	    if( n < 0 ) goto 95
          end if
          n = iscand(line(ioff:),d,nval)
          if( n == 0 ) cycle                    !empty line
          if( nv == 0 ) nv = n
          if( nv /= n ) goto 98
	  i = i + 1
          if( nmax == 0 ) cycle
	  if( i > nmax ) goto 97
          times(i) = atime
          values(i,1:nval) = d(1:nval)
        end do

	nval = nv
        nmax = i
	rewind(iu)

        return
   95   continue
        write(6,*) 'error reading time: ',trim(line)
        stop 'error stop read_timeseries: error reading time'
   96   continue
        write(6,*) 'error parsing time string: ',trim(aline)
        stop 'error stop read_timeseries: error parsing time'
   97   continue
        write(6,*) 'i,nmax: ',i,nmax
        stop 'error stop read_timeseries: i>nmax'
   98   continue
        write(6,*) 'n,nv: ',n,nv
        stop 'error stop read_timeseries: n/=nv'
   99   continue
	write(6,*) 'read error in line ',i
        stop 'error stop read_timeseries: read error'
        end

!*******************************************************************
!*******************************************************************
!*******************************************************************
! testing routines
!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine test_read_ts

	implicit none

	integer iu,nmax,nval,i
	double precision, allocatable :: times(:)
	real, allocatable :: values(:,:)

	write(77,*) 0, 50, 110
	write(77,*) 3600, 50, 120
	write(77,*) 
	write(77,*) 7200, 70, 130
	write(77,*) 8000
	write(77,*) 10800, 80, 140
	close(77)

	call test_read_file(77,4,2)

	write(88,*) '2013-07-01::01:30:00', 200
	write(88,*) '2013-07-01::02:30:00', 300
	write(88,*) '2013-07-01::03:30:00', 400
	write(88,*) '2013-07-01::04:30:00', 500
	close(88)

	call test_read_file(88,4,1)

	end

!*******************************************************************

	subroutine test_read_file(iu,nmax0,nval0)

	implicit none

	integer iu,nmax0,nval0

	integer nmax,nval,i
	double precision, allocatable :: times(:)
	real, allocatable :: values(:,:)

	nmax=0
	nval=0
        call read_timeseries(iu,nmax,nval,times,values)
	write(6,*) 'nmax0,nval0: ',nmax0,nval0
	write(6,*) 'nmax,nval: ',nmax,nval
	if( nmax /= nmax0 .or. nval /= nval0 ) stop 'nmax/nval mismatch'
	allocate(times(nmax),values(nmax,nval))
        call read_timeseries(iu,nmax,nval,times,values)
	do i=1,nmax
	  write(6,*) times(i),values(i,:)
	end do
	
	end

!*******************************************************************

!	program main_test_read_ts
!	call test_read_ts
!	end

!*******************************************************************
!
! gfortran -cpp rd_ts.f90 ../utils/generic/dts.f90 ../utils/generic/file.f90 ../utils/generic/convert.f90
!
!*******************************************************************
