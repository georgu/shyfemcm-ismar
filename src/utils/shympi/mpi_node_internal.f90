
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015-2019  Georg Umgiesser
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

! mpi routines
!
! contents :
!
! revision log :
!
! 24.11.2015	ggu	project started
! 22.06.2016	ggu	added sum option to shympi_reduce
! 07.12.2017	ggu	changed VERS_7_5_40
! 24.01.2018	ggu	changed VERS_7_5_41
! 10.04.2018	ggu	code to exchange arrays
! 19.04.2018	ggu	changed VERS_7_5_45
! 26.04.2018	ggu	changed VERS_7_5_46
! 11.05.2018	ggu	bug fix in exchange arrays for zeta levels
! 16.02.2019	ggu	changed VERS_7_5_60
! 22.04.2021	ggu	bug fix in shympi_allgather_*()
! 02.04.2022	ggu	new routines shympi_rectify_internal_*()
! 03.04.2022	ggu	new routine shympi_bcast_d_internal()
! 06.04.2022	ggu	new routines for handling double precision
! 10.04.2022	ggu	bug fix in shympi_bcast_d_internal() - val was real
! 01.06.2022	ggu	new routine shympi_gather_d_internal()
! 09.10.2022	ggu	rectify 3d arrays with nlv+1 (nextra)
! 27.03.2023	ggu	new routines shympi_receive_internal_*()
! 13.04.2023	ggu	introduced bnode, belem (distinguish calls to node/elem)
! 27.03.2024	ggu	introduced aux array to make dims equal in gather
! 05.04.2024	ggu	changes in shympi_exchange_internal_r()
! 06.09.2024    lrp     nuopc-compliant
! 21.11.2024    ggu     change in shympi_abort_internal() -> hangs for ever
! 30.11.2024    ggu     new routine shympi_icomment()
! 04.12.2024    ggu     new routines shympi_gather_internal_r/i()
! 05.12.2024    ggu     renamed module from shympi_aux to shympi_internal
! 07.12.2024    ggu     big changes: check if arrays have same length for gather
!
!******************************************************************

!==================================================================
        module shympi_internal
!==================================================================

! this module is only used inside this file

        use mpi

        implicit none

        !include 'mpif.h'

	logical, parameter :: bcomment = .false.	!write comment
	integer, parameter :: iuc = 999          	!unit for comment
	character*4, parameter :: indent = '    '	!indent for comment

	logical, parameter :: bpdebug = .false.		!write debug messages
	integer, save :: my_p_unit_base = 800	!base unit for debug output

!	---------------------------------
!	next variables are set internally
!	---------------------------------

	integer, save :: n_p_threads = 1	!total number of threads
	integer, save :: my_p_id = 0		!id of this thread
	integer, save :: status_p_size = 0	!size of status array
	logical, save :: bpmaster = .true.	!is this the master?
	logical, save :: bpmpi = .false.	!use mpi? (threads > 1)
	logical, save :: bpoperall = .true.	!do all gather/reduce?
	integer, save :: my_p_unit = 0		!unit for debug output

        INTERFACE shympi_icomment
        MODULE PROCEDURE   shympi_icomment_only &
     &                   , shympi_icomment_logical
        END INTERFACE

!==================================================================
        contains
!==================================================================

	subroutine shympi_icomment_only(text)

	character*(*), intent(in) :: text

	if( bcomment ) then
	  if( bpmaster ) then
	    write(iuc,*) indent,trim(text)
	  end if
	end if

	end subroutine shympi_icomment_only

!**********************************

	subroutine shympi_icomment_logical(text,blogical)

	character*(*), intent(in) :: text
	logical blogical

	if( bcomment ) then
	  if( bpmaster ) then
	    write(iuc,*) indent,trim(text),' ',blogical
	  end if
	end if

	end subroutine shympi_icomment_logical

!******************************************************************

	subroutine shympi_internal_reduce_what(what,mpi_what)

! returns what to reduce

	character*(*), intent(in) :: what
	integer, intent(out) :: mpi_what

         if( what == 'min' ) then
	  mpi_what = MPI_MIN
         else if( what == 'max' ) then
	  mpi_what = MPI_MAX
         else if( what == 'sum' ) then
	  mpi_what = MPI_SUM
	else
	  write(6,*) 'what = ',trim(what)
	  stop 'error stop shympi_internal_reduce_what: unknown what'
	end if

	end subroutine shympi_internal_reduce_what

!==================================================================
        end module shympi_internal
!==================================================================

!******************************************************************

	subroutine shympi_operate_all_internal(boperall)

	use shympi_internal

	implicit none

	logical boperall

	bpoperall = boperall

	end subroutine shympi_operate_all_internal

!******************************************************************

	subroutine shympi_error(routine,what,ierr)

	use shympi_internal

	implicit none

	character*(*) routine,what
	integer ierr

	integer eslen,iserr
	character*80 estring

	if( ierr /= 0 ) then
	  eslen = 0
	  iserr = 0
	  estring = ' '
	  !call MPI_ERROR_STRING(ierr,estring,eslen,iserr)
	  if( eslen > 80 .or. iserr /= 0 ) then
	    estring = 'no description for error'
	  end if
	  write(6,*) 'error in routine ',trim(routine)
	  write(6,*) 'ierr = ',ierr,' while doing ',trim(what)
	  write(6,*) 'error: ',trim(estring)
	  stop 'error stop shympi_error'
	end if

	end subroutine shympi_error

!******************************************************************

	subroutine shympi_init_internal(my_id,n_threads,binit)

	use shympi_internal

	implicit none

	integer my_id,n_threads
	logical :: binit

	integer ierr,iberr
	integer required,provided

	required = MPI_THREAD_MULTIPLE
	required = MPI_THREAD_SERIALIZED

	ierr = 0
        if( binit ) call MPI_INIT_THREAD( required, provided, ierr )
	!write(6,*) 'thread safety: ',required, provided
	!write(6,*) 'initializing MPI: ',ierr

	call MPI_BARRIER( MPI_COMM_WORLD, iberr)
	call shympi_error('shympi_init_internal','init',ierr)
        call MPI_COMM_RANK( MPI_COMM_WORLD, my_id, ierr )
	call shympi_error('shympi_init_internal','rank',ierr)
        call MPI_COMM_SIZE( MPI_COMM_WORLD, n_threads, ierr )
	call shympi_error('shympi_init_internal','size',ierr)
	call MPI_BARRIER( MPI_COMM_WORLD, ierr)
	call shympi_error('shympi_init_internal','barrier',ierr)

	n_p_threads = n_threads
	my_p_id = my_id
	bpmaster = ( my_p_id == 0 )
	bpmpi = ( n_p_threads > 1 )

	my_p_unit = my_p_unit_base + my_p_id

	end subroutine shympi_init_internal

!******************************************************************

	subroutine shympi_barrier_internal

	use shympi_internal

	implicit none

	integer ierr

	call MPI_BARRIER( MPI_COMM_WORLD, ierr)
	call shympi_icomment('shympi_barrier_internal')

	end subroutine shympi_barrier_internal

!******************************************************************

        subroutine shympi_abort_internal(ierr_code)

	use shympi_internal

        implicit none

        integer ierr,ierr_code

	!call MPI_BARRIER( MPI_COMM_WORLD, ierr)
	!call MPI_FINALIZE(ierr_code)		!this will hang for ever
	call MPI_ABORT(MPI_COMM_WORLD,ierr_code,ierr)

        end subroutine shympi_abort_internal

!******************************************************************

        subroutine shympi_finalize_internal

	use shympi_internal

        implicit none

        integer ierr

	call MPI_BARRIER( MPI_COMM_WORLD, ierr)
	call MPI_FINALIZE(ierr)

        end subroutine shympi_finalize_internal

!******************************************************************

        function shympi_wtime_internal()

	use shympi_internal

        implicit none

	double precision shympi_wtime_internal

	shympi_wtime_internal = MPI_WTIME()
	call shympi_icomment('shympi_wtime_internal')

        end function shympi_wtime_internal

!******************************************************************

        subroutine shympi_get_status_size_internal(size)

	use shympi_internal

        implicit none

        integer size

	size = MPI_STATUS_SIZE
	status_p_size = MPI_STATUS_SIZE
	call shympi_icomment('shympi_get_status_size_internal')

        end subroutine shympi_get_status_size_internal

!******************************************************************

	subroutine shympi_syncronize_internal

	use shympi_internal

	implicit none

	integer ierr

	flush(6)
	call MPI_BARRIER( MPI_COMM_WORLD, ierr)
	call shympi_icomment('shympi_syncronize_internal')

	end subroutine shympi_syncronize_internal

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_syncronize_initial

	use shympi_internal

        implicit none

	integer my_id,nt
	integer root,ierr,i
	integer count
	integer local(1)
	integer, allocatable :: buf(:)

        call MPI_COMM_RANK( MPI_COMM_WORLD, my_id, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, nt, ierr )

	allocate(buf(nt))

	root = my_id
	count = 1
	root = 0

	if( my_id == root ) then
	  do i=1,nt
	    buf(i) = i
	  end do
	end if

	call MPI_SCATTER (buf,count,MPI_INTEGER &
     &			,local,count,MPI_INTEGER &
     &			,root,MPI_COMM_WORLD,ierr)

	local = local * 2
	root = 0

	call MPI_GATHER (local,count,MPI_INTEGER &
     &			,buf,count,MPI_INTEGER &
     &			,root,MPI_COMM_WORLD,ierr)

	call MPI_BARRIER( MPI_COMM_WORLD, ierr)

	deallocate(buf)

	call shympi_icomment('shympi_syncronize_initial')

	end subroutine shympi_syncronize_initial

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_receive_internal_i(id_from,id_to &
     &						,n,val_in,val_out)

	use shympi_internal
	!use shympi

	implicit none

	integer id_from,id_to
	integer n
	integer val_in(n)
	integer val_out(n)

	integer tag,ir,id
	integer ierr
	integer nb
	integer status(status_p_size,2*n_p_threads)
	integer request(2*n_p_threads)

        tag=151
	ir = 0

!ccgguccc!$OMP CRITICAL

	if( my_p_id == id_to ) then
	  ir = ir + 1
	  id = id_from
          call MPI_Irecv(val_out,n,MPI_INTEGER,id &
     &	          ,tag,MPI_COMM_WORLD,request(ir),ierr)
	end if

	if( my_p_id == id_from ) then
	  ir = ir + 1
	  id = id_to
          call MPI_Isend(val_in,n,MPI_INTEGER,id &
     &	          ,tag,MPI_COMM_WORLD,request(ir),ierr)
	end if

        call MPI_WaitAll(ir,request,status,ierr)

!ccgguccc!$OMP END CRITICAL

	call shympi_icomment('shympi_receive_internal_i')

	end subroutine shympi_receive_internal_i

!******************************************************************

	subroutine shympi_receive_internal_r(id_from,id_to &
     &						,n,val_in,val_out)

	use shympi_internal
	!use shympi

	implicit none

	integer id_from,id_to
	integer n
	real val_in(n)
	real val_out(n)

	integer tag,ir,id,iu
	integer ierr
	integer nb
	integer status(status_p_size,2*n_p_threads)
	integer request(2*n_p_threads)

        tag=152
	ir = 0
	ierr = 0

!ccgguccc!$OMP CRITICAL

	if( my_p_id == id_to ) then
	  ir = ir + 1
	  id = id_from
          call MPI_Irecv(val_out,n,MPI_REAL,id &
     &	          ,tag,MPI_COMM_WORLD,request(ir),ierr)
	end if
	if( ierr /= 0 ) write(6,*) 'internal error 1: ',ierr

	if( my_p_id == id_from ) then
	  ir = ir + 1
	  id = id_to
          call MPI_Isend(val_in,n,MPI_REAL,id &
     &	          ,tag,MPI_COMM_WORLD,request(ir),ierr)
	end if
	if( ierr /= 0 ) write(6,*) 'internal error 2: ',ierr

	if( ir > 0 ) then
          call MPI_WaitAll(ir,request,status,ierr)
	  if( ierr /= 0 ) write(6,*) 'internal error 3: ',ierr
	end if

!ccgguccc!$OMP END CRITICAL

	if( ierr > 0 ) then
	  stop 'error stop shympi_receive_internal: exchange error'
	end if

	call shympi_icomment('shympi_receive_internal_r')

	end subroutine shympi_receive_internal_r

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_exchange_internal_i(belem,n0,nlvddi,n,il &
     &						,g_in,g_out,val)

	use shympi_internal
	use shympi

	implicit none

	logical belem
	integer n0,nlvddi,n
	integer il(n)
	integer g_in(n_ghost_max,n_ghost_areas)
	integer g_out(n_ghost_max,n_ghost_areas)
	integer val(n0:nlvddi,n)

	integer tag,ir,ia,id
	integer i,k,nc,ierr
	integer nb
	integer iout,iin
	integer nbs(2,n_ghost_areas)
	integer status(status_p_size,2*n_p_threads)
	integer request(2*n_p_threads)
	integer, allocatable :: buffer_in(:,:)
	integer, allocatable :: buffer_out(:,:)

        tag=121
	ir = 0

	iout = 2
	iin = 3
	if( belem ) then
	  iout = 4
	  iin = 5
	end if

	nb = (nlvddi-n0+1) * n_ghost_max
	call shympi_alloc_buffer(nb)
	allocate(buffer_in(nb,n_ghost_areas))
	allocate(buffer_out(nb,n_ghost_areas))

	do ia=1,n_ghost_areas
	  nc = ghost_areas(iout,ia)
	  call count_buffer(n0,nlvddi,n,nc,il,g_out(:,ia),nb)
	  nbs(1,ia) = nb
	  nc = ghost_areas(iin,ia)
	  call count_buffer(n0,nlvddi,n,nc,il,g_in(:,ia),nb)
	  nbs(2,ia) = nb
	end do

!ccgguccc!$OMP CRITICAL

	do ia=1,n_ghost_areas
	  ir = ir + 1
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(iout,ia)
	  nb = nbs(1,ia)
          call MPI_Irecv(buffer_out(:,ia),nb,MPI_INTEGER,id &
     &	          ,tag,MPI_COMM_WORLD,request(ir),ierr)
	end do

	do ia=1,n_ghost_areas
	  ir = ir + 1
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(iin,ia)
	  nb = nbs(2,ia)
	  call to_buffer_i(n0,nlvddi,n,nc,il &
     &		,g_in(:,ia),val,nb,buffer_in(:,ia))
          call MPI_Isend(buffer_in(:,ia),nb,MPI_INTEGER,id &
     &	          ,tag,MPI_COMM_WORLD,request(ir),ierr)
	end do

        call MPI_WaitAll(ir,request,status,ierr)

!ccgguccc!$OMP END CRITICAL

	do ia=1,n_ghost_areas
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(iout,ia)
	  nb = nbs(1,ia)
	  call from_buffer_i(n0,nlvddi,n,nc,il &
     &		,g_out(:,ia),val,nb,buffer_out(:,ia))
	end do

	call shympi_icomment('shympi_exchange_internal_i')

	end subroutine shympi_exchange_internal_i

!******************************************************************

	subroutine shympi_exchange_internal_r(belem,n0,nlvddi,n,il &
     &						,g_in,g_out,val)

	use shympi_internal
	use shympi

	implicit none

	logical belem
	integer n0,nlvddi,n
	integer il(n)
	integer g_in(n_ghost_max,n_ghost_areas)
	integer g_out(n_ghost_max,n_ghost_areas)
	real val(n0:nlvddi,n)

	integer tag,ir,ia,id
	integer i,k,nc,ierr
	integer nb,nvert
	integer iaux
	integer iout,iin
	integer nbs(2,n_ghost_areas)
	integer status(status_size,2*n_threads)
	integer request(2*n_threads)
	real, allocatable :: buffer_in(:,:)
	real, allocatable :: buffer_out(:,:)
	logical bw

	bw = .false.

        tag=122
	ir = 0

	iout = 2
	iin = 3
	if( belem ) then
	  iout = 4
	  iin = 5
	end if

	nvert = nlv_global
	if( nvert == 0 ) then
	  nvert = 1000
	  if( nlvddi > nvert ) then
	    write(6,*) nlvddi,nvert,nlv_global
	    stop 'error stop shympi_exchange_internal_r: nlvddi>nvert'
	  end if
	  write(6,*) 'increasing nvert ',nvert,my_id
	end if

	nb = (nvert-n0+1) * n_ghost_max_global
	call shympi_alloc_buffer(nb)
	allocate(buffer_in(nb,n_ghost_areas))
	allocate(buffer_out(nb,n_ghost_areas))
	buffer_in = 11111.
	buffer_out = 11111.

	do ia=1,n_ghost_areas
	  nc = ghost_areas(iout,ia)
	  nbs(1,ia) = nb
	  nc = ghost_areas(iin,ia)
	  nbs(2,ia) = nb
	  if( nbs(1,ia) /= nbs(2,ia) ) then
	    nb = (nlvddi-n0+1) * n_ghost_max
	    write(6,*) nb,nbs(1,ia),nbs(2,ia),my_id
	    stop 'error stop nbs'
	  end if
	end do

!ccgguccc!$OMP CRITICAL

	do ia=1,n_ghost_areas
	  ir = ir + 1
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(iout,ia)
	  nb = nbs(1,ia)
	if(bw) write(6,*) 'irec ',belem,ia,nc,nb,id,my_id
          call MPI_Irecv(buffer_out(:,ia),nb,MPI_REAL,id &
     &	          ,tag,MPI_COMM_WORLD,request(ir),ierr)
	end do

	do ia=1,n_ghost_areas
	  ir = ir + 1
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(iin,ia)
	  call to_buffer_r(n0,nlvddi,n,nc,il &
     &		,g_in(:,ia),val,nb,buffer_in(:,ia))
	  nb = nbs(2,ia)
	if(bw) write(6,*) 'isend ',belem,ia,nc,nb,id,my_id
	if(bw) write(6,*) 'buffer before: ',nc,buffer_in(nc,ia),my_id
          call MPI_Isend(buffer_in(:,ia),nb,MPI_REAL,id &
     &	          ,tag,MPI_COMM_WORLD,request(ir),ierr)
	end do

        call MPI_WaitAll(ir,request,status,ierr)

!ccgguccc!$OMP END CRITICAL

	do ia=1,n_ghost_areas
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(iout,ia)
	  nb = nbs(1,ia)
	if(bw) write(6,*) 'copy ',belem,ia,nc,nb,id,my_id
	if(bw) write(6,*) 'buffer after: ',nc,buffer_out(nc,ia),my_id
	  call from_buffer_r(n0,nlvddi,n,nc,il &
     &		,g_out(:,ia),val,nb,buffer_out(:,ia))
	end do

	call shympi_icomment('shympi_exchange_internal_r')

	end subroutine shympi_exchange_internal_r

!******************************************************************

	subroutine shympi_exchange_internal_d(belem,n0,nlvddi,n,il &
     &						,g_in,g_out,val)

	use shympi_internal
	use shympi
	
	logical belem
	integer n0,nlvddi,n
	integer il(n)
	integer g_in(n_ghost_max,n_ghost_areas)
	integer g_out(n_ghost_max,n_ghost_areas)
	double precision val(n0:nlvddi,n)

        integer tag,ir,ia,id
        integer i,k,nc,ierr
        integer nb
        integer iout,iin
        integer nbs(2,n_ghost_areas)
        integer status(status_size,2*n_threads)
        integer request(2*n_threads)
        double precision, allocatable :: buffer_in(:,:)
        double precision, allocatable :: buffer_out(:,:)

        tag=123
        ir = 0

        iout = 2
        iin = 3
        if( belem ) then
          iout = 4
          iin = 5
        end if

        nb = (nlvddi-n0+1) * n_ghost_max
        call shympi_alloc_buffer(nb)
        allocate(buffer_in(nb,n_ghost_areas))
        allocate(buffer_out(nb,n_ghost_areas))

        do ia=1,n_ghost_areas
          nc = ghost_areas(iout,ia)
          call count_buffer(n0,nlvddi,n,nc,il,g_out(:,ia),nb)
          nbs(1,ia) = nb
          nc = ghost_areas(iin,ia)
          call count_buffer(n0,nlvddi,n,nc,il,g_in(:,ia),nb)
          nbs(2,ia) = nb
        end do

!ccgguccc!$OMP CRITICAL

        do ia=1,n_ghost_areas
          ir = ir + 1
          id = ghost_areas(1,ia)
          nc = ghost_areas(iout,ia)
          nb = nbs(1,ia)
          call MPI_Irecv(buffer_out(:,ia),nb,MPI_DOUBLE,id &
     &            ,tag,MPI_COMM_WORLD,request(ir),ierr)
        end do

        do ia=1,n_ghost_areas
          ir = ir + 1
          id = ghost_areas(1,ia)
          nc = ghost_areas(iin,ia)
          nb = nbs(2,ia)
          call to_buffer_d(n0,nlvddi,n,nc,il &
     &          ,g_in(:,ia),val,nb,buffer_in(:,ia))
          call MPI_Isend(buffer_in(:,ia),nb,MPI_DOUBLE,id &
     &            ,tag,MPI_COMM_WORLD,request(ir),ierr)
        end do

        call MPI_WaitAll(ir,request,status,ierr)

!ccgguccc!$OMP END CRITICAL

        do ia=1,n_ghost_areas
          id = ghost_areas(1,ia)
          nc = ghost_areas(iout,ia)
          nb = nbs(1,ia)
          call from_buffer_d(n0,nlvddi,n,nc,il &
     &          ,g_out(:,ia),val,nb,buffer_out(:,ia))
        end do

	call shympi_icomment('shympi_exchange_internal_d')

	end subroutine shympi_exchange_internal_d

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine shympi_allgather_internal_i(ni,no,val,vals)

	use shympi_internal

	implicit none

	integer ni,no
        integer val(ni)
        integer vals(no,n_p_threads)

        integer ierr,nn
	integer, allocatable :: aux(:)

	if( ni /= no ) then
	  write(6,*) 'ni /= no... ',ni,no
	  !commenting next statement creates mpi error and backtrace
	  stop 'error stop shympi_allgather_internal_i: ni /= no'
	end if

	allocate(aux(no))
	aux = 0.
	aux(1:ni) = val(1:ni)
	nn = no

	if( bpmpi ) then
          call MPI_ALLGATHER (aux,nn,MPI_INTEGER &
     &                  ,vals,no,MPI_INTEGER &
     &                  ,MPI_COMM_WORLD,ierr)
	  call shympi_error('shympi_allgather_internal_i' &
     &			,'allgather',ierr)
	else
	  vals(1:ni,1) = val(:)
	end if

	call shympi_icomment('shympi_allgather_internal_i')

        end subroutine shympi_allgather_internal_i

!******************************************************************

        subroutine shympi_allgather_internal_r(ni,no,val,vals)

	use shympi_internal

	implicit none

	integer ni,no
        real val(ni)
        real vals(no,n_p_threads)

        integer ierr,nn
        real, allocatable :: aux(:)

	if( ni /= no ) then
	  write(6,*) 'ni /= no... ',ni,no
	  !commenting next statement creates mpi error and backtrace
	  stop 'error stop shympi_allgather_internal_r: ni /= no'
	end if

	allocate(aux(no))
	aux = 0.
	aux(1:ni) = val(1:ni)
	nn = no

	if( bpmpi ) then
          call MPI_ALLGATHER (aux,nn,MPI_REAL &
     &                  ,vals,no,MPI_REAL &
     &                  ,MPI_COMM_WORLD,ierr)
	  call shympi_error('shympi_allgather_internal_r' &
     &			,'allgather',ierr)
	else
	  vals(1:ni,1) = val(:)
	end if

	call shympi_icomment('shympi_allgather_internal_r')

        end subroutine shympi_allgather_internal_r

!******************************************************************

        subroutine shympi_allgather_internal_d(ni,no,val,vals)

	use shympi_internal

	implicit none

	integer ni,no
        double precision val(ni)
        double precision vals(no,n_p_threads)

        integer ierr,nn
        double precision, allocatable :: aux(:)

	if( ni /= no ) then
	  write(6,*) 'ni /= no... ',ni,no
	  !commenting next statement creates mpi error and backtrace
	  stop 'error stop shympi_allgather_internal_d: ni /= no'
	end if

	allocate(aux(no))
	aux = 0.
	aux(1:ni) = val(1:ni)
	nn = no

	if( bpmpi ) then
          call MPI_ALLGATHER (aux,nn,MPI_DOUBLE &
     &                  ,vals,no,MPI_DOUBLE &
     &                  ,MPI_COMM_WORLD,ierr)
	  call shympi_error('shympi_allgather_internal_d' &
     &			,'allgather',ierr)
	else
	  vals(1:ni,1) = val(:)
	end if

	call shympi_icomment('shympi_allgather_internal_d')

        end subroutine shympi_allgather_internal_d

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine shympi_gather_internal_i(ni,no,val,vals)

	use shympi_internal

	implicit none

	integer ni,no
        integer val(ni)
        integer vals(no,n_p_threads)

        integer ierr

	if( ni /= no ) then
	  write(6,*) 'ni /= no... ',ni,no
	  !commenting next statement creates mpi error and backtrace
	  stop 'error stop shympi_gather_internal_i: ni /= no'
	end if

	if( bpmpi ) then
          call MPI_GATHER (val,ni,MPI_INTEGER &
     &                  ,vals,no,MPI_INTEGER &
     &			,0 &
     &                  ,MPI_COMM_WORLD,ierr)
	  call shympi_error('shympi_gather_internal_i' &
     &			,'gather',ierr)
	else
	  vals(1:ni,1) = val(:)
	end if

	call shympi_icomment('shympi_gather_internal_i')

        end subroutine shympi_gather_internal_i

!******************************************************************

        subroutine shympi_gather_internal_r(ni,no,val,vals)

	use shympi_internal

	implicit none

	integer ni,no
        real val(ni)
        real vals(no,n_p_threads)

        integer ierr

	if( ni /= no ) then
	  write(6,*) 'ni /= no... ',ni,no
	  !commenting next statement creates mpi error and backtrace
	  stop 'error stop shympi_gather_internal_r: ni /= no'
	end if

	if( bpmpi ) then
          call MPI_GATHER (val,ni,MPI_REAL &
     &                  ,vals,no,MPI_REAL &
     &			,0 &
     &                  ,MPI_COMM_WORLD,ierr)
	  call shympi_error('shympi_gather_internal_r' &
     &			,'gather',ierr)
	else
	  vals(1:ni,1) = val(:)
	end if

	call shympi_icomment('shympi_gather_internal_r')

        end subroutine shympi_gather_internal_r

!******************************************************************

        subroutine shympi_gather_internal_d(ni,no,val,vals)

	use shympi_internal

	implicit none

	integer ni,no
        double precision val(ni)
        double precision vals(no,n_p_threads)

        integer ierr

	if( ni /= no ) then
	  write(6,*) 'ni /= no... ',ni,no
	  !commenting next statement creates mpi error and backtrace
	  stop 'error stop shympi_gather_internal_d: ni /= no'
	end if

	if( bpmpi ) then
          call MPI_GATHER (val,ni,MPI_DOUBLE &
     &                  ,vals,no,MPI_DOUBLE &
     &			,0 &
     &                  ,MPI_COMM_WORLD,ierr)
	  call shympi_error('shympi_gather_internal_d' &
     &			,'gather',ierr)
	else
	  vals(1:ni,1) = val(:)
	end if

	call shympi_icomment('shympi_gather_internal_d')

        end subroutine shympi_gather_internal_d

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine shympi_bcast_internal_i(n,val)

	use shympi_internal

	implicit none

	integer n
        integer val(n)

        integer ierr

	if( bpmpi ) then
          call MPI_BCAST(val,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	  call shympi_error('shympi_bcast_internal_i','bcast',ierr)
	end if

	call shympi_icomment('shympi_bcast_internal_i')

        end subroutine shympi_bcast_internal_i

!******************************************************************

        subroutine shympi_bcast_internal_r(n,val)

	use shympi_internal

	implicit none

	integer n
        real val(n)

        integer ierr

	if( bpmpi ) then
          call MPI_BCAST(val,n,MPI_REAL,0,MPI_COMM_WORLD,ierr)
	  call shympi_error('shympi_bcast_internal_r','bcast',ierr)
	end if

	call shympi_icomment('shympi_bcast_internal_r')

        end subroutine shympi_bcast_internal_r

!******************************************************************

        subroutine shympi_bcast_internal_d(n,val)

	use shympi_internal

	implicit none

	integer n
        double precision val(n)

        integer ierr

	if( bpmpi ) then
          call MPI_BCAST(val,n,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
	  call shympi_error('shympi_bcast_internal_d','bcast',ierr)
	end if

	call shympi_icomment('shympi_bcast_internal_d')

        end subroutine shympi_bcast_internal_d

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_reduce_internal_d(what,val)

	use shympi_internal

	implicit none

	character*(*) what
	double precision val

        integer ierr
	integer mpi_what
	double precision valin(1)
	double precision valout(1)

	ierr = 0
	if( bpmpi ) then
	  valin(1) = val
	  call shympi_internal_reduce_what(what,mpi_what)
	  call MPI_ALLREDUCE(valin,valout,1,MPI_DOUBLE,mpi_what &
     &				,MPI_COMM_WORLD,ierr)
	  val = valout(1)
	end if

	call shympi_error('shympi_reduce_internal_d','reduce',ierr)
	call shympi_icomment('shympi_reduce_internal_d')

	end subroutine shympi_reduce_internal_d

!******************************************************************

	subroutine shympi_reduce_internal_r(what,val)

	use shympi_internal

	implicit none

	character*(*) what
	real val

        integer ierr
	integer mpi_what
	real valin(1)
	real valout(1)

	ierr = 0
	if( bpmpi ) then
	  valin(1) = val
	  call shympi_internal_reduce_what(what,mpi_what)
	  call MPI_ALLREDUCE(valin,valout,1,MPI_REAL,mpi_what &
     &				,MPI_COMM_WORLD,ierr)
	  val = valout(1)
	end if

	call shympi_error('shympi_reduce_internal_r','reduce',ierr)
	call shympi_icomment('shympi_reduce_internal_r')

	end subroutine shympi_reduce_internal_r

!******************************************************************

	subroutine shympi_reduce_internal_i(what,val)

	use shympi_internal

	implicit none

	character*(*) what
	integer val

        integer ierr
	integer mpi_what
	integer valin(1)
	integer valout(1)

	ierr = 0
	if( bpmpi ) then
	  valin(1) = val
          call shympi_internal_reduce_what(what,mpi_what)
          call MPI_ALLREDUCE(valin,valout,1,MPI_INTEGER,mpi_what &
     &                          ,MPI_COMM_WORLD,ierr)
          val = valout(1)
	end if

	call shympi_error('shympi_reduce_internal_i','reduce',ierr)
	call shympi_icomment('shympi_reduce_internal_i')

	end subroutine shympi_reduce_internal_i

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_reduce_array_internal_i(what,n,vals)

	use shympi_internal

	implicit none

	character*(*) what
	integer n
	integer vals(n)

        integer ierr
	integer mpi_what
	integer valout(n)

	ierr = 0
	if( bpmpi ) then
          call shympi_internal_reduce_what(what,mpi_what)
	  if( bpoperall ) then
            call MPI_ALLREDUCE(vals,valout,n,MPI_INTEGER,mpi_what &
     &                          ,MPI_COMM_WORLD,ierr)
	  else
            call MPI_REDUCE(vals,valout,n,MPI_INTEGER,mpi_what,0 &
     &                          ,MPI_COMM_WORLD,ierr)
	  end if
          vals = valout
	end if

	call shympi_error('shympi_reduce_array_internal_i','reduce',ierr)
	call shympi_icomment('shympi_reduce_array_internal_i',bpoperall)

	end subroutine shympi_reduce_array_internal_i

!******************************************************************

	subroutine shympi_reduce_array_internal_r(what,n,vals)

	use shympi_internal

	implicit none

	character*(*) what
	integer n
	real vals(n)

        integer ierr
	integer mpi_what
	real valout(n)

	ierr = 0
	if( bpmpi ) then
          call shympi_internal_reduce_what(what,mpi_what)
	  if( bpoperall ) then
            call MPI_ALLREDUCE(vals,valout,n,MPI_REAL,mpi_what &
     &                          ,MPI_COMM_WORLD,ierr)
	  else
            call MPI_REDUCE(vals,valout,n,MPI_REAL,mpi_what,0 &
     &                          ,MPI_COMM_WORLD,ierr)
	  end if
          vals = valout
	end if

	call shympi_error('shympi_reduce_array_internal_r','reduce',ierr)
	call shympi_icomment('shympi_reduce_array_internal_r',bpoperall)

	end subroutine shympi_reduce_array_internal_r

!******************************************************************

	subroutine shympi_reduce_array_internal_d(what,n,vals)

	use shympi_internal

	implicit none

	character*(*) what
	integer n
	double precision vals(n)

        integer ierr
	integer mpi_what
	double precision valout(n)

	ierr = 0
	if( bpmpi ) then
          call shympi_internal_reduce_what(what,mpi_what)
	  if( bpoperall ) then
            call MPI_ALLREDUCE(vals,valout,n,MPI_DOUBLE,mpi_what &
     &                          ,MPI_COMM_WORLD,ierr)
	  else
            call MPI_REDUCE(vals,valout,n,MPI_DOUBLE,mpi_what,0 &
     &                          ,MPI_COMM_WORLD,ierr)
	  end if
          vals = valout
	end if

	call shympi_error('shympi_reduce_array_internal_d','reduce',ierr)
	call shympi_icomment('shympi_reduce_array_internal_d',bpoperall)

	end subroutine shympi_reduce_array_internal_d

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_getvals_internal_r(kind,nlvddi,n &
     &						,val_in,val_out)

	use shympi_internal

	implicit none

	integer kind(2)
	integer nlvddi,n
	real val_in(nlvddi,n)
	real val_out(nlvddi)

	integer id,k,nb,lmax
	integer ierr

	id = kind(2) - 1
	k = kind(1)
	lmax = nlvddi
	nb = lmax

	if( my_p_id == id ) then
	  val_out(1:lmax) = val_in(1:lmax,k)
	end if

        call MPI_BCAST(val_out,nb,MPI_REAL,id &
     &	          ,MPI_COMM_WORLD,ierr)

	call shympi_icomment('shympi_getvals_internal_r')

	end subroutine shympi_getvals_internal_r

!******************************************************************

	subroutine shympi_getvals_internal_i(kind,nlvddi,n &
     &						,val_in,val_out)

	use shympi_internal

	implicit none

	integer kind(2)
	integer nlvddi,n
	integer val_in(nlvddi,n)
	integer val_out(nlvddi)

	integer id,k,nb,lmax
	integer ierr

	id = kind(2) - 1
	k = kind(1)
	lmax = nlvddi
	nb = lmax

	if( my_p_id == id ) then
	  val_out(1:lmax) = val_in(1:lmax,k)
	end if

        call MPI_BCAST(val_out,nb,MPI_INTEGER,id &
     &	          ,MPI_COMM_WORLD,ierr)

	call shympi_icomment('shympi_getvals_internal_i')

	end subroutine shympi_getvals_internal_i

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_get_array_internal_r(nlvddi,n &
     &						,val_in,val_out)

	use shympi_internal
	use shympi

	implicit none

	integer kind(2)
	integer nlvddi,n
	real val_in(nlvddi,*)
	real val_out(nlvddi,n)

	logical bnode,belem
	integer i,ir,ns,nb,tag,id
	integer ierr
	integer ip(0:n_threads)
	integer req(2*n_threads)
	integer status(status_size,2*n_threads)

        tag=131
	ir = 0

	bnode = ( n == nkn_global )
	belem = ( n == nel_global )

	if( bnode ) then
	  ip = nkn_cum_domains
	else if( belem ) then
	  ip = nel_cum_domains
	else
	  write(6,*) 'n,nkn_global,nel_global: ',n,nkn_global,nel_global
	  call shympi_stop('error stop shympi_get_array_internal_i:'// &
     &				' size of out array')
	end if

!ccgguccc!$OMP CRITICAL

	if( my_id == 0 ) then
	  do i=2,n_threads
	    id = i - 1
	    ir = ir + 1
	    ns = nlvddi*ip(i-1) + 1
	    nb = nlvddi*(ip(i) - ip(i-1))
            call MPI_Irecv(val_out(1,ns),nb,MPI_REAL,id &
     &	          ,tag,MPI_COMM_WORLD,req(ir),ierr)
	  end do
	  nb = ip(1)
	  val_out(:,1:nb) = val_in(:,1:nb)
	else
	    i = my_id + 1
	    ir = ir + 1
	    nb = nlvddi*(ip(i) - ip(i-1))
            call MPI_Isend(val_in,nb,MPI_REAL,0 &
     &	          ,tag,MPI_COMM_WORLD,req(ir),ierr)
	end if

        call MPI_WaitAll(ir,req,status,ierr)

!ccgguccc!$OMP END CRITICAL

	call shympi_icomment('shympi_get_array_internal_r')

	end subroutine shympi_get_array_internal_r

!******************************************************************

	subroutine shympi_get_array_internal_i(nlvddi,n &
     &						,val_in,val_out)

	use shympi_internal
	use shympi

	implicit none

	integer nlvddi,n
	integer val_in(nlvddi,*)
	integer val_out(nlvddi,n)

	logical bnode,belem
	integer i,ir,ns,nb,tag,id
	integer ierr
	integer ip(0:n_threads)
	integer req(2*n_threads)
	integer status(status_size,2*n_threads)

        tag=132
	ir = 0

	bnode = ( n == nkn_global )
	belem = ( n == nel_global )

	if( bnode ) then
	  ip = nkn_cum_domains
	else if( belem ) then
	  ip = nel_cum_domains
	else
	  write(6,*) 'n,nkn_global,nel_global: ',n,nkn_global,nel_global
	  call shympi_stop('error stop shympi_get_array_internal_i:'// &
     &				' size of out array')
	end if

!ccgguccc!$OMP CRITICAL

	if( my_id == 0 ) then
	  do i=2,n_threads
	    id = i - 1
	    ir = ir + 1
	    ns = nlvddi*ip(i-1) + 1
	    nb = nlvddi*(ip(i) - ip(i-1))
            call MPI_Irecv(val_out(1,ns),nb,MPI_INTEGER,id &
     &	          ,tag,MPI_COMM_WORLD,req(ir),ierr)
	  end do
	  nb = ip(1)
	  val_out(:,1:nb) = val_in(:,1:nb)
	else
	    i = my_id + 1
	    ir = ir + 1
	    nb = nlvddi*(ip(i) - ip(i-1))
            call MPI_Isend(val_in,nb,MPI_INTEGER,0 &
     &	          ,tag,MPI_COMM_WORLD,req(ir),ierr)
	end if

        call MPI_WaitAll(ir,req,status,ierr)

!ccgguccc!$OMP END CRITICAL

	call shympi_icomment('shympi_get_array_internal_i')

	end subroutine shympi_get_array_internal_i

!******************************************************************
!******************************************************************
!******************************************************************
! next rectify routines are useless and can be deleted
!******************************************************************
!******************************************************************
!******************************************************************

        subroutine shympi_rectify_internal_r0(nv,nh,vals)

	use shympi

        integer nv,nh
        real vals(nh*nv,n_threads)

	integer ia,ip,it,nlv,nextra
        real, allocatable :: vaux(:)

        allocate(vaux(nh*nv))

	nextra = 0
        if( nv == nlv_global+1 ) nextra = 1

        do ia=1,n_threads
          nlv = nlv_domains(ia) + nextra
          vaux(:) = vals(:,ia)
          vals(:,ia) = 0.
          ip = 0
          it = 0
          do n=1,nh
            vals(it+1:it+nlv,ia) = vaux(ip+1:ip+nlv)
	    ip = ip + nlv
	    it = it + nv
          end do
        end do

	end subroutine shympi_rectify_internal_r0

!******************************************************************

        subroutine shympi_rectify_internal_i0(nv,nh,vals)

	use shympi

        integer nv,nh
        integer vals(nh*nv,n_threads)

	integer ia,ip,it,nlv,nextra
        integer, allocatable :: vaux(:)

        allocate(vaux(nh*nv))

	nextra = 0
        if( nv == nlv_global+1 ) nextra = 1

        do ia=1,n_threads
          nlv = nlv_domains(ia) + nextra
          vaux(:) = vals(:,ia)
          vals(:,ia) = 0.
          ip = 0
          it = 0
          do n=1,nh
            vals(it+1:it+nlv,ia) = vaux(ip+1:ip+nlv)
	    ip = ip + nlv
	    it = it + nv
          end do
        end do

	end subroutine shympi_rectify_internal_i0

!******************************************************************

        subroutine shympi_rectify_internal_d0(nv,nh,vals)

	use shympi

        integer nv,nh
        double precision vals(nh*nv,n_threads)

	integer ia,ip,it,nlv,nextra
        double precision, allocatable :: vaux(:)

        allocate(vaux(nh*nv))

	nextra = 0
        if( nv == nlv_global+1 ) nextra = 1

        do ia=1,n_threads
          nlv = nlv_domains(ia) + nextra
          vaux(:) = vals(:,ia)
          vals(:,ia) = 0.
          ip = 0
          it = 0
          do n=1,nh
            vals(it+1:it+nlv,ia) = vaux(ip+1:ip+nlv)
	    ip = ip + nlv
	    it = it + nv
          end do
        end do

	end subroutine shympi_rectify_internal_d0

!******************************************************************

