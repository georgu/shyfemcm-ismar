
!==========================================================================
        module mod_dbg_flx
!==========================================================================

        implicit none

        logical, save :: bsilent = .false.      !be silent
        logical, save :: bquiet = .false.       !be quiet
        logical, save :: bverbose = .false.     !be verbose
        logical, save :: bcheck = .true.        !check for differences
        logical, save :: bstop = .true.         !stop on error
        logical, save :: bnostop = .false.      !do not stop on differences
        logical, save :: bnodiff = .true.       !do not show differences
        logical, save :: bsummary = .true.      !only show summary
        logical, save :: bbalance = .true.      !balance time step
        integer, save :: maxdiff = 0.           !max difference allowed

!==========================================================================
        end module mod_dbg_flx
!==========================================================================

        program dbg_flx

        use clo
        use mod_dbg_flx

        implicit none

        integer nc,ierr

        call dbg_flx_init

        nc = clo_number_of_files()

        if( nc == 0 ) then
          call clo_usage
        else if( nc == 1 ) then
          call read_dbg_flx_file
        else if( nc == 2 ) then
          call compare_dbg_flx_files(ierr)
        else
          write(6,*) 'nc = ',nc
          stop 'error stop dbg_flx: wrong number of files'
        end if

        if( ierr > 0 ) then
          if( ierr == 99 ) ierr = 100   !terrible hack - FIXME
          call exit(ierr)
        else
          call exit(99)
        end if

        end

!**************************************************************************

        subroutine read_dbg_flx_file

! reads one file and outputs info

        use clo
        use mod_dbg_flx

        implicit none

        integer nc
        integer ntime,nrec
        integer nvers,nsect,kfluxm,idtflx
        integer nlmax,nvar,ivar,iv,ierr
	integer is
	integer lmax,l
	integer nh,nv,nt
	integer iu,iu1
        integer ios
        double precision dtime,atime0,atime
        character*60 name_one,text
        character*80 title,femver

	integer, save, allocatable :: kflux(:)
	integer, save, allocatable :: nlayers(:)
	real, save, allocatable :: fluxes(:,:,:)
	character*80, save, allocatable :: strings(:)

!--------------------------------------------------
! open file(s)
!--------------------------------------------------

        call clo_get_file(1,name_one)

	iu1 = 1
        open(iu1,file=name_one,status='old',form='unformatted',iostat=ios)

        if( ios /= 0 ) then
          write(6,*) 'no such file: ',trim(name_one)
          stop 'error opening file'
        end if

        if( .not. bquiet ) write(6,*) 'opened file: ',trim(name_one)

	call flx_is_flx_file(iu1,nvers)
	if( nvers == 0 ) then
          write(6,*) 'file is not a flx file: ',trim(name_one)
          stop 'error opening file'
	end if

!--------------------------------------------------
! read header(s)
!--------------------------------------------------

	iu = iu1
        call flx_read_header(iu1,nvers,nsect,kfluxm,idtflx &
     &                                  ,nlmax,nvar,ierr)
	if( ierr /= 0 ) goto 99

	allocate(kflux(kfluxm))
	allocate(nlayers(nsect))
	allocate(strings(nsect))
	allocate(fluxes(0:nlmax,3,nsect))

        call flx_read_header2(iu1,nvers,nsect,kfluxm &
     &                          ,kflux,nlayers &
     &                          ,atime0,title,femver,strings,ierr)
	if( ierr /= 0 ) goto 98

!--------------------------------------------------
! output info
!--------------------------------------------------

	if( bverbose ) then
	  write(6,*) 'total sections:      ',nsect
	  write(6,*) 'total nodes:         ',kfluxm
	  write(6,*) 'output time step:    ',idtflx
	  write(6,*) 'maximum layers:      ',nlmax
	  write(6,*) 'number of variables: ',nvar
	  write(6,*) '    section     nlayers     description'
	  do is=1,nsect
	    write(6,*) is,nlayers(is),trim(strings(is))
	  end do
	end if

!--------------------------------------------------
! start of time loop
!--------------------------------------------------

        ntime = 0
	iv = 0

        do while(.true.)

	  iu = iu1
          call flx_read_record(iu1,nvers,atime &
     &                  ,nlmax,nsect,ivar &
     &                  ,nlayers,fluxes,ierr)

	  if( ierr < 0 ) exit
	  if( ierr > 0 ) goto 97

	  iv = iv + 1
          if( mod(iv,nvar) == 1 ) ntime = ntime + 1
          if( .not. bquiet ) write(6,*) 'time = ',atime,ntime,ivar
	  if( bverbose ) then
	    do is=1,nsect
	      lmax = nlayers(is)
	      write(6,'(2i8)') is,lmax
	      do l=1,lmax
	        write(6,*) l,fluxes(l,:,is)
	      end do
	    end do
	  end if

        end do

!--------------------------------------------------
! end of time loop
!--------------------------------------------------

        if( .not. bsilent ) write(6,*) 'time records read: ',ntime

!--------------------------------------------------
! end of routine
!--------------------------------------------------

	return
   97	continue
	write(6,*) 'on unit = ',iu
	write(6,*) 'error reading record: ',ierr
	stop 'error stop read_dbg_flx_file: error reading record'
   98	continue
	write(6,*) 'on unit = ',iu
	write(6,*) 'error reading header2: ',ierr
	stop 'error stop read_dbg_flx_file: error reading header2'
   99	continue
	write(6,*) 'on unit = ',iu
	write(6,*) 'error reading header: ',ierr
	stop 'error stop read_dbg_flx_file: error reading header'
        end

!**************************************************************************

        subroutine compare_dbg_flx_files(ierr)

! reads one file and outputs info

        use clo
        use mod_dbg_flx

        implicit none

	integer ierr

        integer nc
        integer ntime,nrec
        integer nvers,nsect,kfluxm,idtflx,nlmax,nvar
        integer nvers2,nsect2,kfluxm2,idtflx2,nlmax2,nvar2
        integer ivar,ivar2,iv
	integer is
	integer lmax,l
	integer nh,nv,nt
	integer iu,iu1,iu2
        integer ios
	double precision eps
        double precision dtime,atime0,atime02,atime,atime2
        character*60 name_one,name_two,text
        character*80 title,title2,femver

	integer, save, allocatable :: kflux(:)
	integer, save, allocatable :: nlayers(:)
	real, save, allocatable :: fluxes(:,:,:)
	character*80, save, allocatable :: strings(:)
	integer, save, allocatable :: kflux2(:)
	integer, save, allocatable :: nlayers2(:)
	real, save, allocatable :: fluxes2(:,:,:)
	character*80, save, allocatable :: strings2(:)

!--------------------------------------------------
! open file(s)
!--------------------------------------------------

        call clo_get_file(1,name_one)

	iu1 = 1
        open(iu1,file=name_one,status='old',form='unformatted',iostat=ios)

        if( ios /= 0 ) then
          write(6,*) 'no such file: ',trim(name_one)
          stop 'error opening file'
        end if

        if( .not. bquiet ) write(6,*) 'opened file: ',trim(name_one)

	call flx_is_flx_file(iu1,nvers)
	if( nvers == 0 ) then
          write(6,*) 'file is not a flx file: ',trim(name_one)
          stop 'error opening file'
	end if

        call clo_get_file(2,name_two)

	iu2 = 2
        open(iu2,file=name_one,status='old',form='unformatted',iostat=ios)

        if( ios /= 0 ) then
          write(6,*) 'no such file: ',trim(name_one)
          stop 'error opening file'
        end if

        if( .not. bquiet ) write(6,*) 'opened file: ',trim(name_one)

	call flx_is_flx_file(iu2,nvers)
	if( nvers == 0 ) then
          write(6,*) 'file is not a flx file: ',trim(name_one)
          stop 'error opening file'
	end if

!--------------------------------------------------
! read header(s)
!--------------------------------------------------

	iu = iu1
        call flx_read_header(iu1,nvers,nsect,kfluxm,idtflx &
     &                                  ,nlmax,nvar,ierr)
	if( ierr /= 0 ) goto 99

	allocate(kflux(kfluxm))
	allocate(nlayers(nsect))
	allocate(strings(nsect))
	allocate(fluxes(0:nlmax,3,nsect))

        call flx_read_header2(iu1,nvers,nsect,kfluxm &
     &                          ,kflux,nlayers &
     &                          ,atime0,title,femver,strings,ierr)
	if( ierr /= 0 ) goto 98

	iu = iu2
        call flx_read_header(iu2,nvers2,nsect2,kfluxm2,idtflx2 &
     &                                  ,nlmax2,nvar2,ierr)
	if( ierr /= 0 ) goto 99

	if( nvers /= nvers2 ) goto 95
	if( nsect /= nsect2 ) goto 95
	if( kfluxm /= kfluxm2 ) goto 95
	if( idtflx /= idtflx2 ) goto 95
	if( nlmax /= nlmax2 ) goto 95
	if( nvar /= nvar2 ) goto 95

	allocate(kflux2(kfluxm2))
	allocate(nlayers2(nsect2))
	allocate(strings2(nsect2))
	allocate(fluxes2(0:nlmax2,3,nsect2))

        call flx_read_header2(iu2,nvers2,nsect2,kfluxm2 &
     &                          ,kflux2,nlayers2 &
     &                          ,atime02,title2,femver,strings2,ierr)
	if( ierr /= 0 ) goto 98

	if( any( kflux /= kflux2 ) ) goto 94
	if( any( nlayers /= nlayers2 ) ) goto 94

!--------------------------------------------------
! output info
!--------------------------------------------------

	if( bverbose ) then
	  write(6,*) 'total sections:      ',nsect
	  write(6,*) 'total nodes:         ',kfluxm
	  write(6,*) 'output time step:    ',idtflx
	  write(6,*) 'maximum layers:      ',nlmax
	  write(6,*) 'number of variables: ',nvar
	  write(6,*) '    section     nlayers     description'
	  do is=1,nsect
	    write(6,*) is,nlayers(is),trim(strings(is))
	  end do
	end if

!--------------------------------------------------
! start of time loop
!--------------------------------------------------

        ntime = 0
	iv = 0
	eps = maxdiff

        do while(.true.)

	  iu = iu1
          call flx_read_record(iu1,nvers,atime &
     &                  ,nlmax,nsect,ivar &
     &                  ,nlayers,fluxes,ierr)

	  if( ierr < 0 ) exit
	  if( ierr > 0 ) goto 97

	  iu = iu2
          call flx_read_record(iu2,nvers,atime2 &
     &                  ,nlmax,nsect2,ivar2 &
     &                  ,nlayers2,fluxes2,ierr)

	  if( ierr < 0 ) exit
	  if( ierr > 0 ) goto 97

	  if( atime /= atime2 ) goto 93
	  if( nsect /= nsect2 ) goto 93
	  if( ivar /= ivar2 ) goto 93
	  if( any( nlayers /= nlayers2 ) ) goto 93
	  
	  iv = iv + 1
          if( mod(iv,nvar) == 1 ) ntime = ntime + 1
          if( .not. bquiet ) write(6,*) 'time = ',atime,ntime,ivar
	  do is=1,nsect
	    lmax = nlayers(is)
	    do l=1,lmax
	      if( any( abs(fluxes(l,:,is)-fluxes(l,:,is)) > eps ) ) then
	        write(6,*) is,l,lmax
	        write(6,*) fluxes(l,:,is)
	        write(6,*) fluxes2(l,:,is)
	      end if
	    end do
	  end do

        end do

!--------------------------------------------------
! end of time loop
!--------------------------------------------------

        if( .not. bsilent ) write(6,*) 'time records read: ',ntime

!--------------------------------------------------
! end of routine
!--------------------------------------------------

	return
   93	continue
	write(6,*) 'incompatible parameters while reading records'
	write(6,*) 'nsect: ',nsect,nsect2
	write(6,*) 'ivar: ',ivar,ivar2
	write(6,*) 'nlayers: ',nlayers
	write(6,*) 'nlayers2: ',nlayers2
	stop 'error stop compare_dbg_flx_files: incompatible records'
   94	continue
	write(6,*) 'incompatible flx arrays:'
	write(6,*) 'kflux: ',kflux
	write(6,*) 'kflux2: ',kflux2
	write(6,*) 'nlayers: ',nlayers
	write(6,*) 'nlayers2: ',nlayers2
	stop 'error stop compare_dbg_flx_files: incompatible arrays'
   95	continue
	write(6,*) 'incompatible flx parameters:'
	write(6,*) 'nvers: ',nvers,nvers2
	write(6,*) 'nsect: ',nsect,nsect2
	write(6,*) 'kfluxm: ',kfluxm,kfluxm2
	write(6,*) 'idtflx: ',idtflx,idtflx2
	write(6,*) 'nlmax: ',nlmax,nlmax2
	write(6,*) 'nvar: ',nvar,nvar2
	stop 'error stop compare_dbg_flx_files: incompatible parameters'
   97	continue
	write(6,*) 'on unit = ',iu
	write(6,*) 'error reading record: ',ierr
	stop 'error stop compare_dbg_flx_files: error reading record'
   98	continue
	write(6,*) 'on unit = ',iu
	write(6,*) 'error reading header2: ',ierr
	stop 'error stop compare_dbg_flx_files: error reading header2'
   99	continue
	write(6,*) 'on unit = ',iu
	write(6,*) 'error reading header: ',ierr
	stop 'error stop compare_dbg_flx_files: error reading header'
        end

!**************************************************************************
!**************************************************************************
!**************************************************************************

        subroutine dbg_flx_init

        use clo
        use mod_dbg_flx

        implicit none

        logical baux
        character*80 version

        version = '1.0'

        call clo_init('dbg_flx_init','flx-file(s)',trim(version))

        call clo_add_info('reads flx files')

        call clo_add_sep('general options:')
        call clo_add_option('silent',.false.,'be silent')
        call clo_add_option('quiet',.false.,'be quiet')
        call clo_add_option('nodiff',.false.,'do not show differences')
        call clo_add_option('verbose',.false.,'be verbose')
        call clo_add_option('nostop',.false.,'do not stop at error')
        call clo_add_option('summary',.false.,'do only summary')
        call clo_add_option('balance',.false.,'balance time records')
        call clo_add_option('maxdiff',0.,'maximum tolerated difference')

        call clo_parse_options

        call clo_get_option('silent',bsilent)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('nodiff',bnodiff)
        call clo_get_option('verbose',bverbose)
        call clo_get_option('nostop',baux)
        call clo_get_option('summary',bsummary)
        call clo_get_option('balance',bbalance)
        call clo_get_option('maxdiff',maxdiff)

        if( baux ) bstop = .false.
        if( bsilent ) bquiet = .true.
        if( bquiet ) bverbose = .false.

        end

!**************************************************************************


