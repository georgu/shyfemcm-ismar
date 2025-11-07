
!=============================================================
	module mod_assimil_admin
!=============================================================

	implicit none

	integer, parameter, private :: ndim = 100
	integer, save, private :: idlast = 0

	type, private :: info
	  integer :: num
	  integer :: nobs
	  integer :: ivar
	  real, allocatable :: corlength(:)
	  real, allocatable :: maxlength(:)
	  real, allocatable :: sigobs(:)
	  real, allocatable :: sigback(:)
	  character*80      :: filecoords
	  character*80      :: fileobs
	  real, allocatable :: xobs(:)
	  real, allocatable :: yobs(:)
	  double precision, allocatable :: times(:)
	  real, allocatable :: zobss(:,:)
	end type info

	type(info), save, target, allocatable :: pinfo(:)

!=============================================================
	contains
!=============================================================

	subroutine assimil_init

	integer, save :: icall = 0

	if( icall /= 0 ) return
	allocate(pinfo(ndim))
	icall = 1

	end

!*************************************************************

	function is_num_in_list(num)

	logical is_num_in_list
	integer num

	integer i

	is_num_in_list = .true.

	do i=1,idlast
	  if( pinfo(i)%num == num ) return
	end do

	is_num_in_list = .false.

	end

!*************************************************************

	function is_var_in_list(ivar)

	logical is_var_in_list
	integer ivar

	integer i

	is_var_in_list = .true.

	do i=1,idlast
	  if( pinfo(i)%ivar == ivar ) return
	end do

	is_var_in_list = .false.

	end

!*************************************************************

	function get_new_id()

	integer get_new_id

	idlast = idlast + 1

	if( idlast > ndim ) then
	  write(6,*) 'too many assimil sections'
	  write(6,*) 'please increase ndim in assimil_admin.f90'
	  stop 'error stop get_new_id: idlast > ndim'
	end if

	get_new_id = idlast

	end

!=============================================================
	end module mod_assimil_admin
!=============================================================

	subroutine rdassimil(num)

	use mod_assimil_admin

	implicit none

	integer num

	call assimil_init

	if( num <= 0 ) then
	  write(6,*) 'section assimil ',num,' not allowed'
	  write(6,*) 'must specify the number of the nudging section'
	  stop 'error stop rdassimil: num <= 0 or not given'
	end if

	if( is_num_in_list(num) ) then
	  write(6,*) 'section assimil ',num,' is already defined'
	  stop 'error stop rdassimil: section not unique'
	end if

	call init_assimil_section(num)

	end

!*************************************************************

	subroutine init_assimil_section(num)

	use mod_assimil_admin
	use para
	use nls

	implicit none

	integer num

	integer ivar,nobs
	integer n,id

	real getpar

        call sctpar('assimil')
        call sctfnm('assimil')

	call addpar('ivar',-1.)
	call addpar('nobs',0.)

	call para_add_array_value('corlength',0.)
	call para_add_array_value('maxlength',0.)
	call para_add_array_value('sigobs',0.)
	call para_add_array_value('sigback',0.)

	call addfnm('filecoords',' ')
	call addfnm('fileobs',' ')

	call nls_read_namelist('assimil')		!reads all parameters

	nobs = nint(getpar('nobs'))
	ivar = nint(getpar('ivar'))
	if( ivar < 0 .or. nobs <= 0 ) goto 99
	if( is_var_in_list(ivar) ) goto 98

	id = get_new_id()	

	allocate(pinfo(id)%corlength(nobs))
	allocate(pinfo(id)%maxlength(nobs))
	allocate(pinfo(id)%sigobs(nobs))
	allocate(pinfo(id)%sigback(nobs))

	pinfo(id)%nobs = nobs
	pinfo(id)%ivar = ivar

	call para_set_array_value('corlength',nobs,pinfo(id)%corlength)
	call para_set_array_value('maxlength',nobs,pinfo(id)%maxlength)
	call para_set_array_value('sigobs',nobs,pinfo(id)%sigobs)
	call para_set_array_value('sigback',nobs,pinfo(id)%sigback)

	call getfnm('filecoords',pinfo(id)%filecoords)
	call getfnm('fileobs',pinfo(id)%fileobs)

	if( any( pinfo(id)%corlength<=0 ) ) goto 97
	if( any( pinfo(id)%sigobs<=0 ) ) goto 97
	if( any( pinfo(id)%sigback<=0 ) ) goto 97
	where( pinfo(id)%maxlength <= 0. )
	  pinfo(id)%maxlength = 3.*pinfo(id)%corlength
	end where
	if( pinfo(id)%filecoords == ' ' ) goto 96
	if( pinfo(id)%fileobs == ' ' ) goto 96

	call delete_section('assimil')

	return
   96	continue
	write(6,*) 'both filecoords and fileobs must be given'
	write(6,*) 'filecoords = ',trim(pinfo(id)%filecoords)
	write(6,*) 'fileobs    = ',trim(pinfo(id)%fileobs)
	stop 'error stop init_assimil_section: file names missing'
   97	continue
	write(6,*) 'some parameters have value 0 which is not allowed'
	write(6,*) 'corlength: ',pinfo(id)%corlength
	write(6,*) 'sigobs: ',pinfo(id)%sigobs
	write(6,*) 'sigback: ',pinfo(id)%sigback
	stop 'error stop init_assimil_section: error in parameters'
   98	continue
	write(6,*) 'ivar has already be specified: ',ivar
	write(6,*) 'can only specify one section with this ivar'
	stop 'error stop init_assimil_section: ivar not unique'
   99	continue
	write(6,*) 'nobs,ivar: ',nobs,ivar
	stop 'error stop init_assimil_section: nobs or ivar not set'
	end

!*************************************************************

