!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
module mod_manage_obs

  use mod_para
  use mod_init_enkf
  use mod_obs_states

  implicit none

  ! Observations
  !
  type(files), allocatable, dimension(:), private ::  ofile

  integer :: islev,isvel,istemp,issalt
  integer :: n_0dlev,n_1dlev,n_2dlev    ! number of obs level (gauge,alt track,field)
  integer :: n_0dvel,n_2dvel		! number of obs vel (gauge,field)
  integer :: n_0dtemp,n_1dtemp,n_2dtemp	! number of obs temp (gauge,profile,field)
  integer :: n_0dsalt,n_1dsalt,n_2dsalt	! number of obs salt (gauge,profile,field)
  integer :: nobs_tot			! total number of good obs to be assimilated
  type(scalar_0d), allocatable :: o0dlev(:),o0dtemp(:),o0dsalt(:)
  type(vector_2d), allocatable :: o2dvel(:)

  ! status of an observation:
  ! 0 = normal obs (assimilated)
  ! 1 = super-observation (assimilated)
  ! 2 = observation merged into a super-observation (not-assimilated)
  ! 3 = observation out of range (not-assimilated)
  ! 4 = observation with a flag value (not-assimilated)

contains

!********************************************************

  subroutine read_obs

  implicit none

  integer ios
  character(len=80) :: line

  integer n,nfile
  integer kinit,kend
  logical linit
  integer nobs
  double precision tobs
  real xobs,yobs,zobs,vobs,stdobs,rho
  integer statobs

  ! flag for the presence of obs
  islev = 0
  isvel = 0
  istemp = 0
  issalt = 0

  ! dimension of the obs types
  n_0dlev = 0
  n_1dlev = 0
  n_2dlev = 0
  n_0dtemp = 0
  n_1dtemp = 0
  n_2dtemp = 0
  n_0dsalt = 0
  n_1dsalt = 0
  n_2dsalt = 0
  n_0dvel = 0
  n_2dvel = 0

  write(*,*) 'Reading observations...'
  !-------------------------------
  ! Read a list of obs files
  !-------------------------------
     n = 1
     open(25,file = obsfile, status = 'old', form = 'formatted', iostat = ios)
     if (ios /= 0) error stop 'read_obs: error opening file'
 88   read(25,*,end=98) line
      n = n + 1
     goto 88
 98  continue
     nfile = n - 1
     rewind(unit=25)

     allocate(ofile(nfile))
     do n = 1,nfile
        read(25,*,err=95) ofile(n)%ty,ofile(n)%name
     end do
     close(25)

  !-------------------------------
  ! Read every file one time to find the number and the type of each obs
  !-------------------------------
  do n = 1,nfile

    select case (trim(ofile(n)%ty))

        case default
             error stop 'read_obs: Unknown file type'

        ! Level timeseries
        !
        case ('0DLEV')
             linit = .true.
             kinit = n_0dlev
             call read_scalar_0d('0DLEV',linit,trim(ofile(n)%name),TEPS,&
                  kinit,kend,tobs,xobs,yobs,zobs,&
                  vobs,stdobs,statobs,rho)
             if (kend > kinit) then
		n_0dlev = n_0dlev + 1
		islev = 1
	     endif

        ! Temperature timeseries
        !
        case ('0DTEM')
             linit = .true.
             kinit = n_0dtemp
             call read_scalar_0d('0DTEM',linit,trim(ofile(n)%name),TEPS,&
                  kinit,kend,tobs,xobs,yobs,zobs,&
                  vobs,stdobs,statobs,rho)
             if (kend > kinit) then
		n_0dtemp = n_0dtemp + 1
                istemp = 1
	     endif

        ! Salinity timeseries
        !
        case ('0DSAL')
             linit = .true.
             kinit = n_0dsalt
             call read_scalar_0d('0DSAL',linit,trim(ofile(n)%name),TEPS,&
                  kinit,kend,tobs,xobs,yobs,zobs,&
                  vobs,stdobs,statobs,rho)
             if (kend > kinit) then
		n_0dsalt = n_0dsalt + 1
                issalt = 1
	     endif

        ! 2d current fem files
        !
        case ('2DVEL')
             n_2dvel = n_2dvel + 1
             isvel = 1

    end select

  end do

  ! check that all the files have the same variable type. It's better to assimilate
  ! only one type every time. Use different times if you have more types (e.g., level, salinity)
  ! Otherwise the matrices should be normalised (not still implemented).
  if ((islev + isvel + istemp + issalt) > 1) then
     write(*,*) 'islev ',islev
     write(*,*) 'isvel ',isvel
     write(*,*) 'istemp ',istemp
     write(*,*) 'issalt ',issalt
     error stop 'Different type of observations. Assimilate them at different times'
  end if 

  !-------------------------------
  ! allocate obs vars
  !-------------------------------
  if (n_0dlev > 0) allocate(o0dlev(n_0dlev))
  if (n_0dtemp > 0) allocate(o0dtemp(n_0dtemp))
  if (n_0dsalt > 0) allocate(o0dsalt(n_0dsalt))
  if (n_2dvel > 0) allocate(o2dvel(n_2dvel))

  !-------------------------------
  ! read a second time and store
  !-------------------------------
  n_2dvel = 0
  kinit = 0
  nobs_tot = 0
  do n = 1,nfile

    select case (trim(ofile(n)%ty))
        case default
             error stop 'read_obs: Unknown file type'
        case ('0DLEV')
             linit = .false.
             call read_scalar_0d('0DLEV',linit,trim(ofile(n)%name),TEPS,&
                  kinit,kend,tobs,xobs,yobs,zobs,&
                  vobs,stdobs,statobs,rho)
             if (kend > kinit) then 
		if (verbose) write(*,*) 'Station n. ',n
                o0dlev(kend)%t = tobs
                o0dlev(kend)%x = xobs
                o0dlev(kend)%y = yobs
                o0dlev(kend)%z = zobs
                o0dlev(kend)%val = vobs
                o0dlev(kend)%std = stdobs
                o0dlev(kend)%stat = statobs
                o0dlev(kend)%id = n
                o0dlev(kend)%rhol = rho

                nobs_tot = nobs_tot + 1
             end if
             kinit = kend
        case ('0DTEM')
             linit = .false.
             call read_scalar_0d('0DTEM',linit,trim(ofile(n)%name),TEPS,&
                  kinit,kend,tobs,xobs,yobs,zobs,&
                  vobs,stdobs,statobs,rho)
             if (kend > kinit) then 
		if (verbose) write(*,*) 'Station n. ',n
                o0dtemp(kend)%t = tobs
                o0dtemp(kend)%x = xobs
                o0dtemp(kend)%y = yobs
                o0dtemp(kend)%z = zobs
                o0dtemp(kend)%val = vobs
                o0dtemp(kend)%std = stdobs
                o0dtemp(kend)%stat = statobs
                o0dtemp(kend)%id = n
                o0dtemp(kend)%rhol = rho

                nobs_tot = nobs_tot + 1
             end if
             kinit = kend

        case ('0DSAL')
             linit = .false.
             call read_scalar_0d('0DSAL',linit,trim(ofile(n)%name),TEPS,&
                  kinit,kend,tobs,xobs,yobs,zobs,&
                  vobs,stdobs,statobs,rho)
             if (kend > kinit) then 
		if (verbose) write(*,*) 'Station n. ',n
                o0dsalt(kend)%t = tobs
                o0dsalt(kend)%x = xobs
                o0dsalt(kend)%y = yobs
                o0dsalt(kend)%z = zobs
                o0dsalt(kend)%val = vobs
                o0dsalt(kend)%std = stdobs
                o0dsalt(kend)%stat = statobs
                o0dsalt(kend)%id = n
                o0dsalt(kend)%rhol = rho

                nobs_tot = nobs_tot + 1
             end if
             kinit = kend

        case ('2DVEL')
             n_2dvel = n_2dvel + 1
             call read_2dvel(trim(ofile(n)%name),n,n_2dvel,TEPS,nobs)
             nobs_tot = nobs_tot + 2 * nobs ! u and v components
    end select

  end do

  !-------------------------------
  ! create super-observations
  !-------------------------------
  call make_super_2dvel
 
  if (nobs_tot < 1) error stop 'No valid observations, stopping.'

  return

 95 write(*,*) 'read_obs error reading file: ',trim(ofile(n)%name)
    stop

  end subroutine read_obs

!******************************************************

  subroutine read_scalar_0d(olabel,linit,filin,eps,kinit,kend,atime_obs,xv,yv,zv,vv,stdvv, &
		  ostatusv,rho)

  use iso8601
  implicit none

  character(len=*),intent(in)  :: olabel
  logical,intent(in)           :: linit
  character(len=*),intent(in)  :: filin
  double precision, intent(in) :: eps
  integer, intent(in)          :: kinit
  integer, intent(out)         :: kend
  double precision, intent(out):: atime_obs
  real, intent(out)            :: xv,yv,zv,vv,stdvv,rho
  integer, intent(out)         :: ostatusv
  integer ios
  real x,y,z,v,stdv
  integer ostatus
  character*80 dstring
  integer nvar
  integer ierr
  integer date, time
  integer k

  xv = -999.
  yv = -999.
  zv = -999.
  vv = -999.
  stdvv = -999.
  ostatusv = -999

  ! read info file for the obs file with coords and std
  !
  open(27,file=trim(filin)//'.info', status = 'old', form = 'formatted', iostat = ios)
  if (ios /= 0) error stop 'read_scalar_0d: error opening info file'
  read(27,*) x,y,z,stdv,rho
  close(27)

  ! read obs file
  !
  open(26,file=trim(filin), status = 'old', form = 'formatted', iostat = ios)
  if (ios /= 0) error stop 'read_scalar_0d: error opening file'

  select case(linit)
  
     ! just look if there is a valid obs inside
     !
     case (.true.)
          k = kinit

 90       read(26,*,end=100) dstring,v
          if (isnan(v)) v = OFLAG

          call string2date(trim(dstring),date,time,ierr)
          if (ierr /= 0) error stop "read_scalar_0d: error reading string"
          call dts_to_abs_time(date,time,atime_obs)

          ! Take only records with times near atime_an
          !
          if (abs(atime_obs - atime_an) < eps) then
             ostatus = 0
             ! check if the obs value is out of range
             !
             call check_obs(olabel,v,v,OFLAG,ostatus)
             k = k + 1
          end if
          goto 90
 100      close(26)
          kend = k

     ! store the obs
     !
     case (.false.)
          k = kinit

 91       read(26,*,end=101) dstring,v
          if (isnan(v)) v = OFLAG

          call string2date(trim(dstring),date,time,ierr)
          if (ierr /= 0) error stop "read_scalar_0d: error reading string"
          call dts_to_abs_time(date,time,atime_obs)

          ostatus = 4
          ! Take only records with times near atime_an
          !
          if (abs(atime_obs - atime_an) < eps) then
             ostatus = 0
             ! check if the obs value is out of range
             !
             call check_obs(olabel,v,v,OFLAG,ostatus)
             k = k + 1
             xv = x
             yv = y
             zv = z
             vv = v
             stdvv = stdv
             ostatusv = ostatus

             goto 101	!take just 1 obs
          end if
          goto 91
 101      close(26)
          kend = k

  end select

  end subroutine read_scalar_0d

!********************************************************

  subroutine read_2dvel(filin,fid,nrec,eps,nobs)

  implicit none

  character(len=*),intent(in)  :: filin
  integer,intent(in)           :: fid
  integer, intent(in)          :: nrec
  double precision,intent(in)  :: eps
  integer, intent(out)          :: nobs

  integer ios
  integer np,iformat,iunit
  integer irec,i,ii,jj
  integer nvers           !version of file format
  integer lmax            !vertical values
  integer nvar            !number of variables to write
  integer ntype           !type of information contained
  double precision tt
  integer datetime(2)     !date and time information
  integer ierr            !return error code
  integer nlvddi
  real*4,allocatable :: hhlv(:)   !vertical structure
  real*4 regpar(7)               !regular array params
  integer nx,ny
  real flag,dx,dy,x0,y0
  integer,allocatable :: ilhkv(:)
  real*4,allocatable :: hd(:)
  real*4,allocatable :: idata(:,:,:)
  character(len=50) :: string
  integer ostatus
  logical bdata
  real x,y,uu,vv,ostd,rho
  integer ix,iy
  double precision atime_obs

  nobs = 0

  ! read info file for the obs file with std
  !
  open(27,file=trim(filin)//'.info', status = 'old', form = 'formatted', iostat = ios)
  if (ios /= 0) error stop 'read_2dvel: error opening file'
  read(27,*) ostd,rho
  close(27)

  ! open fem file
  !
  np = 0
  write(*,*) 'Opening velocity fem-file: ',trim(filin)
  call fem_file_read_open(trim(filin),np,iformat,iunit)

  irec = 0
  bdata = .false.
  do
    irec = irec + 1

    ! Reads headers
    !
    call fem_file_read_params(iformat,iunit,tt,nvers,np,lmax,nvar,ntype,datetime,ierr)
    if (ierr < 0) exit

    call dts_convert_to_atime(datetime,tt,atime_obs)

    allocate(hhlv(lmax))
    nlvddi = lmax
    call fem_file_read_2header(iformat,iunit,ntype,lmax,hhlv,regpar,ierr)

    nx = nint(regpar(1))
    ny = nint(regpar(2))
    x0 = regpar(3)
    y0 = regpar(4)
    dx = regpar(5)
    dy = regpar(6)
    flag = regpar(7)

    if (flag /= OFLAG) error stop 'read_2dvel: bad flag'

    if (abs(atime_obs - atime_an) > eps) then
       do i=1,nvar
          call fem_file_skip_data(iformat,iunit,nvers,np,lmax,string,ierr)
          if (ierr /= 0) error stop 'read_2dvel: error reading file'
       end do
       cycle
    end if

    allocate(ilhkv(np),hd(np),idata(1,nx,ny))
    allocate(o2dvel(nrec)%x(nx,ny),o2dvel(nrec)%y(nx,ny),&
             o2dvel(nrec)%u(nx,ny),o2dvel(nrec)%v(nx,ny),&
             o2dvel(nrec)%std(nx,ny),o2dvel(nrec)%stat(nx,ny))

    do i = 1,nvar
       idata = flag
       call fem_file_read_data(iformat,iunit,nvers,np,lmax,string,ilhkv,hd,nlvddi,&
                               idata,ierr)

       select case (i)
              case default
                   error stop 'read_2dvel: too many variables'
              case (1)
                   o2dvel(nrec)%u = idata(1,:,:)
              case (2)
                   o2dvel(nrec)%v = idata(1,:,:)
       end select
    end do
    deallocate(ilhkv,hd,idata)

    ! find and assign the coords
    !
    do ii = 1,nx
    do jj = 1,ny
       o2dvel(nrec)%x(ii,jj) = x0 + dx * (ii-1)	! x coords
       o2dvel(nrec)%y(ii,jj) = y0 + dy * (jj-1)	! y coords
    end do
    end do
    o2dvel(nrec)%z = 0.		! depth

    ! check observation values
    !
    do ix = 1,nx
    do iy = 1,ny
       uu = o2dvel(nrec)%u(ix,iy)
       vv = o2dvel(nrec)%v(ix,iy)

       call check_obs('2DVEL',uu,vv,flag,ostatus)
       o2dvel(nrec)%stat = ostatus

       if (ostatus == 0) then
          bdata = .true. !at least 1 rec
          nobs = nobs + 1
       end if

    end do
    end do

    ! last variables to assign
    !
    o2dvel(nrec)%nx = nx	! n of x points
    o2dvel(nrec)%ny = ny	! n od y points
    o2dvel(nrec)%std = ostd	! standard dev
    o2dvel(nrec)%id = fid	! id of the original file
    o2dvel(nrec)%stat = ostatus ! status of the obs

    ! exit after finding 1 field
    !
    exit

  end do
  close(iunit)

  if (.not.bdata) error stop 'read_2dvel: 2dvel file without valid data'

  end subroutine read_2dvel

!********************************************************

  subroutine check_obs(ty,v1,v2,flag,stat)

  implicit none

  character(len=*), intent(in) :: ty
  real, intent(in) :: v1,v2
  real, intent(in) :: flag
  integer, intent(out) :: stat

  real vmin,vmax

  stat = 0

  if (nint(v1-flag) == 0 .or. nint(v2-flag) == 0) then
     stat = 4
     return
  end if

  if (trim(ty) == '0DLEV') then
     vmin = SSH_MIN
     vmax = SSH_MAX
  else if (trim(ty) == '0DTEM') then
     vmin = TEM_MIN
     vmax = TEM_MAX
  else if (trim(ty) == '0DSAL') then
     vmin = SAL_MIN
     vmax = SAL_MAX
  else if (trim(ty) == '2DVEL') then
     vmin = VEL_MIN
     vmax = VEL_MAX
  else
     write(*,*) 'Observation not implemented yet'
     error stop
  end if

  if (v1 <= vmin .or. v1 >= vmax) then
     stat = 3
  end if
  if (v2 <= vmin .or. v2 >= vmax) then
     stat = 3
  end if

  end subroutine check_obs

!********************************************************

  subroutine make_super_2dvel

  implicit none

  integer nobs,nx,ny
  real,allocatable :: x(:),y(:),v1(:),v2(:)
  integer,allocatable :: stat(:)
  integer n,ii,jj


  if (n_2dvel < 1) return

  do n = 1,n_2dvel
     nx = o2dvel(n)%nx
     ny = o2dvel(n)%ny
     nobs = nx * ny

    allocate(x(nobs),y(nobs),v1(nobs),v2(nobs),stat(nobs))

    do ii = 1,nx
    do jj = 1,ny
     x(ii*jj) = o2dvel(n)%x(ii,jj)
     y(ii*jj) = o2dvel(n)%y(ii,jj)
     stat(ii*jj) = o2dvel(n)%stat(ii,jj)
     v1(ii*jj) = o2dvel(n)%u(ii,jj)
     v2(ii*jj) = o2dvel(n)%v(ii,jj)
    end do
    end do

    call superobs_horiz_el(nobs,x,y,stat,v1,v2)

    do ii = 1,nx
    do jj = 1,ny
     o2dvel(n)%stat(ii,jj) = stat(ii*jj)
     o2dvel(n)%u(ii,jj) = v1(ii*jj)
     o2dvel(n)%v(ii,jj) = v2(ii*jj)
    end do
    end do

    deallocate(x,y,v1,v2,stat)

  end do

  end subroutine make_super_2dvel

!********************************************************

  ! Compute the innovation vectors and change the std if
  ! it is too low
  !
  subroutine check_spread(d,stdv,mval,mvalm)

  use levels
  use mod_ens_state

  implicit none

  real, intent(in) :: d
  real, intent(in) :: mvalm,mval(nrens)
  real, intent(inout) :: stdv
  integer ne

  real ens_std,stdv_new

  ! if lower than zero do not modify the std
  !
  if (KSTD <= 0) return

  ! Compute the ensemble spread (std)
  !
  ens_std = 0.0d0
  do ne = 1,nrens
     ens_std = ens_std + (mval(ne) - mvalm)**2
  end do
  ens_std = sqrt(ens_std/float(nrens-1))

  ! Formula from: Sakov, 2012 (topaz)
  !
  if (abs(d) > KSTD * ens_std) then
     stdv_new = sqrt( sqrt( (ens_std**2 + stdv**2)**2 +&
                 (1/KSTD * ens_std * d)**2 ) - ens_std**2 )

     if (verbose)&
     write(*,'(a18,2f8.4)') ' changing obs std ',&
             stdv,stdv_new

     stdv = stdv_new
  end if

  end subroutine check_spread


!********************************************************

 subroutine check_spread_speed(d1,d2,stdv,um,vm,umm,vmm)

  use mod_ens_state

  implicit none

  real, intent(in) :: d1,d2
  real, intent(inout) :: stdv
  real, intent(in) :: um(nrens),vm(nrens),umm,vmm
  real cs(nrens),csm
  real inn
  integer ne

  real ens_std,stdv_new

  ! if lower than zero do not modify the std
  !
  if (KSTD <= 0) return

  cs  = sqrt(um**2 + vm**2)
  csm  = sqrt(umm**2 + vmm**2)
  ! Compute the ensemble spread (std)
  !
  ens_std = 0.0d0
  do ne = 1,nrens
     ens_std = ens_std + (cs(ne) - csm)**2
  end do
  ens_std = sqrt(ens_std/float(nrens-1))
  inn = sqrt(d1**2 + d2**2)

  ! Formula from: Sakov, 2012 (topaz)
  !
  if (abs(inn) > KSTD * ens_std) then
     stdv_new = sqrt( sqrt( (ens_std**2 + stdv**2)**2 +&
                 (1/KSTD * ens_std * inn)**2 ) - ens_std**2 )

     if (verbose)&
     write(*,'(a18,2f8.4)') ' changing obs std ',&
             stdv,stdv_new

     stdv = stdv_new
  end if

 end subroutine check_spread_speed


end module mod_manage_obs
