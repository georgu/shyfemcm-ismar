module mod_manage_obs

  use mod_para
  use mod_init_enkf
  use mod_obs_states

  implicit none

  ! Observations
  !
  type(files), allocatable, dimension(:), private ::  ofile

  integer :: n_0dlev,n_2dvel		! number of obs files of different type
  integer :: nobs_tot			! total number of good obs to be assimilated
  type(levels), allocatable :: o0dlev(:)
  type(currentf), allocatable :: o2dvel(:)

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

  write(*,*) 'reading observations...'
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
        read(25,*,err=95) ofile(n)
     end do
     close(25)

  !-------------------------------
  ! Read every file one time to find the number and the type of each obs
  !-------------------------------
  n_0dlev = 0
  n_2dvel = 0
  do n = 1,nfile

    select case (trim(ofile(n)%ty))

        case default
             error stop 'read_obs: Unknown file type'

        ! Level timeseries
        !
        case ('0DLEV')
             linit = .true.
             kinit = n_0dlev
             call read_level(linit,trim(ofile(n)%name),TEPS,&
                             ofile(n)%id,kinit,kend,nobs)
             if (kend > kinit) n_0dlev = n_0dlev + 1

        ! 2d current fem files
        !
        case ('2DVEL')
             n_2dvel = n_2dvel + 1

    end select

  end do


  !-------------------------------
  ! allocate obs vars
  !-------------------------------
  if (n_0dlev > 0) allocate(o0dlev(n_0dlev))
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
             call read_level(linit,trim(ofile(n)%name),TEPS,&
                             ofile(n)%id,kinit,kend,nobs)
             kinit = kend
             nobs_tot = nobs_tot + nobs
        case ('2DVEL')
             n_2dvel = n_2dvel + 1
             call read_2dvel(trim(ofile(n)%name),ofile(n)%id,n_2dvel,TEPS,nobs)
             nobs_tot = nobs_tot + 2 * nobs ! u and v components
    end select

  end do

  !-------------------------------
  ! create super-observations
  !-------------------------------
  call make_super_2dvel

  return

 95 write(*,*) 'read_obs error reading file: ',trim(ofile(n)%name)
    stop

  end subroutine read_obs

!******************************************************

  subroutine read_level(linit,filin,eps,id,kinit,kend,nobs)

  use iso8601
  implicit none

  logical,intent(in)           :: linit
  character(len=*),intent(in)  :: filin
  double precision, intent(in) :: eps
  integer, intent(in)	       :: id
  integer, intent(in)          :: kinit
  integer, intent(out)         :: kend
  integer, intent(out)         :: nobs
  integer ios
  character*80 dstring
  integer nvar
  integer ierr
  integer date, time
  double precision oatime
  real x, y, v, stdv
  integer k
  integer ostatus

  nobs = 0

  ! read info file for the obs file with coords and std
  !
  open(27,file=trim(filin)//'.info', status = 'old', form = 'formatted', iostat = ios)
  if (ios /= 0) error stop 'read_level: error opening info file'
  read(27,*) x,y,stdv
  close(27)

  ! read obs file
  !
  open(26,file=trim(filin), status = 'old', form = 'formatted', iostat = ios)
  if (ios /= 0) error stop 'read_level: error opening file'

  select case(linit)
  
     ! just look if there is a valid obs inside
     !
     case (.true.)
          k = kinit

 90       read(26,*,end=100) dstring

          call string2date(trim(dstring),date,time,ierr)
          if (ierr /= 0) error stop "read_level: error reading string"
          call dts_to_abs_time(date,time,oatime)

          ! Take only records with times near atime
          !
          if (abs(oatime - atime) < eps) then
             ostatus = 0
             ! check if the obs value is out of range
             !
             call check_obs('0DLEV',v,v,OFLAG,ostatus)
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

          call string2date(trim(dstring),date,time,ierr)
          if (ierr /= 0) error stop "read_level: error reading string"
          call dts_to_abs_time(date,time,oatime)

          ostatus = 4
          ! Take only records with times near atime
          !
          if (abs(oatime - atime) < eps) then
             ostatus = 0
             ! check if the obs value is out of range
             !
             call check_obs('0DLEV',v,v,OFLAG,ostatus)
             k = k + 1
             o0dlev(k)%t = oatime
             o0dlev(k)%x = x
             o0dlev(k)%y = y
             o0dlev(k)%val = v
             o0dlev(k)%std = stdv
             o0dlev(k)%status = ostatus
             o0dlev(k)%id = id

             nobs = nobs + 1

             goto 101	!take just 1 lev obs
          end if
          goto 91
 101      close(26)
          kend = k

  end select

  end subroutine read_level

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
  real*4,allocatable :: data(:,:),dataens(:,:)
  character(len=50) :: string
  integer ostatus
  logical bdata
  real x,y,uu,vv,ostd
  integer ix,iy
  double precision oatime

  nobs = 0

  ! read info file for the obs file with std
  !
  open(27,file=trim(filin)//'.info', status = 'old', form = 'formatted', iostat = ios)
  if (ios /= 0) error stop 'read_2dvel: error opening file'
  read(27,*) ostd
  close(27)

  ! open fem file
  !
  np = 0
  call fem_file_read_open(trim(filin),np,iformat,iunit)

  irec = 0
  bdata = .false.
  do
    irec = irec + 1

    ! Reads headers
    !
    call fem_file_read_params(iformat,iunit,tt,nvers,np,lmax,nvar,ntype,datetime,ierr)
    if (ierr < 0) exit

    call dts_convert_to_atime(datetime,tt,oatime)

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

    if (abs(oatime - atime) > eps) then
       do i=1,nvar
          call fem_file_skip_data(iformat,iunit,nvers,np,lmax,string,ierr)
          if (ierr /= 0) error stop 'read_2dvel: error reading file'
       end do
       cycle
    end if

    allocate(ilhkv(np),hd(np),data(nx,ny))
    allocate(o2dvel(nrec)%x(nx,ny),o2dvel(nrec)%y(nx,ny),&
             o2dvel(nrec)%u(nx,ny),o2dvel(nrec)%v(nx,ny),&
             o2dvel(nrec)%std(nx,ny),o2dvel(nrec)%status(nx,ny))

    do i = 1,nvar
       call fem_file_read_data(iformat,iunit,nvers,np,lmax,string,ilhkv,hd,nlvddi,&
                               data,ierr)

       select case (i)
              case default
                   error stop 'read_2dvel: too many variables'
              case (1)
                   o2dvel(nrec)%u = data
              case (2)
                   o2dvel(nrec)%v = data
       end select
    end do
    deallocate(ilhkv,hd,data)

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
       o2dvel(nrec)%status = ostatus

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
    o2dvel(nrec)%status = ostatus ! status of the obs

    ! exit after finding 1 field
    !
    exit

  end do
  close(iunit)

  if (.not.bdata) error stop 'read_2dvel: 2dvel file without valid data'

  end subroutine read_2dvel

!********************************************************

  subroutine check_obs(ty,v1,v2,flag,status)

  implicit none

  character(len=*), intent(in) :: ty
  real, intent(in) :: v1,v2
  real, intent(in) :: flag
  integer, intent(out) :: status

  real vmin,vmax,v

  status = 0

  if (v1 == flag .or. v2 == flag) then
     status = 4
     return
  end if

  if (trim(ty) == '0DLEV') then
     vmin = SSH_MIN
     vmax = SSH_MAX
     v = v1
  else if (trim(ty) == '2DVEL') then
     vmin = VEL_MIN
     vmax = VEL_MAX
     v = sqrt(v1**2 + v2**2)
  else
     write(*,*) 'Observation not implemented yet'
     error stop
  end if

  if (v < vmin .or. v > vmax) then
     status = 3
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
     stat(ii*jj) = o2dvel(n)%status(ii,jj)
     v1(ii*jj) = o2dvel(n)%u(ii,jj)
     v2(ii*jj) = o2dvel(n)%v(ii,jj)
    end do
    end do

    call superobs_horiz_el(nobs,x,y,stat,v1,v2)

    do ii = 1,nx
    do jj = 1,ny
     o2dvel(n)%status(ii,jj) = stat(ii*jj)
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
  subroutine check_obs_inn(ty,x,y,z,v1,v2,stdv,d1,d2,status)

  use levels
  use mod_ens_state

  implicit none

  character(len=*), intent(in) :: ty
  real, intent(in) :: x,y,z
  integer, intent(in) :: status
  real, intent(in) :: v1,v2
  real, intent(inout) :: stdv
  real, intent(out) :: d1,d2

  real vmod,vmod_ens(nrens)
  real v
  integer ie,ik,ne
  real h_1st_layer
  real ens_std,stdv_new
  real inn

  ! Exit if it is a bad obs, not to be assimilated
  !
  if (status /= 0) return

  ! Find the nearest node/element
  !
  call find_el_node(x,y,ie,ik)

  select case (trim(ty))

    case default
     write(*,*) 'Observation not implemented yet'
     error stop

    case ('0DLEV')

     vmod = Am%z(ik)
     vmod_ens = A(:)%z(ik)
     v = v1
     d1 = v - vmod
     d2 = v - vmod
     inn = d1

    case ('2DVEL')

     if (size(hlv) <= 1)&
     error stop 'check_obs_inn: a 3D sim is necessary to assimilate surface currents'

     vmod = sqrt(Am%u(1,ie)**2 + Am%v(1,ie)**2)
     vmod_ens = sqrt(A(:)%u(1,ie)**2 + A(:)%v(1,ie)**2)
     h_1st_layer = hlv(1) + Am%z(ik)

     v = sqrt(v1**2 + v2**2)
     d1 = v1*h_1st_layer - Am%u(1,ie)
     d2 = v2*h_1st_layer - Am%v(1,ie)
     !d1 = v1 - Am%u(1,ie)
     !d2 = v2 - Am%v(1,ie)
     inn = sqrt(d1**2 + d2**2)

  end select

  ! if lower than zero do not modify the std
  !
  if (KSTD <= 0) return

  ! Compute the ensemble spread (std)
  !
  ens_std = 0.0d0
  do ne = 1,nrens
     ens_std = ens_std + (vmod_ens(ne) - vmod)**2
  end do
  ens_std = sqrt(ens_std/float(nrens-1))

  ! Formula from: Sakov, 2012 (topaz)
  !
  if (abs(inn) > KSTD * ens_std) then
     stdv_new = sqrt( sqrt( (ens_std**2 + stdv**2)**2 +&
                 (1/KSTD * ens_std * inn)**2 ) - ens_std**2 )

     if (verbose)&
     write(*,'(a18,2x,a5,2f8.4,2x,2f8.4)') 'changing obs std:',&
             trim(ty),inn,ens_std,stdv,stdv_new

     stdv = stdv_new
  end if

  end subroutine check_obs_inn

end module mod_manage_obs
