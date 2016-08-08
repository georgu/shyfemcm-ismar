  module mod_enKF

  implicit none

  character (len=80), save :: basfile, obsfile, rstfile

  integer, save :: nens, na, sdim

  ! Observations
  double precision, save :: tobs
  integer, save :: nobs
  character(len=6), save, allocatable :: tyobs(:)
  real, save, allocatable :: xobs(:), yobs(:), zobs(:)
  real, save, allocatable :: vobs(:), stdvobs(:)	!Values and standard deviations

  real, save, allocatable :: A(:,:)		! Matrix holding the states
  real, save, allocatable :: Am(:)		! Array with the mean state
  real, save, allocatable :: Am2d(:,:)		! Like Am with one dummy dim for the save routines

  real, save, allocatable :: D(:,:)		! matrix holding perturbed measurments
  real, save, allocatable :: E(:,:)		! matrix holding perturbations (mode=?3)

  real, save, allocatable :: R(:,:)		! Obs error cov matrix
  real, save, allocatable :: S(:,:)		! matrix holding HA`
  real, save, allocatable :: innov(:)		! innovation vector holding d-H*mean(A)


  contains


!********************************************************
!********************************************************
!********************************************************

  subroutine read_info

  implicit none

  integer n

  open(20, file='analysis.info', status='old')

  read(20,*) nens	! number of ens members
  read(20,*) na		! analysis step
  read(20,*) basfile	! name of bas file (no extension)
  read(20,*) tobs	! current time of the observations
  read(20,*) obsfile	! name of obs file list

  close(20)

  end subroutine read_info

!********************************************************

  subroutine push_state(ncol)

! Push a model state into the matrix A in the column ncol.
! TODO: add wet and dry and conz.
  use basin

  use mod_conz
  use mod_geom_dynamic
  use mod_ts
  use mod_hydro_vel
  use mod_hydro
  use levels, only : nlvdi,nlv
  use mod_restart

  implicit none

  integer, intent(in) :: ncol

  integer dimbr,dimbc,dimcon
  integer kend,kinit

  integer nl,ie

  ! Compute dimensions
  dimbr = nkn + 2 * ( nel * nlv )	! barotropic
  dimbc = 2 * ( nkn * nlv )		! temp and salt
  dimcon = iconz * ( nkn * nlv )	! tracers

  sdim = dimbr + dimcon
  if( ibarcl_rst.gt.0 ) sdim = sdim + dimbc

  if( .not. allocated(A) ) allocate(A(sdim,nens))
  if( .not. allocated(Am) ) allocate(Am(sdim))
 
  ! Fill with znv
  kinit = 1
  call pushA(kinit,nkn,sdim,nens,ncol,znv,kend,A)

!  ! Fill with zenv
!  do ie = 1,nel
!     kinit = kend + 1
!     call pushA(kinit,3,sdim,nens,ncol,zenv(:,ie),kend,A)
!  end do

  ! Fill with transports
  do nl = 1,nlv
     kinit = kend + 1
     call pushA(kinit,nel,sdim,nens,ncol,utlnv(nl,:),kend,A)
     kinit = kend + 1
     call pushA(kinit,nel,sdim,nens,ncol,vtlnv(nl,:),kend,A)
  end do

  ! Fill with temperature and salinity
  if( ibarcl_rst.gt.0 ) then
    do nl = 1,nlv
       kinit = kend + 1
       call pushA(kinit,nkn,sdim,nens,ncol,saltv(nl,:),kend,A)
       kinit = kend + 1
       call pushA(kinit,nkn,sdim,nens,ncol,tempv(nl,:),kend,A)
    end do
  end if
 
  if( iconz.ne.0 ) stop 'Concentration not still implemented.'
  !if( iconz.ne.0 ) TODO

  if( kend.ne.sdim ) then
    write(*,*) 'Wrong dimensions: ',kend,sdim
    stop 'push_state: Error'
  end if

  end subroutine push_state

!********************************************************

  subroutine pull_state(ncol)

! Pulls a model state stored in the column ncol of the matrix A and
! writes it in the standard shyfem variables, ready to be written in
! a restart file.

  use basin

  use mod_conz
  use mod_geom_dynamic
  use mod_ts
  use mod_hydro_vel
  use mod_hydro
  use levels, only : nlvdi,nlv
  use mod_restart

  implicit none
  integer, intent(in) :: ncol

  integer kinit,kend
  integer k
  integer nl,ie,ii

  kinit = 1
  call pullA(kinit,nkn,sdim,nens,ncol,znv,kend,A)

  ! Compute zenv (only without dry elements) TODO
  do ie = 1,nel
     do ii = 1,3
        k = nen3v(ii,ie)
        zenv(ii,ie) = znv(k)
     end do
  end do

!  do ie = 1,nel
!     kinit = kend + 1
!     call pullA(kinit,3,sdim,nens,ncol,zenv(:,ie),kend,A)
!  end do

  do nl = 1,nlv
     kinit = kend + 1
     call pullA(kinit,nel,sdim,nens,ncol,utlnv(nl,:),kend,A)
     kinit = kend + 1
     call pullA(kinit,nel,sdim,nens,ncol,vtlnv(nl,:),kend,A)
  end do

  ! Fill with temperature and salinity
  if( ibarcl_rst.gt.0 ) then
    do nl = 1,nlv
       kinit = kend + 1
       call pullA(kinit,nkn,sdim,nens,ncol,saltv(nl,:),kend,A)
       kinit = kend + 1
       call pullA(kinit,nkn,sdim,nens,ncol,tempv(nl,:),kend,A)
    end do
  end if
 
  !if( iconz.eq. 1 ) TODO

  if( kend.ne.sdim ) then
    write(*,*) 'Wrong dimensions: ',kend,sdim
    stop 'pull_state: Error'
  end if

  end subroutine pull_state

!********************************************************

  subroutine pull_av_state
  use basin

  use mod_conz
  use mod_geom_dynamic
  use mod_ts
  use mod_hydro_vel
  use mod_hydro
  use levels, only : nlvdi,nlv
  use mod_restart
  implicit none
  integer kinit,kend
  integer k
  integer nl,ie,ii
  integer :: ncol = 1

  if( .not. allocated(Am2d)) allocate(Am2d(sdim,1))

  Am2d(:,1) = Am

  kinit = 1
  call pullA(kinit,nkn,sdim,1,ncol,znv,kend,Am2d)

  ! Compute zenv (only without dry elements) TODO
  do ie = 1,nel
     do ii = 1,3
        k = nen3v(ii,ie)
        zenv(ii,ie) = znv(k)
     end do
  end do

!  do ie = 1,nel
!     kinit = kend + 1
!     call pullA(kinit,3,sdim,1,ncol,zenv(:,ie),kend,Am2d)
!  end do

  do nl = 1,nlv
     kinit = kend + 1
     call pullA(kinit,nel,sdim,1,ncol,utlnv(nl,:),kend,Am2d)
     kinit = kend + 1
     call pullA(kinit,nel,sdim,1,ncol,vtlnv(nl,:),kend,Am2d)
  end do

  ! Fill with temperature and salinity
  if( ibarcl_rst.gt.0 ) then
    do nl = 1,nlv
       kinit = kend + 1
       call pullA(kinit,nkn,sdim,1,ncol,saltv(nl,:),kend,Am2d)
       kinit = kend + 1
       call pullA(kinit,nkn,sdim,1,ncol,tempv(nl,:),kend,Am2d)
    end do
  end if
 
  !if( iconz.eq. 1 ) TODO

  if( kend.ne.sdim ) then
    write(*,*) 'Wrong dimensions: ',kend,sdim
    stop 'pull_state: Error'
  end if

  deallocate(Am2d)

  end subroutine pull_av_state

!********************************************************

  subroutine read_obs
  implicit none

  integer ios
  character(len=80) :: line
  character(len=80), allocatable :: ofile(:)
  integer, allocatable :: nrec(:)

  integer :: nmax = 10000
  character(len=6), allocatable :: tyobs_a(:)
  real, allocatable :: xobs_a(:), yobs_a(:), zobs_a(:)
  real, allocatable :: vobs_a(:), stdvobs_a(:)
  character(len=6) :: ty
  real x, y, z, v, stdv

  double precision :: eps = 300.	! 300 seconds
  double precision tt
  integer n,nfile,k

     write(*,*) 'Observation file list: ',trim(obsfile)

!    Reads the obs list
     n = 1
     open(25,file = obsfile, status = 'old', form = 'formatted', iostat = ios)
     if( ios.ne.0 ) stop 'read_obs: error opening file'
 88   read(25,*,end=98) line
      n = n + 1
     goto 88
 98  continue
     nfile = n - 1
     rewind(unit=25)

     allocate(ofile(nfile),nrec(nfile))

     do n = 1,nfile
        read(25,*,err=95) ofile(n)
     end do
     close(25)

! Reads every file
  allocate(tyobs_a(nmax), xobs_a(nmax), yobs_a(nmax), zobs_a(nmax),  &
           vobs_a(nmax), stdvobs_a(nmax))

  k = 0
  do n = 1,nfile

     open(26,file=ofile(n), status = 'old', form = 'formatted', iostat = ios)
     if( ios.ne.0 ) stop 'read_obs: error opening file'

 89   read(26,*,end=99) tt, ty, x, y, z, v, stdv

      ! Takes only records with times near tobs
      if( abs(tt - tobs) .lt. eps ) then
        if( trim(ty).ne.'level' ) stop 'Observation type not still implemented.'
        !write(*,*) 'Observation found: ',tt,trim(ty),x,y,z,v,stdv
        k = k + 1
        tyobs_a(k) = ty
        xobs_a(k) = x
        yobs_a(k) = y
        zobs_a(k) = z
        vobs_a(k) = v
        stdvobs_a(k) = stdv
      end if

     goto 89

 99  close(26)

  end do
  ! total numer of records
  nobs = k


  ! allocates the global vars and deallocate local ones
  allocate(tyobs(nobs), xobs(nobs), yobs(nobs), zobs(nobs),  &
           vobs(nobs), stdvobs(nobs))
  tyobs = tyobs_a(1:nobs)
  xobs = xobs_a(1:nobs)
  yobs = yobs_a(1:nobs)
  zobs = zobs_a(1:nobs)
  vobs = vobs_a(1:nobs)
  stdvobs = stdvobs_a(1:nobs)
  deallocate(tyobs_a, xobs_a, yobs_a, zobs_a, vobs_a, stdvobs_a)

  return

 95 write(*,*) 'read_obs error reading file: ',trim(ofile(n))
    stop

  end subroutine read_obs

!********************************************************

  subroutine make_D_E_R

! R (only used if mode=?1 or ?2) (no for low-rank sq root)

  implicit none
  real rand_v(nens)
  integer n,ne

  allocate(D(nobs,nens),E(nobs,nens),R(nobs,nobs))

  R = 0.	!Observations are indipendent
  do n = 1,nobs
     if( trim(tyobs(n)).ne.'level' ) stop 'Observation type not still implemented.'
     ! Makes a random vector
     call random2(rand_v,nens)
     do ne = 1,nens
        rand_v(ne) = (stdvobs(n) * rand_v(ne)) + vobs(n)
     end do
     ! TODO: add more observation types and different perturbation methods

     D(n,:) = rand_v
     !E(n,:) = rand_v - vobs(n)
     E(n,:) = rand_v - sum(rand_v)/nens		! E must have mean 0

     R(n,n) = stdvobs(n)**2
  end do

  end subroutine make_D_E_R

!********************************************************

  subroutine make_S_innov

  use regular
  implicit none

  integer iel, ikn
  integer iA
  integer n

  allocate (S(nobs,nens),innov(nobs))
  
  do n = 1,nobs

     ! Finds the nearest element
     call find_element(xobs(n),yobs(n),iel)

     if( trim(tyobs(n)).eq.'level' ) then
        ! Finds the nearest node
        call find_node(xobs(n),yobs(n),iel,ikn)
        iA = ikn
     else
        stop 'Observation type not still implemented.'
     end if

     S(n,:) = A(iA,:) - Am(iA)
     innov(n) = vobs(n) - Am(iA)

     write(*,*) 'Observation: ',n,trim(tyobs(n)),vobs(n),Am(iA),innov(n)
     
  end do

  end subroutine make_S_innov

  end module mod_enKF

!********************************************************
!********************************************************
!********************************************************


!********************************************************

  subroutine read_basin(basfile)
 
  use basin

  implicit none
  character(len=*),intent(in) :: basfile

  open(21, file=basfile, status='old', form='unformatted')
  call basin_read_by_unit(21)
  close(21)

  end subroutine read_basin

!********************************************************

  subroutine rst_read(ne,nan,tt)

  use basin
  use mod_conz
  use mod_geom_dynamic
  use mod_ts
  use mod_hydro_vel
  use mod_hydro
  use levels, only : nlvdi,nlv
  use mod_restart

  implicit none

  integer, intent(in) :: ne, nan
  double precision, intent(in) :: tt

  character(len=3) :: nrel,nal
  integer io
  integer it,nvers,nrec,iflag,ierr
  double precision atime
  integer date,time
  character(len=16) rstname

  integer, save :: nkn0,nel0,nlv0
  integer, save :: icall = 0

  call num2str(ne,nrel)
  call num2str(nan,nal)

  rstname='an'//nal//'_'//'en'//nrel//'b.rst'

!---- first call
  if( icall.eq.0 ) then

    ! checks dimensions and time
    open(25,file=rstname,status='old',form='unformatted',iostat=io)
    if( io.ne.0 ) stop 'rst_read: Error opening file'
    call rst_skip_record(25,atime,it,nvers,nrec,nkn0,nel0,nlv0,iflag,ierr)
    close(25)

    nkn = nkn0
    nel = nel0
    nlv = nlv0	!To read a restart file this var needs to be set

    ! init shyfem variables
    call mod_geom_dynamic_init(nkn0,nel0)
    call mod_hydro_init(nkn0,nel0,nlv0)
    call mod_hydro_vel_init(nkn0,nel0,nlv0)
    call mod_ts_init(nkn0,nlv0)  
    if( iconz_rst.gt.0 ) call mod_conz_init(iconz_rst,nkn0,nlv0)

    icall = icall + 1

  end if
!---- end first call

  open(24,file=trim(rstname),status='old',form='unformatted',iostat=io)
  if( io.ne.0 ) stop 'rst_read: Error opening file'

  call rst_read_record(atime,it,24,ierr)

  close(24)

  if(( nkn.ne.nkn0 ) .or. ( nel.ne.nel0 ) .or. ( nlv.ne.nlv0 )) then
       write(*,*) 'Different dimensions in the restart files'
       write(*,*) nkn0,nel0,nlv0
       write(*,*) nkn,nel,nlv
       stop
  end if

  if( it.ne.nint(tt) ) stop 'Error in rst file time'

  end subroutine rst_read

!********************************************************

  subroutine rst_write(ne,nan,tt)

  use mod_restart

  implicit none

  integer, intent(in) :: ne, nan
  double precision, intent(in) :: tt

  character(len=3) :: nrel,nal
  character(len=1) :: stype
  integer it
  character(len=16) rstname
  double precision ddate,dtime

  stype = 'a'

  if( ne.eq.-1 ) then
    write(*,*) 'Writing average background state...'
    nrel='avr'
    stype='b'
  elseif( ne.eq.-2 ) then
    write(*,*) 'Writing average analysis state...'
    nrel='avr'
  else
    call num2str(ne,nrel)
  end if

  call num2str(nan,nal)

  rstname = 'an'//nal//'_'//'en'//nrel//stype//'.rst'

  ! adds parameters
  call addpar('ibarcl',float(ibarcl_rst))
  call addpar('iconz',float(iconz_rst))
  call addpar('ibfm',0.)
  ddate = date_rst
  dtime = time_rst
  call daddpar('date',ddate)
  call daddpar('time',dtime)

  open(34,file=rstname,form='unformatted')
  it = nint(tt)
  call rst_write_record(it,34)
  close(34)

  end subroutine rst_write

!********************************************************

  subroutine pushA(kinit,dimv,dimA1,dimA2,ncol,v,kend,A)

  implicit none

  integer, intent(in) :: kinit,dimv,dimA1,dimA2,ncol
  real, intent(in) :: v(dimv)
  integer, intent(out) :: kend
  real, intent(inout) :: A(dimA1,dimA2)

  kend = kinit + dimv - 1
  if( kend.gt.dimA1 ) stop 'pushA: Dimension error'

  A(kinit:kend,ncol) = v

  end subroutine pushA

!********************************************************

  subroutine pullA(kinit,dimv,dimA1,dimA2,ncol,v,kend,A)

  implicit none

  integer, intent(in) :: kinit,dimv,dimA1,dimA2,ncol
  real, intent(out) :: v(dimv)
  integer, intent(out) :: kend
  real, intent(in) :: A(dimA1,dimA2)

  kend = kinit + dimv - 1
  if( kend.gt.dimA1 ) stop 'pullA: Dimension error'

  v = A(kinit:kend,ncol)	

  end subroutine pullA

!********************************************************

  subroutine average_mat(M,Mmean,ni,nj)

  implicit none
  integer, intent(in) :: ni, nj
  real, intent(in) :: M(ni,nj)
  real, intent(out) :: Mmean(ni)
  integer i

  do i = 1,ni
    Mmean(i) = sum(M(i,1:nj)) / nj
  end do

  end subroutine average_mat

!********************************************************

  subroutine random2(work1,n)
! Returns a vector of random values N(variance=1,mean=0)
! From Evensen's code 
   implicit none
   integer, intent(in) :: n
   real,   intent(out) :: work1(n)
   real,   allocatable :: work2(:)
   real, parameter :: pi=3.14159253589

   allocate (work2(n))

   call random_number(work1)
   call random_number(work2)
   work1= sqrt(-2.0*log(work1))*cos(2.0*pi*work2)

   deallocate(work2)

  end subroutine random2

!********************************************************

  subroutine find_node(x,y,ie,ik)
  use shyfile
  use basin
  implicit none

  real, intent(in) :: x,y
  integer, intent(in) :: ie
  integer, intent(out) :: ik
  
  integer i,ikmin
  real d,dmin

  dmin=10e10
  do i = 1,3
     ik = nen3v(i,ie)
     d = ( x - xgv(ik) )**2 + ( y - ygv(ik) )**2
     if( d.le.dmin ) then
       dmin = d
       ikmin = ik
     end if
  end do
  ik = ikmin

  end subroutine find_node

!********************************************************

  subroutine num2str(num,str)
  implicit none
  integer, intent(in) :: num
  character(len=3), intent(out) :: str

  if( (num.ge.0).and.(num.lt.10) ) then
    write(str,'(a2,i1)') '00',num
  elseif( (num.ge.10).and.(num.lt.100) ) then
    write(str,'(a1,i2)') '0',num
  elseif( (num.ge.100).and.(num.lt.1000) ) then
    write(str,'(i3)') num
  else
    stop 'num2str: num out of range'
  end if

  end subroutine num2str
