  module mod_enKF

  implicit none

  integer, save :: sdate, stime
  double precision, save :: tobs
  character (len=80), save :: basfile,obsfile,rstfile

  integer, save :: nens,sdim

  ! Observations
  integer, save :: nobs
  character(len=6), save, allocatable :: tyobs(:)
  real, save, allocatable :: xobs(:), yobs(:), zobs(:)
  real, save, allocatable :: vobs(:), stdvobs(:)	!Values and standard deviations

  real, save, allocatable :: A(:,:)		! Matrix holding the states
  real, save, allocatable :: Am(:)		! Array with the mean state

  real, save, allocatable :: D(:,:)		! matrix holding perturbed measurments
  real, save, allocatable :: E(:,:)		! matrix holding perturbations (mode=?3)

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

  read(20,*) sdate	! date at time 0
  read(20,*) stime	! time at time 0
  read(20,*) tobs	! time in sec from time 0
  read(20,*) nens	! number of ens members
  read(20,*) basfile	! name of bas file (no extension)
  read(20,*) rstfile	! name of restart files. No extension, just basename,
			! then the file must finish with "_enDDD.rst" with DDD 
			! three digits denoting the n of ens member.
  read(20,*) obsfile	! name of obs file

  close(20)

  end subroutine read_info

!********************************************************

  subroutine read_basin
 
  use basin
  implicit none

  open(21, file=trim(basfile)//'.bas', status='old', form='unformatted')
  call sp13_get_par(21,nkn,nel,ngr,mbw)
  call sp13rr(21,nkn,nel)
  close(21)

  end subroutine read_basin

!********************************************************

  subroutine rst_read(ne)

  implicit none

  integer, intent(in) :: ne	! number of ens member

  character(len=3) :: nel
  integer io
  double precision atime
  integer ierr

  if( (ne.ge.0).and.(ne.lt.10) ) then
    write(nel,*) '00',ne
  elseif( (ne.ge.10).and.(ne.lt.100) ) then
    write(nel,*) '0',ne
  elseif( (ne.ge.100).and.(ne.lt.1000) ) then
    write(nel,*) ne
  else
    stop 'rst_read: ne out of range'
  end if

  open(24,file=rstfile//'_'//nel//'.rst',status='old',form='unformatted',iostat='io')

  if( io.ne.0 ) stop 'rst_read: Error opening file'

  call dts_to_abs_time(sdate, stime, atime)
  call rst_read_restart_file(atime + tobs, iunit, ierr)

  close(24)

  end subroutine rst_read

!********************************************************

  subroutine push_state(ne)

! The state is given by znv, utlnv, vtlnv, saltv, tempv, conzv
! I think that after the analysis step the variables: iwegv, zenv, 
! hm3v, rhov must be updated before saving the restarts

  use basin

  use mod_conz
  use mod_geom_dynamic
  use mod_ts
  use mod_hydro_vel
  use mod_hydro
  use levels, only : nlvdi,nlv
  use mod_restart

  implicit none

  integer, intent(in) :: ne

  integer dimbr,dimbc,dimcon
  integer kend,kinit

  integer nl

  ! Compute dimensions
  dimbr = nkn + ( 2 * nel ) * nlv 	! barotropic
  dimbc = 2 * nkn * nlv			! temp and salt
  dimcon = iconz * ( nkn * nlv )	! tracers

  sdim = dimbr + dimcon
  if( ibarcl_rst.gt.0 ) sdim = sdim + dimbc

  if(.not.(allocated(A))) allocate(A(sdim,nens))
 
  ! Fill with znv
  call pushA(1,nkn,sdim,nens,ne,znv,kend,A)

  ! Fill with transports
  do nl = 1,nlv
     kinit = kend + 1
     call pushA(kinit,nel,sdim,nens,ne,utlnv(nl,:),kend,A)
     kinit = kend + 1
     call pushA(kinit,nel,sdim,nens,ne,vtlnv(nl,:),kend,A)
  end do

  ! Fill with temperature and salinity
  if( ibarcl_rst.gt.0 ) then
    do nl = 1,nlv
       kinit = kend + 1
       call pushA(kinit,nkn,sdim,nens,ne,saltv(nl,:),kend,A)
       kinit = kend + 1
       call pushA(kinit,nkn,sdim,nens,ne,tempv(nl,:),kend,A)
    end do
  end if
 
  !if( iconz.eq. 1 ) TODO

  if( kend.ne.sdim ) then
    write(*,*) 'Wrong dimensions: ',kend,sdim
    stop 'push_state: Error'
  end if

  end subroutine push_state

!********************************************************

  subroutine pushA(kinit,dimv,dimA1,dimA2,ne,v,kend,A)

  implicit none

  integer, intent(in) :: kinit,dimv,dimA1,dimA2,ne
  real, intent(in) :: v(dimv)
  integer, intent(out) :: kend
  real, intent(inout) :: A(dimA1,dimA2)

  kend = kinit + dimv - 1
  if( kend.gt.dimA1 ) stop 'pushA: Dimension error'

  A(kinit:kend,ne) = v

  end subroutine pushA

!********************************************************

  subroutine pullA(kinit,dimv,dimA1,dimA2,ne,v,kend,A)

  implicit none

  integer, intent(in) :: kinit,dimv,dimA1,dimA2,ne
  real, intent(out) :: v(dimv)
  integer, intent(out) :: kend
  real, intent(in) :: A(dimA1,dimA2)

  kend = kinit + dimv - 1
  if( kend.gt.dimA1 ) stop 'pullA: Dimension error'

  v = A(kinit:kend,ne)

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

  subroutine read_obs
  implicit none

  integer ios

  double precision tt
  integer n

  write(*,*) 'Observation file: ',obsfile

  open(25,file = obsfile, status = 'old', form = 'formatted', iostat = ios)
  if( ios.ne.0 ) stop 'read_obs: error opening file'

  read(25,*) tt,nobs
  if( tt.ne.tobs ) stop 'read_obs: error in times'

  allocate(tyobs(nobs), xobs(nobs), yobs(nobs), zobs(nobs), vobs(nobs), stdvobs(nobs))

  do n = 1,nobs
     read(25,*) tyobs(n), xobs(n), yobs(n), zobs(n), vobs(n), stdvobs(n)
     if( trim(tyobs(n)).ne.'level' ) stop 'Observation type not still implemented.'
  end do

  close(25)

  end subroutine read_obs

!********************************************************

  subroutine make_D_E

  use m_random2
  implicit none
  real rand_v(nens)
  integer n

  allocate(D(nobs,nens),E(nobs,nens))

  do n = 1,nobs
     if( trim(tyobs(n)).ne.'level' ) stop 'Observation type not still implemented.'
     ! Makes a random vector
     call random2(rand_v,nens)
     rand_v = (stdvobs(n) * rand_v) + vobs(n)

     D(n,:) = rand_v
     E(n,:) = rand_v - vobs(n)
  end do

  end subroutine make_D_E

!********************************************************

  subroutine make_S_innov

  use regular
  implicit none

  integer iel
  integer n

  allocate (S(nobs,nens),innov(nobs))
  
  do n = 1,nobs

     if( trim(tyobs(n)).ne.'level' ) stop 'Observation type not still implemented.'
     ! Finds the nearest element
     call find_element(xobs(n),yobs(n),iel)
     ! Finds the nearest node
     call find_node(xobs(n),yobs(n),iel,ikn)

     ! This is true only because levels are the first stored
     S(n,:) = A(ikn,:)
     innov(n) = vobs(n) - Am(ikn)
     
  end do

  end subroutine make_S_innov

!********************************************************

  subroutine find_node(x,y,ie,ik)
  use shyfile
  implicit none

  real x,y
  integer ie,ik
  real k,d(3)

  do i = 1,3
     k = nen3v(i,iel)
     d(i) = ( x - xgv(k) )**2 + ( y - ygv(k) )**2
  end do

  i = minloc(d)
  ik = nen3v(i,iel)

  end subroutine find_node

!********************************************************

  subroutine pull_state(ne)

  use basin

  use mod_conz
  use mod_geom_dynamic
  use mod_ts
  use mod_hydro_vel
  use mod_hydro
  use levels, only : nlvdi,nlv
  use mod_restart

  implicit none
  integer ne

  integer kinit,kend
  integer nl

  call pullA(1,nkn,sdim,nens,ne,znv,kend,A)
  do nl = 1,nlv
     kinit = kend + 1
     call pullA(kinit,nel,sdim,nens,ne,utlnv(nl,:),kend,A)
     kinit = kend + 1
     call pullA(kinit,nel,sdim,nens,ne,vtlnv(nl,:),kend,A)
  end do

  ! Fill with temperature and salinity
  if( ibarcl_rst.gt.0 ) then
    do nl = 1,nlv
       kinit = kend + 1
       call pullA(kinit,nkn,sdim,nens,ne,saltv(nl,:),kend,A)
       kinit = kend + 1
       call pushA(kinit,nkn,sdim,nens,ne,tempv(nl,:),kend,A)
    end do
  end if
 
  !if( iconz.eq. 1 ) TODO

  if( kend.ne.sdim ) then
    write(*,*) 'Wrong dimensions: ',kend,sdim
    stop 'pull_state: Error'
  end if

  end subroutine pull_state

!********************************************************

  subroutine rst_write(ne)

  use basin

  use mod_conz
  use mod_geom_dynamic
  use mod_ts
  use mod_hydro_vel
  use mod_hydro
  use levels, only : nlvdi,nlv
  use mod_restart

  implicit none
  integer ne

  TODO


  end subroutine rst_write

  end module mod_enKF
