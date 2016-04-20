  module mod_enKF


  implicit none

  integer sdate, stime
  double precision tobs
  integer nrens, xdim, nanl
  character (len=80) :: basfile
  character (len=80) :: baseobs
  character (len=80), allocatable :: rstfile(:)

  real, allocatable :: A(:,:)
  real, allocatable :: Am(:)

  ! Observations
  integer nobs
  character(len=6), allocatable :: tyobs(:)
  real, allocatable :: xobs(:), yobs(:), zobs(:)
  real, allocatable :: vobs(:), stdvobs(:)	!Values and standard deviations


  contains


!********************************************************

  subroutine read_info

  implicit none

  integer n

  open(20, file='analysis.info', status='old')
  read(20,*) sdate
  read(20,*) stime
  read(20,*) tobs
  read(20,*) basfile
  read(20,*) nrens
  read(20,*) nanl
  read(20,*) baseobs

  allocate(rstfile(nrens))
  do n = 1,nrens
     read(20,*) rstfile(n)
  end do

  close(20)

  end subroutine read_info

!********************************************************

  subroutine read_basin
 
  use basin
  implicit none

  open(21, file=trim(basfile), status='old', form='unformatted')
  call sp13_get_par(21,nkn,nel,ngr,mbw)
  call sp13rr(21,nkn,nel)
  close(21)

  end subroutine read_basin

!********************************************************

  subroutine rst_read(iunit,filin)

  implicit none

  integer iunit
  character(len=*), intent(in) :: filin
  double precision atime
  integer ierr

  open(iunit, file=filin, status='old', form='unformatted')
  call dts_to_abs_time(sdate, stime, atime)
  call rst_read_restart_file(atime + tobs, iunit, ierr)
  close(iunit)

  end subroutine rst_read


!********************************************************

  subroutine store_state(ne)

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

  xdim = dimbr + dimcon
  if( ibarcl_rst.gt.0 ) xdim = xdim + dimbc

  if(.not.(allocated(A))) allocate(A(xdim,nrens))
 
  ! Fill with znv
  call fillA(1,nkn,nrens,xdim,ne,znv,kend,A)

  ! Fill with transports
  do nl = 1,nlv
     kinit = kend
     call fillA(kinit,nel,nrens,xdim,ne,utlnv,kend,A)
     kinit = kend
     call fillA(kinit,nel,nrens,xdim,ne,vtlnv,kend,A)
  end do

  ! Fill with temperature and salinity
  if( ibarcl_rst.gt.0 ) then
    do nl = 1,nlv
       kinit = kend
       call fillA(kinit,nkn,nrens,xdim,ne,saltv,kend,A)
       kinit = kend
       call fillA(kinit,nkn,nrens,xdim,ne,tempv,kend,A)
    end do
  end if
 
  !if( iconz.eq. 1 ) TODO

  if( kend.ne.xdim ) then
    write(*,*) 'Wrong dimensions: ',kend,xdim
    stop 'store_state: Error'
  end if

  end subroutine store_state


!********************************************************

  subroutine fillA(kinit,dimv,dimA1,dimA2,ne,v,kend,A)

  implicit none

  integer, intent(in) :: kinit,dimv,dimA1,dimA2,ne
  real, intent(in) :: v(dimv)
  integer, intent(out) :: kend
  real, intent(inout) :: A(dimA1,dimA2)

  kend = kinit + dimv - 1
  if( kend.gt.dimA2 ) stop 'fillA: Dimension error'

  A(kinit:kend,ne) = v

  end subroutine fillA

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

  subroutine read_obs(tobs,nanl,baseobs)
  implicit none
  double precision, intent(in) :: tobs
  integer, intent(in) :: nanl
  character(len=*), intent(in) :: baseobs

  character(len=80) :: filename
  character(len=3) :: nanll
  integer ios

  double precision tt
  

  write(nanll,'(i0.3)') nanl

  filename = trim(baseobs)//'_an'//nanll//'.obs'
  write(*,*) 'Observation file: ',filename

  open(25,file = filename, status = 'old', form = 'formatted', iostat = ios)
  if( ios.ne.0 ) stop 'read_obs: error opening file'

  read(25,*) tt,nobs
  if( tt.ne.tobs ) stop 'read_obs: error in times'

  allocate(tyobs(nobs), xobs(nobs), yobs(nobs), zobs(nobs), vobs(nobs), stdvobs(nobs))

  do n = 1,nobs
     read(25,*) tyobs(n), xobs(n), yobs(n), zobs(n), vobs(n), stdvobs(n)
  end do

  close(25)

  end subroutine read_obs

  end module mod_enKF
