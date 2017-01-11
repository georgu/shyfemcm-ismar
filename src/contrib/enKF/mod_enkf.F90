  module mod_enkf

  use mod_states

  implicit none

  character (len=80), save :: basfile, obsfile, rstfile

  integer, save :: nrens, na

  ! Observations
  double precision, save :: tobs
  integer, save :: nobs
  character(len=6), save, allocatable :: tyobs(:)
  real, save, allocatable :: xobs(:), yobs(:), zobs(:)
  real, save, allocatable :: vobs(:), stdvobs(:)	!Values and standard deviations

  type(states), allocatable, save  :: A(:) 		! ens states
  type(states), save  :: Am		! mean state

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

  read(20,*) nrens	! number of ens members
  read(20,*) na		! analysis step
  read(20,*) basfile	! name of bas file (no extension)
  read(20,*) tobs	! current time of the observations
  read(20,*) obsfile	! name of obs file list
  
  ! Allocates the type A to store the ens states
  allocate(A(nrens))

  close(20)

  end subroutine read_info

!********************************************************

  subroutine read_basin

  use basin
  implicit none

  open(21, file=basfile, status='old', form='unformatted')
  call basin_read_by_unit(21)
  close(21)

  if( ( nkn.ne.nnkn ).or.( nel.ne.nnel) ) stop "read_basin: dim error"

  end subroutine read_basin


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

  subroutine average_mat(rstw)

  use m_write_ensemble
  implicit none
  integer ne
  integer rstw

  Am = 0.
  do ne = 1,nrens
    Am = Am + A(ne)
  end do
  Am = states_real_mult(Am,1./nrens)

  if(rstw == -1) call rst_write(-1,na,tobs)
  if(rstw == -2) call rst_write(-2,na,tobs)

  end subroutine average_mat


!********************************************************

  subroutine make_D_E_R

! R (only used if mode=?1 or ?2) (no for low-rank sq root)

  implicit none
  real rand_v(nrens)
  integer n,ne

  allocate(D(nobs,nrens),E(nobs,nrens),R(nobs,nobs))

  R = 0.	!Observations are indipendent
  do n = 1,nobs
     if( trim(tyobs(n)).ne.'level' ) stop 'Observation type not still implemented.'
     ! Makes a random vector
     call random2(rand_v,nrens)
     do ne = 1,nrens
        rand_v(ne) = (stdvobs(n) * rand_v(ne)) + vobs(n)
     end do
     ! TODO: add more observation types and different perturbation methods

     D(n,:) = rand_v
     !E(n,:) = rand_v - vobs(n)
     E(n,:) = rand_v - sum(rand_v)/nrens		! E must have mean 0

     R(n,n) = stdvobs(n)**2
  end do

  end subroutine make_D_E_R

!********************************************************

  subroutine make_S_innov

  use regular
  implicit none

  integer iel, ikn
  integer n
  double precision av_mod,inn

  allocate (S(nobs,nrens),innov(nobs))
  
  do n = 1,nobs

     ! Finds the nearest element
     call find_element(xobs(n),yobs(n),iel)

     if( trim(tyobs(n)).eq.'level' ) then
        ! Finds the nearest node
        call find_node(xobs(n),yobs(n),iel,ikn)
     else
        stop 'Observation type not still implemented.'
     end if

     av_mod = Am%z(ikn)
     S(n,:) = A(:)%z(ikn) - av_mod
     inn = vobs(n) - av_mod
     
     call check_innov_val(inn,tyobs(n))

     innov(n) = inn

     write(*,*) 'Observation: ',n,trim(tyobs(n)),vobs(n),av_mod,innov(n)
     
  end do

  end subroutine make_S_innov

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

  end module mod_enKF


