!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
module mod_enkf

  use mod_para
  use mod_manage_obs
  use mod_ens_state

  implicit none

  real, save, allocatable :: D(:,:)		! matrix holding perturbed measurments
  real, save, allocatable :: D1(:,:)		! perturbed innovation vectors
  real, save, allocatable :: E(:,:)		! matrix holding perturbations (mode=?3)
  real, save, allocatable :: R(:,:)		! Obs error cov matrix

  real, save, allocatable :: S(:,:)		! matrix holding HA' mod perturbations
  real, save, allocatable :: innov(:)		! innovation vector holding d-H*mean(Abk)
  real, save, allocatable :: HA(:,:)		! matrix HA, ens model on obs space


contains


!********************************************************

! R (only used if mode=?1 or ?2) (not with low rank R: (N-1) R = EE')
!
  subroutine make_matrices

  implicit none

  integer nook

  write(*,*) 'Building the assimilation arrays...'

  allocate(D(nobs_tot,nrens),E(nobs_tot,nrens),&
           R(nobs_tot,nobs_tot))
  allocate(S(nobs_tot,nrens),innov(nobs_tot))
  allocate(HA(nobs_tot,nrens),D1(nobs_tot,nrens))

  R = 0.
  nook = 0
  ! 0D Levels
  !
  if (n_0dlev > 0) then
     write(*,*) 'Assimilation of sea level'
     call fill_scalar_0d('0DLEV',n_0dlev,nook,o0dlev)
  end if

  ! 0D Temperature
  !
  if (n_0dtemp > 0) then
     write(*,*) 'Assimilation of temperature'
     call fill_scalar_0d('0DTEM',n_0dtemp,nook,o0dtemp)
  end if

  ! 0D Salinity
  !
  if (n_0dsalt > 0) then
     write(*,*) 'Assimilation of salinity'
     call fill_scalar_0d('0DSAL',n_0dsalt,nook,o0dsalt)
  end if

  ! 2D currents
  !
  if (n_2dvel > 0) then
     write(*,*) 'Assimilation of velocities'
     call fill_scurrents(n_2dvel,nook)
  end if

  end subroutine make_matrices

!********************************************************

  subroutine fill_scalar_0d(olabel,nfile,nook,ostate)

  implicit none

  character(len=5), intent(in) :: olabel
  integer, intent(in) :: nfile
  type(scalar_0d), intent(inout) :: ostate(nfile)
  integer, intent(inout) :: nook

  integer nf,ne
  real x,y
  integer iemin,kmin
  real oval,stdv
  real inn1,inn2,mval(nrens),mvalm
  real pvec(nrens)
  character(len=5) :: nal

  do nf = 1,nfile 

     ! create a white/red noise random vector with mean 0 and std 1
     !
     call make_0Dpert(olabel,nrens,nanal,ostate(nf)%id,pvec,atime_an,TTAU_0D)

     ! next if the observation is not good
     !
     if (ostate(nf)%stat > 1) cycle

     nook = nook + 1

     x = ostate(nf)%x
     y = ostate(nf)%y
     call find_el_node(x,y,iemin,kmin)
     if (verbose) write(*,*) 'Internal node nearest to obs: ',kmin

     ! compute the observation errors R
     !
     R(nook,nook) = ostate(nf)%std**2

     ! compute the model perturbed values, S = HA' and HA
     ! Remember for enKF: Aa = Af + Abk' [HA']^t [ U L^-1 U^t ] D' and D' = D-HA
     !
     select case (olabel)
       case ('0DLEV')
	mval = Abk(:)%z(kmin)
	mvalm = Abk_m%z(kmin)
       case ('0DTEM')
	mval = Abk(:)%t(1,kmin)
        mvalm = Abk_m%t(1,kmin)
       case ('0DSAL')
	mval = Abk(:)%s(1,kmin)
        mvalm = Abk_m%s(1,kmin)
     end select

     S(nook,:) = mval - mvalm
     HA(nook,:) = mval

     oval = ostate(nf)%val
     inn1 = oval - mvalm
     innov(nook) = inn1
     stdv = ostate(nf)%std

     call check_spread(inn1,stdv,mval,mvalm)

     if (verbose)&
     write(*,'(a25,2x,i4,3f8.4)') 'nobs, vobs, vmod, innov:',&
              nf,oval,mvalm,inn1
 
     ! compute the perturbations E, the perturbed observations D
     ! and the innovation vectors D1
     !
     E(nook,:) = stdv * pvec
     D(nook,:) = oval + (stdv * pvec)
     D1(nook,:) = D(nook,:) - HA(nook,:)
 
  end do

  end subroutine fill_scalar_0d

!********************************************************

  subroutine fill_scurrents(nfile,nook)

  use levels
  implicit none

  integer, intent(in) :: nfile
  integer, intent(inout) :: nook
  integer nf,ix,iy,ne
  real pvec1(nrens),pvec2(nrens)
  real x,y
  integer iemin,kmin
  real uu,vv,stdv
  real inn1,inn2
  real mvalu(nrens),mvalv(nrens),mvalum,mvalvm
  real h_1st_layer

  if (size(hlv) <= 1)&
     error stop 'fill_scurrents: a 3D sim is necessary to assimilate surface currents'

  do nf = 1,nfile 
     
    ! create a white/red noise random vector with mean 0 and std 1
    !
    call make_0Dpert('2DVEL',nrens,nanal,o2dvel(nf)%id,pvec1,atime_an,TTAU_2D)
    call make_0Dpert('2DVEL',nrens,nanal,o2dvel(nf)%id,pvec2,atime_an,TTAU_2D)

    do iy = 1,o2dvel(nf)%ny
    do ix = 1,o2dvel(nf)%nx

     ! next if the observation is not good
     !
     if (o2dvel(nf)%stat(ix,iy) > 1) cycle

     nook = nook + 2	! 2 obs, u and v components

     x = o2dvel(nf)%x(ix,iy)
     y = o2dvel(nf)%y(ix,iy)
     call find_el_node(x,y,iemin,kmin)
     if (verbose) write(*,*) 'Internal element nearest to obs: ',iemin

     h_1st_layer = hlv(1) + Abk_m%z(kmin)

     ! compute the observation errors R
     !
     R(nook-1,nook-1) = o2dvel(nf)%std(ix,iy)**2
     R(nook,nook) = o2dvel(nf)%std(ix,iy)**2

     ! compute the model perturbed values, S = HA' and HA
     ! Remember for enKF: Aa = Af + Abk' [HA']^t [ U L^-1 U^t ] D' and D' = D-HA
     !
     mvalu = Abk(:)%u(1,iemin)
     mvalv = Abk(:)%v(1,iemin)
     mvalum = Abk_m%u(1,iemin)
     mvalvm = Abk_m%v(1,iemin)
     S(nook-1,:) = mvalu - mvalum
     S(nook,:) = mvalv - mvalvm
     HA(nook-1,:) = mvalu
     HA(nook,:) = mvalv

     ! compute the innovation vector
     !
     uu = o2dvel(nf)%u(ix,iy) * h_1st_layer
     vv = o2dvel(nf)%v(ix,iy) * h_1st_layer
     stdv = o2dvel(nf)%std(ix,iy)

     inn1 = uu - mvalum
     inn2 = vv - mvalvm
     innov(nook-1) = inn1
     innov(nook) = inn2

     call check_spread_speed(inn1,inn2,stdv,mvalu,mvalv,mvalum,mvalvm)

     ! compute the perturbations E, the perturbed observations D
     ! and the innovation vectors D1
     !
     ! u component
     !
     E(nook-1,:) = stdv * pvec1
     D(nook-1,:) = uu + (stdv * pvec1)
     D1(nook-1,:) = D(nook-1,:) - HA(nook-1,:)
     ! v component
     !
     E(nook,:) = stdv * pvec2
     D(nook,:) = vv + (stdv * pvec2)
     D1(nook,:) = D(nook,:) - HA(nook,:)
 
    end do
    end do

  end do

  end subroutine fill_scurrents

end module mod_enkf
