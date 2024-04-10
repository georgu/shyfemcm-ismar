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

  integer n,ix,iy

  write(*,*) 'Building the assimilation arrays...'

  nobs_ok = 0
  if (islev /= 0) then
     do n = 1,n_0dlev
        if (o0dlev(n)%stat == 0) nobs_ok = nobs_ok + 1
     end do 
  else if (issalt /= 0) then
     do n = 1,n_0dsalt
        if (o0dsalt(n)%stat == 0) nobs_ok = nobs_ok + 1
     end do 
  else if (istemp /= 0) then
     do n = 1,n_0dtemp
        if (o0dtemp(n)%stat == 0) nobs_ok = nobs_ok + 1
     end do 		
  else if (isvel /= 0) then
     do n = 1,n_2dvel
       do iy = 1,o2dvel(n)%ny
       do ix = 1,o2dvel(n)%nx
         if (o2dvel(n)%stat(ix,iy) == 0) nobs_ok = nobs_ok + 2
       end do
       end do
     end do
  else
     error stop 'Observation type not valid.'
  end if

  allocate(D(nobs_ok,nrens),E(nobs_ok,nrens),&
           R(nobs_ok,nobs_ok))
  allocate(S(nobs_ok,nrens),innov(nobs_ok))
  allocate(HA(nobs_ok,nrens),D1(nobs_ok,nrens))

  R = 0.
  ! 0D Levels
  !
  if (n_0dlev > 0) then
     write(*,*) 'Assimilation of sea level'
     call fill_scalar_0d('0DLEV',n_0dlev,o0dlev)
  end if

  ! 0D Temperature
  !
  if (n_0dtemp > 0) then
     write(*,*) 'Assimilation of temperature'
     call fill_scalar_0d('0DTEM',n_0dtemp,o0dtemp)
  end if

  ! 0D Salinity
  !
  if (n_0dsalt > 0) then
     write(*,*) 'Assimilation of salinity'
     call fill_scalar_0d('0DSAL',n_0dsalt,o0dsalt)
  end if

  ! 2D currents
  !
  if (n_2dvel > 0) then
     write(*,*) 'Assimilation of velocities'
     call fill_scurrents(n_2dvel)
  end if

  end subroutine make_matrices

!********************************************************

  subroutine fill_scalar_0d(olabel,nfile,ostate)

  implicit none

  character(len=5), intent(in) :: olabel
  integer, intent(in) :: nfile
  type(scalar_0d), intent(inout) :: ostate(nfile)
  integer :: nook

  integer nf,ne
  real x,y
  integer iemin,kmin
  real oval,stdv,stdm
  real inn1,inn2,mval(nrens),mvalm,sinn
  real pvec(nrens)
  character(len=5) :: nal
  real maxinn

  nook = 0
  do nf = 1,nfile 

     ! create a white/red noise random vector with mean 0 and std 1
     !
     call make_0Dpert(olabel,nrens,nanal,ostate(nf)%id,pvec,atime_an,TTAU_0D)

     ! next if the observation is not good
     !
     if (ostate(nf)%stat /= 0) cycle

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
	do ne = 1,nrens
	   mval(ne) = Abk(ne)%z(kmin)
	end do
	mvalm = Abk_m%z(kmin)
       case ('0DTEM')
	do ne = 1,nrens
	   mval(ne) = Abk(ne)%t(1,kmin)
	end do
        mvalm = Abk_m%t(1,kmin)
       case ('0DSAL')
	do ne = 1,nrens
	   mval(ne) = Abk(ne)%s(1,kmin)
	end do
        mvalm = Abk_m%s(1,kmin)
     end select

     stdm = sqrt( sum(mval**2)/nrens - (sum(mval)/nrens)**2 )

     S(nook,:) = mval - mvalm
     HA(nook,:) = mval

     oval = ostate(nf)%val
     stdv = ostate(nf)%std

     inn1 = oval - mvalm


     ! check innovation value
     maxinn = sqrt(inn_alpha*(stdv**2+stdm**2))
     sinn = sign(1.,inn1)
     if (inn1**2 > maxinn**2) then
        if (verbose) write(*,*) 'Innovation too high (inn, max-inn): ',inn1,sinn*maxinn
	inn1 = sinn*maxinn
     end if

     innov(nook) = inn1

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

  if (nook /= nobs_ok) error stop 'The nubmer of good observations is wrong.'

  end subroutine fill_scalar_0d

!********************************************************

  subroutine fill_scurrents(nfile)

  use levels
  implicit none

  integer, intent(in) :: nfile
  integer :: nook
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

  nook = 0
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
     do ne = 1,nrens
        mvalu(ne) = Abk(ne)%u(1,iemin)
        mvalv(ne) = Abk(ne)%v(1,iemin)
     end do
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

  if (nook /= nobs_ok) error stop 'The nubmer of good observations is wrong.'

  end subroutine fill_scurrents

end module mod_enkf
