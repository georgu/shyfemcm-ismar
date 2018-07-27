!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
module mod_enkf

  use mod_para
  use mod_manage_obs
  use mod_ens_state

  implicit none

  real, save, allocatable :: D(:,:)		! matrix holding perturbed measurments
  real, save, allocatable :: D1(:,:)		! Innovation vectors
  real, save, allocatable :: E(:,:)		! matrix holding perturbations (mode=?3)
  real, save, allocatable :: R(:,:)		! Obs error cov matrix

  real, save, allocatable :: S(:,:)		! matrix holding HA' mod perturbations
  real, save, allocatable :: innov(:)		! innovation vector holding d-H*mean(A)
  real, save, allocatable :: HA(:,:)		! matrix HA, ens model on obs space


contains


!********************************************************

! R (only used if mode=?1 or ?2) (not with low rank R: (N-1) R = EE')
!
  subroutine make_matrices

  implicit none

  integer nook

  write(*,*) 'preparing the observations for the assimilation...'

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
		!olabel,linit,filin,eps,kinit,kend,oatime,xv,yv,zv,vv,stdvv,ostatusv

  implicit none

  character(len=5), intent(in) :: olabel
  integer, intent(in) :: nfile
  type(scalar_0d), intent(inout) :: ostate(nfile)
  integer, intent(inout) :: nook

  integer nf,ne
  real x,y
  integer iemin,kmin
  real oval,stdv
  integer ostatus
  real inn1,inn2
  real pvec(nrens)
  character(len=3) :: nal,nfil

  do nf = 1,nfile 

     ! create a white/red noise random vector with mean 0 and std 1
     !
     call make_0Dpert(olabel,nrens,nanal,ostate(nf)%id,pvec,atime,TTAU_0D)

     ! next if the observation is not good
     !
     if (ostate(nf)%status > 1) cycle

     nook = nook + 1

     x = ostate(nf)%x
     y = ostate(nf)%y
     call find_el_node(x,y,iemin,kmin)

     ! compute the observation errors R
     !
     R(nook,nook) = ostate(nf)%std**2

     ! compute the model perturbed values, S = HA' and HA
     ! Remember for enKF: Aa = Af + A' [HA']^t [ U L^-1 U^t ] D' and D' = D-HA
     !
     select case (olabel)

       case ('0DLEV')
        do ne = 1,nrens
           S(nook,ne) = A(ne)%z(kmin) - Am%z(kmin)
           HA(nook,ne) = A(ne)%z(kmin)
        end do

       case ('0DTEM')
        if (ostate(nf)%z /= 0) error stop 'fill_scalar_0d: deep temperature not implemented yet'
        do ne = 1,nrens
           S(nook,ne) = A(ne)%t(1,kmin) - Am%t(1,kmin)
           HA(nook,ne) = A(ne)%t(1,kmin)
        end do

       case ('0DSAL')
        if (ostate(nf)%z /= 0) error stop 'fill_scalar_0d: deep salinity not implemented yet'
        do ne = 1,nrens
           S(nook,ne) = A(ne)%s(1,kmin) - Am%s(1,kmin)
           HA(nook,ne) = A(ne)%s(1,kmin)
        end do

     end select

     ! check the obs std and compute the innovation vector
     !
     oval = ostate(nf)%val
     ostatus = ostate(nf)%status
     stdv = ostate(nf)%std
     call check_obs_inn(olabel,x,y,0.,oval,oval,stdv,inn1,inn2,ostatus)
     ostate(nf)%std = stdv
     innov(nook) = inn1

     if (verbose)&
     write(*,'(a25,2x,i4,3f8.4)') 'nobs, vobs, vmod, innov:',&
              nf,ostate(nf)%val,Am%z(kmin),inn1
 
     ! compute the perturbations E, the perturbed observations D
     ! and the innovation vectors D1
     !
     E(nook,:) = ostate(nf)%std * pvec
     D(nook,:) = ostate(nf)%val + (ostate(nf)%std * pvec)
     D1(nook,:) = D(nook,:) - HA(nook,:)
 
  end do

  end subroutine fill_scalar_0d

!********************************************************

  subroutine fill_scurrents(nfile,nook)

  implicit none

  integer, intent(in) :: nfile
  integer, intent(inout) :: nook
  integer nf,ix,iy,ne
  real pvec1(nrens),pvec2(nrens)
  real x,y
  integer iemin,kmin
  real uu,vv,stdv
  integer ostatus
  real inn1,inn2

  do nf = 1,nfile 
     
    ! create a white/red noise random vector with mean 0 and std 1
    !
    call make_0Dpert('u',nrens,nanal,o2dvel(nf)%id,pvec1,atime,TTAU_2D)
    call make_0Dpert('v',nrens,nanal,o2dvel(nf)%id,pvec2,atime,TTAU_2D)

    do iy = 1,o2dvel(nf)%ny
    do ix = 1,o2dvel(nf)%nx

     ! next if the observation is not good
     !
     if (o2dvel(nf)%status(ix,iy) > 1) cycle

     nook = nook + 2	! 2 obs, u and v components

     x = o2dvel(nf)%x(ix,iy)
     y = o2dvel(nf)%y(ix,iy)
     call find_el_node(x,y,iemin,kmin)

     ! compute the observation errors R
     !
     R(nook-1,nook-1) = o2dvel(nf)%std(ix,iy)**2
     R(nook,nook) = o2dvel(nf)%std(ix,iy)**2

     ! compute the model perturbed values, S = HA' and HA
     ! Remember for enKF: Aa = Af + A' [HA']^t [ U L^-1 U^t ] D' and D' = D-HA
     !
     do ne = 1,nrens
        S(nook-1,ne) = A(ne)%u(1,iemin) - Am%u(1,iemin)
        S(nook,ne) = A(ne)%v(1,iemin) - Am%v(1,iemin)
        HA(nook-1,ne) = A(ne)%u(1,iemin)
        HA(nook,ne) = A(ne)%v(1,iemin)
     end do

     ! compute the innovation vector
     !
     ostatus = o2dvel(nf)%status(ix,iy)
     uu = o2dvel(nf)%u(ix,iy)
     vv = o2dvel(nf)%v(ix,iy)
     stdv = o2dvel(nf)%std(ix,iy)
     call check_obs_inn('2DVEL',x,y,0.,uu,vv,stdv,inn1,inn2,ostatus)
     o2dvel(nf)%std(ix,iy) = stdv

     innov(nook-1) = inn1
     innov(nook) = inn2
 
     ! compute the perturbations E, the perturbed observations D
     ! and the innovation vectors D1
     !
     ! u component
     !
     E(nook-1,:) = o2dvel(nf)%std(ix,iy) * pvec1
     D(nook-1,:) = o2dvel(nf)%u(ix,iy) + (o2dvel(nf)%std(ix,iy) * pvec1)
     D1(nook-1,:) = D(nook-1,:) - HA(nook-1,:)
     ! v component
     !
     E(nook,:) = o2dvel(nf)%std(ix,iy) * pvec2
     D(nook,:) = o2dvel(nf)%v(ix,iy) + (o2dvel(nf)%std(ix,iy) * pvec2)
     D1(nook,:) = D(nook,:) - HA(nook,:)
 
    end do
    end do

  end do

  end subroutine fill_scurrents

end module mod_enkf
