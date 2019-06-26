!
! Copyright (C) 2019, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
! Local analysis subroutines
!

!*************************************************************

  subroutine local_analysis

  use mod_mod_states	  
  use mod_para
  use mod_enkf
  use mod_ens_state
  use mod_manage_obs
  use basin

  implicit none

  integer k,ne,i
  integer nk_l,ne_l
  integer no,nno
  real dist,w
  real xe,ye

  integer, parameter ::  lkdim = 2*nnlv+1
  integer, parameter ::  lnedim = 2*nnlv
  integer :: nobs_tot_k,nobs_tot_e

  real, allocatable :: xobs(:),yobs(:)
  real, allocatable :: Ak_loc(:,:),Ak_an(:,:),Ak_bk(:,:)
  real, allocatable :: Ane_loc(:,:),Ane_an(:,:),Ane_bk(:,:)

  integer, allocatable :: ido(:)
  real, allocatable :: wo(:)

  real, allocatable :: innovl(:),D1l(:,:),Sl(:,:),El(:,:),Rl(:,:)

  real, parameter :: eps_la = 1e-4   !Minimum value to limit the local analysis window
  
  allocate(xobs(nobs_tot),yobs(nobs_tot))
  allocate(Ak_loc(lkdim,nrens),Ak_an(lkdim,nrens),Ak_bk(lkdim,nrens))
  allocate(Ane_loc(lnedim,nrens),Ane_an(lnedim,nrens),Ane_bk(lnedim,nrens))
  allocate(ido(nobs_tot),wo(nobs_tot))

  ! find the coordinates of the measurements
  ! note that in read_obs only obs of the same type
  ! are allowed
  do no = 1,nobs_tot
     if (islev /= 0) then
	  if (no <= n_0dlev) then
		  xobs(no) = o0dlev(no)%x
		  yobs(no) = o0dlev(no)%y
	  else if ((no > n_0dlev).and.(no <= n_0dlev+n_1dlev)) then
	          write(*,*) 'Do it!'
	          stop
		  !xobs(no) = o1dlev(no)%x
		  !yobs(no) = o1dlev(no)%y
	  else if (no > n_0dlev+n_1dlev) then
	          write(*,*) 'Do it!'
	          stop
		  !xobs(no) = o2dlev(no)%x
		  !yobs(no) = o3dlev(no)%y
	  end if
     else if (istemp /= 0) then
	  write(*,*) 'Do it!'
	  stop
     else if (issalt /= 0) then
	  write(*,*) 'Do it!'
	  stop
     else if (isvel /= 0) then
	  write(*,*) 'Do it!'
	  stop
     end if

  end do


  !----------------------nodes-------------------------
  nk_l = 0 	! number of nodes corrected
  update_randrot = .true.
  do k = 1,nnkn

     call type_to_kmat(Ak_bk,Ak_an,k,lkdim)

     Ak_loc = Ak_an	! set the local analysis equal to the global

     ido = 0
     wo = 0
     nno = 0	! number of local observations
     do no = 1,nobs_tot

	dist = sqrt( (xgv(k)-xobs(no))**2 + (ygv(k)-yobs(no))**2 )
	call find_weight(rho_loc,dist,w)

	if ( w > eps_la ) then
           nno = nno + 1
           ido(nno) = no
           wo(nno) = w
	end if

     end do

     nobs_tot_k = nno
     if (nobs_tot_k > 0) then
	     
        allocate(innovl(nobs_tot_k),D1l(nobs_tot_k,nrens),Sl(nobs_tot_k,nrens))
        allocate(El(nobs_tot_k,nrens),Rl(nobs_tot_k,nobs_tot_k))

        do no = 1,nobs_tot_k

	     ! tapering of the innovations and of the ens anomalies
	     innovl(no) = innov(ido(no)) * wo(no)
	     D1l(no,:) = D1(ido(no),:) * wo(no)
	     Sl(no,:) = S(ido(no),:) * wo(no)

	     ! no tapering of the obs covariance
	     El(no,:) = E(ido(no),:)
	     Rl(no,no) = R(ido(no),ido(no))

	end do

	nk_l = nk_l + 1
	Ak_loc = Ak_bk	! set equal to the state before the analysis

	call analysis(Ak_loc,Rl,El,Sl,D1l,innovl,lkdim,nrens,nobs_tot_k,.false.,&
		           truncation,rmode,update_randrot)

	update_randrot = .false.	! true just the first time. Beware in OMP

        !call save_X5('local',atime_an)

        deallocate(innovl,D1l,Sl,El,Rl)

	! This is to check the position
	!write(120,*) '1',nk_l,'7',xgv(k),ygv(k)

     end if

     Ak_an = Ak_an + (Ak_loc - Ak_an)

     call kmat_to_type(Ak_an,k,lkdim)

  end do
  write(*,*) 'Number of nodes with local analysis: ',nk_l
  !----------------------end nodes-------------------------

  !----------------------elements-------------------------
  ne_l = 0 	! number of nodes corrected
  do ne = 1,nnel

     call type_to_emat(Ane_bk,Ane_an,ne,lnedim)

     Ane_loc = Ane_an	! set the local analysis equal to the global

     ido = 0
     wo = 0
     nno = 0	! number of local observations
     do no = 1,nobs_tot

	xe =   0.
	ye =   0.
	do i = 1,3
	   k = nen3v(i,ne)
	   xe = xe + xgv(k) 
	   ye = ye + ygv(k) 
	end do
	xe = xe/3.
	ye = ye/3.
	dist = sqrt( (xe-xobs(no))**2 + (ye-yobs(no))**2 )

	call find_weight(rho_loc,dist,w)

	if ( w > eps_la ) then
           nno = nno + 1
           ido(nno) = no
           wo(nno) = w
	end if

     end do

     nobs_tot_e = nno
     if (nobs_tot_e > 0) then

        allocate(innovl(nobs_tot_e),D1l(nobs_tot_e,nrens),Sl(nobs_tot_e,nrens))
        allocate(El(nobs_tot_e,nrens),Rl(nobs_tot_e,nobs_tot_e))

        do no = 1,nobs_tot_e

	     ! tapering of the innovations and of the ens anomalies
	     innovl(no) = innov(ido(no)) * wo(no)
	     D1l(no,:) = D1(ido(no),:) * wo(no)
	     Sl(no,:) = S(ido(no),:) * wo(no)

	     ! no tapering of the obs covariance
	     El(no,:) = E(ido(no),:)
	     Rl(no,no) = R(ido(no),ido(no))

	end do

	ne_l = ne_l + 1
	Ane_loc = Ane_bk	! set equal to the state before the analysis

	call analysis(Ane_loc,Rl,El,Sl,D1l,innovl,lnedim,nrens,nobs_tot_e,.false.,&
	           truncation,rmode,update_randrot)

        !call save_X5('local',atime_an)

        deallocate(innovl,D1l,Sl,El,Rl)

	! This is to check the position
	!write(121,*) '1',ne_l,'8',xe,ye

     end if

     Ane_an = Ane_an + (Ane_loc - Ane_an)

     call emat_to_type(Ane_an,ne,lnedim)

  end do
  write(*,*) 'Number of elements with local analysis: ',ne_l
  !----------------------end elements-------------------------

  end subroutine local_analysis


!*************************************************************

  subroutine find_weight(rho,dst,w)
  implicit none
  real, intent(in) :: rho,dst
  real, intent(out) :: w
  real s

  s = dst/rho

  ! Taper function of Gaspari-Cohn
  if ((s >=0) .and. (s<1)) then
     w = 1 - (5./3. * s**2) + (5./8. * s**3)  + (1./2. * s**4) - (1./4. * s**5)
  else if ((s>=1) .and. (s<2)) then
     w = - (2./3. * s**-1) + 4 - (5. * s) + (5./3. * s**2) + (5./8. * s**3)& 
         - (1./2. * s**4) + (1./12. * s**5)
  else
     w = 0.
  end if

  end subroutine find_weight


!*************************************************************

  subroutine type_to_kmat(Ak_bk,Ak_an,k,kdim)
  use mod_ens_state
  implicit none
  integer, intent(in) :: k,kdim
  real, intent(out) :: Ak_bk(kdim,nrens),Ak_an(kdim,nrens)
  integer n

  do n = 1,nrens
    Ak_bk(1,n) = Abk(n)%z(k)
    Ak_bk(2:nnlv+1,n) = Abk(n)%t(:,k)
    Ak_bk(nnlv+2:2*nnlv+1,n) = Abk(n)%s(:,k)
  end do

  do n = 1,nrens
    Ak_an(1,n) = Aan(n)%z(k)
    Ak_an(2:nnlv+1,n) = Aan(n)%t(:,k)
    Ak_an(nnlv+2:2*nnlv+1,n) = Aan(n)%s(:,k)
  end do

  end subroutine type_to_kmat

!*************************************************************

  subroutine type_to_emat(Ane_bk,Ane_an,ne,nedim)
  use mod_ens_state
  implicit none
  integer, intent(in) :: ne,nedim
  real, intent(out) :: Ane_bk(nedim,nrens),Ane_an(nedim,nrens)
  integer n

  do n = 1,nrens
    Ane_bk(1:nnlv,n) = Abk(n)%u(:,ne)
    Ane_bk(nnlv+1:2*nnlv,n) = Abk(n)%v(:,ne)
  end do

  do n = 1,nrens
    Ane_an(1:nnlv,n) = Aan(n)%u(:,ne)
    Ane_an(nnlv+1:2*nnlv,n) = Aan(n)%v(:,ne)
  end do

  end subroutine type_to_emat

!*************************************************************

  subroutine kmat_to_type(Ak_an,k,kdim)
  use mod_ens_state
  implicit none
  integer, intent(in) :: k,kdim
  real, intent(in) :: Ak_an(kdim,nrens)
  integer n

  do n = 1,nrens
    Aan(n)%z(k) = Ak_an(1,n)
    Aan(n)%t(:,k) = Ak_an(2:nnlv+1,n)
    Aan(n)%s(:,k) = Ak_an(nnlv+2:2*nnlv+1,n)
  end do

  end subroutine kmat_to_type

!*************************************************************

  subroutine emat_to_type(Ane_an,ne,nedim)
  use mod_ens_state
  implicit none
  integer, intent(in) :: ne,nedim
  real, intent(in) :: Ane_an(nedim,nrens)
  integer n
  
  do n = 1,nrens
    Aan(n)%u(:,ne) = Ane_an(1:nnlv,n)
    Aan(n)%v(:,ne) = Ane_an(nnlv+1:2*nnlv,n)
  end do

  end subroutine emat_to_type
