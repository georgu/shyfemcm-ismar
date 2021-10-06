!
! Copyright (C) 2019, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
! Local analysis subroutines
!

!*************************************************************

  subroutine local_analysis

  use mod_ens_state
  use mod_manage_obs
  use mod_restart , only : ibarcl_rst

  implicit none

  integer k,ne,i
  integer nk_l,ne_l
  integer no,nno
  real dist,w
  real xe,ye

  logical :: local_upd_rand = .true.

  integer ::  lkdim 
  integer ::  lnedim 

  real, allocatable :: xobs(:),yobs(:),rho_loc(:)
  real, allocatable :: Ak_an(:,:),Ak_bk(:,:)
  real, allocatable :: Ane_an(:,:),Ane_bk(:,:)

  if (ibarcl_rst == 0) then
	  lkdim = 1
	  lnedim = 2*nnlv
  else
	  lkdim = 1 + 2*nnlv
	  lnedim = 2*nnlv
  end if
  
  allocate(xobs(nobs_tot),yobs(nobs_tot),rho_loc(nobs_tot))
  allocate(Ak_an(lkdim,nrens),Ak_bk(lkdim,nrens))
  allocate(Ane_an(lnedim,nrens),Ane_bk(lnedim,nrens))

  ! find the coordinates of the measurements
  ! note that in read_obs only obs of the same type
  ! are allowed, otherwise this is wrong.
  if (verbose) write(*,*) 'Observations in local analysis:'
  do no = 1,nobs_tot
     if (islev /= 0) then
	  if (no <= n_0dlev) then		!e.g. from timeseries
		  xobs(no) = o0dlev(no)%x
		  yobs(no) = o0dlev(no)%y
		  rho_loc(no) = o0dlev(no)%rhol
	  else if ((no > n_0dlev).and. &
		  (no <= n_0dlev+n_1dlev)) then	!e.g. altimeter track
	          write(*,*) 'Do it!'
	          stop
		  !xobs(no) = o1dlev(no)%x
		  !yobs(no) = o1dlev(no)%y
	  else if (no > n_0dlev+n_1dlev) then	!e.g. altimeter map
	          write(*,*) 'Do it!'
	          stop
		  !xobs(no) = o2dlev(no)%x
		  !yobs(no) = o2dlev(no)%y
	  end if
     else if (istemp /= 0) then
	  if (no <= n_0dtemp) then			!e.g. from timeseries
		  xobs(no) = o0dtemp(no)%x
		  yobs(no) = o0dtemp(no)%y
		  rho_loc(no) = o0dtemp(no)%rhol
	  else if ((no > n_0dtemp).and. &
		  (no <= n_0dtemp+n_1dtemp)) then	!e.g. from profiles
	          write(*,*) 'Do it!'
	          stop
		  !xobs(no) = o1dtemp(no)%x
		  !yobs(no) = o1dtemp(no)%y
	  else if (no > n_0dtemp+n_1dtemp) then		!e.g. from sst maps
	          write(*,*) 'Do it!'
	          stop
		  !xobs(no) = o2dtemp(no)%x
		  !yobs(no) = o2dtemp(no)%y
	  end if
     else if (issalt /= 0) then
	  if (no <= n_0dsalt) then			!e.g. from timeseries
		  xobs(no) = o0dsalt(no)%x
		  yobs(no) = o0dsalt(no)%y
		  rho_loc(no) = o0dsalt(no)%rhol
	  else if ((no > n_0dsalt).and. &		!e.g. from profiles
		  (no <= n_0dsalt+n_1dsalt)) then
	          write(*,*) 'Do it!'
	          stop
		  !xobs(no) = o1dsalt(no)%x
		  !yobs(no) = o1dsalt(no)%y
	  else if (no > n_0dsalt+n_1dsalt) then		!e.g. from surface maps
	          write(*,*) 'Do it!'
	          stop
		  !xobs(no) = o2dsalt(no)%x
		  !yobs(no) = o2dsalt(no)%y
	  end if
     else if (isvel /= 0) then
	  write(*,*) 'Do it!'
	  stop
     end if

     if (verbose) write(*,*) 'x,y,rho = ',xobs(no),yobs(no),rho_loc(no)

  end do

  !----------------------nodes-------------------------
  nk_l = 0 	! number of nodes corrected
!$OMP PARALLEL PRIVATE(k,Ak_bk,Ak_an),SHARED(ibarcl_rst,lkdim,nrens,nobs_tot,local_upd_rand,xobs,yobs,rho_loc)
!$OMP DO
  do k = 1,nnkn

     call type_to_kmat(ibarcl_rst,Ak_bk,k,lkdim,nrens)

     call locan_k(k,lkdim,nrens,nobs_tot,xobs,yobs,rho_loc,Ak_bk,local_upd_rand,nk_l,Ak_an)

     call kmat_to_type(ibarcl_rst,Ak_an,k,lkdim,nrens)

  end do
!$OMP ENDDO
!$OMP END PARALLEL
  write(*,*) 'Number of nodes with local analysis: ',nk_l
  !----------------------end nodes-------------------------

  !----------------------elements-------------------------
  ne_l = 0 	! number of elements corrected
!$OMP PARALLEL PRIVATE(ne,Ane_bk,Ane_an),SHARED(ibarcl_rst,lnedim,nrens,nobs_tot,local_upd_rand,xobs,yobs,rho_loc)
!$OMP DO
  do ne = 1,nnel

     call type_to_emat(Ane_bk,ne,lnedim,nrens)

     call locan_e(ne,lnedim,nrens,nobs_tot,xobs,yobs,rho_loc,Ane_bk,local_upd_rand,ne_l,Ane_an)

     call emat_to_type(Ane_an,ne,lnedim,nrens)

  end do
!$OMP ENDDO
!$OMP END PARALLEL 
  write(*,*) 'Number of elements with local analysis: ',ne_l
  !----------------------end elements-------------------------

  end subroutine local_analysis


!*************************************************************

     subroutine locan_k(nk,lkdim,nren,no_tot,xo,yo,rhoo,Ak_bk,local_upd_rand,nk_l,Ak_an)
     use basin
     use mod_enkf
     use mod_para
     implicit none
     integer, intent(in) :: nk,lkdim,nren,no_tot
     real, intent(in) :: Ak_bk(lkdim,nren)
     real, intent(in) :: xo(no_tot),yo(no_tot),rhoo(no_tot)
     logical, intent(inout) :: local_upd_rand
     integer, intent(inout) :: nk_l
     real, intent(out) :: Ak_an(lkdim,nren)

     real :: Ak_loc(lkdim,nren)
     integer :: no,nno,no_k
     integer, allocatable :: ido(:)
     real, allocatable :: wo(:)

     real dist,w
     real, allocatable :: innovl(:),D1l(:,:),Sl(:,:),El(:,:),Rl(:,:)
     
     real, parameter :: eps_la = 1e-4   !Minimum value to limit the local analysis window

     integer, save :: icall = 0

     allocate(ido(no_tot),wo(no_tot))

     Ak_loc = Ak_bk	! set the local analysis equal to the global

     ido = 0
     wo = 0
     nno = 0	! number of local observations
     do no = 1,no_tot

	dist = sqrt( (xgv(nk)-xo(no))**2 + (ygv(nk)-yo(no))**2 )
	call find_weight(rhoo(no),dist,w)

	if ( w > eps_la ) then
           nno = nno + 1
           ido(nno) = no
           wo(nno) = w
	end if

     end do

     no_k = nno
     if (no_k > 0) then
	     
        allocate(innovl(no_k),D1l(no_k,nren),Sl(no_k,nren))
        allocate(El(no_k,nren),Rl(no_k,no_k))
	!write(*,*) 'Analysis of node: ',nk

        do no = 1,no_k

	     ! tapering of the innovations and of the ens anomalies
	     innovl(no) = innov(ido(no)) * wo(no)
	     D1l(no,:) = D1(ido(no),:) * wo(no)
	     Sl(no,:) = S(ido(no),:) * wo(no)

	     ! no tapering of the obs covariance
	     El(no,:) = E(ido(no),:)
	     Rl(no,no) = R(ido(no),ido(no))

	end do

	call analysis(Ak_loc,Rl,El,Sl,D1l,innovl,lkdim,nren,no_k,.false.,&
		           truncation,rmode,lrandrot,local_upd_rand,lsymsqrt,&
			   inflate,infmult)

        !call save_X5('local',atime_an)

        deallocate(innovl,D1l,Sl,El,Rl)

	! This is to check the position
	!write(120,*) '1',nk_l,'7',xgv(nk),ygv(nk)

!$OMP CRITICAL
	! true just for the first call to analysis
	if (icall == 0) then
		write(*,*) 'Setting the update rand rot to false...'
		local_upd_rand = .false.
		icall = 1
	end if

        nk_l = nk_l + 1
!$OMP END CRITICAL

     end if

     Ak_an = Ak_bk + (Ak_loc - Ak_bk)

     end subroutine locan_k


!*************************************************************

     subroutine locan_e(ne,lnedim,nren,no_tot,xo,yo,rhoo,Ane_bk,local_upd_rand,ne_l,Ane_an)
     use basin
     use mod_enkf
     use mod_para
     implicit none
     integer, intent(in) :: ne,lnedim,nren,no_tot
     real, intent(in) :: Ane_bk(lnedim,nren)
     real, intent(in) :: xo(no_tot),yo(no_tot),rhoo(no_tot)
     logical, intent(in) :: local_upd_rand
     integer, intent(inout) :: ne_l
     real, intent(out) :: Ane_an(lnedim,nren)

     integer  :: i,k
     real :: xe,ye
     real :: Ane_loc(lnedim,nren)
     integer :: no,nno,no_e
     integer, allocatable :: ido(:)
     real, allocatable :: wo(:)

     real dist,w
     real, allocatable :: innovl(:),D1l(:,:),Sl(:,:),El(:,:),Rl(:,:)
     
     real, parameter :: eps_la = 1e-4   !Minimum value to limit the local analysis window

     integer, save :: icall = 0

     allocate(ido(no_tot),wo(no_tot))

     Ane_loc = Ane_bk	! set the local analysis equal to the global

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
	dist = sqrt( (xe-xo(no))**2 + (ye-yo(no))**2 )

	call find_weight(rhoo(no),dist,w)

	if ( w > eps_la ) then
           nno = nno + 1
           ido(nno) = no
           wo(nno) = w
	end if

     end do

     no_e = nno
     if (no_e > 0) then

        allocate(innovl(no_e),D1l(no_e,nrens),Sl(no_e,nrens))
        allocate(El(no_e,nrens),Rl(no_e,no_e))
	!write(*,*) 'Analysis of element: ',ne

        do no = 1,no_e

	     ! tapering of the innovations and of the ens anomalies
	     innovl(no) = innov(ido(no)) * wo(no)
	     D1l(no,:) = D1(ido(no),:) * wo(no)
	     Sl(no,:) = S(ido(no),:) * wo(no)

	     ! no tapering of the obs covariance
	     El(no,:) = E(ido(no),:)
	     Rl(no,no) = R(ido(no),ido(no))

	end do

	call analysis(Ane_loc,Rl,El,Sl,D1l,innovl,lnedim,nrens,no_e,.false.,&
	           truncation,rmode,lrandrot,local_upd_rand,lsymsqrt,&
		   inflate,infmult)

        !call save_X5('local',atime_an)

        deallocate(innovl,D1l,Sl,El,Rl)

	! This is to check the position
	!write(121,*) '1',ne_l,'8',xe,ye

!$OMP CRITICAL
	ne_l = ne_l + 1
!$OMP END CRITICAL

     end if

     Ane_an = Ane_bk + (Ane_loc - Ane_bk)


     end subroutine locan_e


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

  subroutine type_to_kmat(ibrcl,Ak_bk,k,kdim,nre)
  use mod_ens_state
  implicit none
  integer, intent(in) :: ibrcl
  integer, intent(in) :: k,kdim,nre
  real, intent(out) :: Ak_bk(kdim,nre)
  integer n

  do n = 1,nre
    Ak_bk(1,n) = Abk(n)%z(k)
  end do

  if (ibrcl > 0) then
     do n = 1,nre
        Ak_bk(2:nnlv+1,n) = Abk(n)%t(:,k)
        Ak_bk(nnlv+2:2*nnlv+1,n) = Abk(n)%s(:,k)
     end do
  end if

  end subroutine type_to_kmat

!*************************************************************

  subroutine type_to_emat(Ane_bk,ne,nedim,nre)
  use mod_ens_state
  implicit none
  integer, intent(in) :: ne,nedim,nre
  real, intent(out) :: Ane_bk(nedim,nre)
  integer n

  do n = 1,nre
    Ane_bk(1:nnlv,n) = Abk(n)%u(:,ne)
    Ane_bk(nnlv+1:2*nnlv,n) = Abk(n)%v(:,ne)
  end do

  end subroutine type_to_emat

!*************************************************************

  subroutine kmat_to_type(ibrcl,Ak_an,k,kdim,nre)
  use mod_ens_state
  implicit none
  integer, intent(in) :: ibrcl
  integer, intent(in) :: k,kdim,nre
  real, intent(in) :: Ak_an(kdim,nre)
  integer n

  do n = 1,nre
    Aan(n)%z(k) = Ak_an(1,n)
  end do

  if (ibrcl > 0) then
     do n = 1,nre
        Aan(n)%t(:,k) = Ak_an(2:nnlv+1,n)
        Aan(n)%s(:,k) = Ak_an(nnlv+2:2*nnlv+1,n)
     end do
  end if

  end subroutine kmat_to_type

!*************************************************************

  subroutine emat_to_type(Ane_an,ne,nedim,nre)
  use mod_ens_state
  implicit none
  integer, intent(in) :: ne,nedim,nre
  real, intent(in) :: Ane_an(nedim,nre)
  integer n
  
  do n = 1,nre
    Aan(n)%u(:,ne) = Ane_an(1:nnlv,n)
    Aan(n)%v(:,ne) = Ane_an(nnlv+1:2*nnlv,n)
  end do

  end subroutine emat_to_type
