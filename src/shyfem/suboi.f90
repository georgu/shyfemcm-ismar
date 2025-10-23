
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003-2004,2009,2014-2015,2018-2019  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

! revision log :
!
! 20.08.2003	ggu	new laplacian interpolation
! 02.09.2003	ggu	some comments, write to .dat file
! 30.10.2003	ggu	subroutine prepare_bc_l included in this file
! 04.03.2004	ggu	writes also number of variables (1)
! 11.03.2009	ggu	bug fix -> declare hev() here
! 10.02.2014	ggu	finished optimal interpolation
! 26.02.2015	ggu	changed VERS_7_1_5
! 26.03.2018	ggu	some parameters are now arrays (sigma,rr)
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.02.2019	ggu	changed VERS_7_5_60
! 22.11.2020	ggu	some more comments
! 20.02.2025	ggu	some refactoring
!
! notes :
!
!****************************************************************

        subroutine opt_intp(nobs,xobs,yobs,zobs,bobs                    &
     &                  ,nback,bback,xback,yback,zback                  &
     &                  ,rl,rlmax,sigma,rr,zanal)

! computes optimal interpolation
!
! normally only zanal is computed and returned
! however, if bobs is .true. (background grid is given) then
! bobs and zback are initialized to average values and returned back
!
! values of zobs == flag indicate no observation available
! in this case the error covariance is increased by a factor of fact
! (flag == -999 and fact == 100.)

	implicit none

	integer,intent(in):: nobs	!size of observations
	real,intent(in):: xobs(nobs)	!x-coordinates of observations
	real,intent(in):: yobs(nobs)	!y-coordinates of observations
	real,intent(in):: zobs(nobs)	!values of observations
	real,intent(inout):: bobs(nobs)	!background values at observation points
	integer,intent(in):: nback	!size of background field
	logical,intent(in):: bback	!background values are given
	real,intent(in):: xback(nback)	!x-coordinates of background
	real,intent(in):: yback(nback)	!y-coordinates of background
	real,intent(inout):: zback(nback)	!values of background
	real,intent(in):: rl(nobs)	!length scale for covariance
	real,intent(in):: rlmax(nobs)	!max radius to be considered
	real,intent(in):: sigma(nobs)	!std of background on obs points
	real,intent(in):: rr(nobs)	!std of observation error matrix
	real,intent(out):: zanal(nback)	!analysis on return

	integer ivec(nobs)		!aux vector (n)
	double precision dobs(nobs)	!increments at obs points (z-h(x^b))
	double precision rvec(nobs)	!aux vector (n)
	double precision rr2(nobs)	!square of observation error (n)
	double precision sigma2(nobs)	!square of background field (n)
	double precision rl2(nobs)	!square of background field (n)
	double precision rlmax2(nobs)	!square of background field (n)
	double precision rmat(nobs,nobs)!aux matrix (nxn)

	double precision ani,ano,anr,cond
	double precision rorig(nobs,nobs)!aux matrix (nxn)
	double precision rind(nobs,nobs)!aux matrix (nxn)
	double precision zobs1(nobs)	!zobserved changed internally

	logical bcheck,bverbose
	integer i,j,ki,kj,k,n,iacu,no
	double precision, parameter :: flag = -999.
	double precision, parameter :: fact = 100.
	double precision rmean
	double precision xi,yi,xj,yj,xk,yk
	double precision dist2,r,acu

!	------------------------------------------
!	set some parameters
!	------------------------------------------

	bverbose = .true.		!writes infomation to terminal
	bverbose = .false.		!writes infomation to terminal
	bcheck = .true.			!computes and writes some checks
	n = nobs

	rl2(:) = rl(:)**2
	rlmax2(:) = rlmax(:)**2
	sigma2(:) = sigma(:)**2
	rr2(:) = rr(:)**2
	zobs1 = zobs			!save zobs to not alter input value

!	------------------------------------------
!	if background not given create it
!	------------------------------------------

	if( .not. bback ) then

	  !------------------------------------------
	  !compute mean of observations - set background to mean
	  !------------------------------------------

	  no = 0
	  rmean = 0.
	  do i=1,n
	    if( zobs(i) /= flag ) then
	      no = no + 1
	      rmean = rmean + zobs(i)
	    else
	      rr2(i) = ( fact*rr(i) )**2
	    end if
	  end do
	  if( no == 0 ) stop 'error stop: no observations left'

	  rmean = rmean/no
	  zback = rmean
	  bobs = rmean

	  where( zobs1 == flag ) zobs1 = rmean

	  if( bverbose ) then
	    write(6,*) 'no background field given... using mean: ',rmean
	  end if

	end if

	!write(6,*) bback,rmean

!	------------------------------------------
!	create observational innovation vector
!	------------------------------------------

	dobs = zobs1 - bobs

!	------------------------------------------
!	set up covariance matrix H P^b H^T
!	------------------------------------------

	if( bverbose ) then
	  write(6,*) 'setting up covariance matrix'
	end if

	rmat = 0.
	do j=1,n
	  if( zobs(j) == flag ) cycle
	  xj = xobs(j)
	  yj = yobs(j)
	  do i=1,n
	    if( zobs(i) == flag ) cycle
	    xi = xobs(i)
	    yi = yobs(i)
	    dist2 = (xi-xj)**2 + (yi-yj)**2
	    r = 0.
	    if( dist2 <= rlmax2(i) ) r = sigma2(i) * exp( -dist2/rl2(i) )
	    rmat(i,j) = r
	  end do
	end do

!	------------------------------------------
!	add observation error matrix
!	------------------------------------------

	do j=1,n
	  rmat(j,j) = rmat(j,j) + rr2(j)
	end do

!	------------------------------------------
!	invert matrix
!	------------------------------------------

	if( bverbose ) then
	  write(6,*) 'inverting covariance matrix'
	end if

	if( bcheck ) rorig = rmat
	call dmatinv(rmat,ivec,rvec,n,n)

	if( bcheck ) then
	  call dmatnorm(ani,rmat,n,n)
	  call dmatnorm(ano,rorig,n,n)
	  cond = ani * ano
          call dmatmult(rorig,rmat,rind,n,n)
          do i=1,n
            rind(i,i) = rind(i,i) - 1.
          end do
          call dmatnorm(anr,rind,n,n)
	  !write(6,'(a,3f12.4,e12.4)') 'cond: ',ano,ani,cond,anr
	end if

!	------------------------------------------
!	multiply inverted matrix with observational innovation
!	------------------------------------------

	do i=1,n
	  if( zobs(i) == flag ) cycle
	  acu = 0.
	  do j=1,n
	    acu = acu + rmat(i,j) * dobs(j)
	  end do
	  rvec(i) = acu
	end do

!	------------------------------------------
!	multiply P^b H^T with vector and add to background to obtain analysis
!	------------------------------------------

	do k=1,nback
	  xk = xback(k)
	  yk = yback(k)
	  iacu = 0
	  acu = 0.
	  do j=1,n
	    if( zobs(j) == flag ) cycle
	    xj = xobs(j)
	    yj = yobs(j)
	    dist2 = (xk-xj)**2 + (yk-yj)**2
	    r = 0.
	    if( dist2 .le. rlmax2(j) ) then
	      iacu = iacu + 1
	      r = sigma2(j) * exp( -dist2/rl2(j) )
	    end if
	    acu = acu + r * rvec(j)
	  end do
	  !if( iacu == 0 ) goto 99
	  zanal(k) = zback(k) + acu
	end do

!	------------------------------------------
!	end of routine
!	------------------------------------------

	return
   99	continue
	write(6,*) j,acu,rl(j),rlmax(j)
	write(6,*) 'radius for rlmax too small: no points inside'
	stop 'error stop opt_intp: rlmax'
	end

!****************************************************************
!****************************************************************
!****************************************************************
! test routines
!****************************************************************
!****************************************************************
!****************************************************************

	subroutine opt_intp_test_1_obs

	implicit none

	integer nobs
	parameter (nobs=1)
	integer nb
	parameter (nb=10)

	integer nb1,nback
	parameter (nb1=nb+1,nback=nb1)

	logical bback
	integer i,j,ib
	real x,y,dx,dy
	real xobs(nobs)
	real yobs(nobs)
	real zobs(nobs)
	real bobs(nobs)
	real rl(nobs)
	real rlmax(nobs)
	real rr(nobs)
	real sigma(nobs)
	real xback(nback)
	real yback(nback)
	real zback(nback)
	real zanal(nback)

	write(6,*) 'running 1 obs test'

	bback = .true.
	rl = 0.1
	rlmax = 2.0
	rr = 0.01
	sigma = 1.

	do i=1,nobs
	  xobs(i) = 0.5
	  yobs(i) = 0.
	  zobs(i) = 1.
	  bobs(i) = 0.
	end do

	dx = 1./nb

	ib = 0
	do i=0,nb
	  ib = ib + 1
	  xback(ib) = i*dx
	  yback(ib) = 0.
	  zback(ib) = 0.
	end do

	write(6,*) nb,nback,ib
	if( ib .ne. nback ) stop 'error stop: ib .ne. nback'

        call opt_intp(nobs,xobs,yobs,zobs,bobs                          &
     &                  ,nback,bback,xback,yback,zback                  &
     &                  ,rl,rlmax,sigma,rr,zanal)

	write(6,*) nobs
	do i=1,nobs
	  write(6,*) xobs(i),yobs(i),zobs(i),bobs(i)
	end do
	write(6,*) nback,bback
	do i=1,nback
	  write(6,*) xback(i),yback(i),zback(i),zanal(i)
	end do

	end

!****************************************************************

	subroutine opt_intp_test_sin_cos

	implicit none

	integer nobs
	parameter (nobs=50)
	integer nb
	parameter (nb=10)

	integer nb1,nback
	parameter (nb1=nb+1,nback=nb1*nb1)

	logical bback
	integer i,j,ib
	real x,y,dx,dy
	double precision rms
	real xobs(nobs)
	real yobs(nobs)
	real zobs(nobs)
	real bobs(nobs)
	real rl(nobs)
	real rlmax(nobs)
	real rr(nobs)
	real sigma(nobs)
	real xback(nback)
	real yback(nback)
	real zback(nback)
	real zsol(nback)
	real zanal(nback)

	write(6,*) 'running sin/cos test'

	bback = .false.
	rl = 0.4
	rlmax = 2.0
	rr = 0.01
	sigma = 1.

	do i=1,nobs
	  !x = rand(0)
	  !y = rand(0)
	  call random_number(x)
	  call random_number(y)
	  xobs(i) = x
	  yobs(i) = y
	  zobs(i) = sin(6.*x) * cos(6.*y)
	end do

	dx = 1./nb
	dy = 1./nb

	ib = 0
	do i=0,nb
	  do j=0,nb
	    ib = ib + 1
	    x = j*dx
	    y = i*dy
	    xback(ib) = x
	    yback(ib) = y
	    zsol(ib) = sin(6.*x) * cos(6.*y)
	  end do
	end do

	write(6,*) nb,nback,ib
	if( ib .ne. nback ) stop 'error stop: ib .ne. nback'

        call opt_intp(nobs,xobs,yobs,zobs,bobs                          &
     &                  ,nback,bback,xback,yback,zback                  &
     &                  ,rl,rlmax,sigma,rr,zanal)

	write(6,*) nobs
	do i=1,nobs
	  write(6,*) xobs(i),yobs(i),zobs(i),bobs(i)
	end do
	write(6,*) nback,bback
	rms = 0.
	do i=1,nback
	  write(6,*) xback(i),yback(i),zsol(i),zanal(i)
	  rms = rms + (zsol(i)-zanal(i))**2
	end do
	rms = sqrt(rms/nback)
	write(6,*) 'rms: ',rms

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine opt_intp_test
	call opt_intp_test_1_obs
	call opt_intp_test_sin_cos
	end

!******************************************************************

!	program opt_intp_main
!	call opt_intp_test
!	end

!******************************************************************

