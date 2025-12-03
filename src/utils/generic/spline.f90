
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017-2019  Georg Umgiesser
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

! spline routines
!
!**********************************************************************

! revision log :
!
! 05.12.2017	ggu	changed VERS_7_5_39
! 24.01.2018	ggu	changed VERS_7_5_41
! 14.02.2019	ggu	changed VERS_7_5_56
! 02.12.2025	ggu	double version, 1st and 2nd derivative
! 03.12.2025	ggu	integrated spline_compute_extremum()

!===============================================================
        module spline
!===============================================================

        implicit none

        private

        INTERFACE spline_init
        MODULE PROCEDURE   spline_init_real &
     &			 , spline_init_double
        END INTERFACE

        INTERFACE spline_eval
        MODULE PROCEDURE   spline_eval_real &
     &			 , spline_eval_double
        END INTERFACE

        INTERFACE spline_1st_deriv
        MODULE PROCEDURE   spline_1st_deriv_real &
     &			 , spline_1st_deriv_double
        END INTERFACE

        INTERFACE spline_2nd_deriv
        MODULE PROCEDURE   spline_2nd_deriv_real &
     &			 , spline_2nd_deriv_double
        END INTERFACE

        INTERFACE spline_compute_extremum
        MODULE PROCEDURE   spline_compute_extremum_real &
     &			 , spline_compute_extremum_double
        END INTERFACE

        public :: spline_init
        public :: spline_eval
        public :: spline_1st_deriv
        public :: spline_2nd_deriv
        public :: spline_compute_extremum

! Usage:
!
!	call spline_init(n,xa,ya,y2a)
!	call spline_eval(n,xa,ya,y2a,x,y)
!	call spline_1st_deriv(n,xa,ya,ya2,x,y1)
!	call spline_2nd_deriv(n,xa,ya,ya2,x,y2)
!	call spline_compute_extremum(n,xa,ya,ya2,j,xe,ierr)

!===============================================================
        contains
!===============================================================

	subroutine spline_init_real(n,xa,ya,y2a)

! spline initialization

	implicit none

	integer n	!length of arrays
	real xa(n)	!array containing ordered x values
	real ya(n)	!array containing y values
	real y2a(n)	!array containing computed second derivs at return

	real, parameter :: high = 1.e+30
	real aux(n)	!aux array

	call spline_NR_real(xa,ya,n,high,high,aux,y2a)

	end

!**********************************************************************

	subroutine spline_eval_real(n,xa,ya,y2a,x,y)

! spline evaluation

	implicit none

	integer n	!length of arrays
	real xa(n)	!array containing ordered x values
	real ya(n)	!array containing y values
	real y2a(n)	!array containing computed second derivs at return
	real x		!x value for which to evaluate spline
	real y		!computed y value at return

	call splint_NR_real(xa,ya,y2a,n,x,y)

	end

!**********************************************************************

	subroutine spline_init_double(n,xa,ya,y2a)

! spline initialization

	implicit none

	integer n		!length of arrays
	double precision xa(n)	!array containing ordered x values
	double precision ya(n)	!array containing y values
	double precision y2a(n)	!array containing computed second dv at return

	double precision, parameter :: high = 1.e+30
	double precision aux(n)	!aux array

	call spline_NR_double(xa,ya,n,high,high,aux,y2a)

	end

!**********************************************************************

	subroutine spline_eval_double(n,xa,ya,y2a,x,y)

! spline evaluation

	implicit none

	integer n		!length of arrays
	double precision xa(n)	!array containing ordered x values
	double precision ya(n)	!array containing y values
	double precision y2a(n)	!array containing computed second dv at return
	double precision x	!x value for which to evaluate spline
	double precision y	!computed y value at return

	call splint_NR_double(xa,ya,y2a,n,x,y)

	end

!**********************************************************************
!**********************************************************************
!**********************************************************************

      SUBROUTINE SPLINE_NR_real(X,Y,N,YP1,YPN,U,Y2)
      implicit none
      integer n
      real X(N),Y(N),Y2(N),U(N)
      real yp1,ypn,qn,un,p,sig
      integer i,k
      IF (YP1.GT..99E30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) &
     &      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      END DO
      IF (YPN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
      END DO
      RETURN
      END

!************************

      SUBROUTINE SPLINT_NR_real(XA,YA,Y2A,N,X,Y)
      implicit none
      integer n
      real XA(N),YA(N),Y2A(N)
      real x,y
      integer klo,khi,k
      real h,a,b
      KLO=1
      KHI=N
      do while (KHI-KLO.GT.1)
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      end do
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) STOP 'Bad XA input.'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+ &
     &      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      RETURN
      END

!************************

      SUBROUTINE HUNT_real(N,XX,X,JLO)
      implicit none
      integer n
      real XX(N)
      real x
      integer jlo,inc,jhi,jm
      LOGICAL ASCND
      ASCND=XX(N).GT.XX(1)
      IF(JLO.LE.0.OR.JLO.GT.N)THEN
        JLO=0
        JHI=N+1
        GO TO 3
      ENDIF
      INC=1
      IF(X.GE.XX(JLO).EQV.ASCND)THEN
1       JHI=JLO+INC
        IF(JHI.GT.N)THEN
          JHI=N+1
        ELSE IF(X.GE.XX(JHI).EQV.ASCND)THEN
          JLO=JHI
          INC=INC+INC
          GO TO 1
        ENDIF
      ELSE
        JHI=JLO
2       JLO=JHI-INC
        IF(JLO.LT.1)THEN
          JLO=0
        ELSE IF(X.LT.XX(JLO).EQV.ASCND)THEN
          JHI=JLO
          INC=INC+INC
          GO TO 2
        ENDIF
      ENDIF
3     IF(JHI-JLO.EQ.1)RETURN
      JM=(JHI+JLO)/2
      IF(X.GT.XX(JM).EQV.ASCND)THEN
        JLO=JM
      ELSE
        JHI=JM
      ENDIF
      GO TO 3
      END

!**********************************************************************

      SUBROUTINE SPLINE_NR_double(X,Y,N,YP1,YPN,U,Y2)
      implicit none
      integer n
      double precision X(N),Y(N),Y2(N),U(N)
      double precision yp1,ypn,qn,un,p,sig
      integer i,k
      IF (YP1.GT..99E30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) &
     &      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      END DO
      IF (YPN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
      END DO
      RETURN
      END

!************************

      SUBROUTINE SPLINT_NR_double(XA,YA,Y2A,N,X,Y)
      implicit none
      integer n
      double precision XA(N),YA(N),Y2A(N)
      double precision x,y
      integer klo,khi,k
      double precision h,a,b
      KLO=1
      KHI=N
      do while (KHI-KLO.GT.1)
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      end do
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) STOP 'Bad XA input.'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+ &
     &      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      RETURN
      END

!************************

      SUBROUTINE HUNT_double(N,XX,X,JLO)
      implicit none
      integer n
      double precision XX(N)
      double precision x
      integer jlo,inc,jhi,jm
      LOGICAL ASCND
      ASCND=XX(N).GT.XX(1)
      IF(JLO.LE.0.OR.JLO.GT.N)THEN
        JLO=0
        JHI=N+1
        GO TO 3
      ENDIF
      INC=1
      IF(X.GE.XX(JLO).EQV.ASCND)THEN
1       JHI=JLO+INC
        IF(JHI.GT.N)THEN
          JHI=N+1
        ELSE IF(X.GE.XX(JHI).EQV.ASCND)THEN
          JLO=JHI
          INC=INC+INC
          GO TO 1
        ENDIF
      ELSE
        JHI=JLO
2       JLO=JHI-INC
        IF(JLO.LT.1)THEN
          JLO=0
        ELSE IF(X.LT.XX(JLO).EQV.ASCND)THEN
          JHI=JLO
          INC=INC+INC
          GO TO 2
        ENDIF
      ENDIF
3     IF(JHI-JLO.EQ.1)RETURN
      JM=(JHI+JLO)/2
      IF(X.GT.XX(JM).EQV.ASCND)THEN
        JLO=JM
      ELSE
        JHI=JM
      ENDIF
      GO TO 3
      END

!**********************************************************************
!**********************************************************************
!**********************************************************************

	subroutine spline_1st_deriv_real(n,xa,ya,ya2,x,y1)

! computes first derivative at x and gives it back in y1

	implicit none

	integer n
	real xa(n),ya(n),ya2(n)
	real x,y1

	integer j
	real dx,dy,a,b

	j = 0
        call hunt_real(n,xa,x,j)

	dy = ya(j+1)-ya(j)
	dx = xa(j+1)-xa(j)

	a = (xa(j+1)-x)/dx
	b = 1.-a

	y1 = dy/dx - (3.*a*a-1.)*dx*ya2(j)/6. + (3.*b*b-1.)*dx*ya2(j+1)/6.

	end

!************************

	subroutine spline_2nd_deriv_real(n,xa,ya,ya2,x,y2)

! computes second derivative at x and gives it back in y2

	implicit none

	integer n
	real xa(n),ya(n),ya2(n)
	real x,y2

	integer j
	real dx,dy,a,b

	j = 0
        call hunt_real(n,xa,x,j)

	a = (xa(j+1)-x)/(xa(j+1)-xa(j))
	b = 1.-a

	y2 = a*ya2(j) + b*ya2(j+1)

	end

!**********************************************************************

	subroutine spline_1st_deriv_double(n,xa,ya,ya2,x,y1)

! computes first derivative at x and gives it back in y1

	implicit none

	integer n
	double precision xa(n),ya(n),ya2(n)
	double precision x,y1

	integer j
	double precision dx,dy,a,b

	j = 0
        call hunt_double(n,xa,x,j)

	dy = ya(j+1)-ya(j)
	dx = xa(j+1)-xa(j)

	a = (xa(j+1)-x)/dx
	b = 1.-a

	y1 = dy/dx - (3.*a*a-1.)*dx*ya2(j)/6. + (3.*b*b-1.)*dx*ya2(j+1)/6.

	end

!************************

	subroutine spline_2nd_deriv_double(n,xa,ya,ya2,x,y2)

! computes second derivative at x and gives it back in y2

	implicit none

	integer n
	double precision xa(n),ya(n),ya2(n)
	double precision x,y2

	integer j
	double precision dx,dy,a,b

	j = 0
        call hunt_double(n,xa,x,j)

	a = (xa(j+1)-x)/(xa(j+1)-xa(j))
	b = 1.-a

	y2 = a*ya2(j) + b*ya2(j+1)

	end

!**********************************************************************

        subroutine spline_compute_extremum_real(n,x,y,y2,j,xe,ierr)

! look for extremum between data j and j+1
!
! ierr:
!
!       -3      j is out of bounds
!       -2      a == 0 (flat, linear first derivative)
!       -1      no extremum found, discriminant < 0
!       0       no extremum found (out of bounds)
!       1       one extremum found, is in xe(1)
!       2       two extrema found, are in xe(1) and xe(2)

        implicit none

        integer n
        real x(n),y(n),y2(n)
        integer j
        real xe(2)
        integer ierr

        double precision dx,dy
        double precision r1,r2,r3,r
        double precision a,b,c
        double precision disc
        double precision x1,x2

        xe = 0

        ierr = -3
        if( j >= n ) return
        if( j <= 0 ) return

        dy = y(j+1) - y(j)
        dx = x(j+1) - x(j)

        r1 = dy/dx
        r2 = dx*y2(j)/6.
        r3 = dx*y2(j+1)/6.
        r = (r1+r2-r3)/3.

        a = (r3-r2)
        b = 2.*r2
        c = r-r2

        ierr = -2
        if( a == 0. ) return

        disc = b*b - 4*a*c
        ierr = -1
        if( disc < 0. ) return

        x1 = (-b + sqrt(disc) ) / (2.*a)
        x2 = (-b - sqrt(disc) ) / (2.*a)

        ierr = 0

        if( x1 >= 0. .and. x1 <= 1. ) then
          ierr = ierr + 1
          xe(ierr) = x1
        end if

        if( x2 >= 0. .and. x2 <= 1. ) then
          ierr = ierr + 1
          xe(ierr) = x2
        end if

        xe = dx*xe + x(j)

        end

!**********************************************************************

        subroutine spline_compute_extremum_double(n,x,y,y2,j,xe,ierr)

! look for extremum between data j and j+1
!
! ierr:
!
!       -3      j is out of bounds
!       -2      a == 0 (flat, linear first derivative)
!       -1      no extremum found, discriminant < 0
!       0       no extremum found (out of bounds)
!       1       one extremum found, is in xe(1)
!       2       two extrema found, are in xe(1) and xe(2)

        implicit none

        integer n
        double precision x(n),y(n),y2(n)
        integer j
        double precision xe(2)
        integer ierr

        double precision dx,dy
        double precision r1,r2,r3,r
        double precision a,b,c
        double precision disc
        double precision x1,x2

        xe = 0

        ierr = -3
        if( j >= n ) return
        if( j <= 0 ) return

        dy = y(j+1) - y(j)
        dx = x(j+1) - x(j)

        r1 = dy/dx
        r2 = dx*y2(j)/6.
        r3 = dx*y2(j+1)/6.
        r = (r1+r2-r3)/3.

        a = (r3-r2)
        b = 2.*r2
        c = r-r2

        ierr = -2
        if( a == 0. ) return

        disc = b*b - 4*a*c
        ierr = -1
        if( disc < 0. ) return

        x1 = (-b + sqrt(disc) ) / (2.*a)
        x2 = (-b - sqrt(disc) ) / (2.*a)

        ierr = 0

        if( x1 >= 0. .and. x1 <= 1. ) then
          ierr = ierr + 1
          xe(ierr) = x1
        end if

        if( x2 >= 0. .and. x2 <= 1. ) then
          ierr = ierr + 1
          xe(ierr) = x2
        end if

        xe = dx*xe + x(j)

        end

!===============================================================
        end module spline
!===============================================================

!**********************************************************************
!**********************************************************************
!**********************************************************************

        subroutine spline_tst

! test spline

	use spline

        implicit none

        integer, parameter :: ndim = 8
        integer, parameter :: inc = 8		!increment for plot

        real xa(ndim), ya(ndim)
        real ya2(ndim)

        integer i,n,nmax
        real x,y
	real xmax,fact

        data xa /0.,1.,2.,4.,6.,7.,9.,12./
        data ya /0.,1.,3.,4.,2.,0.,2.,5./

        n = ndim
	xmax = xa(ndim)
	nmax = xmax * inc
	fact = 1./inc

        do i=1,n
          write(70,'(2f12.4)') xa(i),ya(i)
        end do

        call spline_init(n,xa,ya,ya2)

        do i=0,nmax
          x = i * fact
          !if( xa .gt. 12. ) xa = 12.
          call spline_eval(n,xa,ya,ya2,x,y)
          write(71,'(2f12.4)') x,y
        end do

	write(6,*) 'spline test successfully finished: fort.70/71'

        end

!***************************************************************

!       program spline_test
!       call spline_tst
!       end

!***************************************************************


