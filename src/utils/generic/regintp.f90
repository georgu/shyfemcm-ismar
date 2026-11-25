
!--------------------------------------------------------------------------
!
!    Copyright (C) 2026  Georg Umgiesser
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

! interpolation from/onto regular grid
!
! revision log :
!
! 15.11.2026    ggu     created from scratch
! 25.11.2026    ggu     eliminated useless routine
! 25.11.2026    ggu     call to get/set_regpar has changed

!****************************************************************

	subroutine bilinear(x,y,xa,ya,za,z)

	implicit none

	real x,y
	real xa(4),ya(4)
	real za(4)
	real z

	end

!****************************************************************

	subroutine obtain_angle(nx,ny,xa,ya,phi,dx,dy)

	implicit none

	integer nx,ny
	real xa(nx,ny),ya(nx,ny)
	real phi,dx,dy

	real dxy21,dxy41
	real xo,yo,xr,yr

	dxy21 = sqrt( (xa(nx,1)-xa(1,1))**2 + (ya(nx,1)-ya(1,1))**2 )
	dxy41 = sqrt( (xa(1,ny)-xa(1,1))**2 + (ya(1,ny)-ya(1,1))**2 )
	
	dx = xa(nx,1)-xa(1,1)
	dy = ya(nx,1)-ya(1,1)

	xo = dx
	yo = dy
	xr = dxy21
	yr = 0.

	phi = atan2(xo*yr-xr*yo,xo*xr+yo*yr)

	end

!****************************************************************

        subroutine fem2reg(nel,nen3v,nkn,afem,xgv,ygv,regpar,nx,ny,areg)

! interpolation of av onto a regular net (nodal values)
!
! afem                  array to be interpolated
! areg                  matrices of interpolated values (u,v,z,h)
! nx,ny                 dimension of matrices
!
! pxareg,pyareg         coordinates of lower left point of matrix
! pxdreg,pydreg         grid size of matrix
! pzlreg                value of z for land points

        !use basin

        implicit none

	integer nel
	integer nen3v(3,nel)
	integer nkn
        real afem(nkn)
        real xgv(nkn)
        real ygv(nkn)
	real regpar(7)
        integer nx,ny
        real areg(nx,ny)

        logical bwater(nel)

        bwater = .true.

        call fem2reg_mask(nel,bwater,nen3v,nkn,afem,xgv,ygv,regpar,nx,ny,areg)

        end

!****************************************************************

        subroutine fem2reg_mask(nel,bwater,nen3v,nkn,afem,xgv,ygv,regpar,nx,ny,areg)

! interpolation of av onto a regular net (nodal values) with mask
!
! bwater                mask for water points
! av                    array to be interpolated
! am                    matrices of interpolated values (u,v,z,h)
! ip,jp                 dimension of matrices
!
! pxareg,pyareg         coordinates of lower left point of matrix
! pxdreg,pydreg         grid size of matrix
! pzlreg                value of z for land points

        implicit none

! arguments
	logical bwater(nel)
	integer nel
	integer nen3v(3,nel)
	integer nkn
        real afem(nkn)
        real xgv(nkn)
        real ygv(nkn)
	real regpar(7)
        integer nx,ny
        real areg(nx,ny)

! parameter
        double precision eps
        parameter ( eps = 1.d-14 )
! local
        real pxareg,pyareg,pxdreg,pydreg,pzlreg
        real x0,y0,dx,dy,flag
        integer i,j,ii,iii,ie,k,kn,iin
        integer imin,imax,jmin,jmax
        integer iwet,nxx,nyy
        double precision x(3),y(3),z(3),a(3),b(3),c(3)
        double precision zh,fh,f,xp,yp
        double precision xmin,xmax,ymin,ymax
! function
        integer intriad

	nxx = nint(regpar(1))
	nyy = nint(regpar(2))
	if( nx /= nxx .or. ny /= nyy ) goto 91
	pxareg = regpar(3)
	pyareg = regpar(4)
	pxdreg = regpar(5)
	pydreg = regpar(6)
	flag = regpar(7)

        areg = flag

        do ie=1,nel
          if( bwater(ie) ) then !wet
            iwet = 0
            do i=1,3
                kn=nen3v(i,ie)
                x(i)=xgv(kn)
                y(i)=ygv(kn)
                z(i)=afem(kn)
                if( z(i) > flag ) iwet = iwet + 1
            end do
            if( iwet .ne. 3 ) cycle

            !f=0.
            do i=1,3
                ii=mod(i,3)+1
                iii=mod(ii,3)+1
                a(i)=x(ii)*y(iii)-x(iii)*y(ii)
                b(i)=y(ii)-y(iii)
                c(i)=x(iii)-x(ii)
                !f=f+a(i)
            end do
            f = c(3)*b(2) - c(2)*b(3)           !bug_f_64bit
            if( f .le. eps ) goto 99

            xmax=max(x(1),x(2),x(3))
            xmin=min(x(1),x(2),x(3))
            ymin=min(y(1),y(2),y(3))
            ymax=max(y(1),y(2),y(3))

            imin=(xmin-pxareg)/pxdreg+1.99
            imax=(xmax-pxareg)/pxdreg+1.01
            jmin=(ymin-pyareg)/pydreg+1.99
            jmax=(ymax-pyareg)/pydreg+1.01

            if(imin.lt.1) imin=1
            if(imax.gt.nx)imax=nx
            if(jmin.lt.1) jmin=1
            if(jmax.gt.ny)jmax=ny

            do i=imin,imax
                do j=jmin,jmax
                    xp=(i-1)*pxdreg+pxareg
                    yp=(j-1)*pydreg+pyareg

                    iin=intriad(x,y,xp,yp)

                    if(iin.ne.0) then
                        zh=0.
                        do k=1,3
                           fh=(a(k)+xp*b(k)+yp*c(k))/f
                           zh=zh+z(k)*fh
                        end do
                        areg(i,j)=zh
                    end if
                end do
            end do
          end if
        end do

        return
   91   continue
	write(6,*) 'nx,nxx: ',nx,nxx
	write(6,*) 'ny,nyy: ',ny,nyy
	stop 'error stop fem2reg_mask: mismatch between dimensions and regpar'
   99   continue
        write(6,*) ie,f
        write(6,*) x
        write(6,*) y
        write(6,*) a
        write(6,*) b
        write(6,*) c
        stop 'error stop fem2reg_mask: area of element'
        end

!****************************************************************

        function intriad(x,y,xp,yp)

! point in triangle or not (double precision version)
!
! x,y           array of coordinates of vertices of triangle
! xp,yp         coordinates of point
! intri         1: point is in triangle  0: point outside (return value)

        implicit none

! arguments
        integer intriad
        double precision x(3),y(3),xp,yp
! local
        integer k1,k2
        double precision x21,y21,xn,yn
        double precision scal,eps
! save
        save eps
        data eps /1.e-13/

        intriad=0

        do k1=1,3
           k2=mod(k1,3)+1
           x21=x(k2)-x(k1)
           y21=y(k2)-y(k1)
           yn = x21
           xn = -y21
           scal=(xp-x(k1))*xn+(yp-y(k1))*yn
           if(scal.lt.0.) return
        end do

        intriad=1        !inside

        end

!****************************************************************

        subroutine create_regpar(dx,dy,xmin,ymin,xmax,ymax,flag,regpar)

        implicit none

        real, intent(in) :: dx,dy
        real, intent(in) :: xmin,ymin,xmax,ymax
        real, intent(in) :: flag
        real, intent(out) :: regpar(7)

        integer nx,ny
        real x0,y0,x1,y1

        if( dx <= 0. .or. dy <= 0. ) goto 96
        if( xmin >= xmax ) goto 96
        if( ymin >= ymax ) goto 96

        x0 = int(xmin/dx)*dx
        y0 = int(ymin/dy)*dy
	if( x0 > xmin ) x0 = x0 - dx
	if( y0 > ymin ) y0 = y0 - dy

        nx = 2 + (xmax - xmin) / dx
        x1 = x0 + (nx-1)*dx
        if( x1 < xmax ) nx = nx + 1
        x1 = x0 + (nx-1)*dx

        ny = 2 + (ymax - ymin) / dy
        y1 = y0 + (ny-1)*dy
        if( y1 < ymax ) ny = ny + 1
        y1 = y0 + (ny-1)*dy

        if( x0 > xmin .or. x1 < xmax ) goto 95
        if( y0 > ymin .or. y1 < ymax ) goto 95

        call set_regpar(regpar,nx,ny,dx,dy,x0,y0,x1,y1,flag)

        return
   95   continue
        write(6,*) x0,xmin,xmax,x1
        write(6,*) y0,ymin,ymax,y1
        stop 'error stop create_reg: internal error (1)'
   96   continue
        write(6,*) 'dx,dy: ',dx,dy
        write(6,*) 'xmin,xmax: ',xmin,xmax
        write(6,*) 'ymin,ymax: ',ymin,ymax
        stop 'error stop create_reg: error in input parameters'
        end

!******************************************************

	subroutine create_rxry(regpar,nx,ny,rxx,ryy)

	implicit none

	real regpar(7)
	integer nx,ny
	real rxx(nx)
	real ryy(ny)

	integer nxx,nyy,ix,iy
	real x0,y0,dx,dy

	nxx = nint(regpar(1))
	nyy = nint(regpar(2))
	if( nx /= nxx .or. ny /= nyy ) goto 91
	x0 = regpar(3)
	y0 = regpar(4)
	dx = regpar(5)
	dy = regpar(6)

	do ix=1,nx
	  rxx(ix) = x0 + (ix-1)*dx
	end do
	  
	do iy=1,ny
	  ryy(iy) = y0 + (iy-1)*dy
	end do
	  
	return
   91   continue
	write(6,*) 'nx,nxx: ',nx,nxx
	write(6,*) 'ny,nyy: ',ny,nyy
	stop 'error stop create_rxry: mismatch between dimensions and regpar'
	end

!******************************************************

        subroutine set_regpar(regpar,nx,ny,dx,dy,x0,y0,x1,y1,flag)

        implicit none

        real regpar(7)
        integer nx,ny
        real dx,dy
        real x0,y0,x1,y1
        real flag

        regpar(1) = nx
        regpar(2) = ny
        regpar(3) = x0
        regpar(4) = y0
        regpar(5) = dx
        regpar(6) = dy
        regpar(7) = flag

	!x1,y1 are not used

        end

!******************************************************

        subroutine get_regpar(regpar,nx,ny,dx,dy,x0,y0,x1,y1,flag)

        implicit none

        real regpar(7)
        integer nx,ny
        real dx,dy
        real x0,y0,x1,y1
        real flag

        nx = nint(regpar(1))
        ny = nint(regpar(2))
        x0 = regpar(3)
        y0 = regpar(4)
        dx = regpar(5)
        dy = regpar(6)
        flag = regpar(7)

        x1 = x0 + (nx-1)*dx
        y1 = y0 + (ny-1)*dy

        end

!******************************************************

        subroutine get_regpar_dim(regpar,x0,y0,x1,y1)

! returns dimension of regular grid

        implicit none

        real regpar(7)
        real x0,y0,x1,y1

        integer nx,ny
        real dx,dy,flag

        nx = nint(regpar(1))
        ny = nint(regpar(2))
        x0 = regpar(3)
        y0 = regpar(4)
        dx = regpar(5)
        dy = regpar(6)

        x1 = x0 + (nx-1)*dx
        y1 = y0 + (ny-1)*dy

        end

!******************************************************

        subroutine print_regpar(regpar)

        implicit none

        real regpar(7)
        integer nx,ny
        real x0,y0,dx,dy,x1,y1
        real flag

        call get_regpar(regpar,nx,ny,dx,dy,x0,y0,x1,y1,flag)

        write(6,'(4x,a,2i16)') 'nx,ny: ',nx,ny
        write(6,'(4x,a,2f16.6)') 'x0,y0: ',x0,y0
        write(6,'(4x,a,2f16.6)') 'x1,y1: ',x1,y1
        write(6,'(4x,a,2f16.6)') 'dx,dy: ',dx,dy
        write(6,'(4x,a,2f16.6)') 'flag : ',flag

        end

!******************************************************

	subroutine print_matrix(nx,ny,areg,flag)

	implicit none

	integer nx,ny
	real areg(nx,ny)
	real flag

	integer ix,iy
	character*256 line
	character*1 c

	do iy=ny,1,-1
	  line = ' '
	  do ix=1,nx
	    if( areg(ix,iy) == flag ) then
	      c = '.'
	    else
	      c = '*'
	    end if
	    line(ix:ix) = c
	  end do
	  write(66,*) trim(line)
	end do

	end

!******************************************************

