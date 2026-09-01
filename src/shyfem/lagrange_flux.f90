
!--------------------------------------------------------------------------
!
!    Copyright (C) 2009-2012,2014-2016,2019  Georg Umgiesser
!    Copyright (C) 2012  Christian Ferrarin
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

! routines handling fluxes
!
! revision log :
!
! 05.02.2009	ggu	copied from other files
! 28.04.2009	ggu	links re-structured
! 23.03.2010	ggu	changed v6.1.1
! 20.10.2011	ggu	new routines for fluxes implemented
! 24.10.2011	ggu	3d routines implemented
! 04.11.2011	ggu	part with bsigma not yet finished
! 16.12.2011	ggu	bug fix - flux2d_aux was integer
! 24.01.2012	ggu	changed VERS_6_1_41
! 23.03.2012	ggu	bug fix - dst was not available in setup_vl_3d
! 30.03.2012	ggu	changed VERS_6_1_51
! 28.08.2012	ggu&ccf	bug fix - ok now for 3d surface tracking
! 19.12.2014	ggu	changed VERS_7_0_10
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 23.04.2015	ggu	changed VERS_7_1_8
! 17.07.2015	ggu	changed VERS_7_1_53
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 24.07.2015	ggu	changed VERS_7_1_82
! 16.11.2015	ggu	changed VERS_7_3_14
! 16.12.2015	ggu	changed VERS_7_3_16
! 19.02.2016	ggu	changed VERS_7_5_2
! 16.02.2019	ggu	changed VERS_7_5_60
!
!****************************************************************

	subroutine lagr_setup_timestep

! initializes length of element sides and fluxes

	implicit none

	call rprs
	call setup_fluxes_2d
	call setup_fluxes_3d	!this overwrites velocities in vel_ie (2d)

	end

!****************************************************************

	subroutine rprs

! initializes length of element sides

	use mod_lagrange
	use basin

	implicit none
	
	real x1,x2,y1,y2,dx,dy,ddx,ddy,d
	integer ie,k,i1,i2,p1,p2

	integer icall
	save icall
	data icall / 0 /

	if( icall .ne. 0 ) return

	icall = 1

	do ie=1,nel
	 do k=1,3
          i1=mod(k,3)+1
          i2=mod(i1,3)+1
	  p1=nen3v(i1,ie)
	  p2=nen3v(i2,ie)
	  x1=xgv(p1)
	  x2=xgv(p2)
	  y1=ygv(p1)
	  y2=ygv(p2)
	  dx=x2-x1
	  dy=y2-y1
	  ddx=dx**2
	  ddy=dy**2
	  d=sqrt(ddx+ddy)
	  dvert(k,ie)=d
	 end do
	end do

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine setup_fluxes_3d

! sets up fluxes in 3d - has to be done every time step

	use mod_rungekutta, only : get_rungekutta_idiag_weights
	use mod_lagrange
	use mod_geom
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer k,ie,ii,i
	integer nn,ne
	integer n
	integer l,lmax,lkmax
	integer itype
	real ai_ll
	real tdif
	
	integer elems(maxlnk)
	real rflux(nlvdi,maxlnk)       !fluxes across finite volume k
	real tflux(nlvdi,maxlnk)       !fluxes across sides of element
	
	integer flxtype

!	--------------------------------------------
!	initialization
!	--------------------------------------------

	call get_rungekutta_idiag_weights(1,ai_ll)
	azlgr = ai_ll		!store in include file

        do ie=1,nel
	  lmax = ilhv(ie)
          do ii=1,3
	    do l=1,lmax
              flux3d(l,ii,ie)=0
	    end do
          end do
        end do

!	--------------------------------------------
!	loop on nodes
!	--------------------------------------------

        do k=1,nkn
          itype=flxtype(k)

	  n = maxlnk
	  call get_elems_around(k,maxlnk,ne,elems)
	  call flx3d_k(k,itype,ai_ll,lkmax,n,rflux,ne,elems)
	  call make_fluxes_3d(k,itype,lkmax,n,rflux,tflux)
  	  call setup_flux3d(k,lkmax,n,tflux,ne,elems)
        end do

!	--------------------------------------------
!	compute velocities
!	--------------------------------------------

        call setup_vl_3d

!	--------------------------------------------
!	end of routine
!	--------------------------------------------

        end

!******************************************************************

	subroutine setup_fluxes_2d

! sets up fluxes - has to be done every time step

	use mod_rungekutta, only : get_rungekutta_idiag_weights
	use mod_lagrange
	use mod_geom
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer k,ie,ii,i
	integer nn,ne
	integer n,j
	integer itype
	real ai_ll
	real tdif
	
	integer elems(maxlnk)
	real rflux(maxlnk)       !fluxes across finite volume k
	real tflux(maxlnk)       !fluxes across sides of element
	real tflux_aux(maxlnk)   !fluxes across sides of element (aux)
	real flux2d_aux(3,nel)
	
	integer flxtype

        do ie=1,nel
          do ii=1,3
            flux2d(ii,ie)=0
            flux2d_aux(ii,ie)=0
          end do
        end do

	call get_rungekutta_idiag_weights(1,ai_ll)

        do k=1,nkn
          itype=flxtype(k)

!	  --------------------------------------------
!	  new way to set-up rflux,tflux
!	  --------------------------------------------

	  n = ngr
	  call get_elems_around(k,maxlnk,ne,elems)
	  call flx2d_k(k,itype,ai_ll,n,rflux,ne,elems)
	  call make_fluxes_2d(k,itype,n,rflux,tflux)
  	  call setup_flux2d(k,n,tflux,flux2d,ne,elems)
        end do

!	--------------------------------------------
!	error check ... delete
!	--------------------------------------------

!	--------------------------------------------
!	compute velocities
!	--------------------------------------------

        call setup_vl_2d

!	--------------------------------------------
!	end of routine
!	--------------------------------------------

        end

!******************************************************************
!******************************************************************
!******************************************************************

  	subroutine setup_flux3d(k,lkmax,n,tflux,ne,elems)

! computes 3d fluxes in element
!
! maybe not existing fluxes are set in flux3d
! we do not use lkmax, but lmax

	use mod_lagrange
	use mod_geom
	use levels

        implicit none

        integer k,lkmax,n
	real tflux(nlvdi,ne)      !fluxes across sides of element
        integer ne                !total number of elements in elems
        integer elems(ne)         !elements around k

        integer i,ie,l,lmax
        integer j,n1,n2

        integer ibhnd
        integer inext

        do i=1,ne
          ie=elems(i)
          n1=ibhnd(k,ie)
          n2=inext(k,ie)
          j=i+1
          if(i.eq.n) j=1
	  lmax = ilhv(ie)
	  do l=1,lmax
            flux3d(l,n2,ie)=flux3d(l,n2,ie)-tflux(l,j)
            flux3d(l,n1,ie)=flux3d(l,n1,ie)+tflux(l,i)      
	  end do
        end do

	end

!******************************************************************

  	subroutine setup_flux2d(k,n,tflux,flux2d_loc,ne,elems)

! computes fluxes in element

	use mod_lagrange
	use mod_geom

        implicit none

        integer k,n
        real tflux(ne)
	real flux2d_loc(3,ne)
        integer ne              !total number of elements in elems
        integer elems(ne)       !elements around k

        integer i,ie
        integer j,n1,n2

        integer ibhnd
        integer inext

        do i=1,ne
          ie=elems(i)
          n1=ibhnd(k,ie)
          n2=inext(k,ie)
          j=i+1
          if(i.eq.n) j=1
          flux2d_loc(n2,ie)=flux2d_loc(n2,ie)-tflux(j)
          flux2d_loc(n1,ie)=flux2d_loc(n1,ie)+tflux(i)      
        end do

	end

!*********************************************************************

        subroutine setup_vl_3d
	
! computes velocities in element

! use first layer thickness which already contains water level
!
! dH = dh + dz		dh undisturbed depth  dz average water level
!
! dz = (z1+z2+z3)/3.
!
! dp = dh + (z1+z2)/2 = dH - dz + (z1+z2)/2 
!		= dH + (3*z1+3*z2-2*z1-2*z2-2*z3)/6
!		= dH + (z1+z2-2*z3)/6
!
! all this has to be revised for sigma layers

	use mod_lagrange
	use mod_layer_thickness
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	logical bsigma
	integer ie,ii,i1,i2
	integer l,lmax
	real flx,dh,zi1,zi2,zi3,dp,ar,dst
	real bfact
	
	call get_bsigma(bsigma)

	bfact = 1.
	if( bback ) bfact = -1.

	do ie=1,nel
	  dh=hdenv(1,ie)		!this already contains water level
	  lmax = ilhv(ie)
	  do l=1,lmax
            do ii=1,3
	      dp = dh
	      dst=dvert(ii,ie)
	      if( .not. bsigma .and. l .eq. 1 ) then	!z-levels, first layer
                i1=mod(ii,3)+1
                i2=mod(i1,3)+1
	        zi1=zenv(i1,ie)
	        zi2=zenv(i2,ie)
	        zi3=zenv(ii,ie)
	        dp=dh+((zi1+zi2-2.*zi3)/6.)	!compute aver depth between 1&2
	      end if
	      ar=dp*dst				!section area
	      flx=flux3d(l,ii,ie)		!flux
	      vel3d_ie(l,ii,ie)=bfact*flx/ar	!velocity
	      !if( l .eq. 1 ) then
	      !  vel_ie(ii,ie)=vel3d_ie(l,ii,ie) !velocity first layer
	      !end if
	    end do
          end do
	  do ii=1,3
	    vel_ie(ii,ie)=vel3d_ie(1,ii,ie)	!velocity first layer
	  end do
	end do
	
	end

!*********************************************************************

        subroutine setup_vl_2d
	
! computes velocities in element

! use first layer thickness which already contains water level
!
! dH = dh + dz		dh undisturbed depth  dz average water level
!
! dz = (z1+z2+z3)/3.
!
! dp = dh + (z1+z2)/2 = dH - dz + (z1+z2)/2 
!		= dH + (3*z1+3*z2-2*z1-2*z2-2*z3)/6
!		= dH + (z1+z2-2*z3)/6
!
! all this has to be revised for sigma layers

	use mod_lagrange
	use mod_layer_thickness
	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer ie,ii,i1,i2
	real flx,dh,zi1,zi2,zi3,dp,ar,dst
	real bfact,eps

	eps=1.e-5
	bfact = 1.
	if( bback ) bfact = -1.

	do ie=1,nel
	 dh=hdenv(1,ie)			!this already contains water level
         do ii=1,3
	  dst=dvert(ii,ie)
          i1=mod(ii,3)+1
          i2=mod(i1,3)+1
	  zi1=zenv(i1,ie)
	  zi2=zenv(i2,ie)
	  zi3=zenv(ii,ie)
	  dp=dh+((zi1+zi2-2.*zi3)/6.)	!compute aver depth between 1&2
	  ar=dp*dst			!section area
	  if(ar.le.eps) ar = eps	!area keep above 0
	  flx=flux2d(ii,ie)		!flux
	  vel_ie(ii,ie)=bfact*flx/ar	!velocity
         end do
	end do
	
	end

!**********************************************************************

