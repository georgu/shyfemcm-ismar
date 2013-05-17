c
c $Id: ousutil.f,v 1.3 2009-09-14 08:20:58 georg Exp $
c
c utilities for OUS files
c
c revision log :
c
c 16.12.2010	ggu	copied from ousextr_gis.f
c 03.06.2011	ggu	some routines transfered to genutil.f
c 08.06.2011	ggu	new routine transp2nodes()
c 10.11.2011    ggu     new routines for hybrid levels
c 02.12.2011    ggu     bug fix for call to get_sigma_info() (missing argument)
c 21.01.2013    ggu     added two new routines comp_vel2d, comp_barotropic
c
c******************************************************************

        subroutine transp2vel(nel,nkn,nlv,nlvdim,hev,zenv,nen3v
     +				,ilhv,hlv,utlnv,vtlnv
     +                          ,uprv,vprv,weight,hl)

c transforms transports at elements to velocities at nodes

        implicit none

	include 'ev.h'

        integer nel
        integer nkn
	integer nlv
        integer nlvdim
        real hev(1)
        real zenv(3,1)
	integer nen3v(3,1)
	integer ilhv(1)
	real hlv(1)
        real utlnv(nlvdim,1)
        real vtlnv(nlvdim,1)
        real uprv(nlvdim,1)
        real vprv(nlvdim,1)
        real weight(nlvdim,1)		!aux variable for weights
	real hl(1)			!aux variable for real level thickness

	logical bsigma,bzeta
        integer ie,ii,k,l,lmax,nsigma,nlvaux
        real hmed,u,v,area
	real hsigma

	bzeta = .true.		!use zeta for depth computation

	call get_sigma_info(nlvaux,nsigma,hsigma)
	if( nlvaux .gt. nlvdim ) stop 'error stop transp2vel: nlvdim'
	bsigma = nsigma .gt. 0

	do k=1,nkn
	  do l=1,nlv
	    weight(l,k) = 0.
	    uprv(l,k) = 0.
	    vprv(l,k) = 0.
	  end do
	end do
	      
        do ie=1,nel

	  area = 12. * ev(10,ie)
	  lmax = ilhv(ie)
	  call get_layer_thickness_e(ie,lmax,bzeta,nsigma,hsigma,hl)

	  do l=1,lmax
	    hmed = hl(l)
	    u = utlnv(l,ie) / hmed
	    v = vtlnv(l,ie) / hmed
	    do ii=1,3
	      k = nen3v(ii,ie)
	      uprv(l,k) = uprv(l,k) + area * u
	      vprv(l,k) = vprv(l,k) + area * v
	      weight(l,k) = weight(l,k) + area
	    end do
	  end do
	end do

	do k=1,nkn
	  do l=1,nlv
	    area = weight(l,k)
	    if( area .gt. 0. ) then
	      uprv(l,k) = uprv(l,k) / area
	      vprv(l,k) = vprv(l,k) / area
	    end if
	  end do
	end do
	      
	end

c******************************************************************

        subroutine transp2nodes(nel,nkn,nlv,nlvdim,hev,zenv,nen3v
     +				,ilhv,hlv,utlnv,vtlnv
     +                          ,utprv,vtprv,weight)

c transforms transports at elements to transports at nodes

        implicit none

        integer nel
        integer nkn
	integer nlv
        integer nlvdim
        real hev(1)
        real zenv(3,1)
	integer nen3v(3,1)
	integer ilhv(1)
	real hlv(1)
        real utlnv(nlvdim,1)
        real vtlnv(nlvdim,1)
        real utprv(nlvdim,1)
        real vtprv(nlvdim,1)
        real weight(nlvdim,1)

        integer ie,ii,k,l,lmax
        real u,v,w

	do k=1,nkn
	  do l=1,nlv
	    weight(l,k) = 0.
	    utprv(l,k) = 0.
	    vtprv(l,k) = 0.
	  end do
	end do
	      
        do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    u = utlnv(l,ie)
	    v = vtlnv(l,ie)
	    do ii=1,3
	      k = nen3v(ii,ie)
	      utprv(l,k) = utprv(l,k) + u
	      vtprv(l,k) = vtprv(l,k) + v
	      weight(l,k) = weight(l,k) + 1.
	    end do
	  end do
	end do

	do k=1,nkn
	  do l=1,nlv
	    w = weight(l,k)
	    if( w .gt. 0. ) then
	      utprv(l,k) = utprv(l,k) / w
	      vtprv(l,k) = vtprv(l,k) / w
	    end if
	  end do
	end do
	      
	end

c***************************************************************

        subroutine comp_vel2d(nel,hev,zenv,ut2v,vt2v,u2v,v2v
     +				,umin,vmin,umax,vmax)

c computes 2D velocities from 2D transports - returns result in u2v,v2v

        implicit none

        integer nel
        real hev(1)
        real zenv(3,1)
        real ut2v(1)
        real vt2v(1)
	real u2v(1), v2v(1)
        real umin,vmin
        real umax,vmax

        integer ie,ii
        real zmed,hmed,u,v

	umin = +1.e+30
	vmin = +1.e+30
	umax = -1.e+30
	vmax = -1.e+30

        do ie=1,nel
          zmed = 0.
          do ii=1,3
            zmed = zmed + zenv(ii,ie)
          end do
          zmed = zmed / 3.
          hmed = hev(ie) + zmed

          u = ut2v(ie) / hmed
          v = vt2v(ie) / hmed

	  u2v(ie) = u
	  v2v(ie) = v

          umin = min(umin,u)
          vmin = min(vmin,v)
          umax = max(umax,u)
          vmax = max(vmax,v)
        end do

        end

c***************************************************************

	subroutine comp_barotropic(nel,nlvdim,ilhv
     +			,utlnv,vtlnv,ut2v,vt2v)

c computes barotropic transport

	implicit none

	integer nel,nlvdim
	integer ilhv(1)
	real utlnv(nlvdim,1)
	real vtlnv(nlvdim,1)
	real ut2v(1)
	real vt2v(1)

	integer ie,l,lmax
	real utot,vtot

	do ie=1,nel
	  lmax = ilhv(ie)
	  utot = 0.
	  vtot = 0.
	  do l=1,lmax
	    utot = utot + utlnv(l,ie)
	    vtot = vtot + vtlnv(l,ie)
	  end do
	  ut2v(ie) = utot
	  vt2v(ie) = vtot
	end do

	end

c***************************************************************

	subroutine compute_volume(nel,zenv,hev,volume)

	implicit none

	include 'param.h'
	include 'ev.h'

	integer nel
	real zenv(3,neldim)
	real hev(neldim)
	real volume

	integer ie,ii
	real zav,area
	double precision vol,voltot,areatot

	voltot = 0.
	areatot = 0.

	do ie=1,nel
	  zav = 0.
	  do ii=1,3
	    zav = zav + zenv(ii,ie)
	  end do
	  area = 12. * ev(10,ie)
	  vol = area * (hev(ie) + zav/3.)
	  voltot = voltot + vol
	  !areatot = areatot + area
	end do

	volume = voltot

	end

c***************************************************************

        subroutine debug_write_node(ks,it,nrec
     +		,nkndim,neldim,nlvdim,nkn,nel,nlv
     +          ,nen3v,zenv,znv,utlnv,vtlnv)

c debug write

        implicit none

	integer ks	!internal node number to output (0 for none)
        integer it,nrec
        integer nkndim,neldim,nlvdim,nkn,nel,nlv
        integer nen3v(3,neldim)
        real znv(nkndim)
        real zenv(3,neldim)
        real utlnv(nlvdim,neldim)
        real vtlnv(nlvdim,neldim)

        integer ie,ii,k,l
        logical bk

	if( ks .le. 0 ) return

        write(66,*) 'time: ',it,nrec
        write(66,*) 'kkk: ',ks,znv(ks)

        do ie=1,nel
          bk = .false.
          do ii=1,3
            k = nen3v(ii,ie)
            if( k .eq. ks ) then
              write(66,*) 'ii: ',ii,ie,zenv(ii,ie)
              bk = .true.
            end if
          end do
          if( bk ) then
          do l=1,nlv
            write(66,*) 'ie: ',ie,l,utlnv(l,ie),vtlnv(l,ie)
          end do
          end if
        end do

        end

c***************************************************************
c***************************************************************
c***************************************************************

        subroutine open_ous_file(name,status,nunit)

        implicit none

        character*(*) name,status
        integer nunit

        integer nb
        character*80 file

        integer ifileo

        call mkname(' ',name,'.ous',file)
        nb = ifileo(0,file,'unform',status)
        if( nb .le. 0 ) stop 'error stop open_ous_file: opening file'

        nunit = nb

        end

c***************************************************************

        subroutine qopen_ous_file(text,status,nunit)

c asks for name and opens ous file

        implicit none

        character*(*) text,status
        integer nunit

        character*80 name

        write(6,*) text
        read(5,'(a)') name
        write(6,*) name

        call open_ous_file(name,status,nunit)

        end

c***************************************************************
c***************************************************************
c***************************************************************

