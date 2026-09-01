
!--------------------------------------------------------------------------
!
!    Copyright (C) 1992,1994,1997-1998,2001,2003-2004  Georg Umgiesser
!    Copyright (C) 2006,2010-2019  Georg Umgiesser
!    Copyright (C) 2012  Debora Bellafiore
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

! routines handling flooding and drying
!
! contents :
!
! subroutine setweg(iweich,iw,hzmin)	sets array iwegv
! subroutine setuvd			sets velocities in dry areas
! subroutine setzev			sets array zenv from znv
! subroutine setznv			sets array znv from zenv
! subroutine zuniq(zv,av)		makes z values unique
!
! revision log :
!
! 01.07.1992	ggu	$$lump  - lumping of matrix
! 05.08.1992	ggu	$$ibtyp3 - implementation of ibtyp=3
! 01.09.1992	ggu	$$eps  - introduction of eps (setuvd)
! 12.01.1994	ggu	$$hzon  - use new variable hzon (setweg)
! 12.01.1994	ggu	$$eps0  - use eps only in last control (setuvd)
! 12.01.1994	ggu	$$99  - do not jump to 99 in loop (setuvd)
! 05.02.1994	ggu	$$azpar - use az to compute velocities (setuvd)
! 04.03.1994	ggu	$$azuvdry - one az too much in formula (setuvd)
! 27.10.1997	ggu	$$isum - better identification of error 99 (setuvd)
! 27.10.1997	ggu	$$dpisum - use double prec. for key values (setuvd)
! 27.03.1998	ggu	eliminated /bnd/, /irv/
! 27.04.1998	ggu	$$NKNEL - do not call nknel in tstlnk
! 08.05.1998	ggu	new routine pntfla -> absolute element index
! 20.05.1998	ggu	hard coded unit 88 substituted (use ifileo)
! 14.07.1998	ggu	$$ibtyp4 - boundary type 4 integrated
! 21.08.1998	ggu	file sublnk splitted into lnk/dry
! 21.08.1998	ggu	routine setczg removed
! 21.08.1998	ggu	xv eliminated
! 21.11.2001	ggu	extra bdebug in setweg
! 21.11.2001	ggu	more debug information in setuvd
! 13.08.2003	ggu	new routine setznv
! 03.09.2004	ggu	setznv: do not stop if znv is not unique (restart)
! 22.09.2004	ggu	debug in setweg()
! 23.03.2006	ggu	changed time step to real
! 23.03.2010	ggu	changed v6.1.1
! 20.05.2011	ggu	different algorithm for element removal (wet&dry)
! 20.05.2011	ggu	new routines set_dry() and set_element_dry(ie)
! 31.05.2011	ggu	changed VERS_6_1_23
! 07.06.2011	ggu	slight changes in wetting and drying
! 14.07.2011	ggu	changed VERS_6_1_27
! 27.01.2012	dbf&ggu	changes for sigma coordinates (bsigma)
! 30.03.2012	ggu	changed VERS_6_1_51
! 12.09.2013	ggu	changed VERS_6_1_67
! 12.12.2014	ggu	changed VERS_7_0_9
! 19.12.2014	ggu	changed VERS_7_0_10
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 05.05.2015	ggu	changed VERS_7_1_10
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 23.09.2015	ggu	changed VERS_7_2_4
! 28.04.2016	ggu	changed VERS_7_5_9
! 05.12.2017	ggu	changed VERS_7_5_39
! 19.04.2018	ggu	changed VERS_7_5_45
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 21.05.2019	ggu	changed VERS_7_5_62
! 06.11.2019	ggu	eliminated femtime
! 23.05.2023	ggu	in setznv() loop over ie_mpi
! 07.06.2023	ggu	new routine compute_dry_elements()
! 24.04.2025	ggu	write out dry elements and nodes
!
!*****************************************************************

	subroutine set_dry

! sets dry elements

	implicit none

	integer iwa

        call setweg(2,iwa)
        call setweg(3,iwa)

	end

!*****************************************************************

	subroutine set_element_dry(ie)

! sets dry elements

	use mod_geom_dynamic

	implicit none

	integer ie

	if( iwegv(ie) .eq. 0 ) then
	  iwegv(ie) = 3
	  iwetv(ie) = 0
	end if
	
	end

!*****************************************************************
!
        subroutine setweg(iweich,iw)
!
! sets array iwegv
!
! iweich 	flag
!	-1: initialize  
!	 0: first call  
!	 1: only take away (hzmin)
!	 2: only take away (hzoff)    
!	 3: only add
! iw		if on return iw != 0 iwegv has changed
!
! 1 is used to check, if iteration has to be repeated
! 2 is just to switch off an element in time, so no iteration
!	has to be repeated
! generally it is true that : 0 <= hzmin <= zmin <= hzoff <= hzon
!
!				in			out
!			   element was inside	   element was outside
!	+---------------------------------------------------------------+
!  h4	|		|	-		|	include		|
!	+- hzon  -------------------------------------------------------+
!  h3	|		|	-		|	-		|
!	+- hzoff -------------------------------------------------------+
!  h2	|		|	exclude		|	-		|
!	+- hzmin -------------------------------------------------------+
!  h1	|		|    exclude, repeat	|	error		|
!	+- zero  -------------------------------------------------------+
!
! iwegv   0:all nodes wet   >0:number of nodes dry -> out of system
!
	use mod_geom_dynamic
	use mod_hydro
	use evgeom
	use basin

        implicit none

! arguments
        integer iweich,iw
! local
        integer ie,ii,iwh,iweg,k,iu
        integer iespec,iwait,iwet
	integer nlv,nsigma
	real hsigma
        real hzg,hzmin,hzoff,hzon,volmin,aomega,zmin
	real hzlim,hztot
	character*10 text
! functions
        real getpar
!	integer ieint,iround,ipint
! aux
	logical debug,bnodry
        logical bdebug,bbdebug
        logical bnewalg,binclude
	logical bnewtime
	logical bsigma

	double precision dtime
	double precision, save :: dtime_old = -1

	debug=.false.
	bnodry=.true.		!drying is not allowd (for debug)
	bnodry=.false.
        bdebug=.true.
        bdebug=.false.

        bnewalg=.false.
        bnewalg=.true.		!new algorithm to include elements

	call get_act_dtime(dtime)
	bnewtime = dtime .ne. dtime_old

	iu = 156
        iespec = 4574
        iespec = 19417
        iespec = -1
	bbdebug = .true.
	bbdebug = .false.
	iwait = 5		!wait so long before including

	call get_sigma_info(nlv,nsigma,hsigma)
	bsigma = nsigma .gt. 0

        hzmin=getpar('hzmin')
        hzoff=getpar('hzoff')
        hzon=getpar('hzon')			!$$hzon
        volmin=getpar('volmin')

        iwh = 0
	iw = 0

	if( bsigma .and. iweich .ge. 1 ) return

        if(iweich.eq.-1) then                   !initialize iwegv
          do ie=1,nel
            iwegv(ie)=0
            iwetv(ie)=0
          end do
        else if(iweich.eq.0) then               !first call
          do ie=1,nel
            aomega=ev(10,ie)
            zmin=hzmin+volmin/(4.*aomega)
            if(zmin.gt.hzoff) zmin=hzmin+0.5*(hzmin+hzoff) !just in case...

            iweg=0
            do ii=1,3
              hzg = zenv(ii,ie) + hm3v(ii,ie)
              if(hzg.lt.hzoff) iweg=iweg+1
              if(hzg.lt.zmin) then
		zenv(ii,ie)=zmin-hm3v(ii,ie)
	        if( bsigma ) znv(nen3v(ii,ie)) = zenv(ii,ie)
	      end if
            end do
	    if( bsigma ) iweg = 0
            iwegv(ie)=iweg
	    if( iweg .eq. 0 ) iwetv(ie) = 1
	    if( iweg > 0 ) iwh = iwh + 1
          end do
        else if(iweich.eq.1.or.iweich.eq.2) then  !only take away
	  if( iweich .eq. 1 ) then
	    hzlim = hzmin
	    text = 'setweg 1: '
	  else
	    hzlim = hzoff
	    text = 'setweg 2: '
	  end if
          do ie=1,nel
	   if( iwegv(ie) .eq. 0 ) then
            debug = ie .eq. iespec
            iweg=0
            do ii=1,3
              hzg = znv(nen3v(ii,ie)) + hm3v(ii,ie)
              if(hzg.lt.hzlim) iweg=iweg+1
              if(debug .or. bdebug.and.hzg.lt.hzlim) then
                      write(iu,*) text,ie,iweg,hzg
                      write(iu,*) znv(nen3v(ii,ie)),hm3v(ii,ie)
              end if
            end do
            if(iweg.gt.0) then		!case h1/h2
	      if( bnodry ) goto 77
              !if(iwegv(ie).gt.0 .and.iweich.eq.1) goto 78	case h1out
              if(iwegv(ie).eq.0) then	!element was inside - case h1/h2in
		iwh=iwh+1
                iwetv(ie)=0
	      end if
              iwegv(ie)=iweg
            end if
	   end if
          end do
        else if(iweich.eq.3) then               !only add
          do ie=1,nel
	   if( iwegv(ie) .gt. 0 ) then
            debug = ie .eq. iespec
            iweg=0
	    hztot = 0.
            do ii=1,3
              hzg = znv(nen3v(ii,ie)) + hm3v(ii,ie)
	      hztot = hztot + hzg
              if(hzg.lt.hzon) iweg=iweg+1 !$$hzon ...(220792)
              if(debug .or. bdebug.and.hzg.lt.hzon) then
		      k = nen3v(ii,ie)
                      write(iu,*) 'setweg 3: ',ie,ii,k,iweg,hzg
                      write(iu,*) znv(nen3v(ii,ie)),hm3v(ii,ie)
              end if
            end do
	    hztot = hztot / 3.
	    !binclude = hztot .ge. hzon
	    binclude = hztot .ge. hzon .and. -iwetv(ie) .gt. iwait
! %%%%%%%%% this is wrong -> test also for iweg > 0	!FIXME
	    if( bnewalg .and. binclude ) then		!new algorithm
	      iweg = 0
              if(bdebug) then
                      write(iu,*) 'setweg 3a: ',ie,iweg,hztot
	      end if
	    end if
            if(iweg.eq.0) then		!case h4
              if(iwegv(ie).gt.0) then	!element was out - case h4out
		iwh=iwh+1
                iwetv(ie)=0
	      end if
              iwegv(ie)=0
            end if
	   end if
          end do
        end if

	do ie=1,nel
	  iweg = iwegv(ie)
	  iwet = iwetv(ie)
	  if( bnewtime .or. iwet .eq. 0 ) then	!new time step or changed
	    if( iweg .eq. 0 ) then
	      iwetv(ie) = iwetv(ie) + 1
	    else
	      iwetv(ie) = iwetv(ie) - 1
	    end if
	  end if
	end do

	dtime_old = dtime
        iw=iwh

	if( bbdebug ) then
	  !write(iu,*) '++++ ',iweich,iw,iwegv(19417),iwegv(19428)
	end if

        return
   77	continue
	write(6,*) 'drying is not allowed...'
	write(6,*) iweich,ie,iwh,iwegv(ie)
	write(6,*) hzmin,hzoff,hzon,volmin
        do ii=1,3
	  k = nen3v(ii,ie)
          hzg = znv(k) + hm3v(ii,ie)
	  write(6,*) ii,k,znv(k),hzg
	  write(6,*) '     ',zenv(ii,ie),hm3v(ii,ie)
	end do
	stop 'error stop setweg'
        end

!****************************************************************
!
        subroutine setuvd(curr_stage,coeff_irk)
!
! sets velocities in dry areas
!
! ie    element
! dt    time step
! hzmin smallest z allowed
! b,c   form functions
!
	use mod_rungekutta, only: n_rkstages,urk_reg,vrk_reg
	use mod_geom_dynamic
	use mod_hydro_baro
	use mod_hydro
	use evgeom
	use levels
	use basin

        implicit none

! arguments
	integer curr_stage
	real coeff_irk(n_rkstages+1)
! local
        integer ie,ii,i1,i2,isum,itot,l,js
        integer jlevel,ilevel
        integer i3,i4,i5,i6,i7,i8
!        real zm,zmed,d1,d2,det
        double precision zm,zmed,d1,d2,det	!$$dpisum
        real adt,axdt,dt,hzmin
        real z(3),b(3),c(3),uo,vo,u,v,zz
        real up(curr_stage-1),vp(curr_stage-1)
	real zn(3)
        integer nnn,itmin
	double precision dtime
	character*20 aline
! save
        real eps
        save eps
        data eps /1.e-5/
! functions
        real getpar
        integer ieext
        logical iseout
        iseout(ie) = iwegv(ie).gt.0
!
        nnn=-802 	!<0 if not needed
!	itmin=1550100	!0 if not needed
	itmin=0		!0 if not needed
!
	call get_timestep(dt)
        hzmin=getpar('hzmin')
!
        i7=0            !mean value is too low for node
          i8=0          !
            i3=0
              i4=0
                i5=0    !two nodes in element with mean too low
                  i6=0  !node is drying with computed velocity

        do ie=1,nel
        if( iseout(ie) ) then
!
	ilevel=ilhv(ie)
	jlevel=jlhv(ie)

	uo = 0.		!use barotropic velocities
	vo = 0.
	do l=jlevel,ilevel
          uo = uo + utlcv(l,ie)
          vo = vo + vtlcv(l,ie)
	end do
	up = 0.
	vp = 0.
	do js=1,curr_stage-1
	  do l=jlevel,ilevel
	    up(js) = up(js) + urk_reg(l,ie,js)
	    vp(js) = vp(js) + vrk_reg(l,ie,js)
	  end do
	end do
	axdt=3.*dt    !$$lump
        adt=1./axdt   !$$lump
!
! z average, set b,c
!
        zm=0.
        do ii=1,3
	  zn(ii) = zenv(ii,ie)
          zm=zm+zn(ii)
          b(ii)=ev(3+ii,ie)
          c(ii)=ev(6+ii,ie)
        end do
        zmed=zm/3.
!
! compute final z value for the nodes
!
        itot=0
        isum=0
!        zm=0.
        do ii=1,3
          if(zmed+hm3v(ii,ie).lt.hzmin) then  !$$eps0
            z(ii)=hzmin-hm3v(ii,ie)
            itot=itot+1
            isum=isum+ii
            i7=i7+1
          else
            z(ii)=zmed
          end if
          zm=zm-z(ii)
        end do
!
        if(itot.eq.0) then
          !ok
        else if(itot.eq.1) then
          i1=mod(isum,3)+1
          i2=mod(i1,3)+1
          if(z(i1)+zm*0.5+hm3v(i1,ie).lt.hzmin) then  !$$eps0
            i8=i8+1
            z(i1)=hzmin-hm3v(i1,ie)
            z(i2)=3.*zmed-z(i1)-z(isum)
!            if(z(i2)+hm3v(i2,ie).lt.hzmin-eps) goto 99 !$$99
          else if(z(i2)+zm*0.5+hm3v(i2,ie).lt.hzmin) then !$$eps0
            i3=i3+1
            z(i2)=hzmin-hm3v(i2,ie)
            z(i1)=3.*zmed-z(i2)-z(isum)
!            if(z(i1)+hm3v(i1,ie).lt.hzmin-eps) goto 99  !$$99
          else
            i4=i4+1
            z(i1)=z(i1)+0.5*zm
            z(i2)=z(i2)+0.5*zm
          end if
        else if(itot.eq.2) then
          i5=i5+1
          i1=6-isum
          z(i1)=z(i1)+zm
!          if(z(i1)+hm3v(i1,ie).lt.hzmin-eps) goto 99  !$$99
        else
          goto 99
        end if
!
! control, leave in any case
!
        isum=-1
        zm=0.
        do ii=1,3
          zm=zm+z(ii)
          if(z(ii)+hm3v(ii,ie).lt.hzmin-eps) goto 99
        end do
        isum=-2		!$$isum
        if(abs(zm-3.*zmed).gt.eps) goto 99
!
! now compute velocities
!
        d1 = (z(1)-zn(1))*adt &
     &    - coeff_irk(curr_stage)*( b(1)*uo + c(1)*vo )
        d2 = (z(2)-zn(2))*adt &
     &    - coeff_irk(curr_stage)*( b(2)*uo + c(2)*vo )
	do js=1,curr_stage-1
          d1 = d1 &
     &      - coeff_irk(js)*( b(1)*up(js) + c(1)*vp(js) )
	  d2 = d2 &
     &	    - coeff_irk(js)*( b(2)*up(js) + c(2)*vp(js) )
	end do
        det=1./(b(1)*c(2)-b(2)*c(1))
!
        u = det * ( c(2)*d1 - c(1)*d2 )
        v = det * (-b(2)*d1 + b(1)*d2 )
!
! with this velocity is there some node drying out ?
!
	isum=-1
        itot=0
        do ii=1,3
          zz = axdt * ( &
     &      b(ii)*(coeff_irk(curr_stage)*uo+coeff_irk(curr_stage+1)*u) +  &
     &      c(ii)*(coeff_irk(curr_stage)*vo+coeff_irk(curr_stage+1)*v) )  &
     &                 + zn(ii)
	  do js=1,curr_stage-1
            zz = zz + axdt * coeff_irk(js) * ( b(ii)*up(js)+ c(ii)*vp(js) )
	  end do
          if(zz+hm3v(ii,ie).lt.hzmin-eps) itot=itot+1 !$$eps
        end do
!
        if(itot.gt.0) then    !node is drying, set next z to above values
	  isum=-2
          i6=i6+1
	  if( coeff_irk(curr_stage+1) == 0. ) then
	    write(6,*) 'drying with explicit scheme is not possible'
	    write(6,*) 'use an implicit run if you have wetting and drying'
	    stop 'error stop setuvd: coeff_irk(curr_stage+1)'
	  end if

          d1 = (z(1)-zn(1))*adt &
     &      - coeff_irk(curr_stage)*( b(1)*uo + c(1)*vo )
          d2 = (z(2)-zn(2))*adt &
     &      - coeff_irk(curr_stage)*( b(2)*uo + c(2)*vo )
	  do js=1,curr_stage-1
	    d1 = d1 - coeff_irk(js)*( b(1)*up(js) + c(1)*vp(js) )
	    d2 = d2 - coeff_irk(js)*( b(2)*up(js) + c(2)*vp(js) )
	  end do
          det=1./( coeff_irk(curr_stage+1) * (b(1)*c(2)-b(2)*c(1)) )
!
! the formula should be   det=1./( a_ll * (b(1)*c(2)-b(2)*c(1)) )
!
          u = det * ( c(2)*d1 - c(1)*d2 )
          v = det * (-b(2)*d1 + b(1)*d2 )
        end if
!
! now set u/v/z
!
        itot=0
        zm=0.          !only for control
        do ii=1,3
          zz = axdt * ( &
     &      b(ii)*(coeff_irk(curr_stage)*uo+coeff_irk(curr_stage+1)*u) +  &
     &      c(ii)*(coeff_irk(curr_stage)*vo+coeff_irk(curr_stage+1)*v) )  &
     &                 + zn(ii)
	  do js=1,curr_stage-1
	    zz = zz + axdt * coeff_irk(js) * ( b(ii)*up(js)+ c(ii)*vp(js) )
	  end do
          if(zz+hm3v(ii,ie).lt.hzmin-eps) itot=itot+1 !$$eps
          zenv(ii,ie)=zz
          zm=zm+zz
        end do
        isum=-3		!$$isum
        if(abs(zm-3.*zmed).gt.eps) goto 99
!
        if(itot.gt.0) goto 97

        unv(ie)=u	!put back to new barotropic velocities
        vnv(ie)=v
!
        end if
        end do

!        write(88,*) i7,i8,i3,i4,i5,i6
!
        return

   99   continue  !use isum to decide which branch has been taken
	call get_act_dtime(dtime)
	call get_timeline(dtime,aline)
        write(6,*) '---------------------'
	write(6,*) 'error log setuvd (z computation)'
        write(6,*) 'time: ',dtime,'  ',aline
        write(6,*) 'element number (int/ext): ',ie,ieext(ie)
        write(6,*) ie,itot,isum
        write(6,*) zmed,zm,hzmin
        write(6,*) (hm3v(ii,ie),ii=1,3)
        write(6,*) (zenv(ii,ie),ii=1,3)
        write(6,*) (z(ii),ii=1,3)
        write(6,*) '---------------------'
	call check_set_unit(6)
	call check_elem(ie)
	call check_nodes_in_elem(ie)
        write(6,*) '---------------------'
        stop 'error stop setuvd : z computation'
   97   continue
	call get_act_dtime(dtime)
	call get_timeline(dtime,aline)
        write(6,*) '---------------------'
	write(6,*) 'error log setuvd (u/v computation)'
        write(6,*) 'time: ',dtime,'  ',aline
        write(6,*) 'element number (int/ext): ',ie,ieext(ie)
        write(6,*) ie,itot,isum
        write(6,*) zmed,zm,hzmin
        write(6,*) uo,vo,u,v
        write(6,*) (hm3v(ii,ie),ii=1,3)
        write(6,*) (zenv(ii,ie),ii=1,3)
        write(6,*) (z(ii),ii=1,3)
        write(6,*) (zn(ii),ii=1,3)
        write(6,*) '---------------------'
	call check_set_unit(6)
	call check_elem(ie)
	call check_nodes_in_elem(ie)
        write(6,*) '---------------------'
        stop 'error stop setuvd : u/v computation'
        end
!
!****************************************************************

        subroutine setzev

! sets array zenv from znv

	use mod_geom_dynamic
	use mod_hydro
	use basin

        implicit none

! local
        integer ie,ii

        do ie=1,nel
          if(iwegv(ie).eq.0) then !element is in system
            do ii=1,3
              zenv(ii,ie)=znv(nen3v(ii,ie))
            end do
          end if
	end do

        end

!****************************************************************

        subroutine setznv

! sets array znv from zenv

	use mod_geom_dynamic
	use mod_hydro
	use evgeom
	use basin
	use shympi
	use mkonst

        implicit none

! local
        integer ie,ii,k,ie_mpi
        integer ntot
	real z,area
	real v1v(nkn),v2v(nkn)

!-------------------------------------------------------------
! initialize znv and counters
!-------------------------------------------------------------

        ntot = 0

	znv = flag
	v1v = 0.
	v2v = 0.

!-------------------------------------------------------------
! set znv and accumulate
!-------------------------------------------------------------

        do ie_mpi=1,nel
	  ie = ip_sort_elem(ie_mpi)
          if( iwegv(ie) .eq. 0 ) then		!element is in system
            do ii=1,3
	      k = nen3v(ii,ie)
              z = zenv(ii,ie)
	      if( znv(k) .ne. flag .and. znv(k) .ne. z ) then   !restart
                ntot = ntot + 1
                write(6,*) 'n,ie,ii,k,z,znv(k) ',ntot,ie,ii,k,z,znv(k)
                z=max(z,znv(k))         !using higher value
              end if
	      znv(k) = z
            end do
	  else					!out of system
	    area = 4. * ev(10,ie)
            do ii=1,3
	      k = nen3v(ii,ie)
              z = zenv(ii,ie)
	      v1v(k) = v1v(k) + z * area
	      v2v(k) = v2v(k) + area
            end do
          end if
	end do

!-------------------------------------------------------------
! compute znv for dry areas
!-------------------------------------------------------------

        !call shympi_comment('shympi_elem: exchange v1v, v2v')
        call shympi_exchange_and_sum_2d_nodes(v1v)
        call shympi_exchange_and_sum_2d_nodes(v2v)

	do k=1,nkn
	  if( znv(k) .eq. flag ) then		!out of system
	    znv(k) = v1v(k) / v2v(k)
	  end if
	end do

	!call shympi_comment('exchanging znv in setznv ')
	call shympi_exchange_2d_node(znv)
	!call shympi_barrier

!-------------------------------------------------------------
! write debug status
!-------------------------------------------------------------

        if( ntot .gt. 0 ) then
          write(6,*) ntot, ' nodes with different z value'
        end if

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	return
   99	continue
	write(6,*) 'ie,ii,k,z,znv(k) ',ie,ii,k,z,znv(k)
	stop 'error stop setznv: nodal value not unique'
        end
!
!****************************************************************
!
        subroutine zuniq(zv,av)
!
! makes z values in dynamic system unique
!
! works only for lumped mass matrix
!
! zv    aux vector for z value
! av    aux vector for weighting factors (areas)
!
	use mod_geom_dynamic
	use mod_hydro
	use evgeom
	use basin

        implicit none
!
! arguments
        real zv(1),av(1)
! local
        integer ie,i,k
        real aomega
! functions
        logical isein
        isein(ie) = iwegv(ie).eq.0

! initialize aux vectors

        do k=1,nkn
          zv(k)=0.
          av(k)=0.
        end do

! accumulate contributions

        do ie=1,nel
          if(isein(ie)) then
            aomega=ev(10,ie)
            do i=1,3
              k=nen3v(i,ie)
              zv(k)=zv(k)+aomega*zenv(i,ie)
              av(k)=av(k)+aomega
            end do
          end if
        end do

! scaling of z values with weighting functions

        do k=1,nkn
          if(av(k).gt.0.) then
            zv(k)=zv(k)/av(k)
          end if
        end do

! write back to original vector zenv

        do ie=1,nel
          if(isein(ie)) then
            do i=1,3
              zenv(i,ie)=zv(nen3v(i,ie))
            end do
          end if
        end do

        return
        end
!
!*****************************************************************

        subroutine compute_dry_elements(iloop)

! computes total number of dry elements

	use mod_geom_dynamic
	use basin
	use shympi
	use mod_info_output

        implicit none

	integer iloop

        integer ie,k,idry,iedry,ikdry
	integer, save :: iuinfo = 0
	real array(3)
	character*80 format

	if( .not. binfo ) return

	idry = 0
        do ie=1,nel_unique
          if( iwegv(ie) /= 0 ) idry = idry + 1
	end do

	iedry = shympi_sum(idry)

	idry = 0
	do k=1,nkn_unique
	  if( inodv(k) == -2 ) idry = idry + 1
	end do

	ikdry = shympi_sum(idry)

        !if( iuinfo .le. 0 ) call getinfo(iuinfo)
        !if(shympi_is_master()) then
        !  write(iuinfo,*) 'dry_elems: ',iloop,idry
        !end if
	format = '(1x,a,3f8.0)'
	array = (/real(iloop),real(iedry),real(ikdry)/)
	!call info_output('wet_and_dry','none',3,array,.false.)
	call info_output('wet_and_dry','none',3,array,.true.,format)

	end

!*****************************************************************

