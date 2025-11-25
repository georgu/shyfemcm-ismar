
!--------------------------------------------------------------------------
!
!    Copyright (C) 2024  Georg Umgiesser
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

! pre-processing routine (intermediate version) - can be deleted
!
! revision log :
!
! 08.11.2024    ggu     last intermediate version

!*******************************************************************

	subroutine handle_shypre_optimization(bwrite,bww,bdebug,bauto,bopti &
     &				,nkn,nel,nen3v,ng,ipv &
     &				,kphv,iarnv &
     &				,xgv,ygv,hkv,mbw,ngr)

	implicit none

	logical bwrite,bww,bdebug,bauto,bopti
	integer nkn,nel
	integer nen3v(3,nel)
	integer ng(nkn)
	integer ipv(nkn)
	integer kphv(nkn),iarnv(nkn)
	real xgv(nkn),ygv(nkn),hkv(nkn)
	integer mbw,ngr

	logical bstop
	integer ngr1
	integer nelddi
	integer, allocatable :: iphv(:)
	integer, allocatable :: iknot(:,:)
	integer, allocatable :: kvert(:,:)

	bstop = .false.
	nelddi = nel

	allocate(iphv(nkn))
	allocate(kvert(2,nkn))

	call ininum(nkn,iphv,kphv)

!--------------------------------------------------------
! setup side index and compute ngr and mbw
!--------------------------------------------------------

        if( bwrite ) write(6,*) ' ...setting up side index'

        call estimate_grade(nkn,nel,nen3v,ng,ngr1)
        allocate(iknot(ngr1,nkn))

        call sidei(nkn,nel,nen3v,ng,iknot,ngr1,ngr,bstop)
        if(bstop) goto 99934
        call check_sidei(nkn,nel,nen3v,ipv,ng,ngr1,bstop)
        if(bstop) goto 99935
        call knscr(nkn,ngr,ngr1,iknot)

        call gtest(bdebug,'bandwidth',nelddi,nkn,nel,nen3v)
        call bandw(nel,nen3v,mbw)
        if( bww ) write(6,*) 'Maximum grade of nodes is ',ngr
        if( bww ) write(6,*) 'Bandwidth is ',mbw

!--------------------------------------------------------
! bandwidth optimization
!--------------------------------------------------------

        if( .not. bopti ) then
	  if( bww ) write(6,*) 'Bandwidth is not optimized'
	  return
	end if

        call gtest(bdebug,'bandwidth 1',nelddi,nkn,nel,nen3v)
        if( bwrite ) write(6,*) ' ...optimizing of band width'
        call bandop(nkn,ngr,ipv,iphv,kphv,ng,iknot,kvert &
     &                  ,bopti,bauto,bww)
        call gtest(bdebug,'bandwidth 2',nelddi,nkn,nel,nen3v)

        if(bopti) then
          call bandex(nkn,nel,nen3v,kphv,ipv,iarnv,xgv,ygv,hkv) !kphv is n-rank
          call gtest(bdebug,'bandwidth 3',nelddi,nkn,nel,nen3v)
          call bandw(nel,nen3v,mbw)
          call gtest(bdebug,'bandwidth 4',nelddi,nkn,nel,nen3v)
          if( bww ) write(6,*) 'Optimized bandwidth is ',mbw
        end if

	return
99934   write(6,*)' (34) error: ngr1 is wrong (internal error'
        stop 'error stop handle_optimization: (34)'
99935   write(6,*)' (35) error: irregular connections'
        stop 'error stop handle_optimization: (35)'
	end

!**********************************************************

        subroutine bandop(nkn,ngr1,ipv,iphv,kphv,ng,iknot,kvert &
     &                          ,bopti,bauto,bwrite)

! optimize band width

        implicit none

        integer nkn,ngr1
        logical bopti,bauto
        integer ipv(nkn)
        integer iphv(nkn),kphv(nkn)
        integer ng(nkn)
        integer iknot(ngr1,nkn)
        integer kvert(2,nkn)
        logical bwrite,bloop

        integer iantw

        if( bauto ) then
          call ininum(nkn,iphv,kphv)
          call optest('before optimization: ',nkn,ipv,iphv,kphv)
          call cmgrade(nkn,ngr1,ipv,iphv,kphv,ng,iknot,1,4,bwrite)
          call optest('after Cuthill McKee: ',nkn,ipv,iphv,kphv)
          call revnum(nkn,iphv,kphv,bwrite)
          call optest('after reversing nodes: ',nkn,ipv,iphv,kphv)
          call rosen(nkn,ngr1,iphv,kphv,ng,iknot,kvert,bwrite)
          call optest('after Rosen: ',nkn,ipv,iphv,kphv)
          return
        end if

        bloop = .true.

        do while( bloop )
          call ininum(nkn,iphv,kphv)

!         call anneal(nkn,ngr1,kphv,ng,iknot,iphv,kvert)

          if( iantw(' Cuthill McKee algorithm ?') .gt. 0 ) then
            call cmv(nkn,ngr1,ipv,iphv,kphv,ng,iknot,bwrite)
          end if

          if( iantw(' Reverse numbering of nodes ?') .gt. 0 ) then
            call revnum(nkn,iphv,kphv,bwrite)
          end if

          if( iantw(' Rosen algorithm ?') .gt. 0 ) then
            call rosen(nkn,ngr1,iphv,kphv,ng,iknot,kvert,bwrite)
          end if

          bloop = iantw(' Repeat optimization of bandwidth ?') .gt. 0
        end do

        bopti = .true.

        end

!**********************************************************

        subroutine revnum(nkn,iphv,kphv,bwrite)

! reverses numbering of nodes

        implicit none

        integer nkn
        integer iphv(nkn), kphv(nkn)
        logical bwrite

        integer i

        if( bwrite ) write(6,*) 'Applying reverse algorithm...'

        do i=1,nkn
          kphv(i) = nkn+1-kphv(i)
          iphv(kphv(i)) = i
        end do

        end

!**********************************************************


        subroutine zernum(nkn,iphv,kphv)

! zeros numbering of nodes

        implicit none

        integer nkn
        integer iphv(nkn), kphv(nkn)

        integer i

        do i=1,nkn
          iphv(i) = 0
          kphv(i) = 0
        end do

        end

!**********************************************************

        subroutine bandex(nkn,nel,nen3v,kphv &
     &                      ,ipv,iarnv,xgv,ygv,hkv)

! exchange nodes after optimization

        implicit none

        integer nkn,nel
        integer nen3v(3,nel)
        integer kphv(nkn)
        integer ipv(nkn)
        integer iarnv(nkn)
        real xgv(nkn),ygv(nkn)
        real hkv(nkn)

        integer ie,ii
        integer neaux(3,nel)

        do ie=1,nel
          do ii=1,3
            neaux(ii,ie)=nen3v(ii,ie)
          end do
        end do

        do ie=1,nel
          do ii=1,3
            nen3v(ii,ie)=kphv(neaux(ii,ie))
          end do
        end do

!       copy arrays with kphv as rank table

        call icopy(nkn,ipv,kphv)
        call icopy(nkn,iarnv,kphv)
        call rcopy(nkn,xgv,kphv)
        call rcopy(nkn,ygv,kphv)
        call rcopy(nkn,hkv,kphv)

        return
        end

!**********************************************************

        subroutine optest(text,nkn,ipv,iphv,kphv)

        implicit none

        character*(*) text
        integer nkn
        integer ipv(nkn)
        integer iphv(nkn)
        integer kphv(nkn)

        integer i

        do i=1,nkn
          if( iphv(i) .gt. nkn ) goto 99
          if( kphv(i) .gt. nkn ) goto 99
          if( iphv(kphv(i)) .ne. i ) goto 99
          if( kphv(iphv(i)) .ne. i ) goto 99
        end do

        return
   99   continue
        write(6,*) '*** Error in optest: '
        write(6,*) text
        write(6,*) 'error in pointers...'
        write(6,*) 'problem is close to following node...'
        write(6,*) 'node (intern/extern): ',i,ipv(i)
        write(6,*) i,nkn
        write(6,*) iphv(i),kphv(i)
        write(6,*) iphv(kphv(i)),kphv(iphv(i))
        write(6,*) 'maybe the domain is not connected...'
        stop 'error stop optest: pointers'
        end

!**********************************************************

        subroutine knscr(nkn,ngr,ngr1,iknot)

! de-scrambles side index for cmv and rosen

        implicit none

        integer nkn,ngr,ngr1
        integer iknot(ngr1*nkn)

        integer k,j

        do k=1,nkn
          do j=1,ngr
            iknot((k-1)*ngr+j)=iknot((k-1)*ngr1+j)
          end do
        end do

        return
        end

!**********************************************************
!**********************************************************
!**********************************************************

        subroutine handle_shypre_partition(bwrite,bnepart,grdpart,nn,ne &
     &				,knrank,ierank)

        use grd
        use basin

        implicit none

	logical bwrite
	logical bnepart
	character*80 grdpart
        integer nn,ne
        integer knrank(nn)
        integer ierank(ne)

        logical bgrd
        integer nnpart,nepart
        integer, allocatable :: area_node(:)
        integer, allocatable :: area_elem(:)

        integer i

        if( bwrite ) write(6,*) 'handle_partition: ',nn,ne

        bgrd = ( grdpart /= ' ' )
        if( .not. bnepart .and. .not. bgrd ) return

        if( bnepart .and. bgrd ) then
          write(6,*) 'only one of -partition and -nepart can be given'
          stop 'error stop handle_partition: options'
        end if

        allocate(area_node(nn))
        allocate(area_elem(ne))

        if( bgrd ) then
          write(6,*) 'reading partitioning file ',trim(grdpart)
          call grd_read(grdpart)
          if( nk_grd /= nn ) goto 99
          if( ne_grd /= ne ) goto 99
          area_node = ianv
          area_elem = iaev
        else
          write(6,*) 'using area code as partitioning information'
          area_node = iarnv
          area_elem = iarv
        end if

        call renumber_partition(nn,area_node,nnpart)
        call renumber_partition(ne,area_elem,nepart)

        if( bgrd ) then
          call icopy(nn,area_node,knrank)
          call icopy(ne,area_elem,ierank)
        end if

        call basin_set_partition(nn,ne,nnpart,nepart,area_node,area_elem)

        return
   99   continue
        write(6,*) nk_grd,nn
        write(6,*) ne_grd,ne
        stop 'error stop handle_partition: incompatibility'
        end

!**********************************************************

        subroutine renumber_partition(n,area,npart)

        implicit none

        integer n
        integer area(n)
        integer npart

        integer i,ia,imax
        integer nmin,nmax
        integer, allocatable :: table_in(:),table_out(:)

        nmin = minval(area)
        nmax = maxval(area)

        if( nmin == nmax ) then         !no partition
          npart = 0
          area = 0
          return
        end if

        if( nmin < 0 .or. nmax < 0 ) then
          write(6,*) nmin,nmax
          stop 'error stop renumber_partition: nmin,nmax'
        end if

        allocate(table_in(0:nmax),table_out(0:nmax))

        table_in = 0
        table_out = 0

        do i=1,n
          ia = area(i)
          table_in(ia) = table_in(ia) + 1
        end do

        imax = -1
        do ia=0,nmax
          if( table_in(ia) > 0 ) then
            imax = imax + 1
            table_out(ia) = imax
          end if
        end do
        npart = imax

        do i=1,n
          ia = area(i)
          area(i) = table_out(ia)
        end do

        end

!**********************************************************


!--------------------------------------------------------------------------
!
!    Copyright (C) 1998,2010,2012,2015,2017,2019  Georg Umgiesser
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

! Rosen bandwidth optimization algorithm
!
! revision log :
!
! 05.06.1998	ggu	avoid write to terminal
! 23.03.2010	ggu	changed v6.1.1
! 30.03.2012	ggu	changed VERS_6_1_51
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_53
! 18.12.2015	ggu	changed VERS_7_3_17
! 14.11.2017	ggu	changed VERS_7_5_36
! 16.02.2019	ggu	changed VERS_7_5_60
!
!*********************************************************************

	subroutine rosen(nkn,ngrddi,iphv,kphv,ng,iknot,kvert,bwrite)
!
! rosen algorithmus
!
! k..		neue knotennummern
! i..		alte knotennummern
! ng		!grad des knotens
! iknot		!kantenverzeichnis des knotens
! 		...iknot(n,k) --> iknot(ngr*(k-1)+n)
!		...( n <= ng(k) )
! iphv		pointer fuer knotennummern : neu --> alt
! kphv		pointer fuer knotennummern : alt --> neu
! ipaar		indexpaare die die maximale bandbreite bestimmen
! nppar   anzahl der indexpaare auf ipaar (max. ndim)
! npaara	gesammtanzahl der indexpaare die die
!		...maximale bandbreite bestimmen
! khil		indexpaare die bei vertauschung die bandbreite
!		...unveraendert lassen
! nhil    anzahl der indexpaare auf khil (max. ndim)
! kvert		indexpaare die bereits miteinander
!		...vertauscht worden sind
! nvert		anzahl der indexpaare auf kvert (max. nkn)
!		...wird auf 0 gesetzt wenn es gelingt zwei
!		...indizes zu vertauschen bei denen die
!		...bandbreite erniedrigt wird
! lvert		=1 ==> es ist gelungen, ein indexpaar zu
!		...vertauschen. dabei ist nicht notwendigerweise
!		...die bandweite erniedrigt worden. das program
!		...wird beendet wenn in do-90-loop kein index-
!		...paar mehr vertauscht werden kann.
!
        integer ndim
        parameter (ndim=10)

	integer nkn,ngrddi
	integer iphv(nkn),kphv(nkn)
	integer ng(nkn),iknot(ngrddi*nkn)
	integer kvert(2,nkn)
	logical bwrite

        integer ipaar(2,ndim),khil(2,ndim)

        integer nhil,nvert,ngr
        integer m,npaar,npaara,lvert
        integer i,ii,kn1,kn2
        integer mh,ms,msloop
        integer ipp,ipa1,ipa2,kk,kmax
        integer kanf,k,mhhmax,mhh,mmin
        integer kmerk,k1,k2,kmin

        logical true
        integer msmax
        save true,msmax
        data true  / .true. /
        data msmax /  100   /
!
	if( bwrite ) write(6,*) 'Applying Rosen algorithm...'

	msloop=0
	nhil=0
	nvert=0
        ms=nkn+1
	ngr=ngrddi
	kmerk = 0
!
        do while (true)
!
! bestimmung der indexpaare die maximale
! ...bandbreite m haben und ablegen auf %%%%%%%%%%%%%%%%%%
! ...ipaar (max. ndim)
!
	m=0
	npaar=0
	npaara=0
	lvert=0
	do i=1,nkn
	kn1=kphv(i)
	do ii=1,ng(i)
	kn2=kphv(iknot(ngr*(i-1)+ii))
        mh=abs(kn2-kn1)
	if(mh.ge.m) then
		if(mh.eq.m) then
            if(npaar.lt.ndim) then
				npaar=npaar+1
				ipaar(1,npaar)=i
				ipaar(2,npaar)=iknot(ngr*(i-1)+ii)
			end if
			npaara=npaara+1
		else
			m=mh
			npaar=1
			npaara=1
			ipaar(1,1)=i
			ipaar(2,1)=iknot(ngr*(i-1)+ii)
		end if
	end if
	end do
	end do

!        write(6,'(a,i4,a,i4,a,i4,a,i4)')
!     +       ' m = '            ,m
!     +      ,' pairs = '        ,npaara
!     +      ,' yet exchanged = ',nvert
!     +      ,' tries = '        ,msloop

! new on may 93 t stop loop

        if(m.lt.ms) then
          ms=m
          msloop=msmax
        else if(msloop.eq.0) then
          npaar=0   !do not enter loop -> stop
        else
          msloop=msloop-1
        end if
!
! schleife ueber paare die erniedrigt werden muessen %%%%%%%%%%%%
!
	do ipp=1,npaar
!
	ipa1=ipaar(1,ipp)
	ipa2=ipaar(2,ipp)
	kn1=kphv(ipa1)			!kleinerer index
	kn2=kphv(ipa2)			!groesserer index
        if(kn1.gt.kn2) then
          call iiswap(kn1,kn2)
          call iiswap(ipa1,ipa2)
        end if
	mh=kn2-kn1			!bandbreite
	if(mh.lt.m) goto 90		!naechstes paar
	if(mh.gt.m) stop 'error stop : m to big'
!
! groesseren index kn2 mit einem kleineren vertauschen
!
	kmax=kn1			!maximale knoten-
	do ii=1,ng(ipa2)		!...nummer bestimmen
          kk=kphv(iknot(ngr*(ipa2-1)+ii)) !...die mit kn2
          if(kk.gt.kmax) kmax=kk    !...verbunden ist
	end do
!
	kanf=kmax-m
        if(kanf.lt.1) kanf=1
	mmin=m				!minimale bandbreite
	do k=kn2-1,kanf,-1		!k ist index der mit
	ia=iphv(k)			!...kn2 vertauscht wird
	mhhmax=max(kmax-k,k-kn1)	!...und bandbreite
	if(mhhmax.gt.mmin) goto 10	!...mhhmax ergibt
	do ii=1,ng(ia)		!bestimmung der bandbreite fuer
	mhh=iabs(kn2-kphv(iknot(ngr*(ia-1)+ii))) !...restl. knoten
	if(mhh.eq.0) mhh=kn2-k	!kn2 war im knotenverzeichnis
	if(mhh.gt.mmin) goto 10
	if(mhh.gt.mhhmax) mhhmax=mhh
	end do
	if(mhhmax.lt.mmin) then		!kleinere bandbreite
		mmin=mhhmax		!...gefunden, merke
		kmerk=k			!...k und m
	else if(mhhmax.eq.m) then	!vertauschung laesst
          if(nhil.lt.ndim) then !...bandbreite
			nhil=nhil+1	!...gleich,
			khil(1,nhil)=k	!...trage in
			khil(2,nhil)=kn2!...nhil ein
		end if
	end if
   10	continue
	end do
!
	if(mmin.lt.m) then		!kleinere bandbreite gefunden
		nvert=0			!...--> vertausche
		call vertau(nkn,kmerk,kn2,iphv,kphv,lvert,nhil)
		goto 90
	end if
!
! kleineren index kn1 mit einem groesseren vertauschen
!
	kmin=kn2			!minimale knoten-
	do ii=1,ng(ipa1)		!...nummer bestimmen
	kk=kphv(iknot(ngr*(ipa1-1)+ii))	!...die mit kn1
	if(kk.lt.kmin) kmin=kk		!...verbunden ist
	end do
!
	kend=kmin+m
        if(kend.gt.nkn) kend=nkn
	mmin=m
	do k=kn1+1,kend		!k ist index der mit
	ia=iphv(k)			!...kn1 vertauscht wird
	mhhmax=max(k-kmin,kn2-k)
	if(mhhmax.gt.mmin) goto 11
	do ii=1,ng(ia)
	mhh=iabs(kn1-kphv(iknot(ngr*(ia-1)+ii)))
	if(mhh.eq.0) mhh=k-kn1	!kn1 war im knotenverzeichnis
	if(mhh.gt.mmin) goto 11
	if(mhh.gt.mhhmax) mhhmax=mhh
	end do
	if(mhhmax.lt.mmin) then
		mmin=mhhmax
		kmerk=k
	else if(mhhmax.eq.m) then
          if(nhil.lt.ndim) then
			nhil=nhil+1
			khil(1,nhil)=kn1
			khil(2,nhil)=k
		end if
	end if
   11	continue
	end do
!
	if(mmin.lt.m) then
		nvert=0
		call vertau(nkn,kn1,kmerk,iphv,kphv,lvert,nhil)
		goto 90
	end if
!
! bandbreite kann nicht erniedrigt werden
!
	if(nvert.ge.nkn) stop 'error stop :dimension kvert'
!
	do i=1,nhil		!ueberpruefen ob es noch
	do ii=1,nvert		!...indizes gibt auf nhil
	if(khil(1,i).eq.kvert(1,ii)) then!...die noch nicht
		if(khil(2,i).eq.kvert(2,ii)) goto 21	!...
	end if			!...vertauscht worden sind
	end do
	k1=khil(1,i)
	k2=khil(2,i)
	nvert=nvert+1
	kvert(1,nvert)=k1
	kvert(2,nvert)=k2
	call vertau(nkn,k1,k2,iphv,kphv,lvert,nhil)
	goto 90
   21	continue
	end do
!
! keine vertauschung moeglich
!
   90	continue
	end do
!
	if(lvert.eq.0) exit	!keine vertauschung mehr moegl.
!
	end do	!do while
!
	if( bwrite ) then
		write(6,*) 'minimal bandwidth found =',m
		write(6,*) 'total number of pairs =',npaara
	end if

	return
	end

!*****************************************************************

	subroutine vertau(nkn,k1,k2,iphv,kphv,lvert,nhil)

! exchanges indices

	implicit none

	integer nkn
	integer k1,k2
	integer iphv(nkn),kphv(nkn)
	integer lvert,nhil

	integer i1,i2

	lvert=1
	nhil=0
	i1=iphv(k1)
	i2=iphv(k2)
	iphv(k1)=i2
	iphv(k2)=i1
	kphv(i1)=k2
	kphv(i2)=k1

	return
	end

!******************************************************************

        subroutine iiswap(a,b)

! swaps two integers

        integer a,b,c

        c=a
        a=b
        b=c

        return
	end

!--------------------------------------------------------------------------
!
!    Copyright (C) 1998,2010-2011,2013,2015,2017-2019  Georg Umgiesser
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

! cuthill-mckee algorithm for optimization of bandwidth
!
! revision log :
!
! 20.03.1998	ggu	reorganized, automatic procedure introduced
! 05.06.1998	ggu	avoid write to terminal
! 23.03.2010	ggu	changed v6.1.1
! 14.04.2011	ggu	changed VERS_6_1_22
! 25.10.2013	ggu	changed VERS_6_1_68
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_53
! 18.12.2015	ggu	changed VERS_7_3_17
! 14.11.2017	ggu	changed VERS_7_5_36
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.02.2019	ggu	changed VERS_7_5_60
!
!*******************************************************************

	subroutine cmv(nkn,ngrddi,ipv,iphv,kphv,ng,iknot,bwrite)

! cuthill-mckee algorithmus
!
! k..		neue knotennummern
! i..		alte knotennummern
! ng		grad des knoten
! iknot		kantenverzeichnis des knotens
!		...iknot(n,k) --> iknot(ngr*(k-1)+n)
!		...( n <= ng(k) )
! iphv		pointer fuer knotennummern : neu --> alt
! kphv		pointer fuer knotennummern : alt --> neu
! ikanf		anfangsknoten (erste stufe)
! ikmer		anfangsknoten fuer minimale bandbreite
! mmin		bandbreite fuer anfangsknoten ikmer
! knum		anzahl der bereits nummerierten knoten
! kanf,kend	neue knotennummern der alten stufe

        implicit none

	integer nkn,ngrddi
	integer ipv(nkn)
	integer iphv(nkn),kphv(nkn)
	integer ng(nkn),iknot(ngrddi*nkn)
	logical bwrite

        integer iwei
        integer jgrmin,jgrmax
        integer knum,m

! get switch iwei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	iwei = 1

	do while( iwei .eq. 1 .or. iwei .eq. 2 )

	  call cmvwhat(iwei)

          if(iwei.eq.1) then				!grades
	    call getgrds(ngrddi,nkn,ng,jgrmin,jgrmax)
	    call cmgrade(nkn,ngrddi,ipv,iphv,kphv,ng,iknot &
     &				,jgrmin,jgrmax,bwrite)
	  else if(iwei.eq.2) then			!give first level
	    call getfstl(nkn,iphv,kphv,knum)
	    if( knum .gt. 0 ) then
              call cmalg(nkn,ngrddi,knum,m,iphv,kphv,ng,iknot)
	      write(6,*) 'node =',ipv(iphv(1)),'   mbw =',m
	    end if
	  else if(iwei.eq.3) then			!return old numbering
	    call ininum(nkn,iphv,kphv)
	  end if

	end do	!do while

	end

!********************************************************************

        subroutine cmalg(nkn,ngrddi,knum,m,iphv,kphv,ng,iknot)

! cm-algorithmus
!
! knum		anzahl der vorgegebenen knoten der ersten stufe
! kanf,kend	neue knotennummern der alten stufe
!		(...knum=kend-kanf+1)
! 		iphv,kphv muss fuer erste stufe schon
!		...besetzt sein (der rest null)
! m		enthaelt am ende gefundene bandbreite

        implicit none

        integer nkn,ngrddi,knum,m
	integer iphv(nkn),kphv(nkn)
	integer ng(nkn),iknot(ngrddi*nkn)

        integer kanf,kend,knold
        integer ngr,ka,inh,in
        integer nggmin,ngg,n,ia,mh

        ngr=ngrddi

	m=0
	kanf=1			!grenzen fuer
	kend=knum		!...erste stufe
        knold=0

        do while(knum.gt.knold)
          knold=knum
          do ka=kanf,kend
            ia=iphv(ka)   !knoten der alten stufe
            inh=1
            do while(inh.ne.0)  !inh=0 ==> zu ia gibt es
              inh=0     !...keinen unnummerierten
              nggmin=ngr+1    !...knoten mehr --> neues ia
              do n=1,ng(ia)   !suche knoten in der neuer
                in=iknot((ia-1)*ngr+n)  !...stufe der mit ia
                if(kphv(in).eq.0) then  !...verbunden und noch
                  ngg=ng(in)  !...nicht nummeriert
                  if(ngg.lt.nggmin) then  !...ist und
                    nggmin=ngg  !...kleinsten
                    inh=in  !...grad besitzt und
                  end if    !...und schreibe
                end if      !...ihn auf inh
              end do
              if(inh.ne.0) then !naechster zu
                knum=knum+1 !...nummerierender
                iphv(knum)=inh  !...knoten ist
                kphv(inh)=knum  !...inh und
                mh=iabs(ka-knum)!...die bandweite
                if(mh.gt.m) m=mh!...ist mh
              end if
            end do  !do while(inh.ne.0)
          end do  !ka
          kanf=kend+1   !neue stufe wird
          kend=knum   !...alte stufe
	end do	!do while

        return
	end

!****************************************************************

	subroutine cmvwhat(iwei)

! asks about what to do (CMV algorithm)

	implicit none

	integer iwei

	integer ianz
	real f(10)

	integer inquire_numbers,iround

	write(6,*) '0   save numbering and exit'
	write(6,*) '1   numbering for grades'
	write(6,*) '2   numbering starting from first level'
	write(6,*) '3   quit with numbering prior to routine call'

	iwei = -1
	do while( iwei .lt. 0 .or. iwei .gt. 3 )
	  iwei = 0
	  ianz = inquire_numbers('give number :',f,1)
	  if( ianz .gt. 0 ) iwei = iround(f(1))
	end do

	end

!****************************************************************

	subroutine getgrds(ngr,nkn,ng,jgrmin,jgrmax)

! gets min/max grades to use for a starting point minimization

	implicit none

	integer ngr,nkn
	integer ng(nkn)
	integer jgrmin,jgrmax

	real f(2)
	integer ianz,jz,n,j
	integer inquire_numbers

	do n=1,ngr
            jz=0
            do j=1,nkn
              if(ng(j).eq.n) jz=jz+1
            end do

            if(jz.ne.0) then
              write(6,*) 'grade =',n,'   number of nodes =',jz
            end if
	end do

	ianz = -1
	do while( ianz .lt. 0 .or. ianz .gt. 2 )
	  ianz=inquire_numbers('give min/max grade (default 1/4) :',f,2)
	  if(ianz.eq.2) then
		jgrmin=f(1)
		jgrmax=f(2)
	  else if(ianz.eq.1) then
		jgrmin=f(1)
		jgrmax=f(1)
	  else if(ianz.eq.0) then
		jgrmin=1
		jgrmax=4
	  else
		write(6,*) 'erroneous input'
	  end if
	end do

	end

!****************************************************************

	subroutine getfstl(nkn,iphv,kphv,knum)

! gets first level interactively

	implicit none

	integer nkn
	integer iphv(nkn), kphv(nkn)
	integer knum

	integer ianz
	integer node,noden
	real f(2)

	integer ipint,inquire_numbers,iround

	knum=0
	call zernum(nkn,iphv,kphv)

	write(6,*)'give nodes of first level (<cr> to end)'

	ianz=inquire_numbers(' : ',f,2)
	do while(ianz.gt.0)
		if(ianz.gt.1) then
		  write(6,*)'only one node a time -> repeat'
		else
		  node = iround(f(1))
		  noden = ipint(node)
		  if(noden.ne.0) then
			knum=knum+1
			iphv(knum)=noden
			kphv(noden)=knum
		  else
			write(6,*)'invalid node number'
		  end if
		end if
		ianz=inquire_numbers(' : ',f,2)
	end do

	if(knum.eq.0) then	!no node
		write(6,*)'no valid node - old numbering reinstalled'
		call ininum(nkn,iphv,kphv)
	end if

	end

!****************************************************************

        subroutine nodnum(nkn,iphv,kphv,node,knum)

! initializes numbering of nodes for first node

        implicit none

        integer nkn
        integer iphv(nkn), kphv(nkn)
	integer node,knum

        integer i

        do i=1,nkn
          iphv(i) = 0
          kphv(i) = 0
        end do

        knum=1
        iphv(knum) = node
        kphv(node) = knum

	end

!****************************************************************

	subroutine cmgrade(nkn,ngrddi,ipv,iphv,kphv,ng,iknot, &
     &				jgrmin,jgrmax,bwrite)

! cuthill-mckee algorithm for grades

        implicit none

	integer nkn,ngrddi
	integer ipv(nkn)
	integer iphv(nkn),kphv(nkn)
	integer ng(nkn),iknot(ngrddi*nkn)
        integer jgrmin,jgrmax
	logical bwrite

        integer i
        integer mmin,ikmer,knum,m

	ikmer=0
	mmin=nkn
	m = 0

	if( bwrite ) write(6,*) 'Applying Cuthill-McKee algorithm...'
	!write(6,*) jgrmin,jgrmax
        !write(6,'(16x,a,7x,a,3x,a)') 'node','grade','bandwidth'

        do i=1,nkn
            if(ng(i).ge.jgrmin.and.ng(i).le.jgrmax) then
              call nodnum(nkn,iphv,kphv,i,knum)
              call cmalg(nkn,ngrddi,knum,m,iphv,kphv,ng,iknot)
	      !write(6,*) i,ipv(i),ng(i),m,mmin
	      call optest('inside cmgrade',nkn,ipv,iphv,kphv)
!              write(6,'(8x,3i12)') ipv(i),ng(i),m
              if(m.lt.mmin) then
                mmin=m
                ikmer=i
              end if
            end if
	end do

	if(ikmer.eq.0) then	!no valid node found
	    write(6,*)'no node for this grade -> '
	    write(6,*)'old numbering reinstalled'
	    call ininum(nkn,iphv,kphv)
        else			!repeat mimimal node
            call nodnum(nkn,iphv,kphv,ikmer,knum)
            call cmalg(nkn,ngrddi,knum,m,iphv,kphv,ng,iknot)
            if( bwrite ) write(6,*) 'node =',ipv(iphv(1)),'   mbw =',m
	end if

	end

!****************************************************************

