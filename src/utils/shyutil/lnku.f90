
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003-2004,2007,2009-2012,2014-2016  Georg Umgiesser
!    Copyright (C) 2019  Georg Umgiesser
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

! topological set up routines
!
! contents :
!
! function kthis(i,ie)	        	gets node at position i in ie
! function knext(k,ie) 	        	gets node after k in ie
! function kbhnd(k,ie)          	gets node before k in ie
! function ithis(k,ie) 	        	gets i of k in ie
! function inext(k,ie) 	        	gets i after k in ie
! function ibhnd(k,ie)	        	gets i before k in ie
!
! function ksthis(i,ie,nen3v)	        gets node at position i in ie
! function ksnext(k,ie,nen3v) 	        gets node after k in ie
! function ksbhnd(k,ie,nen3v)           gets node before k in ie
! function isthis(k,ie,nen3v) 	        gets position of k in ie
! function isnext(k,ie,nen3v) 	        gets position after k in ie
! function isbhnd(k,ie,nen3v)	        gets position before k in ie
!
! subroutine link_fill(n)       returns filling of linkv
!
! subroutine get_elem_linkp(k,ipf,ipl)	gets pointer to elements around k
! subroutine get_node_linkp(k,ipf,ipl)	gets pointer to nodes around k
! subroutine get_elem_links(k,n,ibase)	gets pointer to elements around k
! subroutine get_node_links(k,n,ibase)	gets pointer to nodes around k
!
! subroutine get_elems_around(k,ndim,n,elems) returns all elems around node k
! subroutine get_nodes_around(k,ndim,n,nodes) returns all nodes around node k
!
! subroutine find_elems_to_segment(k1,k2,ie1,ie2) finds elements to segment
!
! notes :
!
! the preferred way to get nodes/elems around node k is through
!	get_elems_around()
!	get_nodes_around()
! 
! revision log :
!
! 01.08.2003	ggu	created from sublnk.f
! 16.03.2004	ggu	new routine node_links() and link_fill()
! 10.11.2007	ggu	new routine line_elems
! 28.08.2009	ggu	routine line_elems renamed to find_elems_to_segment
! 28.08.2009	ggu	new routines get_elems_around, get_nodes_around
! 09.09.2009	ggu	bug fix in find_elems_to_segment (BUGip2)
! 23.03.2010	ggu	changed v6.1.1
! 20.10.2011	ggu	check dimension in set_elem_links(), set_node_links()
! 16.12.2011	ggu	in lnk_elems at boundary set last value to 0
! 24.01.2012	ggu	changed VERS_6_1_41
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_2
! 19.01.2015	ggu	changed VERS_7_1_3
! 05.05.2015	ggu	changed VERS_7_1_10
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 02.12.2015	ggu	lnk_elems and lnk_nodes eliminated
! 16.12.2015	ggu	changed VERS_7_3_16
! 28.04.2016	ggu	changed VERS_7_5_9
! 16.02.2019	ggu	changed VERS_7_5_60
! 29.09.2020	ggu	added dimension check in get_nodes/elems_around()
! 22.04.2022	ggu	new locator functions using passed nen3v
! 07.11.2024	ggu	new framework for connection
! 07.10.2025	ggu	bug fix in find_elems_to_segment()
!
!****************************************************************
!****************************************************************
!****************************************************************
! next are locator functions using global nen3v from basin
!****************************************************************
!****************************************************************
!****************************************************************

        function kthis(i,ie)

! gets node at position i in ie

	use basin

        implicit none

	integer kthis
        integer i			!position
        integer ie			!element

	kthis = nen3v(i,ie)

	end

!****************************************************************

        function knext(k,ie)

! gets node after k in ie

	use basin

        implicit none

	integer knext
        integer k			!actual node
        integer ie			!element

        integer i

        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            knext=nen3v(mod(i,3)+1,ie)
            return
          end if
        end do

        knext=0

        end

!****************************************************************

        function kbhnd(k,ie)

! gets node before k in ie

	use basin

        implicit none

	integer kbhnd
        integer k			!actual node
        integer ie			!element

        integer i

        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            kbhnd=nen3v(mod(i+1,3)+1,ie)
            return
          end if
        end do

        kbhnd=0

        end

!****************************************************************

        function ithis(k,ie)

! gets position of k in ie (might also use lenkiiv for this)

	use basin

        implicit none

	integer ithis
        integer k			!actual node
        integer ie			!element

        integer i

        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            ithis=i
            return
          end if
        end do

        ithis=0

        end

!****************************************************************

        function inext(k,ie)

! gets position after k in ie

	use basin

        implicit none

	integer inext
        integer k			!actual node
        integer ie			!element

        integer i

        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            inext=mod(i,3)+1
            return
          end if
        end do

        inext=0

        end

!****************************************************************

        function ibhnd(k,ie)

! gets position before k in ie

	use basin

        implicit none

	integer ibhnd
        integer k			!actual node
        integer ie			!element

        integer i

        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            ibhnd=mod(i+1,3)+1
            return
          end if
        end do

        ibhnd=0

        end

!****************************************************************
!****************************************************************
!****************************************************************
! next are locator functions using passed in nen3v
!****************************************************************
!****************************************************************
!****************************************************************

        function ksthis(i,ie,nen3v)

! gets node at position i in ie

        implicit none

	integer ksthis
        integer i			!position
        integer ie			!element
        integer nen3v(3,*)		!node index

	ksthis = nen3v(i,ie)

	end

!****************************************************************

        function ksnext(k,ie,nen3v)

! gets node after k in ie

        implicit none

	integer ksnext
        integer k			!actual node
        integer ie			!element
        integer nen3v(3,*)		!node index

        integer i

        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            ksnext=nen3v(mod(i,3)+1,ie)
            return
          end if
        end do

        ksnext=0

        end

!****************************************************************

        function ksbhnd(k,ie,nen3v)

! gets node before k in ie

        implicit none

	integer ksbhnd
        integer k			!actual node
        integer ie			!element
        integer nen3v(3,*)		!node index

        integer i

        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            ksbhnd=nen3v(mod(i+1,3)+1,ie)
            return
          end if
        end do

        ksbhnd=0

        end

!****************************************************************

        function isthis(k,ie,nen3v)

! gets position of k in ie (might also use lenkiiv for this)

        implicit none

	integer isthis
        integer k			!actual node
        integer ie			!element
        integer nen3v(3,*)		!node index

        integer i

        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            isthis=i
            return
          end if
        end do

        isthis=0

        end

!****************************************************************


        function isnext(k,ie,nen3v)

! gets position after k in ie

        implicit none

	integer isnext
        integer k			!actual node
        integer ie			!element
        integer nen3v(3,*)		!node index

        integer i

        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            isnext=mod(i,3)+1
            return
          end if
        end do

        isnext=0

        end

!****************************************************************

        function isbhnd(k,ie,nen3v)

! gets position before k in ie

        implicit none

	integer isbhnd
        integer k			!actual node
        integer ie			!element
        integer nen3v(3,*)		!node index

        integer i

        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            isbhnd=mod(i+1,3)+1
            return
          end if
        end do

        isbhnd=0

        end

!****************************************************************
!****************************************************************
!****************************************************************

        subroutine link_fill0(n)

! returns filling of linkv

	use mod_geom
	use basin, only : nkn,nel,ngr,mbw

        implicit none

! arguments
        integer n       !filling of linkv (return)

        n = ilinkv(nkn+1)

        end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine get_elem_linkp(k,ipf,ipl)

! gets pointer to first and last element around k
!
! to loop over the neibor elements, use similar:
!
!       call get_elem_linkp(k,ipf,ipl)
!       do ip=ipf,ipl
!         ien = lenkv(ip)          !ien is number of neibor element
!       end do

	use mod_geom

	implicit none

	integer k,ipf,ipl

        ipf = ilinkv(k)+1
        ipl = ilinkv(k+1)

	if( lenkv(ipl) .eq. 0 ) ipl = ipl - 1	!FIXME

	end

!****************************************************************

	subroutine get_node_linkp(k,ipf,ipl)

! gets pointer to first and last node around k
!
! to loop over the neibor nodes, use similar:
!
!       call get_node_linkp(k,ipf,ipl)
!       do ip=ipf,ipl
!         kn = linkv(ip)          !kn is number of neibor node
!       end do

	use mod_geom

	implicit none

	integer k,ipf,ipl

        ipf = ilinkv(k)+1
        ipl = ilinkv(k+1)

	end

!****************************************************************

	subroutine get_elem_links(k,n,ibase)

! gets pointer and total number of elements around k
!
! to loop over the neibor elements, use similar:
!
!       call get_elem_links(k,n,ibase)
!       do i=1,n
!         ien = lenkv(ibase+i)          !ien is number of neibor element
!       end do

	use mod_geom

	implicit none

	integer k,n,ibase

	n = ilinkv(k+1)-ilinkv(k)
	ibase = ilinkv(k)

	if( lenkv(ibase+n) .eq. 0 ) n = n - 1

        end

!****************************************************************

	subroutine get_node_links(k,n,ibase)

! gets pointer and total number of nodes around k
!
! to loop over the neibor nodes, use similar:
!
!       call get_node_links(k,n,ibase)
!       do i=1,n
!         kn = linkv(ibase+i)          !kn is number of neibor node
!       end do

	use mod_geom

	implicit none

	integer k,n,ibase

	n = ilinkv(k+1)-ilinkv(k)
	ibase = ilinkv(k)

        end

!****************************************************************

        subroutine get_elems_around(k,ndim,n,elems)

! returns all elems around node k

	use mod_geom

        implicit none

        integer, intent(in) :: k               !central node
        integer, intent(in) :: ndim            !dimension of elems()
        integer, intent(out) :: n              !total number of elems around k
        integer, intent(out) :: elems(ndim)    !elems around k

	integer i,ibase
	integer aux(size(elems))

	n = elist(0,k)
	if( n > ndim ) stop 'error stop get_elems_around: n>ndim'
	elems(1:n) = elist(1:n,k)
	if( .not. buselink ) return
	
	n = ilinkv(k+1)-ilinkv(k)
	ibase = ilinkv(k)

	if( lenkv(ibase+n) .eq. 0 ) n = n - 1

	if( n > ndim ) stop 'error stop get_elems_around: n>ndim'

	do i=1,n
	  aux(i) = lenkv(ibase+i)
	end do

	if( any( aux(1:n) /= elems(1:n) ) ) then
	  write(6,*) 'k = ',k
	  write(6,*) elist(0,k),n
	  write(6,*) elist(1:n,k)
	  write(6,*) aux(1:n)
	  stop 'error stop get_elems_around: icompatibility'
	end if

	end

!****************************************************************

        subroutine get_nodes_around(k,ndim,n,nodes)

! returns all nodes around node k

	use mod_geom

        implicit none

        integer, intent(in) :: k               !central node
        integer, intent(in) :: ndim            !dimension of elems()
        integer, intent(out) :: n              !total number of elems around k
        integer, intent(out) :: nodes(ndim)    !nodes around k

	integer i,ibase
	integer aux(size(nodes))

	n = nlist(0,k)
	if( n > ndim ) stop 'error stop get_nodes_around: n>ndim'
	nodes(1:n) = nlist(1:n,k)
	if( .not. buselink ) return
	
	n = ilinkv(k+1)-ilinkv(k)
	ibase = ilinkv(k)

	if( n > ndim ) stop 'error stop get_nodes_around: n>ndim'

	do i=1,n
	  aux(i) = linkv(ibase+i)
	end do

	if( any( aux(1:n) /= nodes(1:n) ) ) then
	  write(6,*) 'k = ',k
	  write(6,*) nlist(0,k),n
	  write(6,*) nlist(1:n,k)
	  write(6,*) aux(1:n)
	  stop 'error stop get_nodes_around: icompatibility'
	end if

	end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine find_elems_to_segment(k1,k2,ie1,ie2)

! finds elements to segment between nodes k1 and k2
!
! returns elements in ie1 and ie2
!
! ie1 is to the left of segment k1-k2, ie2 to the right
! if boundary segment only one ie is set, the other is zero
! if no such segment, both ie are zero

	use basin
	use mod_geom

	implicit none

        integer, intent(in)  :: k1,k2
        integer, intent(out) :: ie1,ie2

        integer k,ipf,ipl,ip,ip2
	integer i,iee1,iee2,n
	integer n1,n2

	ie1 = 0
	ie2 = 0
	
	n1 = nlist(0,k1)
	n2 = nlist(0,k2)
	do i=1,n1
	  k = nlist(i,k1)
	  if( k .eq. k2 ) then
	    ie1 = elist(i,k2)
	    if( i == 1 ) then
	      ie2 = elist(n2,k2)	!last element - if border is 0 (bug fix)
	    else
	      ie2 = elist(i-1,k2)	!previous element
	    end if
	    exit
	  end if
	end do

	if( ie1 == 0 .and. ie2 /= 0 ) then	!node with more boundaries
	  ie1 = ie2
	  ie2 = 0
	end if

	if( .not. buselink ) return

	k = k1
        ipf=ilinkv(k)+1
        ipl=ilinkv(k+1)

	iee1 = 0
	iee2 = 0

	do ip=ipf,ipl
	  k = linkv(ip)
	  if( k .eq. k2 ) then
	    iee1 = lenkv(ip)
	    if( ip .eq. ipf ) then
		ip2 = ipl		!this sets it to 0
	    else
		ip2 = ip - 1		!previous element	!BUGip2
	    end if
	    iee2 = lenkv(ip2)
	    exit
	  end if
	end do

	if( ie1 /= iee1 .or. ie2 /= iee2 ) then
	  write(6,*) ie1,iee1,ie2,iee2
	  stop 'error stop find_elems_to_segment: icompatibility'
	end if

	if( .not. any(k1==nen3v(:,ie1)) ) goto 99
	if( .not. any(k2==nen3v(:,ie1)) ) goto 99
	if( ie2 == 0 ) return
	if( .not. any(k1==nen3v(:,ie2)) ) goto 99
	if( .not. any(k2==nen3v(:,ie2)) ) goto 99

	return
   99	continue
	write(6,*) k1,k2
	write(6,*) ie1,ie2
	write(6,*) nen3v(:,ie1)
	if( ie2 /= 0 ) write(6,*) nen3v(:,ie2)
	stop 'error stop find_elems_to_segment: nodes not found'
	end

!****************************************************************

