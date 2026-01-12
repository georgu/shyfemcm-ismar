
!--------------------------------------------------------------------------
!
!    Copyright (C) 2000,2011  Georg Umgiesser
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

! hash routines
!
! to change for new type of info and key :
!
! probably only changes to make in hashin internal routine
!
! flagf, flagd	must be different from any possible value in key
! ndim          must be greater than about twice the elements
!			to be inserted and must be a prime
! keyt		is key table that should be ok for any use (integer)
! infot		is the information to be inserted  -> change the
!			structure of this array if needed.
!			you must change also routine copyi
!			that actually performs the copy of
!			info and infot. the other functions
!			probably may be left unchanged since
!			they only pass info to hashin
! hash		the hash function -> first statement in hashin
!
! revision log :
!
! 30.06.2000	ggu	error check for hash table half full
! 02.12.2011	ggu	initialize hashin to zero
! 12.01.2026    ggu     large arrays transformed to allocatable
!
!******************************************************************

! integer function insert(key,info) inserts element into hash table
! integer function retriv(key,info) retrieves element from hash table
! integer function delete(key,info) deletes element from hash table
! integer function visit(key,info)  visits hash table
! subroutine reset                  resets hash table for visiting
! subroutine init                   initializes hash table
! subroutine copyi(infos,infod)     copies from source (s) to destination (d)
! subroutine test                   tests hash routines
!
! arguments :
!
! integer function insert(key,info) in: key,info
! integer function retriv(key,info) in: key  out: info
! integer function delete(key,info) in: key  out: info
! integer function visit(key,info)  out: key,info
!
! return values :
!
! insert	1: inserted  0: already there  -1: error
! retriv	1: found     0: not found      -1: error
! delete	1: deleted   0: not there      -1: error
! visit		1: found     0: finish

!******************************************************************

	integer function insert(key,info)

! inserts element into hash table

	implicit none

	integer key
	integer info
	integer hashin

	insert=hashin(key,info,1)

	return
	end

!******************************************************************

	integer function retriv(key,info)

! retrieves element from hash table

	implicit none

	integer key
	integer info
	integer hashin

	retriv=hashin(key,info,2)

	return
	end

!******************************************************************

	integer function delete(key,info)

! deletes element from hash table

	implicit none

	integer key
	integer info
	integer hashin

	delete=hashin(key,info,3)

	return
	end

!******************************************************************

	integer function visit(key,info)

! visits hash table

	implicit none

	integer key
	integer info
	integer hashin

	visit=hashin(key,info,4)

	return
	end

!******************************************************************

	subroutine reset

! resets hash table for visiting

	implicit none

	integer key
	integer info
	integer hashin
	integer dummy

	dummy=hashin(key,info,5)

	return
	end
!******************************************************************

	subroutine init

! initializes hash table

	implicit none

	integer key
	integer info
	integer hashin
	integer dummy

	dummy=hashin(key,info,0)

	return
	end

!******************************************************************

	subroutine copyi(infos,infod)

! copies from source (infos) to destination (infod)

	implicit none

	integer infos,infod

	infod=infos

	return
	end

!******************************************************************
!	internal routines
!******************************************************************

	integer function hashin(key,info,mode)

! internal hash table function

! mode :
!		0	initialize
!		1	insert
!		2	retrieve
!		3	delete
!		4	visit
!		5	reset for visiting
!
! return :
!		-1	error
!		0	not found, nothing more, already there
!		1	ok
! flags :
!		flagf	flags free entry
!		flagd	flags deleted entry (may be reused)
!
!	-> both flags must be different from
!		any possible key used

	implicit none

	integer ndim,ndim2,flagf,flagd
! using a prime number for ndim is a save choice
!	parameter (ndim=7957)
!	parameter (ndim=21973)
!	parameter (ndim=54563)
	parameter (ndim=99991)
	parameter (ndim2=(ndim-1)/2)
	parameter (flagf=0,flagd=-1)

	integer key
	integer mode
	integer info

	integer, allocatable :: keyt(:)
	integer, allocatable :: infot(:)

	integer i,ipos,start
	integer itotal

	integer lookup,nfree

	save ipos,keyt,infot,itotal
	data ipos /1/
	data itotal /0/

	allocate(keyt(ndim))
	allocate(infot(ndim))

	hashin=0

	start=mod(key,ndim)+1

	if(mode.eq.1) then 	!-------------------------- insert
	  if( itotal .gt. ndim2 ) goto 95
	  i=lookup(key,keyt,start,flagf,ndim)
	  if(i.eq.0) then	!not there -> insert
	    i=nfree(keyt,start,flagf,flagd,ndim)
	  else if(i.gt.0) then	!key already there
	    i=0
	  end if
	  if(i.gt.0) then	!insert
	    keyt(i)=key
	    call copyi(info,infot(i))
	    itotal = itotal + 1
	    hashin=1
	  else
	    hashin=i
	  end if
	else if(mode.eq.2) then	!-------------------------- retrieve
	  i=lookup(key,keyt,start,flagf,ndim)
	  if(i.gt.0) then	!found
	    call copyi(infot(i),info)
	    hashin=1
	  else
	    hashin=i
	  end if
	else if(mode.eq.3) then	!-------------------------- delete
	  i=lookup(key,keyt,start,flagf,ndim)
	  if(i.gt.0) then	!found
	    keyt(i)=flagd
	    itotal = itotal - 1
	    hashin=1
	  else
	    hashin=i
	  end if
	else if(mode.eq.4) then	!-------------------------- visit
	  do while(ipos.le.ndim.and. &
     &		(keyt(ipos).eq.flagf.or.keyt(ipos).eq.flagd))
	    ipos=ipos+1
	  end do
	  if(ipos.le.ndim) then
	    key=keyt(ipos)
	    call copyi(infot(ipos),info)
	    ipos=ipos+1
	    hashin=1
	  else
	    hashin=0
	  end if
	else if(mode.eq.5) then	!-------------------------- reset
	  ipos=1
	else if(mode.eq.0) then	!-------------------------- init
	  ipos=1
	  do i=1,ndim
	    keyt(i)=flagf
	  end do
	end if

	return
   95	continue
	write(6,*) 'hash table more than half full... ',itotal
	write(6,*) 'hash table dimension is = ',ndim
	stop 'error stop hashin: dimension of hash table'
	end

!***********************************************************

	integer function lookup(key,keyt,start,flag,ndim)

! looks up index for key in keyt()

! starts search at index start
! flag flags free index
!
! return value :
!			i	index in keyt 
!			0	not found
!			-1	error

	implicit none

	integer key,start
	integer keyt(1)
	integer flag,ndim

	integer icount,inc
	integer i
	integer ndim2

	icount=0
	inc=1
	i=start
	ndim2 = ndim / 2

	lookup=-999
	do while(lookup.eq.-999)
	  if(keyt(i).eq.flag) then
	    lookup=0
	  else if(keyt(i).eq.key) then
	    lookup=i
	  else if(icount.gt.ndim2) then
	    lookup=-1
	  end if
	  i=i+inc
	  if( i .gt. ndim ) i = mod(i,ndim)
	  inc=inc+2
	  icount=icount+1
	end do

	return
	end

!***********************************************************

	integer function nfree(keyt,start,flagf,flagd,ndim)

! looks for first free index in keyt()

! starts search at index start
! flagf, flagd flag free or deleted index
!
! return value :
!			i	index in keyt 
!			-1	error

	implicit none

	integer start
	integer keyt(1)
	integer flagf,flagd,ndim

	integer icount,inc
	integer i
	integer ndim2

	icount=0
	inc=1
	i=start
	ndim2 = ndim / 2

	nfree=-999
	do while(nfree.eq.-999)
	  if(keyt(i).eq.flagf) then
	    nfree=i
	  else if(keyt(i).eq.flagd) then
	    nfree=i
	  else if(icount.gt.ndim2) then
	    nfree=-1
	  end if
	  i=i+inc
	  if( i .gt. ndim ) i = mod(i,ndim)
	  inc=inc+2
	  icount=icount+1
	end do

	return
	end

!*************************************************************

	subroutine test

! tests hash routines

	implicit none

	integer mode
	integer key,info
	integer i

	integer insert,retriv,delete,visit

	call init

	mode=1
	do while(mode.ne.0)
	  write(6,*) '0 exit  1 insert  2 retriv  3 delete  4 visit'
	  read(5,'(i10)') mode

	  if(mode.eq.1) then
		write(6,*) 'key,info :'
		read(5,*) key,info
		i=insert(key,info)
		write(6,*) i,key,info
	  else if(mode.eq.2) then
		write(6,*) 'key :'
		read(5,*) key
		i=retriv(key,info)
		write(6,*) i,key,info
	  else if(mode.eq.3) then
		write(6,*) 'key :'
		read(5,*) key
		i=delete(key,info)
		write(6,*) i,key,info
	  else if(mode.eq.4) then
		call reset
		i=visit(key,info)
		do while(i.gt.0)
		  write(6,*) i,key,info
		  i=visit(key,info)
		end do
	  end if
	end do

	end

!********************************************************************

!	program test
!	call test
!	end
