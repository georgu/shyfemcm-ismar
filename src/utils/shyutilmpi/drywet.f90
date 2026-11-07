
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

!------------------------------------------------------------------
!
! routines for encrytption of iweg in to iwng (nodal variable)
!
! revision log :
!
! 07.11.2026    ggu     started from scratch
!
! notes :
!
! module mod_connect must be initialized
!
!------------------------------------------------------------------

!******************************************************************

	subroutine drywet_encrypt(iweg,iwng)

! encrypt iweg and puts results into iwng

	use basin
	use mod_connect
	use shympi

	implicit none

	integer iweg(nel)
	integer iwng(nkn)

	integer k,ie,j,ne,bits,idry
	integer, allocatable :: nlist(:,:)
	integer, allocatable :: elist(:,:)

	allocate(nlist(0:ngr,nkn),elist(0:ngr,nkn))
	call connect_get_lists(nkn,ngr,nlist,elist)

	do k=1,nkn
	  ne = elist(0,k)
	  bits = 0
	  idry = 0
	  do j=1,ne
	    ie = elist(j,k)
	    if( iweg(ie) > 0 ) then
	      bits = ibset(bits,j)
	      idry = idry + 1
	    end if
	  end do
	  if( idry == ne ) bits = ibset(bits,0)	!all elems around k are dry
	  iwng(k) = bits
	end do

	call shympi_exchange_2d_node(iwng)

	end

!******************************************************************

	subroutine drywet_decrypt(iweg,iwng)

! decrypt iwng and puts results into iweg

	use basin
	use mod_connect
	use shympi

	implicit none

	integer iweg(nel)
	integer iwng(nkn)

	integer k,ie,j,ne,bits
	integer, allocatable :: nlist(:,:)
	integer, allocatable :: elist(:,:)

	allocate(nlist(0:ngr,nkn),elist(0:ngr,nkn))
	call connect_get_lists(nkn,ngr,nlist,elist)

	iweg = 0

	do k=1,nkn
	  ne = elist(0,k)
	  bits = iwng(k)
	  do j=1,ne
	    if( btest(bits,j) ) then
	      ie = elist(j,k)
	      iweg(ie) = iweg(ie) + 1
	    end if
	  end do
	end do

	call shympi_exchange_2d_elem(iweg)

	do ie=1,nel
	  if( iweg(ie) == 0 .or. iweg(ie) == 3 ) then
	    !ok
	  else
	    write(6,*) ie,iweg(ie)
	    write(6,*) 'iweg must be 0 or 3'
	    stop 'error stop drywet_decrypt: internal error (1)'
	  end if
	end do

	end

!******************************************************************

	subroutine drywet_drynode(iwng,idrynode)

! sets idrynode to 1 if node completely dry (no wet elelents attached)

	use basin
	use mod_connect

	implicit none

	integer iwng(nkn)
	integer idrynode(nkn)

	integer k

	idrynode = 0

	do k=1,nkn
	  if( btest(iwng(k),0) ) idrynode(k) = 1
	end do

	end

!******************************************************************

	subroutine drywet_test_encryption(iweg)

! test encryption of iwng

	use basin
	use mod_connect

	implicit none

	integer iweg(nel)

	integer ic
	integer, allocatable :: iwng(:)
	integer, allocatable :: ieaux(:)

	allocate(iwng(nkn),ieaux(nel))

	call drywet_encrypt(iweg,iwng)
	call drywet_decrypt(ieaux,iwng)

	if( any(iweg/=ieaux) ) then
	  ic = count(iweg/=ieaux)
	  write(6,*) 'there are errors in encryption drywet: ',ic
	  stop 'error stop drywet_test_encryption: internal error (2)'
	end if

	end

!******************************************************************

