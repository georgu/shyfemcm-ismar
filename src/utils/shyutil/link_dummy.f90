
!--------------------------------------------------------------------------
!
!    Copyright (C) 2025  Georg Umgiesser
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

! dummy links
!
! revision log :
!
! 15.06.2025    ggu     updated
!
!--------------------------------------------------------------
!
! dummy program not used anymore with new connection framework
!
! update_ielt() should be transfered to lnku.f90
!
!--------------------------------------------------------------

!****************************************************************

        subroutine make_links_old(nkn,nel,nen3v)
	implicit none
	integer nkn,nel,nen3v(3,nel)
	end

        subroutine make_links(nkn,nel,nen3v,ibound,kerr)
	implicit none
	integer nkn,nel,nen3v(3,nel),ibound(nkn),kerr
	end

        subroutine checklenk(nkn,ilinkv,lenkv,nen3v)
	implicit none
	integer nkn,ilinkv(*),lenkv(*),nen3v(3,*)
	end

        subroutine checklink(nkn,ilinkv,linkv)
	implicit none
	integer nkn,ilinkv(*),linkv(*)
	end

        subroutine checkkant(nkn,kantv)
	implicit none
	integer nkn,kantv(2,nkn)
	end

        subroutine checkielt(nel,ieltv,nen3v)
	implicit none
	integer nel,ieltv(3,nel),nen3v(3,nel)
	end

	subroutine link_set_write(bset)
	implicit none
	logical bset
	end

	subroutine link_set_stop(bset)
	implicit none
	logical bset
	end

!****************************************************************

        subroutine update_ielt(nkn,nel,ibound,ieltv,nen3v)

! updates vector ieltv with open boundary nodes

        implicit none

        integer nkn
        integer nel
        integer ibound(nkn)       ! >0 => open boundary node
        integer ieltv(3,nel)
        integer nen3v(3,nel)

        integer k,ie,ii,i
        integer ksnext,isbhnd,ksthis

        do ie=1,nel
          do ii=1,3
            k = ksthis(ii,ie,nen3v)
            if( ibound(k) .gt. 0 .and. ibound(ksnext(k,ie,nen3v)) .gt. 0 ) then
              i = isbhnd(k,ie,nen3v)
              ieltv(i,ie) = -1
            end if
          end do
        end do

        end

!****************************************************************

