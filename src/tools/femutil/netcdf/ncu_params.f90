
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

! routines to store variables for output
!
! revision log :
!
! 12.11.2026    ggu     written from scratch

!===================================================================
        module mod_ncuparams
!===================================================================

        logical, save :: bdebug = .false.
        logical, save :: bverbose = .false.
        logical, save :: bquiet = .false.
        logical, save :: bsilent = .false.

!===================================================================
        end module mod_ncuparams
!===================================================================

        subroutine set_ncu_params(bdebug_l,bverbose_l,bquiet_l,bsilent_l)

        use mod_ncuparams

        implicit none

        logical bdebug_l,bverbose_l,bquiet_l,bsilent_l

        bdebug = bdebug_l
        bverbose = bverbose_l
        bquiet = bquiet_l
        bsilent = bsilent_l

        end

!*******************************************************************

        subroutine get_ncu_params(bdebug_l,bverbose_l,bquiet_l,bsilent_l)

        use mod_ncuparams

        implicit none

        logical bdebug_l,bverbose_l,bquiet_l,bsilent_l

        bdebug_l = bdebug
        bverbose_l = bverbose
        bquiet_l = bquiet
        bsilent_l = bsilent

        end

!*******************************************************************

