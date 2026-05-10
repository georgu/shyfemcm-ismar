
!--------------------------------------------------------------------------
!
!    Copyright (C) 2013,2015,2019-2020  Georg Umgiesser
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

! revision log :
!
! 28.04.2026	lrp	written from scratch
! 04.05.2026	ggu	new routine get_rungekutta_crk()

!==========================================================================
	module mod_rungekutta
!==========================================================================

	implicit none

        real, allocatable, save :: a_erk(:,:)	    !explicit A matrix
        real, allocatable, save :: a_irk(:,:)       !implicit A matrix
        real, allocatable, save :: a_srk(:,:)       !stiffly implicit A matrix
        real, allocatable, save :: b_erk(:)	    !explicit b vector
        real, allocatable, save :: b_irk(:)	    !implicit b vector
        real, allocatable, save :: b_srk(:)	    !stiffly implicit b vector
        real, allocatable, save :: c_rk(:)	    !implicit/explicit c vectors:
                                                    !they are equal for additive runge-kutta

	integer, save :: n_rkstages = 0		    !number of (non-trivial) stages

        double precision, allocatable, save :: uverk_reg(:,:,:) !register vector for explicit momentum
        double precision, allocatable, save :: uvirk_reg(:,:,:) !register vector for implicit momentum
        double precision, allocatable, save :: uvsrk_reg(:,:,:) !register vector for stiffly implicit momentum
        double precision, allocatable, save :: urk_reg(:,:,:)   !register vector for u-transport
        double precision, allocatable, save :: vrk_reg(:,:,:)   !register vector for v-transport
        double precision, allocatable, save :: wrk_reg(:,:,:)   !register vector for w-velocity

!==========================================================================
        contains
!==========================================================================

	subroutine mod_rungekutta_init(nkn,nel,nlv)

        integer nkn, nel, nlv

	integer rk_triplet
	double precision am,at,az,af,av,gamma
	real getpar

	!The ImEx triplet $(s, \sigma, p)$ identifies a scheme where:
	!1/ $s$ is the number of non-trivial stages of the implicit scheme,
	!2/ $\sigma$ is the number of non-trivial stages of the explicit scheme
	!3/ $p$ is the combined order of accuracy.
	rk_triplet = nint(getpar('rkscheme'))

        if( rk_triplet .ne. 111 .and.  &
     &      rk_triplet .ne. 222 ) then
          write(6,*) 'runge-kutta triplet: ', rk_triplet
          stop 'error stop mod_rungekutta_init: incompatible params'
        end if

	if( n_rkstages > 0 ) then
          deallocate(a_erk)
          deallocate(a_irk)
          deallocate(a_srk)
          deallocate(b_erk)
          deallocate(b_irk)
          deallocate(b_srk)
          deallocate(c_rk)

          deallocate(uverk_reg)
          deallocate(uvirk_reg)
          deallocate(uvsrk_reg)

          deallocate(urk_reg)
          deallocate(vrk_reg)
          deallocate(wrk_reg)
        end if

	!The classical ARK(1,1,1) method of Ascher,Ruuth&Spiteri:
	!It is composed of the pair theta−method and Forward Euler.
	!The theta-method uses a coefficient $am$ read from file.
	!Notice that for vertical viscous terms uses a different coefficient
	! $at$, so that an L-stable Backward Euler is used, by default with
	!weight $at=1$.
	if (rk_triplet == 111) then
	  am=getpar('ampar')
	  at=getpar('atpar')

	  az = getpar('azpar')
	  if( az .ne. am ) then
            write(6,*) 'You are using an obsolete configuration with:'
            write(6,*) 'azpar: ', az
            write(6,*) 'ampar: ', am
            write(6,*) 'ampar and azpar must be equal.'
            stop 'error stop mod_rungekutta_init: incompatible params'
          end if
          af = getpar('afpar')
	  if( af .ne. am ) then
            write(6,*) 'You are using an obsolete configuration with:'
            write(6,*) 'afpar: ', af
            write(6,*) 'ampar: ', am
            write(6,*) 'afpar and azpar must be equal.'
            stop 'error stop mod_rungekutta_init: incompatible params'
          end if
          av = getpar('avpar')
	  if( av .ne. 0. ) then
            write(6,*) 'You are using an obsolete configuration with:'
            write(6,*) 'avpar: ', av
            write(6,*) 'avpar must be zero.'
            stop 'error stop mod_rungekutta_init: incompatible params'
          end if

          n_rkstages = 1

          allocate (a_erk(n_rkstages,n_rkstages))
          allocate (a_irk(n_rkstages,n_rkstages+1))
          allocate (a_srk(n_rkstages,n_rkstages+1))
          allocate (b_erk(n_rkstages+1))
          allocate (b_irk(n_rkstages+1))
          allocate (b_srk(n_rkstages+1))
          allocate (c_rk(n_rkstages))

	  a_erk(1,1) = 1.

	  a_irk(1,1) = 1.-am
	  a_irk(1,2) = am

	  a_srk(1,1) = 1.-at
	  a_srk(1,2) = at

	  b_erk(1)  = 1.
	  b_erk(2)  = 0.

	  b_irk(1)  = 1.-am
	  b_irk(2)  = am

	  b_srk(1)  = 1.-at
	  b_srk(2)  = at

	  c_rk(1)  = 1.

	!The classical ARK(2,2,2) method of Ascher,Ruuth&Spiteri:
	!It is composed of the pair second order DIRK22 and its
	!explicit companion scheme ERK22. The same implicit scheme
	!is used for all the stiff terms.
	else if (rk_triplet == 222) then
	  gamma=getpar('gapar')
	  if( gamma <= 0. .or. gamma >= 1. ) then
            write(6,*) 'You are using th rkscheme=222 with:'
            write(6,*) 'gapar: ', gamma
	    write(6,*) 'this parameter is not allowed'
	    write(6,*) 'Please use'
	    write(6,*) '  gapar>0 and ga<1'
            stop 'error stop mod_rungekutta_init: incompatible params'
          end if

          n_rkstages = 2

          allocate (a_erk(n_rkstages,n_rkstages))
          allocate (a_irk(n_rkstages,n_rkstages+1))
          allocate (a_srk(n_rkstages,n_rkstages+1))
          allocate (b_erk(n_rkstages+1))
          allocate (b_irk(n_rkstages+1))
          allocate (b_srk(n_rkstages+1))
          allocate (c_rk(n_rkstages))

	  a_erk(1,1) = gamma
	  a_erk(1,2) = 0.
	  a_erk(2,1) = 1.-1./(2.*gamma)
	  a_erk(2,2) = 1./(2.*gamma)

	  a_irk(1,1) = 0.
	  a_irk(1,2) = gamma
	  a_irk(1,3) = 0.
	  a_irk(2,1) = 0.
	  a_irk(2,2) = 1.-gamma
	  a_irk(2,3) = gamma

	  a_srk(1,1) = 0.
	  a_srk(1,2) = gamma
	  a_srk(1,3) = 0.
	  a_srk(2,1) = 0.
	  a_srk(2,2) = 1.-gamma
	  a_srk(2,3) = gamma

	  b_erk(1)  = 1.-1./(2.*gamma)
	  b_erk(2)  = 1./(2.*gamma)
	  b_erk(3)  = 0.

	  b_irk(1)  = 0.
	  b_irk(2)  = 1.-gamma
	  b_irk(3)  = gamma

	  b_srk(1)  = 0.
	  b_srk(2)  = 1.-gamma
	  b_srk(3)  = gamma

	  c_rk(1)  = gamma
	  c_rk(2)  = 1.
	end if

        if ( n_rkstages > 1 ) then
          allocate(uverk_reg(2*nlv,nel,n_rkstages-1))
          allocate(uvirk_reg(2*nlv,nel,n_rkstages-1))
          allocate(uvsrk_reg(2*nlv,nel,n_rkstages-1))

          allocate(urk_reg(nlv,nel,n_rkstages-1))
          allocate(vrk_reg(nlv,nel,n_rkstages-1))
          allocate(wrk_reg(0:nlv,nkn,n_rkstages-1))

          uverk_reg = 0.
          uvirk_reg = 0.
          uvsrk_reg = 0.

          urk_reg = 0.
          vrk_reg = 0.
          wrk_reg = 0.
        end if


        end subroutine mod_rungekutta_init

!**************************************************************************

	subroutine get_rungekutta_weights(kstage,crk,aerk,airk,asrk)

! returns c_rk for desired stage

	integer, intent(in) :: kstage
	real, intent(out) :: crk
	real, dimension(n_rkstages), intent(out) :: aerk
	real, dimension(n_rkstages+1), intent(out) :: airk
	real, dimension(n_rkstages+1), intent(out) :: asrk

	crk  = c_rk(kstage)
	aerk = a_erk(kstage,:)
	airk = a_irk(kstage,:)
	asrk = a_srk(kstage,:)

	end subroutine get_rungekutta_weights

!==========================================================================
        end module mod_rungekutta
!==========================================================================

!**************************************************************************

!**************************************************************************

