
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

!Butcher tableaux
        real, allocatable, save :: a_erk(:,:)	    		!explicit A matrix
        real, allocatable, save :: a_irk(:,:)       		!implicit A matrix
        real, allocatable, save :: a_srk(:,:)     		!stiffly implicit A matrix
        real, allocatable, save :: b_erk(:)	    		!explicit b vector
        real, allocatable, save :: b_irk(:)	    		!implicit b vector
        real, allocatable, save :: b_srk(:)	    		!stiffly implicit b vector
        real, allocatable, save :: c_rk(:)	    		!implicit/explicit c vectors:
                                                    		!they are equal for all im/ex schemes in additive runge-kutta
!specific Butcher tableaux for concentrations
!vertical concentration advection is discretized by means of additional weights in the semi-implicit method.
!For high-order ImEx schemes only an explicit discretization is available, by now. 
        real, allocatable, save :: a_crk(:,:)       		!generic A matrix for conc. vert adv (can be im or ex)
        real, allocatable, save :: b_crk(:)	    		!generic b vector for conc. vert adv (can be im or ex)

	integer, save :: n_rkstages = 0		    		!number of (non-trivial) stages
	integer, save :: rkscheme = 0		    		!scheme triplet

        double precision, allocatable, save :: uverk_reg(:,:,:) !register vector for explicit momentum rhs
        double precision, allocatable, save :: uvirk_reg(:,:,:) !register vector for implicit momentum rhs
        double precision, allocatable, save :: uvsrk_reg(:,:,:) !register vector for stiffly implicit momentum rhs
        double precision, allocatable, save :: urk_reg(:,:,:)   !register vector for u-transport
        double precision, allocatable, save :: vrk_reg(:,:,:)   !register vector for v-transport
        double precision, allocatable, save :: wrk_reg(:,:,:)   !register vector for w-velocity

        double precision, allocatable, save :: terk_reg(:,:,:) 	!register vector for explicit temperature rhs
        double precision, allocatable, save :: tcrk_reg(:,:,:) 	!register vector for implicit temperature rhs
        double precision, allocatable, save :: tsrk_reg(:,:,:) 	!register vector for stiffly implicit temperature rhs
        double precision, allocatable, save :: serk_reg(:,:,:) 	!register vector for explicit salinity rhs
        double precision, allocatable, save :: scrk_reg(:,:,:) 	!register vector for implicit salinity rhs
        double precision, allocatable, save :: ssrk_reg(:,:,:) 	!register vector for stiffly implicit salinity rhs
        double precision, allocatable, save :: cerk_reg(:,:,:,:)!register vector for explicit concentration rhs
        double precision, allocatable, save :: ccrk_reg(:,:,:,:)!register vector for implicit concentration rhs
        double precision, allocatable, save :: csrk_reg(:,:,:,:)!register vector for stiffly implicit concentration rhs

!==========================================================================
        contains
!==========================================================================

	subroutine mod_rungekutta_init(nkn,nel,nlv)

        integer nkn, nel, nlv

	double precision am,at,az,af,av,ad,aa,gamma,chi,a32
	real getpar
	integer itemp,isalt,iconz

	!The ImEx triplet $(s, \sigma, p)$ identifies a scheme where:
	!1/ $s$ is the number of non-trivial stages of the implicit scheme,
	!2/ $\sigma$ is the number of non-trivial stages of the explicit scheme
	!3/ $p$ is the combined order of accuracy.
	rkscheme = nint(getpar('rkscheme'))

	itemp = nint(getpar('itemp'))
	isalt = nint(getpar('isalt'))
	iconz = nint(getpar('iconz'))

        if( rkscheme .ne. 111  .and.  &
     &      rkscheme .ne. 222  .and.  &
     &      rkscheme .ne. 2221 .and.  &
     &      rkscheme .ne. 232  .and.  &
     &      rkscheme .ne. 33 ) then
          write(6,*) 'runge-kutta triplet: ', rkscheme
          stop 'error stop mod_rungekutta_init: incompatible params'
        end if

	if( n_rkstages > 0 ) then
          deallocate(a_erk)
          deallocate(a_irk)
          deallocate(a_srk)
          deallocate(a_crk)
          deallocate(b_erk)
          deallocate(b_irk)
          deallocate(b_srk)
          deallocate(b_crk)
          deallocate(c_rk)

          deallocate(uverk_reg)
          deallocate(uvirk_reg)
          deallocate(uvsrk_reg)

          deallocate(urk_reg)
          deallocate(vrk_reg)
          deallocate(wrk_reg)

	  if (itemp == 1) then
            deallocate(terk_reg)
            deallocate(tcrk_reg)
            deallocate(tsrk_reg)
	  end if
	  if (isalt == 1) then
            deallocate(serk_reg)
            deallocate(scrk_reg)
            deallocate(ssrk_reg)
	  end if
	  if (iconz == 1) then
            deallocate(cerk_reg)
            deallocate(ccrk_reg)
            deallocate(csrk_reg)
	  end if
        end if

	!The classical ARK(1,1,1) method of Ascher,Ruuth\&Spiteri:
	!It is composed of the pair theta−method and Forward Euler.
	!The theta-method uses a coefficient $am$ read from file.
	!Notice that for vertical viscous terms uses a different coefficient
	! $at$, so that an L-stable Backward Euler is used, by default with
	!weight $at=1$.
	if (rkscheme == 111) then
	  am=getpar('ampar')
	  at=getpar('atpar')
	  aa=getpar('aapar')

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
          ad = getpar('adpar')
	  if( ad .ne. at ) then
            write(6,*) 'You are using an obsolete configuration with:'
            write(6,*) 'adpar: ', ad
            write(6,*) 'atpar: ', at
            write(6,*) 'adpar and atpar must be equal.'
            stop 'error stop mod_rungekutta_init: incompatible params'
          end if

          n_rkstages = 1

          allocate (a_erk(n_rkstages,n_rkstages))
          allocate (a_irk(n_rkstages,n_rkstages+1))
          allocate (a_srk(n_rkstages,n_rkstages+1))
          allocate (a_crk(n_rkstages,n_rkstages+1))
          allocate (b_erk(n_rkstages))
          allocate (b_irk(n_rkstages+1))
          allocate (b_srk(n_rkstages+1))
          allocate (b_crk(n_rkstages+1))
          allocate (c_rk(n_rkstages))

	  a_erk(1,1) = 1.

	  a_irk(1,1) = 1.-am
	  a_irk(1,2) = am

	  a_srk(1,1) = 1.-at
	  a_srk(1,2) = at

	  a_crk(1,1) = 1.-aa
	  a_crk(1,2) = aa

	  b_erk(1)  = 1.

	  b_irk(1)  = 1.-am
	  b_irk(2)  = am

	  b_srk(1)  = 1.-at
	  b_srk(2)  = at

	  b_crk(1)  = 1.-aa
	  b_crk(2)  = aa

	  c_rk(1)  = 1.

	!The classical ARK(2,2,2) method of Ascher,Ruuth\&Spiteri:
	!It is composed of the pair second order DIRK22 and its
	!explicit companion scheme ERK22. The same implicit scheme
	!is used for all the stiff terms.
	else if (rkscheme == 222) then
	  gamma=getpar('gapar')
	  if( gamma <= 0. .or. gamma >= 1. ) then
            write(6,*) 'You are using the rkscheme=222 with:'
            write(6,*) 'gapar: ', gamma
	    write(6,*) 'this parameter is not allowed'
	    write(6,*) 'Please use'
	    write(6,*) '  gapar>0 and gapar<1'
            stop 'error stop mod_rungekutta_init: incompatible params'
          end if

          n_rkstages = 2

          allocate (a_erk(n_rkstages,n_rkstages))
          allocate (a_irk(n_rkstages,n_rkstages+1))
          allocate (a_srk(n_rkstages,n_rkstages+1))
          allocate (a_crk(n_rkstages,n_rkstages+1))
          allocate (b_erk(n_rkstages))
          allocate (b_irk(n_rkstages+1))
          allocate (b_srk(n_rkstages+1))
          allocate (b_crk(n_rkstages+1))
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

	  a_srk = a_irk
	  a_crk = 0.
	  a_crk(:,1:n_rkstages) = a_erk

	  b_erk(1)  = 1.-1./(2.*gamma)
	  b_erk(2)  = 1./(2.*gamma)

	  b_irk(1)  = 0.
	  b_irk(2)  = 1.-gamma
	  b_irk(3)  = gamma

	  b_srk  = b_irk
	  b_crk  = 0.
	  b_crk(1:n_rkstages) = b_erk

	  c_rk(1)  = gamma
	  c_rk(2)  = 1.

	!Another ARK(2,2,2) of (Noelle.2014). The implicit scheme is
	!composed of a first stage with Heun’s method followed by a
	!second stage of the Cranck-Nicholson scheme and the explicit
	!method is the midpoint scheme.
	else if (rkscheme == 2221) then
          n_rkstages = 2

          allocate (a_erk(n_rkstages,n_rkstages))
          allocate (a_irk(n_rkstages,n_rkstages+1))
          allocate (a_srk(n_rkstages,n_rkstages+1))
          allocate (a_crk(n_rkstages,n_rkstages+1))
          allocate (b_erk(n_rkstages))
          allocate (b_irk(n_rkstages+1))
          allocate (b_srk(n_rkstages+1))
          allocate (b_crk(n_rkstages+1))
          allocate (c_rk(n_rkstages))

	  a_erk(1,1) = 0.5
	  a_erk(1,2) = 0.
	  a_erk(2,1) = 0.
	  a_erk(2,2) = 1.

	  a_irk(1,1) = 0.
	  a_irk(1,2) = 0.5
	  a_irk(1,3) = 0.
	  a_irk(2,1) = 0.5
	  a_irk(2,2) = 0.
	  a_irk(2,3) = 0.5

	  a_srk = a_irk
	  a_crk = 0.
	  a_crk(:,1:n_rkstages) = a_erk

	  b_erk(1)  = 0.
	  b_erk(2)  = 1.

	  b_irk(1)  = 0.5
	  b_irk(2)  = 0.
	  b_irk(3)  = 0.5

	  b_srk  = b_irk
	  b_crk  = 0.
	  b_crk(1:n_rkstages) = b_erk

	  c_rk(1)  = 0.5
	  c_rk(2)  = 1.

	!The ARK(2,3,2) where the implicit scheme is TR-BDF2 introduced
	!in (Bank,1985) and analysed in (Hosea,1996) and the explicit scheme
	!is designed to match the coupling and order conditions (Giraldo,2013).
	!This scheme preserve invariants and it is L-stable.
	else if (rkscheme == 232) then
	  chi=getpar('chipar')
	  if( chi <= 0. .or. chi >= 1. ) then
            write(6,*) 'You are using th rkscheme=232 with:'
            write(6,*) 'chipar: ', chi
	    write(6,*) 'this parameter is not allowed'
	    write(6,*) 'Please use'
	    write(6,*) '  chipar>0 and chipar<1'
            stop 'error stop mod_rungekutta_init: incompatible params'
          end if
	  a32=getpar('a32par')
	  if( a32 <= 0. .or. a32 >= 1. ) then
            write(6,*) 'You are using th rkscheme=232 with:'
            write(6,*) 'a32par: ', a32
	    write(6,*) 'this parameter is not allowed'
	    write(6,*) 'Please use'
	    write(6,*) '  a32par>0 and a32par<1'
            stop 'error stop mod_rungekutta_init: incompatible params'
          end if

          n_rkstages = 3

          allocate (a_erk(n_rkstages,n_rkstages))
          allocate (a_irk(n_rkstages,n_rkstages+1))
          allocate (a_srk(n_rkstages,n_rkstages+1))
          allocate (a_crk(n_rkstages,n_rkstages+1))
          allocate (b_erk(n_rkstages))
          allocate (b_irk(n_rkstages+1))
          allocate (b_srk(n_rkstages+1))
          allocate (b_crk(n_rkstages+1))
          allocate (c_rk(n_rkstages))

	  a_erk(1,1) = 2.*chi
	  a_erk(1,2) = 0.
	  a_erk(1,3) = 0.
	  a_erk(2,1) = 1.-a32
	  a_erk(2,2) = a32
	  a_erk(2,3) = 0.
	  a_erk(3,1) = 1.-(1.-2.*chi)/(4.*chi)-chi
	  a_erk(3,2) = (1.-2.*chi)/(4.*chi)
	  a_erk(3,3) = chi

	  a_irk(1,1) = chi
	  a_irk(1,2) = chi
	  a_irk(1,3) = 0.
	  a_irk(1,4) = 0.
	  a_irk(2,1) = 1.-(1.-2.*chi)/(4.*chi)-chi
	  a_irk(2,2) = (1.-2.*chi)/(4.*chi)
	  a_irk(2,3) = chi
	  a_irk(2,4) = 0.
	  a_irk(3,1) = 1.-(1.-2.*chi)/(4.*chi)-chi
	  a_irk(3,2) = (1.-2.*chi)/(4.*chi)
	  a_irk(3,3) = chi
	  a_irk(3,4) = 0.

	  a_srk = a_irk
	  a_crk = 0.
	  a_crk(:,1:n_rkstages) = a_erk

	  b_erk(1)  = 1.-(1.-2.*chi)/(4.*chi)-chi
	  b_erk(2)  = (1.-2.*chi)/(4.*chi)
	  b_erk(3)  = chi

	  b_irk(1)  = 1.-(1.-2.*chi)/(4.*chi)-chi
	  b_irk(2)  = (1.-2.*chi)/(4.*chi)
	  b_irk(3)  = chi
	  b_irk(4)  = 0.

	  b_srk  = b_irk
	  b_crk  = 0.
	  b_crk(1:n_rkstages) = b_erk

	  c_rk(1)  = 2*chi
	  c_rk(2)  = 1.
	  c_rk(3)  = 1.

	!The third order Strong Stability Preserving SSP scheme.
	!The implicit scheme is coded as the explicit one with zero
	!element on the diagonal of the Butcher Tableux. Being fully
	!explicit this scheme is highly inefficient and it is used
	!only for testing or computing reference solutions. It is
	!not documented in the manual.
	else if (rkscheme == 33) then
          n_rkstages = 3

          allocate (a_erk(n_rkstages,n_rkstages))
          allocate (a_irk(n_rkstages,n_rkstages+1))
          allocate (a_srk(n_rkstages,n_rkstages+1))
          allocate (a_crk(n_rkstages,n_rkstages+1))
          allocate (b_erk(n_rkstages))
          allocate (b_irk(n_rkstages+1))
          allocate (b_srk(n_rkstages+1))
          allocate (b_crk(n_rkstages+1))
          allocate (c_rk(n_rkstages))

	  a_erk(1,1) = 1.
	  a_erk(1,2) = 0.
	  a_erk(1,3) = 0.
	  a_erk(2,1) = 1./4.
	  a_erk(2,2) = 1./4.
	  a_erk(2,3) = 0.
	  a_erk(3,1) = 1./6.
	  a_erk(3,2) = 1./6.
	  a_erk(3,3) = 2./3.

	  a_irk(1,1) = 1.
	  a_irk(1,2) = 0.
	  a_irk(1,3) = 0.
	  a_irk(1,4) = 0.
	  a_irk(2,1) = 1./4.
	  a_irk(2,2) = 1./4.
	  a_irk(2,3) = 0.
	  a_irk(2,4) = 0.
	  a_irk(3,1) = 1./6.
	  a_irk(3,2) = 1./6.
	  a_irk(3,3) = 2./3.
	  a_irk(3,4) = 0.

	  a_srk = a_irk
	  a_crk = 0.
	  a_crk(:,1:n_rkstages) = a_erk

	  b_erk(1)  = 1./6.
	  b_erk(2)  = 1./6.
	  b_erk(3)  = 2./3.

	  b_irk(1)  = 1./6.
	  b_irk(2)  = 1./6.
	  b_irk(3)  = 2./3.
	  b_irk(4)  = 0.

	  b_srk  = b_irk
	  b_crk  = 0.
	  b_crk(1:n_rkstages) = b_erk

	  c_rk(1)  = 1.
	  c_rk(2)  = 1./2.
	  c_rk(3)  = 1.

	end if

	!Register vectors to store stage right-hand sides
	!for the momentum and stage transport/velocity components
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

	  if (itemp == 1) then
            allocate(terk_reg(nlv,nkn,n_rkstages-1))
            allocate(tcrk_reg(nlv,nkn,n_rkstages-1))
            allocate(tsrk_reg(nlv,nkn,n_rkstages-1))

            terk_reg = 0.
            tcrk_reg = 0.
            tsrk_reg = 0.
	  end if
	  if (isalt == 1) then
            allocate(serk_reg(nlv,nkn,n_rkstages-1))
            allocate(scrk_reg(nlv,nkn,n_rkstages-1))
            allocate(ssrk_reg(nlv,nkn,n_rkstages-1))

            serk_reg = 0.
            scrk_reg = 0.
            ssrk_reg = 0.
	  end if
	  if (iconz > 0) then !iconz index runs slowest
            allocate(cerk_reg(nlv,nkn,n_rkstages-1,iconz))
            allocate(ccrk_reg(nlv,nkn,n_rkstages-1,iconz))
            allocate(csrk_reg(nlv,nkn,n_rkstages-1,iconz))

            cerk_reg = 0.
            ccrk_reg = 0.
            csrk_reg = 0.
	  end if
        end if


        end subroutine mod_rungekutta_init

!**************************************************************************

	subroutine get_rungekutta_weights(curr_stage,crk,coeff_erk, &
     &    coeff_irk,coeff_srk)

! returns Butcher Tableaux
! c_rk and a_rk for inner stages
! c_rk and b_rk for final stage

	integer, intent(in) :: curr_stage
	real, intent(out) :: crk
	real, dimension(n_rkstages), intent(out) :: coeff_erk
	real, dimension(n_rkstages+1), intent(out) :: coeff_irk
	real, dimension(n_rkstages+1), intent(out) :: coeff_srk

	crk  = c_rk(curr_stage)
	if ( curr_stage < n_rkstages ) then
	  coeff_erk = a_erk(curr_stage,:)
	  coeff_irk = a_irk(curr_stage,:)
	  coeff_srk = a_srk(curr_stage,:)
	else
	  coeff_erk = b_erk(:)
	  coeff_irk = b_irk(:)
	  coeff_srk = b_srk(:)
	end if

	end subroutine get_rungekutta_weights

!**************************************************************************

	subroutine get_rungekutta_weights_tracer(curr_stage,coeff_crk)

! returns additional Butcher Tableaux for tracers/concentrations
! a_rk for inner stages
! b_rk for final stage

	integer, intent(in) :: curr_stage
	real, dimension(n_rkstages+1), intent(out) :: coeff_crk

	if ( curr_stage < n_rkstages ) then
	  coeff_crk = a_crk(curr_stage,:)
	else
	  coeff_crk = b_crk(:)
	end if

	end subroutine get_rungekutta_weights_tracer

!==========================================================================
        end module mod_rungekutta
!==========================================================================

!**************************************************************************

!**************************************************************************

