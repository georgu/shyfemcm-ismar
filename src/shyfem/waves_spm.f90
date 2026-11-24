
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

! parametric wave model
!
! revision log :
!
! 30.05.2026    ggu     extracted from wave file

	subroutine wave_spm(wis,wid,fet,dep,waeh,waep,waed)

!	parametric wave model

	implicit none

	real wis	! wind speed
	real wid	! wind direction
	real fet	! fetch
	real dep	! depth
	real waeh	! wave height (return)
	real waep	! wave period (return)
	real waed	! wave direction (return)

        real ah1,ah2,ah3,eh1,eh2,eh3,eh4
        real at1,at2,at3,et1,et2,et3,et4
!------------------------------------------------------ Hurdle and Stive
!        parameter(ah1=0.25,ah2=0.6,eh1=0.75)
!        parameter(eh2=0.5,ah3=4.3e-05,eh3=1.,eh4=2.)
!        parameter(at1=8.3,at2=0.76,et1=0.375)
!        parameter(et2=1./3.,at3=4.1e-05,et3=1.,et4=3.)
!------------------------------------------------------ SPM
        parameter(ah1=0.283,ah2=0.53,eh1=3./4.)
        parameter(eh2=1.,ah3=0.00565,eh3=1./2.,eh4=1.)
        parameter(at1=7.54,at2=0.833,et1=3./8.)
        parameter(et2=1.,at3=0.0379,et3=1./3.,et4=1.)

	real, parameter :: g = 9.81             !gravity acceleration [m2/s]

	real gh,gx,hg
	real auxh,auxh1,auxt,auxt1
	real hbr

          if( wis > 0. ) then
            gh = (g*dep)/(wis**2.)
            gx = (g*fet)/(wis**2.)
            hg = dep / (g*wis**2.)
            auxh = ah2*gh**eh1
            auxh1 = ah2*hg**eh1
            auxt = at2*gh**et1
            auxt1 = ah2*gx**eh1

            waeh = (tanh(auxh))**eh4
            waeh = (ah3*gx**eh3) / waeh
            waeh = (tanh(waeh))**eh2
            waeh = ah1*tanh(auxh)*waeh
            waeh = waeh * wis**2 / g

            waep = (tanh(auxt))**et4
            waep = (at3*gx**et3) / waep
            waep = (tanh(waep))**et2
            waep = at1*tanh(auxt)*waep
            waep = waep * wis / g
          else
            waeh = 0.
            waep = 0.
          end if

          waed = wid

	  hbr = 0.5 * dep		!limiting wave height
	  if( waeh > hbr ) waeh = hbr

	end

!**********************************************************************

!          waeh(ie) = 0.283 * tanh(0.530*(gh**(3./4.)))*
!     %            tanh((0.00565*(gx**0.5))/
!     %            (tanh(0.530*(gh**(3./4.)))))*((wis**2)/g)
!
!          waep(ie) = 7.54*tanh(0.833*(gh**(3./8.)))*
!     %            tanh((0.0379*(gx**(1./3.)))/
!     %            (tanh(0.833*(gh**(3./8.)))))*(wis/g)
!
!         -----------------------------------------------------------------
!         method of hurdle and stive
!         -----------------------------------------------------------------

!          waeh(ie) = (tanh(auxh))**eh4
!          waeh(ie) = (ah3*gx**eh3) / waeh(ie)
!          waeh(ie) = (tanh(waeh(ie)))**eh2
!          waeh(ie) = ah1*tanh(auxh)*waeh(ie)
!          waeh(ie) = waeh(ie) * wis**2 / g

!          waep(ie) = (tanh(auxt1))**et4
!          waep(ie) = (at3*gx**et3) / waep(ie)
!          waep(ie) = (tanh(waep(ie)))**et2
!          waep(ie) = at1*tanh(auxt)*waep(ie)
!          waep(ie) = waep(ie) * wis / g

!**********************************************************************

	subroutine test_spm

	implicit none

	integer, parameter :: ifmax = 20
	integer, parameter :: ismax = 20
	integer, parameter :: idmax = 20
	integer, parameter :: iffact = 5000
	integer, parameter :: isfact = 1
	integer, parameter :: idfact = 3

	integer is,if,id
	integer iu
	real wis,wid,fet,dep
	real wh,wp,wd
	real hmax

	real aux(idmax)

	wid = 0.
	hmax = 0.

	do if=1,ifmax
	  fet = iffact * if
	  do is=1,ismax
	    wis = isfact * is
	    do id=1,idmax
	      dep = idfact * id
	      call wave_spm(wis,wid,fet,dep,wh,wp,wd)
	      aux(id) = wh
	      hmax = max(hmax,wh)
	    end do
	    iu = 100 + if
	    write(iu,*) is,aux(:)
	  end do
	end do

	write(6,*) 'maximum fetch: ',ifmax*iffact
	write(6,*) 'maximum depth: ',idmax*idfact
	write(6,*) 'maximum speed: ',ismax*isfact
	write(6,*) 'maximum wave height: ',hmax

	end

!**********************************************************************

!	program main_waves_spm
!	call test_spm
!	end program

!**********************************************************************

