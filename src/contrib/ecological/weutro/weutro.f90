
!--------------------------------------------------------------------------
!
!    Copyright (C) 2000,2002-2004,2006,2010,2012,2015-2016  Georg Umgiesser
!    Copyright (C) 2018-2019  Georg Umgiesser
!    Copyright (C) 2002-2004,2013-2014,2016  Donata Melaku Canu
!    Copyright (C) 2003  Isabella Scroccaro
!    Copyright (C) 2006  Michol Ghezzo
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

! weutro - EUTRO box model
!
! contents :
!
! revision log :
!
! 01.01.2000	ggu	original version until April 2002
! 19.04.2002	ggu&dmc	BUGFIX1, print_init, changes in ZOOP
! 20.06.2003	ggu&dmc	new routine for sediments
! 19.12.2003	ggu	new routines param_taranto and param_venezia
! 19.12.2003	isa	new routine rintens_tar()
! 18.02.2004	dmc	dtl non viene passato nel reattore dtl=300/86400
! 26.02.2004	dmc	BUGFIX 2D 3D for depth > 1 (???) (LIGHTFIX)
! 03.03.2004	ggu	in rlim new var btaranto (FIX)
! 04.03.2004	ggu	changes from donata integrated
! 10.06.2004	ggu	new routine param_read for Michol
! 09.08.2004	ggu	marked changes from Donata with GGU
! 17.08.2004	ggu	new routines for debug (eutro_replay)
! 21.08.2004	ggu	rearangments, renaming ($LGA)
! 24.08.2004	ggu	all changes incorporated (see check history 1-11)
! 30.08.2004	ggu	cleanup of settopseg, setbotseg
! 15.02.2006	ggu&mcg	some comments inserted for denitrification (SK18D,SK14D)
! 23.03.2010	ggu	changed v6.1.1
! 08.10.2010	ggu	changed VERS_6_1_13
! 16.03.2012	ggu	dummy restart routines added
! 24.09.2013	dmc	insert direct call to qflux input file from shyfem
! 28.09.2014	dmc	PNO3G1 cancelled
! 10.07.2015	ggu	changed VERS_7_1_50
! 01.04.2016	ggu	changed VERS_7_5_7
! 07.06.2016	ggu	changed VERS_7_5_12
! 17.06.2016	dmc	deleted rlux, link for shyfem 7_5_13
! 27.06.2016	ggu	changed VERS_7_5_16
! 31.08.2018	ggu	changed VERS_7_5_49
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 16.10.2025	ggu	read parameters from namelist - handle_param_init()
!
! notes :
!
! 1
!       debug routines
!       nutlim
!       in steele temp1 = ke and not ke*depth
!       zood ?
!       wsedim, nempar (?)
!
! 2
!       changes signed with GGU
!       ddin1, ddin2, prod, cons
!       loicz()
!       wsedim
!
! 3
!       no denit (?)
!
! 5
!       denit set again
!
! 6
!       new array cold
!       eutro_check
!       param_read
!       wsedim -> cold
!
! 7
!       changes to eutro_check
!       deleted rdtemp()
!       eutro_replay, read/write_debug_param
!       ample re-formatting
!
! 8
!      new default of SOD1D
!      deleted rlim(), rewritten steele
!      changed rintens()
!      deleted rintens_tar()
!
! 9
!       renamed param_init, param_print
!       in source new rlux in arguments
!          set rlghts, call rintens()
!       cleaned up reaeration part
!       consider light attenuation in ditoro()
!       new weutro_test as a test drive for weutro
!       revised steele, rintens
!       moved luxlen, it2n to weutro_light.f
!       moved loicz, wsedim, tempcoef to weutro_sedim.f
!       new comments, reformatting
!
! 10
!       no changes
!
! 11
!       use id, rlux internally (error_check etc.)
!       new routines get/set_internal_param
!       new routine weutro_debug
!       eliminated ITYPE and segtype()
!       call BENTFLUX only if botseg is true
!       topseg used for reaeration
!       debug -> wdebug : weutro_debug(.true.), weutro_debug(.false.)
!
!
!
! State variables used: (Wasp)
!
! nh3           71      1
! no3           72      2
! opo4          73      3
! phyto         74      4
! cbod          75      5
! do            76      6
! on            77      7
! op            78      8
! zoo           79      9
!
!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************

	subroutine handle_param_init

	implicit none

	character*80, save :: wparam_in = 'weutro.nml'
	character*80, save :: wparam_out = 'weutro_out.txt'

	logical bexist
	integer iu

	write(6,*) 'initializing parameters for weutro'

	call param_init		!read hard coded paramer list

	inquire(file=wparam_in,exist=bexist)
	if( bexist ) then
	  write(6,*) 'reading weutro parameter file: ',trim(wparam_in)
	  call read_weutro_param
	end if

	write(6,*) 'writing weutro parameters: ',trim(wparam_out)
	open(iu,file=wparam_out,status='unknown',form='formatted')
	call write_weutro_param(iu)
	close(iu)

	end

!********************************************************************

      subroutine param_init

      implicit none
      INCLUDE 'weutro.h'

	real iavpar

!! initialization of parameters

!       BYPASS OPTIONS FOR SYSTEMS 1-9. 1=BYPASS 0=SIMULATE
      

        PI = 3.14159

!!--------------------------------------------------------------
!!       please do not change anything in this subroutine below this point
!!       all changes to parameters should be done in a custom subroutine
!!       see param_venezia or param_taranto for example
!!--------------------------------------------------------------

!!       start parameter values

!--------------------------------------------------------------
!       general parameters
!--------------------------------------------------------------

!!      !if graztype=0  use wasp formulation
!!      !Set initial value of variable zoo=0
!!      !NUTLIM Default=0  ! GGU new value (was 1) ?
!!      !if CCHLX is variable in the segments, introduce each value else

        graztype=1   !type of grazing, simulate=1  forcing=0
        NUTLIM=1.    !nutrient limitation 0=minimum, 1=multiplicative.
        CCHL=50      !carbon to chlorophyll ratio mg carbon/mb chl a range 20-50
	CCHLX(1) = CCHL

!--------------------------------------------------------------
!       subroutine phyto
!--------------------------------------------------------------

!------ growth and respiration

!       K1RT=1.045      !wasp orig      
!       K1RC=0.125      !wasp orig
!       K1RT=1.068      !adapting the curve to dejak model
!       K1RC=0.096      !adapting the curve to dejak model
!       K1C=2.          !GGU new value (was 1.5) ?
!       KMPHYT=1.

        K1RT=1.068      !respiration temperature coefficient unitless
        K1RC=0.096      !respiration rate constant endogenous term day-1
        K1T=1.068       !growth rate temperature coeff WHEN SYSBY(4)=1
        K1C=2.          !maximum growth rate day-1 WHEN SYSBY(4)=1
        KMPHYT=0.       !half saturation constant for phytop. (eutroint)

!------ decomposition

        KPZDC=0.02      !decomp rate constant in sediment at 20 C, per day
        KPZDT=1.08      !decomposition temperaure coefficient in sediment day-1

!------ nutrients limitation

!       KMNG1=0.025     !for standard model application
!!                      !use a large KMNG1 (when KMNG1=0 PNH3G1=1.0)
!       KMNG1=0.05      !from dejak model Venice lagoon
!!      F(10,10,1)      !spatially variable dissolved fraction of inorganic P
!!      attenzione: nel benthic layer asume valori diversi (0.045-0.001)

        KMNG1=0.05      !N half saturn const for phyto growth microg N/L
        KMPG1=0.01      !P half saturn const for phyto growth microg P/L 
        FOPO4=0.9       !spatially variable dissolved fraction of inorganic P

!------ grazing - zooplankton as a forcing

!       K1G=0.08            
!       K1D=0.02        !wasp orig, 0.12 from dejak model Venice Lagoon
!!      if fraction is variable in segments, input each value, else:

        K1G=0.          !K1G=grazing rate on Phytop. per unit Zoop.L/mgC-day
        K1D=0.12        !unknown
        ZOO=0.7         !ZOO=time funct or const, zooplankton population mg C/L
        ZOOSG(1)=1.	!fraction of herbivorous zooplankto

!------ grazing - zooplankton as a variable

!       KGRZ=1.44       !max grazing rate day-1
!       KGRZ=2.         !prova 12sett
!       EFF=0.5         !zoo-phyto digestion efficiency

        KGRZ=2.         !grazing rate constant day-1
        KPHYZ=0.5       !half saturation constant for phytop in grazing
        EFF=0.7         !grazing efficiency for zooplankton growth
        KDZ=0.192       !death rate for zooplankton (with excrection)

!--------------------------------------------------------------
!       subroutine phosphorus (organop, inorganp, P for phosphorous)
!--------------------------------------------------------------

        PCRB=0.025   !P/C ratio mgP/mgC
        FOP=0.5      !fraction of dead/respired phyto recycled to organic P [1]
        K58T=1.08    !temp coeff: mineralization of dissolved org P unitless
        K58C=0.22    !rate for mineralization of dissolved organic P day-1 
        KOPDC=0.0004 !rate for mineralization  of organic P in sediment day-1
        KOPDT=1.08   !temp coeff for mineralization  of org P in sed unitless
        KPO4=1.0     !half saturation constant for P microgP/L

!--------------------------------------------------------------
!       subroutine nitrogen (organicn, ammonia, nitrate)
!--------------------------------------------------------------

!       K1320C=0.1      !0.09-0.13 day-1 wasp orig, dejak model 0.05 

        NCRB=0.115      !N/C ratio mg N/mg C
        FON=0.5         !fraction of organic Nitrogen
        K1320C=0.05     !rate for nitrification
        K1320T=1.08     !temp coefficient for nitrification unitless
        K140C=0.09      !rate for denitrification day-1
        K140T=1.045     !temp coefficient for denitrification unitless
        KNIT=2.         !half saturation constant for nitrification mgO2/L 
        KNO3=0.1        !half saturation constant for denitrificatiion mgO2/L 
        K1013C=0.075    !rate for mineralization dissolved organic nitrgen day-1
        K1013T=1.08     !temp coefficient for min. dissolved organic N unitless
        KONDC=0.0004    !rate for mineralization of N in sediment day-1     
        KONDT=1.08      !temp coeff for mineralization of N in sediment unitless

!--------------------------------------------------------------
!       subroutine CBOD
!--------------------------------------------------------------

        OCRB=2.66         !O/C ratio mg O2/mg C (=32./12.)
        KDC=0.18          !oxidation rate in water column day-1 (0.16-0.21)
        KDT=1.047         !temp coeff for oxidation in water column unitless
        KDSC=0.000        !oxidation rate in sediment at 20 C day-1
        KDST=1.08         !temp coeff for oxidation in sedim at 20 C unitless
        KBOD=0.5          !half saturation constant for oxidation mg N/L

!--------------------------------------------------------------
!       subroutine dissoxyg
!--------------------------------------------------------------

!!        REARSG(i)=0     !Segment specific reareation rate constant
!!                        !REARSG  used when rear is not calc from wind or hydro

        WIND=3.           !wind speed m/s
        AIRTMP=22.        !air temperature C
        WTYPE=3.          !water body size 1=small 2=medim 3=large
        XICECVR=1.        !ice cover, default: 1 for no ice

        K2=0 !reareation rate at 20 C day-1, 4.1-4.7 ,k2=0 use kawind or kahydra
        SOD1D(1)=0.0      !sediment oxygen demand for segment g/m2-day 0.2-4.0
        SODTA(1)=1.08     !temp coeff for sediment ox. demand unitless

!--------------------------------------------------------------
!       light limitation
!--------------------------------------------------------------

!       IS2=50000.      !optimum light intensity for Steele lux/h
!       IS2=1200000.    !optimum light intensity for Steele lux/day
!!      IS2= 20         !prova W/m2/timestep 300 sec!

        LGHTSW=3        !LIGHT SWITCH: 1=Di Toro 2=Smith  3=Steele
        IS2=1200000.    !optimum light intensity for Steele lux/day

!--------------------------------------------------------------
!       subroutine ditoro
!--------------------------------------------------------------

        FDAY=0.5        !fraction of day that is daylight unitless
        IS1=300.        !saturation light intensity for phyto used langleys/day
        IS1X(1) = IS1

!--------------------------------------------------------------
!       subroutine smith
!--------------------------------------------------------------

!       DTDAY=0.5       !?????? verificare, serve a modulare il seno: day
!       NEWDAY=1        !unknown???? if G.E. 1 calcolo smith

        PHIMX=720.      !mg C fixed per mole of light absorbed mg C/mole photon
        XKC=0.017       !chl a, phytop. self light attenuation m2/mg chl a
        ITOT=500.       !incident solar radiation ly/day
	iavpar = 0.9    !fraction of light used unitless
        iav = iavpar*itot/FDAY
        DTDAY=0.5       !unknown
        NEWDAY=1        !unknown

!--------------------------------------------------------------
!       light attenuation
!--------------------------------------------------------------

!! you can have up to 5 KE time functions here default is 1

        KE(1)=1.        !time var extin coeff funct unitless  (ditoro, smith)
        KEFN(1)=1       !IKE=KEFN=pointer designating time variable
        KESG(1)=1.5     !non algal light att. m-1 (SKE=KESG(ISEG)) range: 0.1-5

      end

!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************


!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************

        subroutine eutro0d(id,t,dt,vol,depth,vel,uws,stp,sal,qss &
     &				,c,loads)

! EUTRO 0-Dimensional

        implicit none

        integer nstate          !total number of state parameters
        parameter( nstate = 9 )

        integer id              !id of box
        real t                  !actual time [day]
        real dt                 !time step [day]
        real vol                !volume [m**3]
        real depth              !depth of box [m]
        real vel                !velocity [m/s]
        real uws                !wind velocity [m/s]
        real stp                !temperature [C]
        real sal                !salinity [psu] == [per mille]
        real qss                !solar radiation [W/m**2]
        real c(nstate)          !state variable [mg/L] == [g/m**3]
        real loads(nstate)      !loading for c [g/(m**3 day)]

        real cold(nstate)       !old state variable (for diagnostic purpose)
        real cds(nstate)        !source term (right hand side) [g/day]


        call source(id,t,vel,uws,stp,sal,vol,depth,qss,c,cds)
        call load0d(cds,loads,vol)
        call euler(nstate,dt,vol,c,cold,cds)
        call eutro_check(id,nstate,t,dt,vel,stp,sal,qss,vol,depth &
     &          ,c,cold,cds,loads)

        end

!********************************************************************

      subroutine load0d(cds,loads,vol)

! integrate loadings

      implicit none

        integer nstate          !total number of state parameters
        parameter( nstate = 9 )

      real cds(nstate)      !source term [g/day]
      real loads(nstate)      !loading for c [g/(m**3 day)]
      real vol            !volume of box [m**3]

      integer i
        
      do i=1,nstate
        cds(i) = cds(i) + vol * loads(i)
      end do

      end

!********************************************************************

      subroutine euler(nstate,dt,vol,c,cold,cds)

! new c is computed
! cold is returned which is just c before call

      implicit none

      integer nstate            !number of state variables
      real dt                  !time step [day]
      real vol            !volume [m**3]
      real c(nstate)            !state variable [mg/L] == [g/m**3]
      real cold(nstate)      !old state variable (return)
      real cds(nstate)      !source term [g/day]

      integer i
      real volold,volnew      !volume [m**3]
      real mass            !mass [g]
      real mder            !derivative [g/day]

      volold = vol
      volnew = vol

      do i=1,nstate
          cold(i) = c(i)
        mass = c(i) * volold
        mder = cds(i)
        c(i) = ( mass + dt * mder ) / volnew
        !if(c(i).lt.0.00001) c(i)=0.00001
      end do

      end
      
!********************************************************************

        subroutine eutro_check(id,nstate,t,dt,vel,stp,sal,qss &
     &          ,vol,depth &
     &          ,c,cold,cds,loads)

! checks result of integration for negatives

        implicit none

      integer id          !id of box
      integer nstate      !number of state variables
      real t              ![day]
      real dt             !time step [day]
      real vel            ![m/s]
      real stp            ![Celsius]
      real sal            ![g/L], e.g., 30.
      real qss            !solar rad from shyfem
      real vol            ![m**3]
      real depth          ![m]
      real c(nstate)      ![mg/L]
      real cold(nstate)   !old state variable (return)
      real cds(nstate)    !source term [g/day]
      real loads(nstate)  !loading for c [g/(m**3 day)]

        integer i
        logical berror

        berror = .false.

        if( vol .lt. 0. ) berror = .true.
        if( depth .lt. 0. ) berror = .true.

       do i=1,nstate
          if( c(i) .lt. 0. ) berror = .true.
        end do

        if( berror ) then
          write(6,*) '*** eutro_check: error in eutro0d'
          call write_debug_param(id,nstate,t,dt &
     &                ,vol,depth,vel,stp,sal,qss,c,cold,cds,loads)
          write(6,*) '--- eutro_check: end of error message'

          write(6,*) 'id=', id
          write(6,*) 't=', t
          write(6,*) 'dt=', dt
          write(6,*) 'vol=', vol
          write(6,*) 'depth=', depth
          write(6,*) 'vel=', vel
          write(6,*) 'stp=', stp
          write(6,*) 'sal=', sal
          write(6,*) 'qss=', qss
          write(6,*) 'c=', c
          write(6,*) 'cold=', cold
          write(6,*) 'cds=', cds
          write(6,*) 'loads=', loads
          write(6,*)
          stop 'error stop eutro_check: negative values'
          do i=1,nstate
            if( c(i) .lt. 0. ) c(i) = 0.00001
          end do
        end if

        end

!********************************************************************

      subroutine eutroini

! initializes eutro routines

      implicit none

      !call param_taranto
      !call param_venezia
      !call param_read        !GGU new for Michol
      !call param_print

      call handle_param_init
      call EUTROINT
      !call write_weutro_param(0)

      end

!********************************************************************

      subroutine param_print

! prints parameters for weutro

      implicit none

      INCLUDE 'weutro.h'

      write (6,*) 'graztype', graztype 
      write (6,*) 'CCHL' ,CCHL
      write (6,*) 'NUTLIM',NUTLIM
      write (6,*)  'growth and respiration'                                
      write (6,*) 'K1RT',K1RT
      write (6,*) 'K1RC',K1RC
      write (6,*) 'K1T',K1T
      write (6,*) 'K1C',K1C
      write (6,*) 'KMPHYT',KMPHYT
      write (6,*) 'decomposition'
      write (6,*) 'KPZDC',KPZDC
      write (6,*) 'KPZDT',KPZDT
      write (6,*) 'nutrients limitation' 
      write (6,*) 'KMNG1', KMNG1
      write (6,*) 'KMPG1', KMPG1
      write (6,*) 'FOPO4',FOPO4
      write (6,*) 'grazing'
      write (6,*) 'graztype',graztype
      write (6,*) 'wasp orig, K1G', K1G
      write (6,*) 'wasp orig, K1D', K1D
      write (6,*) 'wasp orig, ZOO', ZOO
      write (6,*)  'wasp orig, ZOOSG',ZOOSG(1) 
      write (6,*) 'KGRZ',KGRZ
      write (6,*) 'KPHYZ',KPHYZ
      write (6,*) 'EFF',EFF
      write (6,*) 'KDZ',KDZ
      write (6,*) 'subroutine organop, inorganp'
      write (6,*) 'PCRB',PCRB
      write (6,*) 'FOP',FOP
      write (6,*) 'K58T',K58T
      write (6,*) 'K58C',K58C
      write (6,*) 'KOPDC',KOPDC
      write (6,*) 'KOPDT',KOPDT
      write (6,*) 'KPO4',KPO4
      write (6,*) 'subroutine organicn, ammonia,nitrate' 
      write (6,*) 'NCRB',NCRB
      write (6,*) 'FON',FON
      write (6,*) 'K1320C',K1320C
      write (6,*) 'K1320T',K1320T
      write (6,*) 'K140C',K140C
      write (6,*) 'K140T',K140T
      write (6,*) 'KNIT',KNIT
      write (6,*) 'KNO3',KNO3
      write (6,*) 'K1013C',K1013C
      write (6,*) 'K1013T',K1013T
      write (6,*) 'KONDC',KONDC
      write (6,*) 'KONDT',KONDT
      write (6,*) 'subroutine CBOD'
      write (6,*) 'OCRB',OCRB
      write (6,*) 'KDC',KDC
      write (6,*) 'KDT',KDT 
      write (6,*) 'KDSC',KDSC
      write (6,*) 'KDST',KDST
      write (6,*) 'KBOD',KBOD
      write (6,*) ' subroutine dissoxyg'
      write (6,*) 'SOD1D(1)',SOD1D(1)
      write (6,*) 'SODTA(1)',SODTA(1)
      write (6,*) ' light limitation'
      write (6,*) 'LGHTSW',LGHTSW
      write (6,*) 'IS2',IS2
      write (6,*) 'PHIMX',PHIMX 
      write (6,*) 'XKC',XKC
      write (6,*) 'ITOT',ITOT
      write (6,*) 'iav',iav
      write (6,*) 'DTDAY',DTDAY

      end

!********************************************************************

      subroutine source(id,t,vel0,uws0,stp0,sal0,vol0,depth0,qss0,c,cds)

! vel      velocity
! stp      temperature
      
      implicit none

        integer nstate          !total number of state parameters
        parameter( nstate = 9 )

      integer id            !id of box
      real t                  ![day]
      real vel0            ![m/s]
      real uws0		   !wind velocity [m/s]
      real stp0            ![Celsius]
      real sal0            ![g/L], e.g., 30.
      real qss0            !solar radiation from shyfem
      real vol0            ![m**3]
      real depth0            ![m]
      real c(nstate)            ![mg/L]            !FIXME
      real cds(nstate)      ![g/day]      !FIXME

! if rlux is > 0, at return it will be replaced with the attenuation
! factor at the bottom of the box -> therefore it can be used
! with the next call in 3D models looping from surface to bottom
! for 2D models the value of rlux should be 1 (value replaced at return)
! or <= 0 which leaves the value untouched
! for the surface box the value of rlux should be 1

        INCLUDE 'weutro.h'
      
      integer i,n

        idbox = id

      iseg = 1
      sedseg = .false.
      ito = 0
      
      daytime = t

! deal with light climate


!      call rintens(daytime,itot,fday,iinst)      !compute instant light

         NH3 = C (1)
         NO3 = C (2)
         OPO4 = C (3)
         PHYT = C (4)
         CBOD = C (5)
         DO = C (6)
         ON = C (7)
         OP = C (8)
       ZOO = C (9)

      stp = stp0
      stp20 = stp-20.
      sal = sal0
      vel = vel0
      wind = uws0
      h = depth0
      vol = vol0

        SA = vol0/depth0

      iinst=qss0
!      write(6,*)iinst

!                        Compute derivatives
!              For PHYT, OP, OPO4, ON, NH3, NO3, CBOD, DO

      do i=1,nstate
        cd(i,iseg) = 0.
      end do

       CALL ZOOP
         CALL PHYTO
         CALL ORGANOP
         CALL INORGANP
         CALL ORGANICN
         CALL AMMONIA
         CALL NITRATE
         CALL CBODSV
         CALL DISSOXYG

         if( botseg ) CALL BENTFLUX

! adesso sono definiti i cd

      do i=1,nstate
        cds(i) = cd(i,iseg)
      end do

! save light climate for next call


      end

!********************************************************************
!********************************************************************
!********************************************************************

      subroutine wmeteo(tempair,windspeed)

! sets meteo parameters
!
! tempair is the air temperature in degrees C
! windspeed is the wind speed in m/s

      implicit none

      include 'weutro.h'

      real tempair
      real windspeed

      airtmp = tempair
      wind   = windspeed

      end

!********************************************************************

      subroutine wturbid(turbid)

! sets turbidity parameter (units [1/m])

      implicit none

      include 'weutro.h'

      real turbid

      !kefn(1) = 1
      !ke(1)   = 1.

      kesg(1) = turbid

      end

!********************************************************************

      subroutine wlight(fracday,itotal)

! sets light parameters
!
! fracday      the fraction of the day that is light [0-1]
! itotal      average incident light intensity during whole day [ly/day]

      implicit none

      include 'weutro.h'

      real fracday
      real itotal

      fday = fracday
      itot = itotal
      itotmp = itotal            !not used

      iav=0.9 *itot/FDAY
!       IAV = 0.9*ITOTmp*RLGHTS (ISEG, 2)/FDAY !is multiplied in ditoro

      end

!********************************************************************

      subroutine icecover(ice)

! sets ice cover parameter
! ice=1: ice covers all the surface; ice=0 no ice cover      

      implicit none

      include 'weutro.h'

      real ice

      xicecvr = 1. - ice

      end

!********************************************************************

      subroutine grazing(zoopl)

! sets zooplankton state variable
!
! set in zoo the C/l value of grazing zooplankton
! may be used to implement measured zoo

      implicit none

      include 'weutro.h'

      real zoopl

      zoo = zoopl

      end

!********************************************************************

      subroutine sedflux(flnh4,flpo4)

! sets sediment fluxes
! fluxes are in [mg/(m**2 day)]

      implicit none

      include 'weutro.h'

      real flnh4, flpo4

      fnh4(1) = flnh4
      fpo4(1) = flpo4

!      write(6,*) 'sedflux : ',flnh4,flpo4

      end

!********************************************************************
!********************************************************************
!********************************************************************

        subroutine settopseg(bstate)

! sets segment as surface or not

        implicit none
        include 'weutro.h'
        logical bstate

        topseg = bstate

        end

!********************************************************************

        subroutine setbotseg(bstate)

! sets segment as bottom or not

        implicit none
        include 'weutro.h'
        logical bstate

        botseg = bstate

        end

!********************************************************************

        subroutine setsedseg(bstate)

! sets segment as sediment or not

        implicit none
        include 'weutro.h'
        logical bstate

        sedseg = bstate

        end

!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************

      SUBROUTINE EUTROINT

!     Initialize Values
!
!     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:53.

      implicit none

!DDD      changed PHIMAX ---> PHIMX

      INCLUDE 'weutro.h'
      real R0MIN
      parameter( R0MIN = 1.e-15 )


      real xarg,phimax
      integer i,j,initb,mxdmp,mxseg,noseg

        IDBOX = 0

      topseg = .true.
      botseg = .true.
        wdebug = .false.

!      DUMMY = 0.
!      CHLA2 = 0.
      GP1 = 0.
      DP1 = 0.
!      GP2 = 0.
!      DP2 = 0.
!      GZ1 = 0.
!      DZ1 = 0.
!      GZ2 = 0.
!      DZ2 = 0.
!      CFOREA = 1.0
!
      noseg=segmax

      DO 1000 J = 0,noseg 
         RLGHTS (J, 1) = 0.0
         RLGHTS (J, 2) = 1.0
 1000 CONTINUE
!
      INITB = 1
      MXDMP = 4
!
!             Check to see if Michalis Menton constants are Zero
!             and readjust values to prevent floating zero divide
!
!RBA--Date: Tuesday, 1 June 1993.  Time: 09:01:20.
      XARG = ABS(NCRB)
      IF (XARG .LT. R0MIN) NCRB = 0.25
      XARG = ABS(PCRB)
      IF (XARG .LT. R0MIN) PCRB = 0.025
!      XARG = ABS(LGHTSW)
!      IF (XARG .LT. R0MIN) LGHTSW = 1.0
      if( LGHTSW .lt. 0 ) LGHTSW = 0
      if( LGHTSW .gt. 3 ) LGHTSW = 3
      
!CSC
      XARG = ABS(KBOD)
      IF (XARG .LT. R0MIN) KBOD = 1.00E-20
!CSC
      XARG = ABS(KNO3)
      IF (XARG .LT. R0MIN) KNO3 = 1.00E-20
!CSC
      XARG = ABS(KPO4)
      IF (XARG .LT. R0MIN) KPO4 = 1.00E-20
!CSC
      XARG = ABS(KNIT)
      IF (XARG .LT. R0MIN) KNIT = 1.00E-20
!CSC
      XARG = ABS(KMNG1)
      IF (XARG .LT. R0MIN) KMNG1 = 1.00E-20
!CSC
      XARG = ABS(KMPG1)
      IF (XARG .LT. R0MIN) KMPG1 = 1.00E-20
!CSC
      XARG = ABS(KMPHYT)
      IF (XARG .LT. R0MIN) KMPHYT = 1.00E-20
!
!CSC
      XARG = ABS(OCRB)
      IF (XARG .LT. R0MIN) OCRB = 32./12.
!CSC
      XARG = ABS(IS1)
      IF (XARG .LT. R0MIN) IS1 = 300.
!CSC
      XARG = ABS(CCHL)
      IF (XARG .LT. R0MIN) CCHL = 30.
!CSC
      XARG = ABS(FON)
      IF (XARG .LT. R0MIN) FON = 1.0
!CSC
      XARG = ABS(FOP)
      IF (XARG .LT. R0MIN) FOP = 1.0
!CSC
      XARG = ABS(PHIMX)
      IF (XARG .LT. R0MIN) PHIMX = 720.
!CSC
      XARG = ABS(XKC)
      IF (XARG .LT. R0MIN) XKC = 0.017
!
!  Check for Zero Temperature Correction Factors and readjust to 1.0
!
!CSC
      XARG = ABS(K1320T)
      IF (XARG .LT. R0MIN) K1320T = 1.0
!CSC
      XARG = ABS(K140T)
      IF (XARG .LT. R0MIN) K140T = 1.0
!CSC
      XARG = ABS(K1T)
      IF (XARG .LT. R0MIN) K1T = 1.0
!CSC
      XARG = ABS(K1RT)
      IF (XARG .LT. R0MIN) K1RT = 1.0
!CSC
      XARG = ABS(KDT)
      IF (XARG .LT. R0MIN) KDT = 1.0
!CSC
      XARG = ABS(K1013T)
      IF (XARG .LT. R0MIN) K1013T = 1.0
!CSC
      XARG = ABS(KONDT)
      IF (XARG .LT. R0MIN) KONDT = 1.0
!CSC
      XARG = ABS(K58T)
      IF (XARG .LT. R0MIN) K58T = 1.0
!CSC
      XARG = ABS(KOPDT)
      IF (XARG .LT. R0MIN) KOPDT = 1.0
!CSC
      XARG = ABS(KPZDT)
      IF (XARG .LT. R0MIN) KPZDT = 1.0
!CSC
      XARG = ABS(KDST)
      IF (XARG .LT. R0MIN) KDST = 1.0

! sediment fluxes

      do i=1,segmax
          FNH4 (i) = 0.
          FPO4 (i) = 0.
      end do

!      write(6,*) segmax,FNH4(1),FPO4(1)

      RETURN
      END

!********************************************************************

      SUBROUTINE AMMONIA

!     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:55.

      implicit none

      real SR13ON,SR13P,SK13P1

      INCLUDE 'weutro.h'
      include 'donata.h'      !GGU (and others)
!
!       *-*-*-*-*  SYSTEM 1 - AMMONIA (NH3-N)  *-*-*-*-*
!
!                        Sources
!               Mineralization of organic nitrogen
!
      SR13ON = SK1013
      denit = 0.
      denit = SK14D
!
!                  Phytoplankton Death
!
      SR13P = NCRB*DPP*(1.0 - FON)
!
!                        Sinks
!                    Algal Uptake
!
      SK13P1 = PNH3G1*NCRB*GP1*PHYT
!
!                      Nitrification
!
      IF (DO .GT. 1.0E-10) THEN
         SK1314 = (K1320C*K1320T**STP20)*NH3*DO/(KNIT + DO)
      ELSE
         SK1314 = 0.0
      END IF
      IF (STP .LT. 7.) SK1314 = 0.0
       ddin1=sk1013+SR13P-SK13P1        !GGU new
!      write(6,*)ddin1

!
!                   Formulate Derivative
!
!      write(6,*) 'ammonia debug :'
!      write(6,*) SR13P,SR13ON,SK13P1,SK1314,VOL,CD (1, ISEG)
!      write(6,*) PNH3G1,NCRB,GP1,PHYT

      if (wdebug) then
      write(99,*)'ammondebug:',SR13P,denit,SR13ON,SK13P1,SK1314
      end if

!      CD (1, ISEG) = (SR13P +denit+ SR13ON - SK13P1 - SK1314)*VOL      !just for tests
      CD (1, ISEG) = (SR13P + SR13ON - SK13P1 - SK1314)*VOL     !ORIG

!      write(6,*) CD (1, ISEG),iseg
!
      RETURN
      END

!********************************************************************

      SUBROUTINE NITRATE

!     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:56.

      implicit none
      INCLUDE 'weutro.h'
      include 'donata.h'
      real SR1413,SK14P1
!
!    *-*-*-*-* SYSTEM 2 - Nitrate (+Nitrite)  (NO3-N+NO2-N)  *-*-*-*-*
!
!                         Sources
!                      Nitrification
!
      SR1413 = SK1314
!
!                          Sinks
!                      Algal Uptake
!
      SK14P1 = (1. - PNH3G1)*NCRB*GP1*PHYT
!
!                     Denitrification
!
      SK14D = (K140C*K140T**STP20)*NO3
      IF (DO .GT. 0) SK14D = SK14D*KNO3/(KNO3 + DO)
!
! for no denitrification set SK14D to 0.

!      IF (SK14D .LT. 1.00E-24) SK14D = 1.00E-24

      ddin2=-SK14P1-SK14D       !GGU new
!      write (6,*) ddin2


!
!                   Formulate Derivative
!
      if (wdebug) then
        write(99,*)'NO2debug:',SR1413,SK14P1,SK14D,DO,NO3
      end if


      CD (2, ISEG) = (SR1413 - SK14P1 - SK14D)*VOL
!
      RETURN
      END

!********************************************************************

      SUBROUTINE INORGANP

!     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:56.

      implicit none
      INCLUDE 'weutro.h'
      include 'donata.h'
      real SR8P,SR8OP,SK8P
!
!   *-*-*-*-*  SYSTEM 3 - Dissolved Inorganic Phosphorous  *-*-*-*-*
!
!                            Sources
!               Mineralization of Organic Phosphorous
!
      SR8OP = SK58
!
!                         Phytoplankton Death
!
      SR8P = PCRB*DPP*(1. - FOP)
!
!                            Sinks
!                        Algal uptake
!
      SK8P = PCRB*GPP

        prod = SR8OP+SR8P       !GGU new
      cons= SK8P              !GGU new
!      write(6,*) 'prod phosph',prod
!
!                   Formulate Derivative
!
      if (wdebug) then
      write(99,*)'inPdebug:',SR8P,SR8OP,SK8P
      end if

      CD (3, ISEG) = (SR8P + SR8OP - SK8P)*VOL
!
      RETURN
      END

!********************************************************************

      SUBROUTINE PHYTO

!     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:57.

      implicit none

      INCLUDE 'weutro.h'
      real R0MIN
      parameter( R0MIN = 1.e-15 )
      REAL XARG,XDIFF,GIT1,CN,DOPO4
      real ttr, ttr1
!
!        *-*-*-*-*-*  System 4 - Phytoplankton  *-*-*-*-*-*
!
!         open (9, file = 'outgp.dat')      !ggu


!      write(6,*) 'phyto, sedseg,iseg : ',SEDSEG,iseg
      IF (SEDSEG) THEN

!      write(6,*) 'phyto,  ctrl: ',K1C,PHYT,K1T,STP20
!      write(6,*) 'phyto, ctrl: ',GITMX1,KEFN(1)
            GP1 = K1C
            GPP = GP1*PHYT
            PNH3G1 = 0.

         ELSE
!      write(6,*) 'phyto, ctrl : ',K1C,PHYT,K1T,STP20
            GITMX1 = K1C*K1T**STP20
            GITMAX = GITMX1
            IKE = KEFN (ISEG)
      
!         write(6,*) 'phyto, ctrl1 : ',GITMAX
!
!               Compute growth rate reduction due to light conditions
!               using either Dick Smith's or Di Toro's formulation
!
            RLIGHT = 1.
            if( LGHTSW .eq. 1 ) call ditoro
            if( LGHTSW .eq. 2 ) call smith
            if( LGHTSW .eq. 3 ) call steele

             GIT1 = RLIGHT*GITMAX
!            write(6,*) 'phyto, ctrl2 : ',GITMAX,GIT1

!                   Compute ammonia preference
            PNH3G1 = 0.0
!
!
            IF (NH3 .GE. 1.0E-5) &
     &         PNH3G1 = NH3*NO3/((KMNG1 + NH3)*(KMNG1 + NO3)) &
     &         + NH3*KMNG1/((NH3 + NO3)*(KMNG1 + NO3))
!
!              Compute Michaelis Limitations
!
            CN = NH3 + NO3
!      write(6,*)KMNG1,PNH3G1,NH3,NO3,CN
            XEMP1 = CN/(KMNG1 + CN)

               DOPO4 = OPO4*FOPO4
            XEMP2 = DOPO4/(KMPG1 + DOPO4)
      
!
!       Compute Growth Rate Reduction due to Nutrient Limitation
!
!CSC
            XARG = ABS(NUTLIM)
            IF (XARG .LT. R0MIN) RNUTR = AMIN1 (XEMP1, XEMP2)

!      write(6,*) '-------------- ggu --------------'
!      write(6,*) NUTLIM,XARG,R0MIN
!      write(6,*) XEMP1,XEMP2,RNUTR
!      write(6,*) '-------------- ggu --------------'
!CSC
            IF (NUTLIM .EQ. 1.) RNUTR = XEMP1*XEMP2
!            write(6,*) RNUTR
            GP1 = RNUTR*GIT1
            GPP = GP1*PHYT
      
!
!       ********************************************
!                     Respiration Rate
!       ********************************************
!
         RESP = K1RC*K1RT**STP20
!      
!
!       ALGAL RESPIRATION + DEATH + GRAZING
!
!         DP1 = RESP + K1D + K1G*ZOO*ZOOSG(ISEG)      !old

         DP1 = RESP + K1D             !FIXMED 

         DPP = DP1*PHYT      !BUGFIX1   LAA
!         DPP = DPP + DP1*PHYT       
!         RESP = RESP*PHYT      !BUGFIX1

!      write(6,*) DP1,RESP,K1D,K1G,ZOO,ZOOSG(iseg)


!         IF (PHTY .GT. 1.0E-6)THEN
      XEMPRC=1.
!doni            XEMPRC = PHYT/(KMPHYT + PHYT)
!         ELSE
!            XEMPRC = 1.0E-6/(KMPHYT + 1.0E-6)
!         ENDIF

!
      END IF      !sedsed
!
!      write(6,*) 'phyto debug ggu'
      if (wdebug) then
      write(99,*)'phytodebug:',GP1,DPP,GRZ
      end if

!        write(9,'(3(f8.4,2x))') RLIGHT,DP1, RESP       !ggu

!      CD (4, ISEG) = (GP1 - DP1-GRZ)*PHYT*VOL      !BUGFIX1
      CD (4, ISEG) = (GPP - DPP - GRZ)*VOL
      RETURN
      END

!********************************************************************

      subroutine zoop

! zooplankton

      implicit none
      include 'weutro.h'
      real gra

!      *-*-*-*- system 9 zooplankton *-*-*-
!
!      Sources: zooplankton growth

! zoosk      Sink term: zooplankton death -->source term for organop organicn cbodsv
!          grazing inefficiency -> source term for  organop organicn cbodsv

      if( graztype.eq.1 ) then
        DPP = 0.
          GRZ=KGRZ*ZOO*PHYT/(PHYT+KPHYZ)
        zoosk = (1-EFF)*GRZ + KDZ*ZOO
        CD (9, ISEG) = (EFF*GRZ - KDZ*ZOO)*VOL
      else      !BUGFIX1
        DPP = K1G*ZOO*PHYT*ZOOSG(ISEG)      !original formulation
        GRZ = 0.
        zoosk = 0. 
        CD (9, ISEG) = 0.
      end if
      !zood=KDZ*ZOO           !GGU ???

      if (wdebug) then
      write(99,*)'zoodebug:',GRZ
      end if

      return
      end

!********************************************************************

      SUBROUTINE CBODSV

!     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:57.

      implicit none
      INCLUDE 'weutro.h'
      real SR18P,SK18D,GRC
!
!                *-*-*-*-* SYSTEM 5 - CBOD  *-*-*-*-*
!                             Sources
!
      IF ( .NOT. SEDSEG) GO TO 1000
      SR18P = OCRB*DPP
      GO TO 1010
!
!                        Phytoplankton 'DEATH'
!
 1000 CONTINUE
      DEATH = K1D*PHYT
      SR18P = OCRB*DEATH
      GRC= ZOOSK * OCRB
!
!                              Sinks
!                            Oxidation
!
 1010 CONTINUE
      IF ( .NOT. SEDSEG) THEN
!         IF (DO .GT. 1.0E-15) THEN
         IF (DO .GT. 1.0E-15 .AND. CBOD .GT. 1.0E-5) THEN
            SK180 = (KDC*KDT**STP20)*((CBOD*DO)/(KBOD + DO))
         ELSE
            SK180 = 0.0
         END IF
      ELSE
         IF (CBOD .GT. 1.0E-15) THEN
            SK180 = (KDSC*KDST**STP20)*CBOD
         ELSE
            SK180 = 0.0
         END IF
      END IF
!
!                         Denitrification
!
! ggu & mcg: next formula is indepenedent of BOD -> insert dependence

       SK18D = (5./4.)*(32./14.)*SK14D
!
!                      Formulate Derivative
      if (wdebug) then
      write(99,*)'cboddebug:',sr18p,grc,sk180,sk18d
      end if
!
      CD (5, ISEG) = (SR18P+GRC - SK180 - SK18D)*VOL
!
      RETURN
      END

!********************************************************************

      SUBROUTINE DISSOXYG

!     Last Revised:  Date: Monday, 26 August 1991.  Time: 10:37:53.

      implicit none
      INCLUDE 'weutro.h'
      real R0MIN
      parameter( R0MIN = 1.e-15 )
      REAL XARG1, XARG2
      REAL K2WIND, K2HYDRA
      real SR190,SR19PA,SR19PB,SK1913,SK1918,SK19S
      real TK,RLNCS,CL
!
!              *-*-*-*-*  SYSTEM 6 - OXYGEN  *-*-*-*-*
!
!                          Sources
!
      K20 = 0.0
      IF ( .NOT. SEDSEG) GO TO 1000
      SR190 = 0.0
      SR19PA = 0.0
      SR19PB = 0.0
      SK19P = 0.0
      GO TO 1010
!
!                         Reaeration
!
 1000 CONTINUE
      IF ( topseg .AND. XICECVR .GT. 0.0) THEN

! if surface element and not completly covered with ice

! elimino rear e rearsg perch\E9 non ci interessano, sono una time function 
! e una segment specific reareation rate, ma noi usiamo la reareation from
!  kawind o kahydra o imponendo un valore della costante di reareazione

         IF (K2 .EQ. 0.0) THEN
            CALL KAWIND (WIND, STP, AIRTMP, H, WTYPE,  K2WIND) 
            CALL KAHYDRA (K2HYDRA)
            IF (K2WIND .GT. K2HYDRA) THEN
               KA = K2WIND * XICECVR
            ELSE
               KA = K2HYDRA * XICECVR
            END IF
         ELSE
            IF (K2 .GT. 0)THEN
                KA = ((K2*1.028**STP20)* XICECVR)
            ELSE
            KA = 0.
            ENDIF
         ENDIF
      ELSE
         KA = 0.0
      END IF

!       Calculate oxygen saturation level for current water
!       temperature; DOSAT is expressed as mg oxygen per liter

      CL = SAL/1.80655
      TK = STP + 273.
      RLNCS = - 139.34411 + (1.575701E05/TK) - (6.642308E07/TK**2) + &
     &   (1.243800E10/TK**3) - (8.621949E11/TK**4) - &
     &   (CL*(3.1929E-02 - (19.428/TK) + (3.8673E03/TK**2)))
!
      CS = EXP (RLNCS)
!
      SR190 = KA*(CS - DO)

!      write(88,*) KA,CS,SR190
!
!                 Evolution by phytoplankton
!          growth of phytoplankton using CO2 and NH3
!
      SR19PA = PNH3G1*GP1*PHYT*32./12.
!
!       Growth of phytoplankton using CO2 and NO3 (2NO3 = 2NH3 + 302)
!
      SR19PB = (1. - PNH3G1)*GP1*PHYT*32.*(1./12. + 1.5*NCRB/14.)

!
      SR19P = SR19PA + SR19PB
!
!      NOTE: SR19P = GPP*(32/12 + (1.5*NCRB/14)*(1-PNH3G1))
!
!                             Sinks
!                       Algal Respiration

!      SK19P = OCRB*RESP                  !BUGFIX1
      SK19P = OCRB*RESP*PHYT
!
!        Nitrification (NH3-N + 2O2 = NO3-N + H2O + H)
!
 1010 CONTINUE
      SK1913 = 64./14.*SK1314
!
!                      Oxidation of CBOD
!
      SK1918 = SK180
!
!             Sediment Oxygen Demand (1-D Networks)
!
      SK19S = SOD1D (ISEG)*SODTA(ISEG)**stp20/H
!
!======================================================================
!                     Formulate Derivative
!
      if (wdebug) then
      write(99,*)'doxdebug:',SR190,SR19PA,SR19PB,SK19P &
     &                        ,SK1913,SK1918,SK19S
      end if

      CD (6, ISEG) = (SR190 + SR19PA + SR19PB - SK19P &
     &   - SK1913 - SK1918 - SK19S)*VOL

!        write(54,'(8(f8.4,2x))') SR190,SR19PA,SR19PB,SK19P,SK1913,SK1918
!     & ,SK19S 
!        write(54,'(3(f8.4,2x))') SR19PA,SR19PB,GP1      !ggu
       
      RETURN
      END

!********************************************************************

      SUBROUTINE ORGANICN

!     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:59.

      implicit none
      INCLUDE 'weutro.h'
      real SR10P,GRN
!
!
!             *-*-*-*-*  SYSTEM 7 Organic Nitrogen      *-*-*-*-*
!
!                          Sources
!             Phytoplankton Respiration and 'DEATH'
!            GRAZING and zooplankton death
!                  GRN=  ZOO --> ORGANICN

!
      SR10P = NCRB*DPP*FON
      GRN   = ZOOSK*NCRB      !BUGFIX1
!
!                         Sinks
!         Mineralization of Dissolved Organic Nitrogen
!
      IF ( .NOT. SEDSEG) SK1013 = (K1013C*K1013T**STP20)*ON*XEMPRC
      IF (SEDSEG) SK1013 = (KONDC*KONDT**STP20)*ON
!
!                   Formulate Derivative
!
      if (wdebug) then
      write(99,*)'orgNdebug:',SR10P,GRN,SK1013
      end if
      CD (7, ISEG) = (SR10P+ GRN - SK1013)*VOL
!
      RETURN
      END

!********************************************************************

      SUBROUTINE ORGANOP

!     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:59.

      implicit none
      real SR5P,GRP

      INCLUDE 'weutro.h'
!
!        *-*-*-*-*-*  SYSTEM 8 Organic Phosphorous  *-*-*-*-*-*
!
!                            Sources
!                 Phytoplankton respiration and 'death'
!               GRAZING and zooplankton death
!                  GRP=  ZOO --> ORGANOP
!
      SR5P = PCRB*DPP*FOP
      GRP  = ZOOSK * PCRB      !BUGFIX1
!
!                             Sinks
!
!           Mineralization of dissolved organic phosphorus and
!           Phytoplankton respiration and 'death'
!
      IF ( .NOT. SEDSEG) SK58 = (K58C*K58T**STP20)*OP*XEMPRC
      IF (SEDSEG) SK58 = (KOPDC*KOPDT**STP20)*OP
!
!                        Formulate Derivative
!
      if (wdebug) then
      write(99,*)'orgPdebug:',SR5P,GRP,SK58
      end if
      CD (8, ISEG) = (SR5P + GRP - SK58)*VOL
!
      RETURN
      END

!********************************************************************
!
! IAVBOT                bottom light                            not used
! IAVBOTX(ISEG)         bottom light (variable in seg)          not used
! daytime               = t
! KESG(ISEG)            given turbidity [1/m]
! KE(IKE)               time function for turbidity
! SKE                   turbidity (given and phyto) [1/m]
!
! RLIGHT                limiting factor
! RLGHTS(ISEG, 1)       = RLIGHT
! RLGHTS(ISEG, 2)       fraction of surface light in segment (0-1)
!                          is 1 in surface layer
! ITO                   pointer to lower segment
! IS1,IS1X(ISEG)        saturation light intensity
! IS2                   optimum light intensity for Steele
! FDAY                  fraction of day that is daylight (unitless)
! IAV                   average incident light intensity during daylight hours
!                          (iav=0.9*itot/fday)
! ITOT                  total incident light intensity during one day [ly/day]
! ITOTMP                ...                                     not used
!
! smith and steele use ITOT
! ditoro uses IAV
!
! summary: we need IAV, ITOT, IINST, FDAY, IS1, IS2, KESG
!
! primary:   ITOT, FDAY
! derived:   IAV, IINST
! constants: IS1, IS2, KESG
!
! set primary, compute derived, use constants
!
! moreover: 
!      RLGHTS(ISEG,2) is attenuation coefficient (enter) [0-1]
!      RLGHTS(ITO,2) is attenuation coefficient at bottom (return) [0-1]
!
!********************************************************************

      SUBROUTINE DITORO

!*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
!        Di Toro et al Light Formulation
!*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
!
! needs IAV, IS1 and FDAY (plus KESG)

      implicit none
      INCLUDE 'weutro.h'
      real CCHL1,TCHLA,KESHD,SKE
      real temp1, temp2, temp3

      CCHL1 = CCHL
      TCHLA = PHYT/CCHL1
      KESHD = (0.0088*1000.*TCHLA + 0.054*(1000.*TCHLA)**0.6667)
      SKE = KESG (ISEG)
      IF (IKE .GT. 0 .AND. IKE .LE. 5) SKE = SKE*KE (IKE)
      SKE = SKE + KESHD
      TEMP1 = SKE*H

!         Get average solar radiation during daylight hours

      TEMP2 = IAV/IS1
      TEMP2 = RLGHTS(ISEG,2)*IAV/IS1      !consider light attenuation - $LGA
      TEMP3 = EXP ( - TEMP1)
      IAVBOT=IAV * TEMP3
      RLGHTS (ITO, 2) = RLGHTS (ISEG, 2)*TEMP3
      RLIGHT = 2.718*FDAY/TEMP1*(EXP ( - TEMP2* &
     &   TEMP3) - EXP ( - TEMP2))
      RLGHTS (ISEG, 1) = RLIGHT

!      write(6,*) 'ditoro debug ggu'
!      write(6,*) H,SKE,TEMP1,TEMP2,TEMP3,IAVBOT
!      write(6,*) FDAY,RLGHTS (ISEG, 1),RLIGHT
!      write(6,*) ISEG,ITO,IS1,IAV, 2.718*FDAY,temp1
!      write(6,*) H
!       write(10,'(6(f8.5,2x))')KESHD,SKE,PHYT,RLIGHT

      RETURN
      END

!********************************************************************

      SUBROUTINE SMITH

!*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
!          Dick Smith variable carbon/chlorophyll ratio
!*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
!
! needs ITOT and FDAY (plus KESG)
! IS1 is not needed (optimal light value is computed)
! RLGHTS (ISEG, 2) is always 1 as implemented now

      implicit none
      INCLUDE 'weutro.h'

      real I0,IMAX,IAVSG,SUM,KESHD,SKE
!      real IAVBOTX(segmax)
      real CCHL1,TCHLA
       real temp1, temp2, temp3
      integer I
!
!            IAV = IAV/1440.
!
!               (Average solar radiation during daylight hours)
!
!/cm2 day*[0.43 visible/total*10,000 cm2/m2*E/52,000]=E visible/m2-day
!             0.083  is FU (E/M2 - Ly or E/10 Kcal)
!
!RBA--Date: Tuesday, 1 June 1993.  Time: 09:05:36.
      CCHL1 = CCHLX(ISEG)
      IS1 = IS1X(ISEG)
      IAVBOT=IAVBOTX(ISEG)
!      CCHL1 = 0.3*0.083*PHIMX*XKC*IAV/(GITMX1*2.718)
!      IF (CCHL1 .LT. 20.) CCHL1 = 20.
      TCHLA = PHYT/CCHL1
      RLIGHT = RLGHTS (ISEG, 1)
!
      IF (NEWDAY .GE. 1) THEN
!
!           Dick Smith formulation integrated every day
!
         KESHD = XKC*1000.*TCHLA      !FIXME
         SKE = KESG (ISEG)
         IF (IKE .GT. 0 .AND. IKE .LE. 5) SKE = SKE*KE (IKE)
         SKE = SKE + KESHD
         TEMP1 = SKE*H
         TEMP2 = 0.083*PHIMX*XKC/(GITMAX*CCHL1*2.718)
         TEMP3 = EXP ( - TEMP1)
         RLGHTS (ITO, 2) = RLGHTS (ISEG, 2)*TEMP3
         IMAX = PI*ITOT*RLGHTS (ISEG, 2)/(2.*FDAY)
         SUM = 0.0
         DO 1000 I = 1, 25
            DTDAY = FLOAT (I - 1)/24.
            IF (DTDAY .GT. FDAY) GO TO 1010
            I0 = IMAX*SIN (PI*DTDAY/FDAY)
            SUM = SUM + 2.7183/TEMP1* &
     &         (EXP ( - TEMP2*I0*TEMP3) - EXP ( - TEMP2*I0))
 1000    CONTINUE
 1010    CONTINUE
         RLIGHT = SUM/24.
         RLGHTS (ISEG, 1) = RLIGHT
!RBA--Date: Tuesday, 1 June 1993.  Time: 09:06:22.
!        Adapt carbon to chlorophyll ratio:
         IAVSG=IAV*(1.0-TEMP3)/TEMP1
         CCHLX(ISEG)=0.3*0.083*PHIMX*XKC*IAVSG/(GITMX1*2.718)
         IF(CCHLX(ISEG).LT.20.0) CCHLX(ISEG)=20.0
         IS1X(ISEG)=1/TEMP2
         IAVBOTX(ISEG)=IAV*TEMP3
      END IF
      RETURN
      END

!***************************************************************

      subroutine steele

! computes light limitation with steele formula
!
! uses instantaneous light intensity IINST, and not IAV or ITOT
!
! needs IINST and IS2 (plus KESG)

      implicit none

      include 'weutro.h'

      real expon
      parameter ( expon = 2.718 )

      real TCHLA,KESHD,SKE
        real temp1,temp2,temp3

        TCHLA = PHYT/CCHL
        KESHD = (0.0088*1000.*TCHLA + 0.054*(1000.*TCHLA)**0.6667)

        SKE = KESG (ISEG)
        IF (IKE .GT. 0 .AND. IKE .LE. 5) SKE = SKE*KE (IKE)
        SKE = SKE + KESHD

!       write(66,*) daytime,itot,fday,iinst,ske

      temp1 = ske * h
      temp1 = ske                                    !$WRONG_FORMULA
      temp2 = iinst / is2
      temp2 = RLGHTS(ISEG,2) * iinst / is2      !light attenuation $LGA
      temp3 = exp ( - temp1 )

        RLGHTS (ITO, 2) = RLGHTS (ISEG, 2)*TEMP3

      RLIGHT = (expon/temp1) * ( exp(-temp2*temp3) - exp(-temp2) )
        RLIGHT = temp2*temp3*exp(1-temp2*temp3)         !2D_FORMULA

        RLGHTS (ISEG, 1) = RLIGHT
!        write(66,*) ,'2', temp1,temp2,temp3,'fine2'
!        write(66,*) ,'3', RLGHTS,'fine3'
!        write(66,*) ,'4', RLIGHT,'fine4'

      end

!***************************************************************

      subroutine rintens(t,itot,fday,iinst)

! computes light intensity during a day given itot and fday
!
! uses formula 5.7 in WASP user manual
!
! needs ITOT (average incident light intensity over one whole day)
! and FDAY (fraction of day length)

      implicit none

      real t            !day [0-365]
      real itot      !average incident light intensity over one whole day
      real fday      !fraction of day length [0-1]
        real iinst      !instantaneous light intensity at time t (return)
      
      real pi
      parameter( pi = 3.14159 )

      real tday,aux

      real t_old,iinst_old
      save t_old,iinst_old
      data t_old,iinst_old /0.,0./      !hope we are not in Antarktic

      if( t .eq. t_old ) then            !if for same daytime we are done
        iinst = iinst_old
        return
      else
        t_old = t
      end if

!      write(6,*) ,'rintens1',iinst

      tday = t - int(t)                  !fraction in day
      if( tday .lt. 0. ) tday = tday + 1.      !negative days
      tday = tday - 0.5                  !maximum at noon

!        intens = itot * (300./86400.)         !LIGHTFIX !FIXME

      aux = pi / fday

      if( abs(tday) .le. fday/2 ) then
        iinst = 0.5 * itot * aux * cos( tday * aux )
      else
        iinst = 0.
      end if

      iinst_old = iinst

!        write(66,*) ,'intens' ,t, tday,itot,iinst
      end

!********************************************************************

      SUBROUTINE KAWIND (WS, TW, TA, depth, WTYPE,  RK)

!*     THIS SUBROUTINE CALCULATES:
!*
!*              RK = REAERATION COEFFICIENT (RK) (M/DAY)*
!*
!*     GIVEN:
!*              WS = Wind Speed (WS) (m/s)
!*              TA = Temperature of the Air  (Degrees C)
!*              TW = Water Temperature (Degrees C)
!*
!* Using the method presented in:
!*
!*           Jour. of Env Eng, Vol. 109, NO.3,PP.731-752,
!*           June 1983, Author: D.J.O'Connor, TITLE: "Wind Effects on
!*           Gas- Liquid Transfer Coefficients"
!*
!*====================================================================
!*
!* THIS SUBROUTINE WAS WRITTEN BY:
!*     SANDRA BIRD
!*     USAE WATERWAYS EXPERIMENT STATION (WES-ES-Q)
!*     VICKSBURG, MISSISSIPPI
!* AND MODIFIED BY JAMES L. MARTIN
!*
!*
!*====================================================================
!*
!*   Parameters used in the model include:
!*
!*        Transitional Shear Velocity - UT(cm/sec)
!*        Critical Shear Velocity - UC (cm/sec)
!*        Vonkarman's Constant (KARMAN)
!*        Equilibrium Roughness - ZE (cm)
!*        1/LAM Is a Reynold's Number
!*        GAM is a a Nondimensional Coefficient Dependent on
!*        Water Body Size (WTYPE).
!*        LAM, GAM, UT, UC and ZE are Dependent on Water Body
!*        Size (See O'Conners Paper for Table of Values).
!*
!*       UT       UC      ZE    LAM     GAM
!*       10.0     11.    .35    3.0     5.          Large Scale
!*       10.0     11.    .25    3.0     6.5         Intermediate
!*        9.      22.    .25   10.     10.          Small Scale
!*
!C******************************************************************

      implicit none

      real R0MIN
      parameter( R0MIN = 1.e-15 )
!      INCLUDE 'weutro.h'
!*
!*  Declarations:
!*

        real WS         !wind speed, m/s
        real TW         !water temperature C
        real TA         !air temperature C
        real depth      !defined depth(segmax) in geometrical
        real WTYPE      !type of water body
        real RK         !reareation term calculated in kawind

      REAL XARG1,XARG2,XWDIFF
      Character*1 cont1
      REAL*4 KARMAN, LAM, KA3, TMPAIR
       real DIFF       !diffusivity in water, from TW CM**2/SE
        real VW         !auxiliar variable to calculate RK=rearsg
        real VA         !viscosity of air, CM**2/SEC  function of TA
        real PA         !density of air g/cm**3
        real PW         !density of water g/cm**3
        real WH         !constant in kawind: =1000
        real SRCD       !constant in kawin: =0.04
        real EF         !auxiliar variable
        real F1         !auxiliar variable
        real F2         !auxiliar variable
        real FP1        !auxiliar variable
        real FP2        !auxiliar variable
        real FP3        !auxiliar variable
        real FP4        !auxiliar variable
        integer N
        real SRCD2      !auxiliar variable
        real ERR        !auxiliar variable
        real  US        !auxiliar variable
        real Z0         !auxiliar variable
        real RK1        !auxiliar variable
        real GAMU       !!auxiliar variable
        real RK2        !!auxiliar variable
        real RK3        !!auxiliar variable

!       !!!!???!!!decidere      quale tenere depth(iseg) o depth. main1 o param?
        real ut !Transitional Shear Velocity, cm/s defined in kawind
        real uc !Critical Shear Velocity, cm/sec  defined in kawind
        real ze !Equilibrium Roughness, cm
        real gam!Nondimensional Coefficient Dependent on water body type in kaw.

      integer IWTYPE

      real*4 cddrag
      COMMON /KAHOLD/ ut, uc, ze, lam, gam, TMPAIR
      data cont1/'$'/
!*
!*   Determine Water Body Type, if WTYPE=0., then default is large
!*   Water Body:
!*
!CSC
!RBA      XARG1 = ABS(WTYPE)
!dd      XWDIFF= WTYPE - 3.0
!dd      XARG2 = ABS(XWDIFF)
      IWTYPE = NINT(WTYPE)
      IF (IWTYPE .eq. 3) then
!RBA      IF (XARG1.LT.R0MIN.OR.XARG2.LT.R0MIN) THEN
!ddd      IF (XARG2.LT.R0MIN) THEN
         UT = 10.0
         UC = 11.0
         ZE = 0.35
         LAM = 3.0
         GAM = 5.0
!CSC
!ddd         XWDIFF = WTYPE - 1.0
!ddd         XARG1 = ABS(XWDIFF)
      ELSE IF (IWTYPE .eq. 1 ) THEN
!ddd         IF (XARG1 .LT. R0MIN) THEN
            UT = 9.0
            UC = 22.0
            ZE = 0.25
            LAM = 10.0
            GAM = 10.0
      ELSE if( WTYPE .eq. 2 ) then
            UT = 10.0
            UC = 11.0
            ZE = 0.25
            LAM = 3.0
            GAM = 6.5
      else
!            write(6,*) 'iwtype: ',iwtype
            stop 'error stop'
      END IF
!*
! CALCULATE DIFFUSIVITY OF OXYGEN IN WATER (DIFF) (CM**2/SEC), VISCOSIT
! OF WATER (VW) (CM**2/SEC),VISCOSITY OF AIR (VA) (CM**2/SEC),DENSITY
! OF WATER (PW) (G/CM**3), DENSITY OF AIR (PA) (G/CM**3)
      DIFF = 4.58E-07*TW + 1.2E-05
!  NOTE: IF OTHER CHEMICALS WERE USED, THEN THEIR DIFFUSIVITIES
!  MAY VARY. FOR EXAMPLE FOR TCDD:   (JLM)
!      DIFF=4.83E-6
!
      VW = 0.0164 - .00024514*TW
      VA = .133 + .0009*TA
      PA = .00129 - .0000040*TA
      PW = 1.00
      WS = WS*100.
      RK = 1.
! USE NEWTON RAPHSON METHOD TO CALCULATE THE SQUARE ROOT OF THE DRAG
! COEFFICIENT
      N = 0
! PARAMETERS USED IN THE MODEL INCLUDE TRANSITIONAL SHEAR VELOCITY - UT(
! CRITICAL SHEAR VELOCITY - UC (CM/SEC); VONKARMAN'S CONSTANT (KARMAN);
! EQUILIBRIUM ROUGHNESS - ZE (CM); 1/LAM IS A REYNOLD'S NUMBER; GAM IS
! NONDIMENSIONAL COEFFICIENT DEPENDENT ON WATER BODY SIZE.  LAM, GAM, UT
! UC AND ZE ARE DEPENDENT ON WATER BODY SIZE
      KARMAN = 0.4
      KA3 = KARMAN**.3333
      WH = 1000.
! MAKE INITIAL GUESS FOR SQUARE ROOT OF THE DRAG COEFFICIENT
      SRCD = 0.04
 1000 CONTINUE
      N = N + 1
! CALCULATE VALUE OF FUNCTION(F2) AND DERIVATIVE OF FUNCTION(FP)
      EF = EXP ( - SRCD*WS/UT)
      F1 = LOG ((WH/ZE) + (WH*LAM/VA)*SRCD*WS*EF)
      F2 = F1 - KARMAN/SRCD
      FP1 = 1./((WH/ZE) + (LAM*WH/VA)*SRCD*WS*EF)
      FP2 = ((WH*LAM)/(VA*UT))*SRCD*WS**2*EF
      FP3 = (WH*LAM/VA)*WS*EF
      FP4 = FP1*(FP2 + FP3) + (KARMAN/(SRCD**2))
! CALCULATE A NEW GUESS FOR SQUARE ROOT OF DRAG AND COMPARE TO
! PREVIOUS GUESS AND LOOP BACK THROUGH N-R WITH NEW GUESS IF
! APPROPRIATE
      SRCD2 = SRCD - F2/FP4
      ERR = ABS (SRCD - SRCD2)
      IF (ERR .GT. 0.0005 .AND. N .LT. 8) THEN
         SRCD = SRCD2
         GO TO 1000
      END IF
      IF (ERR .GT. 0.005 .AND. N .GE. 8) GO TO 1010
      CDDRAG = SRCD**2
      US = SRCD*WS
      Z0 = 1./((1./ZE) + LAM*US*EXP ( - US/UT)/VA)
      WS = WS/100.
      IF (WS .LT. 6.0) GO TO 1020
      IF (WS .GE. 6.0 .AND. WS .LE. 20.0) GO TO 1030
      IF (WS .GT. 20.0) GO TO 1040
! CALC 1050S VALUES FOR WINDSPEEDS LESS THAN 6.0 M/SEC
 1020 CONTINUE
      RK1 = (DIFF/VW)**0.666667*SRCD*(PA/PW)**0.5
      RK = RK1*KA3*WS/GAM
      RK = RK*3600.*24.
      GO TO 1050
! CALC 1050S VALUES FOR WINDSPEED GREATER THAN 20 M/S
 1040 CONTINUE
!RBA--Date: Tuesday, 1 June 1993.  Time: 09:15:47.
!      RK = (DIFF*PA*VA*US/(0.1*PW*VW))**0.5
      RK = (DIFF*PA*VA*US/(KARMAN*ZE*PW*VW))**0.5
      RK = RK*3600.*24./100.
      GO TO 1050
! CALC 1050S VALUES FOR WINDSPEED BETWEEN 6 AND 20 M/S
 1030 CONTINUE
      GAMU = GAM*US*EXP ( - US/UC + 1.)/UC
      RK1 = (DIFF/VW)**.6667*KA3*(PA/PW)**0.5*US/GAMU
      RK2 = (DIFF*US*PA*VA/(KARMAN*Z0*PW*VW))**0.5
      RK3 = (1./RK1) + (1./RK2)
      RK = 1./RK3
      RK = RK*3600.*24./100.
      GO TO 1050
 1050 continue
      GO TO 1060
 1010 CONTINUE
      WRITE (6, 6000)
 6000 FORMAT(5X,'SOLUTION DID NOT CONVERGE')
 1060 CONTINUE
      rk = rk/depth
!ddd
      RETURN
      END

!********************************************************************

      SUBROUTINE KAHYDRA (K2HYDRA) 

!     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:58.

      implicit none
      INCLUDE 'weutro.h'
      REAL K2HYDRA
      real CFOREA,AVDEPE,AVVELE,REAK,EXPREV,EXPRED,TRANDP,DIF
!
!
!                      Calculate Oxygen Reaeration
!
      CFOREA = 1.0
      AVDEPE = H
      AVVELE = ABS(VEL)
!
!
!         Calculate reaeration coefficient for free-flowing reach
!         Calculate reaeration coefficient as a power function of
!         average hydraulic depth and velocity; determine exponents
!         to depth and velocity terms and assign value to REAK
!
      IF (AVDEPE .GT. 0.61) GO TO 1000
!
!          Use Owen's formulation for reaeration
!
      REAK = 5.349
      EXPREV = 0.67
      EXPRED = - 1.85
      GO TO 1010
!
 1000 CONTINUE
!
!       Calculate transition depth; transition depth determines
!       which method of calculation is used given the current
!       velocity
!
      IF (AVVELE .GE. 0.518) GO TO 1020
      TRANDP = 0.0
      GO TO 1030
 1020 CONTINUE
      TRANDP = 4.411*(AVVELE**2.9135)
 1030 CONTINUE
!
      DIF = AVDEPE - TRANDP
      IF (DIF .GT. 0.0) GO TO 1040
!
!                 Use Churchill's formulation for reaeration
!
      REAK = 5.049
      EXPREV = .969
      EXPRED = - 1.673
      GO TO 1050
!
 1040 CONTINUE
!
!                 Use O'Connor-Dobbins formulation for reaeration
!
      REAK = 3.93
      EXPREV = 0.5
      EXPRED = - 1.5
!
 1050 CONTINUE
!
 1010 CONTINUE
!
!                               Calculate reaeration coefficient
!
      K20 = REAK*(AVVELE**EXPREV)*(AVDEPE**EXPRED)
      IF (K20 .GT. 24.) K20 = 24.
      K2HYDRA = K20*1.028**STP20
 1060 CONTINUE
      RETURN
      END

!****************************************************

        function oxysat(temp,salt)

! computes oxygen level at saturation

        implicit none

        real oxysat
        real temp,salt

        oxysat = (14.6244-0.367134*temp+4.4972E-3*temp**2 &
     &     -0.0966*salt+0.00205*salt*temp+2.739E-4*salt**2)
!     &    /32.

        end

!****************************************************

      SUBROUTINE BENTFLUX

!     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:54.

      implicit none

      INCLUDE 'weutro.h'

      REAL AUX, FLUXP, FLUXN
!
!              Benthic ammonium and phosphate fluxes
!
       AUX = SA * 0.001   !conversion to g/day, fluxes are [mg/(m**2 day)]

         FLUXN = FNH4 (ISEG)*AUX
         CD (1, ISEG) = CD (1, ISEG) + FLUXN
         FLUXP = FPO4 (ISEG)*AUX
         CD (3, ISEG) = CD (3, ISEG) + FLUXP

!      write(77,*) FLUXN,ISEG,AUX,FNH4(ISEG)      !ggu
!
      RETURN
      END

!********************************************************************

        subroutine param_read

        implicit none
        INCLUDE 'weutro.h'

        integer ndim
        parameter(ndim=10)
        integer n_params
        real v_params(ndim)

        integer i

        open(1,file='michol.dat',status='old',form='formatted')
        read(1,*) n_params
        if( n_params .gt. ndim ) stop 'error stop: ndim - n_params'
        read(1,*) (v_params(i),i=1,n_params)
        close(1)

        K1C = v_params(1)
        K1D = v_params(2)

        write(6,*) 'using params in weutro from file michol.dat'

        end


!********************************************************************

      subroutine param_venezia

      implicit none
      INCLUDE 'weutro.h'

!        IS2 = 63       !optimum for Steele input Kj/m2/day, (or Kj/m2/300sec?dmk 20/9/2013)timestep 300 sec
        IS2 = 20        !optimum for Steele input Watt*h/mq/300sec (coherently with the previous line-dmk 20/9/2013)FixME
        IS2 = 200       !optimum for Steele input Watt*h/mq/300sec (coherently with the previous line-dmk 20/9/2013)FixME
        KESG(1)=0.85     !m-1 range: 0.1-5
!        KGRZ=0.8        !run07KGRZ= 1.2 primi 2 run 2013 mgl (inclusi nuovi_input)
        KGRZ=0.80
        KDZ = 0.150       ! 0.168 run08
!        K1C = 2.        !primi 2 run 2013
!        K1C = 2.5       ! test15 e test16 (squentin) 
!         K1C = 1.025      ! test17 (squentin)
         K1C = 0.45       ! test19 (squintin) originale 0.9

!        K58C = 0.44     ! test17 (squentin) increase(*2) remineralization of OP
!        k58C = 0.66     ! test19 (squintin) increase 10% remineralization of P
         k58C = 0.60     ! ultimo valore, provare 0.20
         KMNG1 = 0.125  ! test22 (squintin) decrease KMNG1(0.072),increase ammonia preference (PNH3G1)
         KMPG1= 0.096
      
        EFF = 0.5
        K1013C=0.037    !day-1 prova 4 agosto (primi 2 run)
!       K1013C =0.075   !test21(squintin) 0.035 descrease 50% remineralization DON (0.075 original)
         FON=1.0         !prova 4 agosto
!        FON=0.5        !prova 22 luglio 2015
         FOP=1.          !prova 4 agosto
!       K140C=0.18      !prova 11 agosto !che sia per questo che va a valori negativi? 24 sett 2013
        K140C=1.6      !increase denitrification at 1 d-1. San quintin value
     
        K1320C=0.010    !rate for nitrification 0.09-0.3
        KNO3=0.020          !mgO2/L       

        SOD1D(1)=0.1   

        write(6,*) 'using params in weutro for Venezia'

        end

!********************************************************************

      subroutine param_taranto

      implicit none
      INCLUDE 'weutro.h'

        K1C=2.88
        K1G=1.
        K1D=0.05
        KGRZ=1.8
        KDZ=0.1
        SOD1D(1)=0.2
        IS2=300/2.06
        KE(1)=0.5
        KESG(1)=1

        write(6,*) 'using params in weutro for Taranto'

        end

!********************************************************************

        subroutine weutro_debug(debug)

        implicit none

        logical debug

        include 'weutro.h'

        wdebug = debug

        end

!********************************************************************

        subroutine eutro_replay

        implicit none

        integer nstate
        parameter (nstate=9)

        integer id
        real t                  ![day]
        real dt                 !time step [day]
        real vol                ![m**3]
        real depth              ![m]
        real vel                ![m/s]
        real uws                ![m/s]
        real stp                ![Celsius]
        real sal                ![g/L], e.g., 30.
        real qss                !solar rad from shyfem
        real c(nstate)          ![mg/L]
        real cold(nstate)       !old state variable (return)
        real cds(nstate)        !source term [g/day]
        real loads(nstate)      !loading for c [g/(m**3 day)]

        integer i

        write(6,*) 'initializing...'
        call eutroini

        write(6,*) 'reading...'
        call read_debug_param(id,nstate,t,dt &
     &                ,vol,depth,vel,stp,sal,qss,c,cold,cds,loads)

        call write_debug_param(id,nstate,t,dt &
     &                ,vol,depth,vel,stp,sal,qss,c,cold,cds,loads)

        do i=1,nstate
          c(i) = cold(i)
        end do

        write(6,*) 'running...'
        call eutro0d(id,t,dt,vol,depth,vel,uws,stp,sal,qss,c,loads)
        write(6,*) 'finished...'

        call write_debug_param(id,nstate,t,dt &
     &                ,vol,depth,vel,stp,sal,qss,c,cold,cds,loads)

        end

!********************************************************************

        subroutine read_debug_param(id,nstate,t,dt &
     &                ,vol,depth,vel,stp,sal,qss,c,cold,cds,loads)

        implicit none

        integer id              !id of box
        integer nstate          !number of state variables
        real t                  ![day]
        real dt                 !time step [day]
        real vol                ![m**3]
        real depth              ![m]
        real vel                ![m/s]
        real stp                ![Celsius]
        real sal                ![g/L], e.g., 30.
        real qss                !solar rad from shyfem
        real c(nstate)          ![mg/L]
        real cold(nstate)       !old state variable (return)
        real cds(nstate)        !source term [g/day]
        real loads(nstate)      !loading for c [g/(m**3 day)]

        real tair,wspeed,turbid,fracday,itotal
        real ice,flnh4,flpo4
        logical tops,bots,seds

        integer ns,i

        read(5,*) id,ns,t,dt
        if( ns .ne. nstate ) stop 'error stop read_debug_param: nstate'

        read(5,*) vel,stp,sal,vol,depth,qss
        read(5,*) tair,wspeed,turbid,fracday,itotal
        read(5,*) ice,flnh4,flpo4
        read(5,*) tops,bots,seds
        read(5,*) (c(i),i=1,nstate)
        read(5,*) (cold(i),i=1,nstate)
        read(5,*) (cds(i),i=1,nstate)
        read(5,*) (loads(i),i=1,nstate)

        call set_internal_param(tair,wspeed,turbid,fracday,itotal &
     &                  ,ice,flnh4,flpo4 &
     &                  ,tops,bots,seds)

        end

!***************************************************************

        subroutine write_debug_param(id,nstate,t,dt &
     &                ,vol,depth,vel,stp,sal,qss,c,cold,cds,loads)

        implicit none

        integer id              !id of box
        integer nstate          !number of state variables
        real t                  ![day]
        real dt                 !time step [day]
        real vol                ![m**3]
        real depth              ![m]
        real vel                ![m/s]
        real stp                ![Celsius]
        real sal                ![g/L], e.g., 30.
        real qss                !solar rad from shyfem
        real c(nstate)          ![mg/L]
        real cold(nstate)       !old state variable (return)
        real cds(nstate)        !source term [g/day]
        real loads(nstate)      !loading for c [g/(m**3 day)]

        real tair,wspeed,turbid,fracday,itotal
        real ice,flnh4,flpo4
        logical tops,bots,seds

        integer i

        call get_internal_param(tair,wspeed,turbid,fracday,itotal &
     &                  ,ice,flnh4,flpo4 &
     &                  ,tops,bots,seds)

        write(6,*) id,nstate,t,dt
        write(6,*) vel,stp,sal,vol,depth,qss
        write(6,*) tair,wspeed,turbid,fracday,itotal
        write(6,*) ice,flnh4,flpo4
        write(6,*) tops,bots,seds
        write(6,*) (c(i),i=1,nstate)
        write(6,*) (cold(i),i=1,nstate)
        write(6,*) (cds(i),i=1,nstate)
        write(6,*) (loads(i),i=1,nstate)

        end

!***************************************************************

        subroutine get_internal_param(tair,wspeed,turbid,fracday,itotal &
     &                  ,ice,flnh4,flpo4 &
     &                  ,tops,bots,seds)

! gets internal parameters (that can be altered through functions)

        implicit none

        include 'weutro.h'

        real tair,wspeed,turbid,fracday,itotal
        real ice,flnh4,flpo4
        logical tops,bots,seds

        tair = airtmp
        wspeed = wind
        turbid = kesg(1)
        fracday = fday
        itotal = itot
        ice = xicecvr
        flnh4 = fnh4(1)
        flpo4 = fpo4(1)
        tops = topseg
        bots = botseg
        seds = sedseg

        end

!***************************************************************

        subroutine set_internal_param(tair,wspeed,turbid,fracday,itotal &
     &                  ,ice,flnh4,flpo4 &
     &                  ,tops,bots,seds)

! sets internal parameters (that can be altered through functions)

        implicit none

        include 'weutro.h'

        real tair,wspeed,turbid,fracday,itotal
        real ice,flnh4,flpo4
        logical tops,bots,seds

        airtmp = tair
        wind = wspeed
        kesg(1) = turbid
        fday = fracday
        itot = itotal
        xicecvr = ice
        fnh4(1) = flnh4
        fpo4(1) = flpo4
        topseg = tops
        botseg = bots
        sedseg = seds

        end

!***************************************************************

!      subroutine weutro_test

! test drives weutro

!      integer nstate
!      parameter (nstate=9)

!       integer iddt,nt,iday,n,it
!        integer id
!        real vol,depth,vel,stp,sal
!        real dt,t
!        real rlux
!        real pi
!        real itot,fday,iinst,imax

!      real c(nstate)
!      real loads(nstate)

!       data c /0.05, 0.4, 0.01, 0.05, 2., 11., 0.2, 0.01, 0.015/
!        data loads /0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

!        id = 1234
!      iddt = 300
!        rlat = 45.              !latitude

!        pi = 4.*atan(1.)
!      vol = 1.e+6
!      depth = 1.
!      vel = 1.
!      stp = 20.
!      sal = 18.

!      nt = 86400 / iddt      !time steps per day
!        if( nt*iddt .ne. 86400 ) then
!            stop 'error stop: iddt must devide 86400'
!      end if
!      dt = iddt/86400.
!      t = 0.
!        it = 0

!       call eutroini
!       call luxlen_init('lux.dat')      

!      do iday=1,365

!          stp = 17. + 13. * cos(2.*pi*(iday-172)/365.)   !simple temp curve
!      write(6,*) 'ma qui ci arrivo?'

!        call luxlen(t,itot,fday)
!          itot = itot * (300./86400.)   !for data Donata
!          call get_radiation(iday,rlat,fday,itot,imax)
!          itot = itot / 8.              !just to have same magnitude

!        call wlight(fday,itot)

!        do n=1,nt
!            it = it + iddt
!          t = it / 86400.
!            rlux = 1.
!            call eutro0d(id,t,dt,vol,depth,vel,stp,sal,qss,c,loads)
!        end do

!          write(6,*) t,stp,c
!          write(75,'(11e14.6)') t,stp,c
!      end do

!      end


!**************************************************************
!**************************************************************
!**************************************************************
! restart files -> fill in real routines
!**************************************************************
!**************************************************************
!**************************************************************

        subroutine write_restart_eco(iunit)
        implicit none
        integer iunit
        integer nstate,nkn,i
        nstate = 0
        nkn = 0
        write(iunit) nstate,nkn
        end
        subroutine skip_restart_eco(iunit)
        implicit none
        integer iunit
        integer nstate,nkn,i
        read(iunit) nstate,nkn
        do i=1,nstate
          read(iunit)
        end do
        end
        subroutine read_restart_eco(iunit)
        implicit none
        integer iunit
        call skip_restart_eco(iunit)
        end

!***************************************************************

!        program main
!        call eutro_replay
!        call weutro_test
!        end

!***************************************************************

