
!--------------------------------------------------------------------------
!
!    Copyright (C) 2005-2006,2008-2011,2015-2019  Christian Ferrarin
!    Copyright (C) 2008,2010,2012-2020  Georg Umgiesser
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

! module definitions for sediments

! revision log :
!
! 01.03.2005	ccf	(sedi3d_f1.f) coming from sedi3d_e3_multi7.f
! ...				new cohesive routine
! 01.03.2005	ccf	(sedi3d_f2.f) add mixing thickness
! 01.03.2005	ccf	(sedi3d_f3.f) merge layer 1 and 2 if bedn(1) < bmix
! 01.03.2005	ccf	(sedi3d_f4.f) update element depth
! 01.03.2005	ccf	(sedi3d_f5.f) get viscosity from new routine
! 01.03.2005	ccf	(sedi3d_f6.f) active layer = bottom roughness heigh
! 01.04.2005	ccf	(sedi3d_f7.f) style changes, create sedt05
! 01.04.2005	ccf	(sedi3d_f8.f) initialize percbd from file
! 01.04.2005	ccf	(sedi3d_f8.f) change rhos per cohesive sed
! 01.06.2005	ccf	(sedi3d_f13.f) adapt 3D, bottom layer and total depth
! ...				bugfix in computing bedn(1,3), dimension in tuek
! 01.06.2005	ccf	(sedi3d_f15.f) change deposition density,
! ...				add consolidation
! 01.06.2005	ccf	(sedi3d_f16.f) change units, bug fix in TCONC1
! 01.06.2005	ccf	(sedi3d_f17.f) bug fix in upedepth
! 01.06.2005	ccf	(sedi3d_f18.f) bug fix in checkbed, 
! ...				adapt to sedtrans05-h5.f
! 01.07.2005	ccf	(sedi3d_f20.f) adapt to sedtrans05-h6.f
! 01.07.2005	ccf	(sedi3d_f21.f) bug fix in updepth (as subdry.f)
! 01.07.2005	ccf	(sedi3d_f22.f) bug fix in bedload
! 01.08.2005	ccf	(sedi3d_f23.f) eliminate ripple variables
! 01.08.2005	ccf	(sedi3d_f24.f) bed slope contribution. adjust updepth
! 01.08.2005	ccf	(sedi3d_f25.f) change deposition characteristics
! 01.09.2005	ccf	(sedi3d_f26.f) bug fix in bedman,
! ...				change boundary for cohesive
! 01.09.2005	ccf	(sedi3d_f27.f) separete erosion/deposition in bedman
! ...				number of layer computed each time
! 01.11.2005	ccf	(sedi3d_f28f) adapt for 3D version
! 01.11.2005	ccf	(sedi3d_f29f) bed slope threshold, bug fix in updepth
! 01.11.2005	ccf	(sedi3d_f30f) bug fix in bedslope and in getmd
! 01.11.2005	ccf	(sedi3d_f31f) bed slope by gradients
! 01.11.2005	ccf	(sedi3d_f32f) last layer not smaller than 0.1 m
! 01.11.2005	ccf	(sedi3d_f33f) compute vertical mixing coeffcients
! 01.01.2006	ccf	(sedi3d_f34f) bug fix in upedepth and blimit
! 01.02.2006	ccf	(sedi3d_f35f) new suspco routine
! 01.02.2006	ccf	(sedi3d_f36f) adapt Van Rijn (C0) in sedtrans05-h8.f
! 01.02.2006	ccf	(sedi3d_f37f) bug fix in checkbed and other things
! 01.05.2006	ccf	(sedi3d_f38f) bug fix in line 1119. Introduce KCOES
! 01.05.2006	ccf	(sedi3d_f39f) bugs nonco. non used edr in cohse.
! ...				no limit percbd. pers(1)>0.1 in line 1120
! 01.05.2006	ccf	(sedi3d_f40f) limit BEDCHA(1,2) to 0.1. 
! ...				bugfix in suspco, better conc in cohse
! 01.06.2006	ccf	(sedi3d_f41f) no transport in case of depth< 0.1m
! 01.06.2006	ccf	(sedi3d_f42f) read constants from file
! 01.06.2006	ccf	(sedi3d_f43f) add limcoh
! 01.06.2006	ccf	(sedi3d_f44f) smooth bed elevation, get_timestep
! 01.07.2006	ccf	(sedi3d_f45f) read angle of repose, 
! ...				limit shear velocity, write bathymetry
! 11.04.2008	ggu&ccf	treatment of boundaries slightly changed
! 16.04.2008	ggu&ccf	bugfix calling lin (must pass double precision 0.)
! 22.04.2008	ggu	advection parallelized
! 27.04.2008	ccf	new check, other changes
! 28.04.2008	ggu	call to nrdpar, nrdnls read with double precision
! 29.04.2008	ccf	bug fix in depbed
! 20.10.2008	ccf	add routine for computing the active layer
! 30.04.2009	ccf	add adjtime for initialization to reset bh
! 03.10.2009	ccf	compute transport only for GD
! 23.03.2010	ggu	changed v6.1.1
! 13.04.2010	ccf	introduce SSC effect on density
! 27.07.2010	ccf	settling velocity for cohesives as function of space
! 08.10.2010	ggu	changed VERS_6_1_13
! 20.06.2011	ccf	load for eros/depos for both cohesive and non-cohesive
! 20.06.2011	ccf	deleted suspco routine
! 30.03.2012	ggu	changed VERS_6_1_51
! 21.06.2012	ggu	changed VERS_6_1_54
! 25.01.2013	ggu	changed VERS_6_1_62
! 13.06.2013	ggu	changed VERS_6_1_65
! 18.06.2014	ggu	changed VERS_6_1_77
! 21.10.2014	ggu	changed VERS_7_0_3
! 05.11.2014	ggu	changed VERS_7_0_5
! 05.12.2014	ggu	changed VERS_7_0_8
! 23.12.2014	ggu	changed VERS_7_0_11
! 09.01.2015	ggu	changed VERS_7_0_12
! 19.01.2015	ccf	ia_out1/2 introduced
! 19.01.2015	ggu	changed VERS_7_1_3
! 10.02.2015	ggu	new read for constants
! 26.02.2015	ggu	changed VERS_7_1_5
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 24.07.2015	ggu	changed VERS_7_1_82
! 30.07.2015	ggu	changed VERS_7_1_83
! 23.09.2015	ggu	changed VERS_7_2_4
! 12.10.2015	ggu	changed VERS_7_3_3
! 22.10.2015	ggu	changed VERS_7_3_8
! 05.11.2015	ggu	changed VERS_7_3_12
! 16.12.2015	ggu	changed VERS_7_3_16
! 18.12.2015	ggu	changed VERS_7_3_17
! 11.03.2016	ggu	changed VERS_7_5_5
! 31.03.2016	ggu&ccf	update to new shyfem version
! 01.04.2016	ggu	changed VERS_7_5_7
! 19.04.2016	ccf	krocks for no erosion-deposition in specific e-types
! 25.05.2016	ggu	changed VERS_7_5_10
! 30.09.2016	ggu	changed VERS_7_5_18
! 10.02.2017	ggu	read in init data from fem files (init_file_fem)
! 13.04.2017	ggu	changed VERS_7_5_25
! 29.08.2017	ccf	read parameters as array
! 02.09.2017	ggu	changed VERS_7_5_31
! 26.09.2017	ggu	changed VERS_7_5_32
! 04.11.2017	ggu	changed VERS_7_5_34
! 14.11.2017	ggu	changed VERS_7_5_36
! 17.11.2017	ggu	readsed split in readsed and initsed
! 05.12.2017	ggu	changed VERS_7_5_39
! 22.02.2018	ccf	adjusted for new sinking velocity
! 03.04.2018	ggu	changed VERS_7_5_43
! 03.02.2019	ggu	adjust z0 etc (GGUZ0)
! 05.02.2019	ggu	upgraded to new version from chris
! 14.02.2019	ggu	set negative conz values to 0
! 15.02.2019	ccf	pass PHIB,PHI100 into nonco (bug)
! 13.03.2019	ggu	changed VERS_7_5_61
! 24.10.2019	ggu	better checking of grainsize percentage in inibed
! 16.02.2020	ggu	femtime eliminated
! 22.09.2020    ggu     correct warnings for PGI compiler
! 21.10.2021    ggu     set rhosa to a reasonable value (was infinite) (GGUBS)
! 21.10.2021    ggu     in call to get_sigma_info() protect nlv
! 23.10.2024    ggu     module definition taken out from sedim_admin.f90
! 12.06.2025    clc     set sedpa to 0 (avoid undefined value)
! 
!****************************************************************************

!==================================================================
      module mod_sediment
!==================================================================

      implicit none

      integer, private, save    :: nkn_sedim = 0
      integer, private, save    :: nlv_sedim = 0

      integer, save	        :: isedi   = 0		!sedi call parameter
      integer, parameter        :: nlbdim  = 10         !max number of bed layer

      character*4, save         :: what = 'sedt'

      double precision, save    :: da_sed(4) = 0
      double precision, save    :: da_ssc(4) = 0

      real, allocatable, save	:: tcn(:,:) 	!total sediment concentration [kg/m3]
      real, allocatable, save	:: tmsus(:)	!total suspended sediment load [kg/s] (> 0 -eros, < 0 depo)
      real, allocatable, save	:: bdens(:,:)	!dry bulk density of bottom sediments [kg/m**3]
      real, allocatable, save	:: bleve(:,:)	!depth below sediment surface of sediments [m]

!==================================================================
      contains
!==================================================================

      subroutine mod_sedim_init(nkn,nlv)

      integer  :: nkn
      integer  :: nlv

      if( nlv == nlv_sedim .and. nkn == nkn_sedim ) return

      if( nlv > 0 .or. nkn > 0 ) then
        if( nlv == 0 .or. nkn == 0 ) then
          write(6,*) 'nlv,nkn: ',nlv,nkn
          stop 'error stop mod_sedim_init: incompatible parameters'
        end if
      end if

      if( nkn_sedim > 0 ) then
        deallocate(tcn)
        deallocate(tmsus)
        deallocate(bdens)
        deallocate(bleve)
      end if

      nlv_sedim = nlv
      nkn_sedim = nkn

      if( nkn == 0 ) return

      allocate(tcn(nlv,nkn))
      allocate(tmsus(nkn))
      allocate(bdens(nlbdim,nkn))
      allocate(bleve(nlbdim,nkn))

      tcn = 0.
      tmsus = 0.
      bdens = 0.
      bleve = 0.

      end subroutine mod_sedim_init

!==================================================================
      end module mod_sediment
!==================================================================

!==================================================================
      module mod_sediment_para
!==================================================================

      implicit none

      double precision, save, dimension(8) :: sedpa = 0 !sediment parameter vector
      real, allocatable, save   :: gsc(:)		!grainsize class read from str
      real, allocatable, save   :: tue(:)		!initial erosion threshold
      real, allocatable, save   :: prin(:,:)		!initial percetage [0,1]

      double precision, save    :: KCOES   = 0.15d0     !Fraction of mud for sediment to be cohesive [0-1]
      double precision, save    :: LIMCOH  = 0.000063d0 !grainsize limit for cohesive sediment [m]
      double precision, save    :: SMOOTH  = 1.0d0      !smoothing factor for morphodynamic [0-1]
      double precision, save    :: ANGREP  = 32.0d0     !angle of repose [rad]
      double precision, save    :: MORPHO  = 1.d0       !Morphological factor
      double precision, save    :: RHOSED  = 2650.d0    !sediment grain density
      double precision, save    :: POROS   = 0.4d0      !bed porosity [0,1]
      double precision, save    :: SURFPOR = 0.6d0      !surficial porosity [0,1]
      integer, save             :: IOPT    = 5          !SEDIMENT TRANSPORT FORMULA OPTION NUMBER
      integer, save             :: NBCC                 !same as NBCONC but equal 0 if no cohesive

      real, save		:: difmol	!Molecolar diffusion coefficient [m**2/s]
      logical, save		:: wsetbio = .true.	!temp dependent settling velocity
      !logical, save		:: wsetbio = .false.	!temp dependent settling velocity

      logical, save :: bbudget = .false.		!write budget file

!==================================================================
      end module mod_sediment_para
!==================================================================

