!--------------------------------------------------------------------------
!
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
!
! routines for coupling WW3 spectral wave model
!
! revision log :
!
! 04.07.2019    ggu     written from scratch
! 24.03.2022    ggu     newly started with Aron
! 05.05.2022    aar     lots of changes from Aron
! 15.10.2022    ggu     shympi_exchange_array substituted with shympi_l2g_array
! 03.12.2024    lrp     ww3 restored
! 15.02.2025    ccf     ww3 coupling updated
! 01.04.2025    ccf     wind partitioning
!
!**************************************************************
! DOCS  START   S_wave_ww3

! SHYFEM has been coupled to the WAVEWATCH III (WW3) wave model.
! The SHYFEM model was coupled to WW3 based on a so called ``hard coupling'', e.g.
! binding the models directly without any so called coupling library such as
! (ESMF, OASIS, PGMCL or others). The benefit is on the hand, minimal memory
! usage, very limited code changes in both models. Top-level approach using the
! coupled model framework based on having a routine for initialization,
! computation and finalization, that allows a neat inclusion of WW3 in any kind
! of flow model.
! 
! \subsubsection{Installation and compilation}
! %The coupled model is disributed within the SHYFEM distribution and the WW3 code
! %is located in the |shyfem/WW3| subfolder. The latest version of the WW3 code is
! %available on GitHub at \url{https://github.com/NOAA-EMC/WW3}.
! 
! In order to compile the coupled SHYFEM-WW3 model, the variable |WW3| in the 
! file |Rules.make| should be set as |WW3 = true| and |WW3DIR| should
! point to the folder containing the |WW3| code (e.g., WW3DIR = \$(HOME)/bin/WW3). 
! |SHYFEM| should be compiled for running in parallel by setting 
! |PARALLEL\_MPI = NODE|.
! 
! Moreover, WW3 also needs additional libraries:
! \begin{itemize}
! \item |netcdf| compiled with the same compiler.  You must set |NETCDF = true|
! and specify the variable |NETCDFDIR| to indicate the directory where the 
! libraries and its include files can be found.
! \item |METIS| and |PARMETIS| compiled with the same compiler. You must set
! |PARTS = PARMETIS| and set the variables |METISDIR| and |PARMETISDIR|
! to indicate the directory where the libraries and its include files can be
! found.
! \end{itemize}
! 
! The compilation of |WW3| can be customized changing a set of model options,
! defined with the variable |WW3CFLAGS| in the file |fem3d/Makefile|. They
! correspond to the parameters listed in the file |switch| needed by the stand 
! alone WW3 installation. For a detailed description of the options see the WW3 
! manual (section 5.9).
! 
! You can now follow the general |SHYFEM| installation instruction to compile the 
! code.
! 
! It is however required to download, compile and install the stand alone |WW3| 
! model from the ERDC GitHub at \url{https://github.com/erdc/WW3/tree/ww3_shyfem}. 
! This step is needed for having |WW3| source code and pre- and post-processing 
! tools (e.g., ww3\_grid). 
! It is important that the so called ``switch'' file of |WW3| (see |WW3| documentation) 
! contains the same parameters listed with variable |WW3CFLAGS| are identical 
! for both modules. The default setup provided with |WW3CFLAGS| is basically 
! identical to the setup of SHOM and Meteo France as well as the USACE based 
! on implicit time stepping on unstructured grids. Basically, there is no need 
! for modification of the ``switch'' file with the given settings but all 
! settings are supported as long the unstructured grid option in WW3 is used.
! 
! For installing and comping the stand alone |WW3| model you can follow
! the following procedure (see the |WW3| manual for more informations):
! \begin{enumerate}
! \item download the ww3\_shyfem branch from the ERDC git repository: 
! git clone --branch ww3\_shyfem https://github.com/erdc/WW3/tree/ww3\_shyfem;
! \item move to the WW3 directory;
! \item set the NetCDF path: with |export WWATCH3\_NETCDF=NC4| and 
! |export NETCDF\_CONFIG=/ditectory-of-netcdf/nc-config|;
! \item setup the model using |./model/bin/w3\_setup /home/model/bin/WW3/model -c 
! <cmplr> -s <swtch>| with |cmplr| the compiler |comp| and |link| files (e.g., shyfem
! for files comp.shyfem and link.shyfem) and |swtch| the switch file (e.g., shyfem 
! for a switch file located in ./model/bin having name switch\_shyfem).
! \item compile |WW3| with |w3\_automake|, or for a few programs with
! |./model/bin/w3\_automake ww3\_grid ww3\_shel ww3\_ounf| with name 
! \item the compiled programms are located in |./model/exe/ww3\_grid|
! \end{enumerate}
! 
! \subsubsection{Coupling description}
! The implementation was done by developing two modules, one in
! |SHYFEM| (ww3.f90) and the other one in |WW3| (w3cplshymfen.F90). The 
! SHYFEM coupling module is connected by the ``use'' statement to the 
! modules of the |WW3| code. In this way one can access all fields from 
! the |WW3| in |SHYFEM| and vice versa. In the |WW3| module all spectral based
! quantities are computed and |SHYFEM| is directly assigning the various 
! arrays, based on global arrays, to the certain domains for each of the 
! model decompositions.
! 
! In terms of |WW3| we have used the highest-level implementation based
! on the |WW3| multigrid driver. This allows basically to couple any kind
! of grid type with |SHYFEM| (rectangular, curvi-linear, SMC and
! unstructured). It offers of course the possibility to run |WW3| based
! on a multigrid SETUP and coupled to |SHYFEM|. Both models use their
! native input files for running the model except that the forcing by
! wind within the coupled model is done via |SHYFEM| and the flow field
! is by definition provided based on the coupling to |SHYFEM|. In this
! way both models are fully compatible to the available documentation.
! 
! The two numerical models (SHYFEM and WW3) should exchange all the
! variables that are needed to simulate the current-wave interaction.
! The following variables sare passed by |SHYFEM| to |WW3|:
! \begin{itemize}
! \item water levels;
! \item three-dimensional water currents;
! \item wind components.
! \end{itemize}
! 
! The following physical quantities have been computed in |WW3| and are 
! available for |SHYFEM|:
! \begin{itemize}
! \item significant wave height;
! \item mean wave period;
! \item peak wave frequency;
! \item main wave direction (the where the wave go);
! \item radiation stress,
! \item Charnock parameter;
! \item Friction velocity;
! \item Bernulli head (J term);
! \item Stokes drift velocities;
! \item Stokes drift volume transports;
! \end{itemize}
! 
! The exchange module of SHYFEM contains all the needed infrastructure
! for initializing, calling and finalizing the wave model run. Two
! subroutines (getvarSHYFEM and getvarWW3), which are called before and
! after the call to the wave model in order to obtain currents, water
! levels and on the other integral wave parameters, radiation
! stresses, Stokes velocities and transports.
! 
! The significant wave height, mean wave period and main wave 
! direction are written by |SHYFEM| in the |.wave.shy| file according
! to the values of the variables |idtwav| and |itmwav| in the |wave|
! section of the |SHYFEM| parameter file (see below).
! 
! \subsubsection{Running the coupled SHYFEM-WW3 model}
! Since the wave output are written by |SHYFEM|, the output are
! set to off in |WW3|, which reduces disk usage and significantly reduces 
! the parallel overhead due to the output part. 
! 
! As we have utilized here the implicit scheme in WW3, which was
! developed by Roland \& Partner and the implicit scheme is well
! validated and unconditionally stable. The choice of the time step
! should be according to the physical time scales of the modelled
! processes. 
! 
! In order to run the coupled model, we suggest the following
! procedure:
! \begin{enumerate}
! \item Setup the |SHYFEM| set-up for the region of interest. 
! |SHYFEM| needs a mesh in the .bas format, a parameter file and 
! all forcing files (see |SHYFEM| documentation). In the |SHYFEM| 
! parameter file the section |\$waves| must be activated with 
! \begin{itemize}
! \item |iwave = 2|: coupling via the radiation stress formulation.
! \item |iwave = 3|: coupling via the vortex force formulation.
! \end{itemize}
!
! The time step for coupling with WW3 |dtwave| must be set equal
! to the |SHYFEM| and |WW3| timesteps (TO BE FIXED).
! |idtwav|, |itmwav| should also be set for determining the 
! time step and start time for writing to the output file |.wave.shy| .
! Preprocess the |SHYFEM| grid with a selected number of domains 
! without bandwidth optimization with |shypre -noopti grid.grd|.
! 
! \item Setup the |WW3| model for the region of interest. In the
! running directoty, |WW3| requires:
! \begin{enumerate}
! \item the file |ww3\_grid.nml| containing the spectrum, run, timesteps
! and grid parameterizations (see |WW3| documentation). This file also
! contains the names of mesh and namelist files.
! \item the namelist file containing the variables for setting
! the tunable parameters for source terms, propagation schemes, and 
! numerics.
! \item the file |ww3\_multi.nml| which handles the time steps and 
! the I/O. Since |WW3| receives flow and wind from |SHYFEM| there are 
! no input fields that need to be prescribed. The output part of the 
! wave model itself is handled by |SHYFEM| as well. However, at least
! one variable (e.g. HS) has to be set in |ALLTYPE%FIELD%LIST| as 
! well as the |ALLDATE%FIELD| dates (START, STOP) and timestep (STRIDE).
! \item the numerical mesh in the GMSH format .msh. The |SHYFEM|
! mesh in .grd format can be easily converted to .msh with the command:
! shybas -msh grid.bas (it creates a bas.msh file).
! \end{enumerate}
! 
! Before running the model, the setup |ww3\_grid| tool (found in the
! stand alone WW3/build/bin directory) needs to be run and the output, 
! which is named per default ``mod\_def.ww3'' needs to be renamed 
! with the extension defined in parameter |MODEL(1)\%NAME = 'med'|
! in |ww3\_multi.nml| (e.g., mod\_def.med in the above mentioned case). 
! The |ww3\_grid| tool must be run every time the namelist or grid files 
! are modified. 
! 
! \item run the coupled model with the |shyfem| command. Since the code
! is optimized to run with MPI on domain decomposition, we suggest to
! run the coupled model in parallel using |mpiexec -np np shyfem
! namelist.str| with np the number of processors.
!
! The wave module writes in the WAV file the following output:
! \begin{itemize}
! \item significant wave height [m], variable 231
! \item mean wave period [s], variable 232
! \item mean wave direction [deg], variable 233
! \end{itemize}
!
! \end{enumerate}
!
! DOCS  END
!===========================================================
	module mod_ww3
!===========================================================

	use shympi
	use shympi_internal
	use mod_meteo
	use meteo_forcing_module, only : iatm
	use mod_waves
	use mod_roughness
	use mod_depth 
	use mod_hydro
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw,xgv,ygv

	use wminitmd
	use w3cplshyfem, ONLY: SXXWW3, SXYWW3, SYYWW3, OUTVARWW3
	use w3cplshyfem, ONLY: CHARNWW3, USTWW3
	use w3cplshyfem, ONLY: TAUWWXWW3, TAUWWYWW3
	use w3cplshyfem, ONLY: USSXWW3, USSYWW3, TUSXWW3, TUSYWW3, BHDWW3
	use w3cplshyfem, ONLY: w3cplparam, w3cplrdstr2d, nvarsww3
	use w3gdatmd, ONLY: nseal
	use w3odatmd, ONLY: IAPROC, NAPROC
	use w3idatmd, ONLY: WX0, WXN, WY0, WYN, CX0, CXN, CY0, CYN, WLEV
	use w3idatmd, ONLY: FLCUR, FLLEV
	use yowNodepool, only: np, npa, iplg, ipgl, np_global
	use yowDatapool, only: rtype, istatus, myrank

	implicit none
	logical, save         :: bww3 = .true.
	integer               :: mpiComm = -99
	integer, allocatable  :: nwild_i(:), nwild_gbi(:)
	real, allocatable     :: nwild_gb(:) 
	integer, allocatable  :: nwild_i_shyfem(:), nwild_gbi_shyfem(:)
	real, allocatable     :: nwild_gb_shyfem(:)

	real, allocatable     :: sxxgl(:)      	   !radiation stress xx
	real, allocatable     :: sxx3dlshyfem(:,:)
	real, allocatable     :: sxygl(:)          !radiation stress xy
	real, allocatable     :: sxy3dlshyfem(:,:) 
	real, allocatable     :: syygl(:)          !radiation stress yy
	real, allocatable     :: syy3dlshyfem(:,:)
	real, allocatable     :: charngl(:)	   !Charnock parameter
	real, allocatable     :: charnshyfem(:)
	real, allocatable     :: ustgl(:)	   !Charnock parameter
	real, allocatable     :: ustshyfem(:)
	real, allocatable     :: tauwwxgl(:)       !Sum of wind input to wave and wave to ocean x
	real, allocatable     :: tauwwxshyfem(:)
	real, allocatable     :: tauwwygl(:)       !Sum of wind input to wave and wave to ocean y
	real, allocatable     :: tauwwyshyfem(:) 
	real, allocatable     :: ussxgl(:)	   !Stokes drift velocity x
	real, allocatable     :: ussxshyfem(:)
	real, allocatable     :: ussygl(:)	   !Stokes drift velocity y
	real, allocatable     :: ussyshyfem(:)
	real, allocatable     :: tusxgl(:)         !Stokes drift transport x
	real, allocatable     :: tusxshyfem(:)
	real, allocatable     :: tusygl(:)	   !Stokes drift transport y
	real, allocatable     :: tusyshyfem(:)
	real, allocatable     :: bhdgl(:)	   !Bernulli head (J term)
	real, allocatable     :: bhdshyfem(:)
	real, allocatable     :: currxgl(:)        !Current velocity x
	real, allocatable     :: currxlshyfem(:)
	real, allocatable     :: currygl(:)	   !Current velocity y
	real, allocatable     :: currylshyfem(:)
	real, allocatable     :: wlvgl(:)	   !Water level
	real, allocatable     :: wlvshyfem(:)
	real, allocatable     :: wxvgl(:)	   !Wind velocity x
	real, allocatable     :: wyvgl(:)	   !Wind velocity y
	real, allocatable     :: outvargl(:,:)	   !Output variables
	real, allocatable     :: outvarshyfem(:,:)
	integer, allocatable  :: ipgl_shyfem(:)
	integer, allocatable  :: iplg_shyfem(:)
	integer, allocatable  :: node2domain(:)

!===========================================================

!***********************************************************
!
!***********************************************************
	contains 

	subroutine ww3_init_internal
	USE w3cplshyfem 

	implicit none

	integer ier, idsi, idso, idss, idst, idse
	integer istat, ip, ierr, iproc, j
	logical lexist, has_output_d
	integer iwave, nvar, id
	character(len=30) :: ifname = 'ww3_multi.nml'

	real getpar

	iwave = nint(getpar('iwave'))
	bww3 = ( iwave > 1 )
	if( .not. bww3 ) return

	inquire(file=ifname,exist=lexist)
	if(.not. lexist) then
		ifname = 'ww3_multi.inp'
	endif

!AR: todo: set proper file handles, same as ww3 and check if they are free in shyfem 
	  idsi = 8
	  idso = 9
	  idss = 6
	  idst = 10
	  idse = 6
!WW3 original handles: 8           9           6          10           6
	  mpiComm = MPI_COMM_WORLD 

	  if ( trim(ifname).eq.'ww3_multi.nml' ) then
		call wminitnml ( idsi, idso, idss, idst, idse, trim(ifname), & 
                          mpiComm, './' )
		else

		call wminit ( idsi, idso, idss, idst, idse, & 
                   'ww3_multi.inp', mpiComm, './')
		endif
!
!AR: Init coupling part of shyfem 
!
		call w3cplinit

		allocate(sxxgl(np_global)); sxxgl = 0.
		allocate(sxx3dlshyfem(nlv_global,nkn)); sxx3dlshyfem = 0.
		allocate(sxygl(np_global)); sxygl = 0.
		allocate(sxy3dlshyfem(nlv_global,nkn)); sxy3dlshyfem = 0.
		allocate(syygl(np_global)); syygl = 0.
		allocate(syy3dlshyfem(nlv_global,nkn)); syy3dlshyfem = 0.
		allocate(outvargl(nvarsww3,np_global)); outvargl = 0.
		allocate(outvarshyfem(nvarsww3, nkn)); outvarshyfem = 0.
		allocate(charngl(np_global)); charngl = 0.
		allocate(charnshyfem(nkn)); charnshyfem = 0.
		allocate(ustgl(np_global)); ustgl = 0.
		allocate(ustshyfem(nkn)); ustshyfem = 0.
		allocate(tauwwxgl(np_global)); tauwwxgl = 0.
		allocate(tauwwxshyfem(nkn)); tauwwxshyfem = 0.
		allocate(tauwwygl(np_global)); tauwwygl = 0.
		allocate(tauwwyshyfem(nkn)); tauwwyshyfem = 0.
		allocate(ussxgl(np_global)); ussxgl = 0.
		allocate(ussxshyfem(nkn)); ussxshyfem = 0.
		allocate(ussygl(np_global)); ussygl = 0.
		allocate(ussyshyfem(nkn)); ussyshyfem = 0.
		allocate(tusxgl(np_global)); tusxgl = 0.
		allocate(tusxshyfem(nkn)); tusxshyfem = 0.
		allocate(tusygl(np_global)); tusygl = 0.
		allocate(tusyshyfem(nkn)); tusyshyfem = 0.
		allocate(bhdgl(np_global)); bhdgl = 0.
		allocate(bhdshyfem(nkn)); bhdshyfem = 0.
		allocate(wxvgl(nkn_global)); wxvgl = 0.
		allocate(wyvgl(nkn_global)); wyvgl = 0.
		allocate(currxlshyfem(nkn)); currxlshyfem = 0.
		allocate(currylshyfem(nkn)); currylshyfem = 0.
		allocate(wlvshyfem(nkn)); wlvshyfem = 0.
		allocate(wlvgl(nkn_global)); wlvgl = 0.
		allocate(currxgl(nkn_global)); currxgl = 0.
		allocate(currygl(nkn_global)); currygl = 0.
		allocate(ipgl_shyfem(nkn_global)); ipgl_shyfem = 0
		allocate(iplg_shyfem(nkn)); iplg_shyfem = 0
		allocate(node2domain(nkn_global)); node2domain = 0
!
! ... get the ghost hallo sums for all nodes ...
!
		ALLOCATE(nwild_i(np_global), nwild_gbi(np_global), stat=istat)
		IF (istat/=0) stop 'allocate error nwild_i,nwild_gbi'
		nwild_i = 0
		nwild_gbi = 0
		DO IP = 1, NSEAL
			nwild_i(IPLG(IP)) = 1
		END DO
		call mpi_reduce(nwild_i,nwild_gbi,NP_GLOBAL,MPI_INT, &
                        MPI_SUM, 0, mpiComm, ierr)

		ALLOCATE(nwild_gb(NP_GLOBAL), stat=istat); nwild_gb = 0
		IF (istat/=0) stop 'allocate error nwild_i,nwild_gbi'
		IF (iaproc.eq.1) THEN
			DO IP = 1, NP_GLOBAL
				nwild_gb(IP)=1./REAL(nwild_gbi(IP))
			END DO
		DO iProc=2,naproc
			CALL MPI_SEND(nwild_gb,np_global,MPI_REAL, &
                          iProc-1, 2037, mpiComm, ierr)
		END DO
	ELSE
		CALL MPI_RECV(nwild_gb,np_global,MPI_REAL, &
                       0, 2037, mpiComm, istatus, ierr)
	END IF
	DEALLOCATE(nwild_i, nwild_gbi)
!
! ... now the same for shyfem ...
!
	ALLOCATE(nwild_i_shyfem(nkn_global), &
                nwild_gbi_shyfem(nkn_global), stat=istat)
	IF (istat/=0) stop 'allocate error nwild_i,nwild_gbi'
		nwild_i_shyfem = 0
		nwild_gbi_shyfem = 0
		DO IP = 1, NKN
			nwild_i_shyfem(ip_int_node(IP)) = 1
		END DO
		call mpi_reduce(nwild_i_shyfem,nwild_gbi_shyfem, &
                        nkn_global,MPI_INT, &
                        MPI_SUM, 0, mpiComm, ierr)

		ALLOCATE(nwild_gb_shyfem(nkn_global), stat=istat)
		nwild_gb_shyfem = 0
		IF (istat/=0) stop 'allocate error nwild_i,nwild_gbi'
		IF (iaproc.eq.1) THEN
			DO IP=1,nkn_GLOBAL
				nwild_gb_shyfem(IP)=1./REAL(nwild_gbi_shyfem(IP))
			END DO
			DO iProc=2,naproc
				CALL MPI_SEND(nwild_gb_shyfem,nkn_global,MPI_REAL, &
                          iProc-1, 2037, mpiComm, ierr)
			END DO
		ELSE
			CALL MPI_RECV(nwild_gb_shyfem,nkn_global,MPI_REAL, &
                       0, 2037, mpiComm, istatus, ierr)
		END IF
		DEALLOCATE(nwild_i_shyfem, nwild_gbi_shyfem)

! now compute some basic mappings for shyfem 
		do ip = 1, nkn_inner
			node2domain(ip_int_node(ip)) = myrank + 1
		enddo 

		j = 1
		do ip = 1, nkn_global
			if ( node2domain(ip) == myrank+1 ) then
				iplg_shyfem(j) = ip
				ipgl_shyfem(IP) = j
			j = j + 1
			endif
		end do
!
! init coupling
		call getvarww3(.true.)
		call getvarshyfem(.true.)

! init output
		nvar = 3
		call init_output_d('itmwav','idtwav',da_wav)
		if( has_output_d(da_wav)) then
			call shyfem_init_scalar_file('wave',nvar,.true.,id)
			da_wav(4) = id
		end if

	end subroutine ww3_init_internal
!***********************************************************
!
!***********************************************************
	subroutine ww3_loop_internal

	use basin
	use mod_hydro
	use mod_meteo
        use WMWAVEMD 

	implicit none

	integer n,nsave, tend(2,1), yy,mm,dd,h,m,s, it
	integer ys(8), ierr 

	double precision atime, dtime
	double precision daux

	if( .not. bww3 ) return

	call get_act_dtime(dtime)
	it = dtime
	call convert_time_d('dtwave',daux)
	idcoup = nint(daux)

	if (mod(it,idcoup) .eq. 0 ) then

		call get_absolute_act_time(atime)
		call dts_from_abs_time_to_ys(atime,ys)

		yy = ys(1)
		mm = ys(2)
		dd = ys(3)
		h  = ys(4)
		m  = ys(5)
		s  = ys(6)
 
		call dts_from_abs_time_to_ys(atime+daux,ys)

		yy = ys(1)
		mm = ys(2)
		dd = ys(3)
		h  = ys(4)
		m  = ys(5)
		s  = ys(6)

! 		Set end time of this advance
		tend(1,1) = 10000*yy + 100*mm + dd
		tend(2,1) = 10000*h  + 100*m  + s

		if( .not. bww3 ) return

		call mpi_barrier(mpiComm,ierr)
		call getvarshyfem(.false.) 
		call mpi_barrier(mpiComm,ierr)
		call wmwave ( tend )
		call mpi_barrier(mpiComm,ierr)
		call getvarww3(.false.)
		call mpi_barrier(mpiComm,ierr)
        endif ! idcoup
        
	end subroutine ww3_loop_internal
!***********************************************************
!
!***********************************************************

	subroutine ww3_finalize_internal

	use basin
        use WMFINLMD

	implicit none

	if( .not. bww3 ) return

! here finalize ww3
        CALL WMFINL

	end subroutine ww3_finalize_internal

!***********************************************************
!
!**************************************************************
	subroutine getvarww3(lfirst) 

	use pkonst

	implicit none

! common
! local
	logical, intent(in) :: lfirst 
	integer k,l, ie
	integer ivar
	integer it
	logical next_output_d
	integer id, nvar, ierr
	double precision dtime
	double precision daux
	real wfact,wparam,wsmin,cdw
	real z0alpha

!------------------------------------------------------

	if (.not. bww3) return

!       -----------------------------------------------
!       Same time step, do read
!       -----------------------------------------------

	call get_act_dtime(dtime)

	it = dtime

!-------------------------------------------------------------
! set coupling time step 
!-------------------------------------------------------------

	call convert_time_d('dtwave',daux)
	idcoup = nint(daux)

	if (mod(it,idcoup) .eq. 0 ) then

!         -----------------------------------------------
!         compute stress and wave characteristics
!         -----------------------------------------------
	  call w3cplrdstr2d
	  call w3cplparam

!         -----------------------------------------------
!         compute global arrays from ww3
!         -----------------------------------------------
!         radiation stress ...
  	  call get_global_array_ww3(sxxww3,sxxgl)
	  call get_global_array_ww3(sxyww3,sxygl)
	  call get_global_array_ww3(syyww3,syygl)

!	  integral wave parameter ... 
	  do ivar = 1, nvarsww3
 	     call get_global_array_ww3(outvarww3(ivar,:),outvargl(ivar,:)) 
	  enddo 

!         other variables ... 
	  call get_global_array_ww3(charnww3,charngl) 
	  call get_global_array_ww3(ustww3,ustgl) 
	  call get_global_array_ww3(tauwwxww3,tauwwxgl) 
	  call get_global_array_ww3(tauwwyww3,tauwwygl) 
	  call get_global_array_ww3(ussxww3,ussxgl) 
	  call get_global_array_ww3(ussyww3,ussygl) 
	  call get_global_array_ww3(tusxww3,tusxgl) 
	  call get_global_array_ww3(tusyww3,tusygl) 
	  call get_global_array_ww3(bhdww3,bhdgl) 

!         -----------------------------------------------
!         compute local arrays for shyfem
!         -----------------------------------------------
!         radiation stress ...
	  sxx3dlshyfem = 0.
	  sxy3dlshyfem = 0.
	  syy3dlshyfem = 0.
  	  call fill_local_array_shyfem(sxx3dlshyfem(1,:),sxxgl(:))
   	  call fill_local_array_shyfem(sxy3dlshyfem(1,:),sxygl(:))
  	  call fill_local_array_shyfem(syy3dlshyfem(1,:),syygl(:))
  	  call mpi_barrier(mpicomm,ierr) 

!         integral wave parameter ... 
	  do ivar = 1, nvarsww3
	    call fill_local_array_shyfem(outvarshyfem(ivar,:),outvargl(ivar,:))
	  enddo
  	  call mpi_barrier(mpicomm,ierr) 
	  waveh = outvarshyfem(1,:)
	  wavep = outvarshyfem(2,:)
	  wavepp = 1./outvarshyfem(3,:)         !convert frequency to period
	  waved = outvarshyfem(4,:)

!         other variables ... 
	  call fill_local_array_shyfem(charnshyfem,charngl) 
	  call fill_local_array_shyfem(ustshyfem,ustgl) 
          call fill_local_array_shyfem(tauwwxshyfem,tauwwxgl)
          call fill_local_array_shyfem(tauwwyshyfem,tauwwygl)
	  call fill_local_array_shyfem(ussxshyfem,ussxgl) 
	  call fill_local_array_shyfem(ussyshyfem,ussygl) 
	  call fill_local_array_shyfem(tusxshyfem,tusxgl) 
	  call fill_local_array_shyfem(tusyshyfem,tusygl) 
	  call fill_local_array_shyfem(bhdshyfem,bhdgl) 
	  call mpi_barrier(mpicomm,ierr) 

!         -----------------------------------------------
!         Compute surface roughness z0s = z0alpha*Hs
!           with z0alpha a calibration parameter varying from 0.5 to 1.3 
!           TOBEDONE: set z0alpha as a str parameter
!         -----------------------------------------------
          z0alpha = 0.5
	  z0sk = z0alpha * waveh

!         -----------------------------------------------
!	  Set Charnock parameter
!         -----------------------------------------------
          charn  = charnshyfem

!      	  -----------------------------------------------
!         Compute residual momentum flux from atm to ocean
!           substract wave-supported stress from wind stress
!      	  -----------------------------------------------
          if ( iatm /= 1 ) then		!Wind stress from WW3
  	    wfact = 1. / rowass		!need to be recomputed here
  	    do k = 1,nkn
              wsmin = max(0.01,metws(k))
              cdw = (ustshyfem(k)/wsmin)**2
              windcd(k) = max(0.001, cdw)
  	      wparam = (1.-metice(k))*wfact*windcd(k)*metws(k)
              tauxnv(k) = wparam*wxv(k)
              tauynv(k) = wparam*wyv(k)
            end do
            tauxnv = tauxnv - tauwwxshyfem*(1.-metice)
            tauynv = tauynv - tauwwyshyfem*(1.-metice)
          else				!Wind stress from WRF
            if (icall_nuopc == 1) then  !partition it only at first call
  	      tauxnv = tauxnv - tauwwxshyfem*(1.-metice)
              tauynv = tauynv - tauwwyshyfem*(1.-metice)
            end if
          end if

!         -----------------------------------------------
!         Compute wave induced forces
!         -----------------------------------------------
  	  if (iwave .eq. 2) then
!         	-----------------------------------------------
!         	Radiation stress formulation
!         	-----------------------------------------------
	  	call diffxy(sxx3dlshyfem,sxy3dlshyfem,syy3dlshyfem &
			,wavefx,wavefy)

  	  elseif (iwave .eq. 3) then
!         	-----------------------------------------------
!         	Vortex force formulation 
!		3D formulation TO BE DONE
!         	-----------------------------------------------
		call wave_vortex(ussxshyfem,ussyshyfem,bhdshyfem &
			,wavefx,wavefy)

  	  endif ! iwave

	end if ! dtcoup

!     -----------------------------------------------
!       Writes output to the file.wav 
!     -----------------------------------------------
	if( next_output_d(da_wav) ) then
		id = nint(da_wav(4))
		call get_act_dtime(dtime)
		call shy_write_scalar_record2d(id,dtime,231,waveh)
		call shy_write_scalar_record2d(id,dtime,232,wavep)
		call shy_write_scalar_record2d(id,dtime,233,waved)
		call shy_sync(id)
	end if

	end subroutine getvarww3
!***********************************************************
!
!***********************************************************
	subroutine getvarshyfem(lfirst) 

	use pkonst

	implicit none
! common
! local
	logical, intent(in) :: lfirst
	integer k,l, ie
	integer it
	double precision dtime
	double precision daux

!------------------------------------------------------

	if (.not. bww3) return
        fllev = .true.

!       -----------------------------------------------
!       Same time step, do read
!       -----------------------------------------------

	call get_act_dtime(dtime)

	it = dtime

!-------------------------------------------------------------
! set coupling time step 
!-------------------------------------------------------------

	call convert_time_d('dtwave',daux)
	idcoup = nint(daux)

	if (mod(it,idcoup) .eq. 0 ) then
!
!         	-----------------------------------------------
!	        compute global arrays from shyfem
!       	-----------------------------------------------
!         	wind 
		if (.not. lfirst) then 
			wx0 = wxn
			wy0 = wyn 
		endif 
		call get_global_array_shyfem(wxv,wxvgl)
		call get_global_array_shyfem(wyv,wyvgl)

		if (lfirst) then 
			wx0(:,1) = wxvgl 
			wy0(:,1) = wyvgl
			wxn(:,1) = wxvgl
			wyn(:,1) = wyvgl
		else
			wxn(:,1) = wxvgl
			wyn(:,1) = wyvgl
		endif 

		if (myrank == 0) then 
			write(*,*) 'currxy flag', flcur
			write(*,*) 'water level flag', fllev
		endif

!         curr 
		do k = 1,nkn
			call getuv(1,k,currxlshyfem(k),currylshyfem(k))
		enddo 

		if (flcur) then 
			if (.not. lfirst) then    
				cx0 = cxn
				cy0 = cyn    
			endif 

			call get_global_array_shyfem(currxlshyfem, currxgl)
			call get_global_array_shyfem(currylshyfem, currygl)

			if (lfirst) then 
				cx0(:,1) = currxgl
				cy0(:,1) = currygl
				cxn(:,1) = currxgl
				cyn(:,1) = currygl
			else
				cxn(:,1) = currxgl
				cyn(:,1) = currygl
			endif 
		endif ! flcur

!         	water level
		if (fllev) then 
			wlvshyfem = znv 
			call get_global_array_shyfem(wlvshyfem,wlvgl)
			wlev(:,1) = wlvgl
		endif

	end if ! dtcoup

	end subroutine getvarshyfem

!***********************************************************
!
!***********************************************************
!differenzation of radiation stresses
!**************************************************************************
! Computes wave forcing terms according to the radiation stress formulation

	subroutine diffxy(SXX3D,SYY3D,SXY3D,wavefx,wavefy)

	use evgeom
	use levels
	use basin

	implicit none

	real, intent(in) :: SXX3D(:,:)       !radiation stress xx
	real, intent(in) :: SYY3D(:,:)       !radiation stress yy
	real, intent(in) :: SXY3D(:,:)       !radiation stress xy

	real,intent(out) ::  wavefx(nlv,nel)      !wave forcing term x
	real,intent(out) ::  wavefy(nlv,nel)      !wave forcing term y

	double precision :: b,c           !x and y derivated form function [1/m]
	integer          :: k,ie,ii,l,ilevel
	real             :: radsx,radsy

	wavefx = 0.
	wavefy = 0.
	do ie = 1,nel
		ilevel = ilhv(ie)
		do l=1,ilevel
			radsx = 0.
			radsy = 0.
			do ii = 1,3
				k = nen3v(ii,ie)
				b = ev(3+ii,ie)
				c = ev(6+ii,ie)
				radsx = radsx -(SXX3D(l,k)*b + SXY3D(l,k)*c)
				radsy = radsy -(SXY3D(l,k)*b + SYY3D(l,k)*c)
			end do
			wavefx(l,ie) = -radsx
			wavefy(l,ie) = -radsy
		enddo
	enddo

	end subroutine diffxy

!**************************************************************************
! Computes wave forcing terms according to the vortex force formulation

        subroutine wave_vortex(stokesx,stokesy,wavejb,wavefx,wavefy)

        use mod_internal
        use mod_depth
        use mod_layer_thickness
        use mod_hydro_vel
        use evgeom
        use levels
        use basin
        use pkonst

        implicit none

! arguments
        real stokesx(nlv,nkn)           !x stokes velocity
        real stokesy(nlv,nkn)           !y stokes velocity
        real wavejb(nkn)                !wave pressure
        real wavefx(nlv,nel)            !x wave forcing term
        real wavefy(nlv,nel)            !y wave forcing term

! common

! local
        real, allocatable :: stokesz(:,:)       !z stokes velocity on node k
        real, allocatable :: hk(:)              !layer tickness on nodes
        real, allocatable :: stxe(:,:)          !x stokes transport on elements
        real, allocatable :: stye(:,:)          !y stokes transport on elements
        real, allocatable :: saux1(:,:)
        real, allocatable :: saux2(:,:)
        real, allocatable :: stokesze(:,:)      !z stokes velocity on elements
        real, allocatable :: uaux(:)
        real, allocatable :: vaux(:)
        integer k,ie,ii,l,ilevel
        real f                          !Coriolis parameter on elements
        real h                          !layer thickness
        real u,v                        !velocities at level l and elements
        real stxk, styk                 !stokes transport on nodes
        real auxx, auxy                 !auxiliary variables
        real jbk                        !integrated wave perssure term
        double precision b,c            !x and y derivated form function [1/m]
        real wavesx,wavesy
        real wuz,wvz                    !z vortex force
        real sz,sz1

!       -----------------------------------------------
!       Initialization
!       -----------------------------------------------
        wavefx = 0.d0
        wavefy = 0.d0

        allocate(stokesz(nlv,nkn))
        allocate(hk(nlv))
        allocate(stxe(nlv,nel))
        allocate(stye(nlv,nel))
        allocate(saux1(nlv,nkn))
        allocate(saux2(nlv,nkn))
        allocate(stokesze(0:nlv,nel))
        allocate(uaux(0:nlv+1))
        allocate(vaux(0:nlv+1))

        stokesze = 0.
        stxe = 0.
        stye = 0.

!       -----------------------------------------------
!       Computes wave forcing terms due to horizontal stokes
!       velocities and wave pressure head
!       -----------------------------------------------
        do ie = 1,nel
          ilevel = ilhv(ie)
          f = fcorv(ie)
          do l = 1,ilevel
            wavesx = 0.
            wavesy = 0.
            auxx = 0.
            auxy = 0.
            h = hdenv(l,ie)
            u = ulnv(l,ie)
            v = vlnv(l,ie)
            do ii = 1,3
              k = nen3v(ii,ie)
              h = hdknv(l,k)
              b = ev(3+ii,ie)
              c = ev(6+ii,ie)
              stxk = stokesx(l,k) * h
              styk = stokesy(l,k) * h
              auxx = auxx + stxk
              auxy = auxy + styk
              jbk = wavejb(k) * h * grav / 3.   !???? is it correct to divide by 3?

              wavesx = wavesx - (u*b*styk - v*c*styk) + b*jbk
              wavesy = wavesy + (u*b*stxk - v*c*stxk) + c*jbk

              !Dutour Sikiric
              !wavesx = wavesx - (u*b*stxk + v*b*styk) + b*jbk
              !wavesy = wavesy - (u*c*stxk + v*c*styk) + c*jbk

            end do
            stxe(l,ie) = auxx / 3.
            stye(l,ie) = auxy / 3.

            wavefx(l,ie) = wavesx - f*stye(l,ie)
            wavefy(l,ie) = wavesy + f*stxe(l,ie)
          end do
        end do
!       -----------------------------------------------
!       Check for nan
!       -----------------------------------------------
        call nantest(nel*nlv,wavefx,'WAVEFX')
        call nantest(nel*nlv,wavefy,'WAVEFy')

!       -----------------------------------------------
!       Computes vertical stokes velocity
!       -----------------------------------------------
        call stokes_vv(saux1,saux2,stxe,stye,stokesz,stokesze)

!       -----------------------------------------------
!       Computes wave forcing terms due to vertical stokes
!       velocity
!       -----------------------------------------------
        do ie = 1,nel
          ilevel = ilhv(ie)

          do l = 1,ilevel
            uaux(l) = ulnv(l,ie)
            vaux(l) = vlnv(l,ie)
          end do
          uaux(0) = 0.
          vaux(0) = 0.
          uaux(ilevel+1) = 0.
          vaux(ilevel+1) = 0.

          do l = 1,ilevel
            sz = stokesze(l,ie)
            sz1 = stokesze(l-1,ie)
            wuz = 0.
            wvz = 0.

            if ( sz1 .lt. 0. ) then
              wuz = wuz - sz1*uaux(l-1)
              wvz = wvz - sz1*vaux(l-1)
            else
              wuz = wuz - sz1*uaux(l)
              wvz = wvz - sz1*vaux(l)
            end if

            if ( sz .gt. 0. ) then
              wuz = wuz + sz*uaux(l+1)
              wvz = wvz + sz*vaux(l+1)
            else
              wuz = wuz + sz*uaux(l)
              wvz = wvz + sz*vaux(l)
            end if

            wavefx(l,ie) = wavefx(l,ie) + wuz           !check stokesz first
            wavefy(l,ie) = wavefy(l,ie) + wvz
          end do
        enddo

        end subroutine wave_vortex

!**************************************************************************
! Computes vertical stokes velocity

        subroutine stokes_vv(vf,va,stxe,stye,auxstz,stokesze)

        use evgeom
        use levels
        use basin

        implicit none

! parameters
! arguments
        real vf(nlv,nkn)                !auxiliary array
        real va(nlv,nkn)                !auxiliary array
        real stxe(nlv,nel)      	!x stokes transport on elements
        real stye(nlv,nel)      	!y stokes transport on elements
        real auxstz(nlv,nkn)    	!z stokes velocity on node k for plot
        real stokesze(0:nlv,nel)        !z stokes velocity on elements

! local
        real, allocatable :: stokesz(:,:)       !z stokes velocity on node k
        logical debug
        integer k,ie,ii,kk,l,lmax
        integer ilevel
        double precision b,c            !x and y derivated form function [1/m]
        real aj,ff,atop,acu
        logical is_zeta_bound,is_boundary_node

! initialize
        allocate(stokesz(0:nlv,nkn))
        vf = 0.
        va = 0.
        stokesz = 0.
        auxstz = 0.

! compute difference of velocities for each layer
        do ie=1,nel
          aj=4.*ev(10,ie)               !area of triangle / 3
          ilevel = ilhv(ie)
          do l=1,ilevel
            do ii=1,3
               kk=nen3v(ii,ie)
               b = ev(ii+3,ie)
               c = ev(ii+6,ie)
               ff = stxe(l,ie)*b + stye(l,ie)*c
               vf(l,kk) = vf(l,kk) + 3. * aj * ff
               va(l,kk) = va(l,kk) + aj
            end do
          end do
        end do
! from vel difference get absolute velocity (w_bottom = 0)
!       -> stokesz(nlv,k) is already in place !
!       -> stokesz(nlv,k) = 0 + stokesz(nlv,k)
! w of bottom of last layer must be 0 ! -> shift everything up
! stokesz(nlv,k) is always 0
!
! dividing stokesz [m**3/s] by area [vv] gives vertical velocity
!
! in vv(l,k) is the area of the upper interface: a(l) = a_i(l-1)
! =>  w(l-1) = flux(l-1) / a_i(l-1)  =>  w(l-1) = flux(l-1) / a(l)

        do k=1,nkn
          lmax = ilhkv(k)
          stokesz(lmax,k) = 0.
          debug = k .eq. 0
          do l=lmax,1,-1
            stokesz(l-1,k) = stokesz(l,k) + vf(l,k)
          end do
          stokesz(0,k) = 0.        ! ensure no flux across surface - is very small
        end do

        do k=1,nkn
          lmax = ilhkv(k)
          debug = k .eq. 0
          do l=2,lmax
            atop = va(l,k)
            if( atop .gt. 0. ) then
              stokesz(l-1,k) = stokesz(l-1,k) / atop
              if( debug ) write(6,*) k,l,atop,stokesz(l-1,k)
            end if
            auxstz(l,k) = stokesz(l,k)
          end do
          auxstz(lmax,k) = stokesz(lmax,k)
        end do

! set w to zero at open boundary nodes (new 14.08.1998)
        do k=1,nkn
            !if( is_zeta_bound(k) ) then
            if( is_boundary_node(k) ) then
              do l=0,nlv
                stokesz(l,k) = 0.
                auxstz(l,k) = 0.
              end do
            end if
        end do

!-----------------------------------------------------------
! convert values to elelemts
!-----------------------------------------------------------
        do ie=1,nel
          lmax = ilhv(ie)
          do l = 1,lmax
            acu = 0.
            do ii=1,3
              k = nen3v(ii,ie)
              acu = acu + stokesz(l,k)
            end do
            stokesze(l,ie) = acu / 3.
          end do
          stokesze(0,ie) = 0.
        end do

        return

        end subroutine stokes_vv

!**************************************************************
!
!**************************************************************
	subroutine get_global_array_ww3(localarray,globalarray)

	implicit none

	real, intent(in)  :: localarray(:)
	real, intent(out) :: globalarray(:)
	real, allocatable :: tmp(:)
	integer           :: ip, ierr

	allocate(tmp(NP_GLOBAL)); tmp = 0.
	globalarray = 0. 

	if (size(globalarray).ne. NP_GLOBAL) then
		write(*,*) size(globalarray), NP_GLOBAL
		stop 'error in grid sizes get_global_array_ww3' 
	endif

	DO IP = 1, NSEAL 
		tmp(iplg(IP)) = localarray(IP)
	END DO

	call mpi_allreduce(tmp,globalarray,NP_GLOBAL, &
                        MPI_REAL,MPI_SUM,mpiComm,ierr)

!         if(myrank==0) then
		do IP = 1, NP_GLOBAL
			globalarray(IP) = globalarray(IP)*nwild_gb(IP)
		enddo !IP
!         endif !myrank

	end subroutine get_global_array_ww3
!**************************************************************
!
!**************************************************************
	subroutine get_global_array_shyfem(localarray,globalarray)

	implicit none

	real, intent(in)  :: localarray(:)
	real, intent(out) :: globalarray(:)

	integer           :: ip, ierr

	if (size(globalarray).ne. nkn_GLOBAL) then
		write(*,*) size(globalarray), nkn_GLOBAL
		stop 'error in grid sizes get_global_array_shyfem'
	endif

	call shympi_l2g_array(localarray,globalarray)

	end subroutine get_global_array_shyfem
!**************************************************************
!
!**************************************************************
	subroutine fill_local_array_ww3(localarray,globalarray)

	implicit none

	real, intent(out)  :: localarray(:)
	real, intent(in)   :: globalarray(:)
	real, allocatable :: tmp(:)

	integer           :: ip, ierr

	allocate(tmp(NP_GLOBAL)); tmp = 0.
	localarray = 0.

	if (size(globalarray).ne. NP_GLOBAL) then
		write(*,*) size(globalarray), NP_GLOBAL
		stop 'error in grid sizes fill_local_array_ww3'
	endif

	DO IP = 1, NP_GLOBAL
		IF (ipgl(ip) .gt. 0) THEN
			localarray(ipgl(ip)) = globalarray(IP)
		END IF
	ENDDO

	end subroutine fill_local_array_ww3
!**************************************************************
!
!**************************************************************
	subroutine fill_local_array_shyfem(localarray,globalarray)

	implicit none

	real, intent(out) :: localarray(:)
	real, intent(in)  :: globalarray(:)

	integer           :: ip, ierr

	localarray = 0.

	if (size(globalarray).ne. nkn_GLOBAL) then
		write(*,*) size(globalarray), nkn_GLOBAL
		stop 'error in grid sizes fill_local_array_shyfem'
	endif

	call shympi_g2l_array(globalarray,localarray)

	end subroutine fill_local_array_shyfem
!**************************************************************************
	end module
!**************************************************************************

        subroutine ww3_init
        use mod_ww3
        call ww3_init_internal
        end subroutine ww3_init

        subroutine ww3_loop
        use mod_ww3
        call ww3_loop_internal
        end subroutine ww3_loop

        subroutine ww3_finalize
        use mod_ww3
        call ww3_finalize_internal
        end subroutine ww3_finalize

!------------------------------------------------------------------------
