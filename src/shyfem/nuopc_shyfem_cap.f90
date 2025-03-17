!==============================================================================
! Earth System Modeling Framework
! Copyright (c) 2002-2023, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

	!! \section{The ocean component: SHYFEM}
	!!
	!! We are cascading into the objects hierarchy. We miss the last step: we need
	!! to specify the methods for the children object. Here we do this for the 
        !! ocean model object.  The ocean component also has preprocessors for fine tuning.
	!! We have added the following preprocessor to selectively decouple the 
        !! ocean component from the rest of the earth system. Disabling all 
	!! the following macros, will result in an ocean component that does not 
	!! advertise any importable/exportable Field. You should use this only if you
	!! want to drive the model independently.
        !! You could also couples the ocean component only through
	!! momentum or heat flux or sst only. These options may be useful for
	!! debugging and testing for example. In general, all the macros must be defined.
        !! Switching on/off the different fields is possible through the each model
        !! configuration file.
#define WITHFIELDS_MOMENTUMFLUX
#define WITHFIELDS_HEATFLUX
#define WITHFIELDS_MASSFLUX
#define WITHFIELDS_SST
#define WITHFIELDS_SEAROUGHNESS

	!! The ocean component is coded into a Fortran module.
        module ocean_shyfem

	!! We call the modules of the libraries. First the SHYFEM library:
	use mod_shyfem
	!! and of course ESMF and NUOPC:
	use ESMF
        use NUOPC
	use NUOPC_Model, &
     &	  modelSS => SetServices

	implicit none

	!! We declare private members and methods. \textsf{SHYFEM\_FieldMetadata}
	!! is an object that was not present in ESMF. We have created it in the cap
	!! layer and for this reason it has the SHYFEM prefix. This object contains
	!! important information to create the ESMF fields and will help us
	!! in factorizing field creation within do loops, instead of doing evreything
	!! manually. The info concerns the field names and mesh staggering
	!! and are initialized here with default values.
	private

	type SHYFEM_Metadata
	  character(4)                :: fieldName = " "
          character(50)               :: longName = " "
	  type(ESMF_MeshLoc)          :: meshloc = ESMF_MESHLOC_NODE
	end type

        type SHYFEM_Mesh
          integer                         :: nkn_ghost
          integer, allocatable            :: id_node_ghost(:)
          integer, allocatable            :: ipv_ghost(:)
          real(ESMF_KIND_R8), allocatable :: xgv_ghost(:)
          real(ESMF_KIND_R8), allocatable :: ygv_ghost(:)
          integer, allocatable            :: nen3v_ghost(:,:)
	  integer, allocatable            :: table_ghostToLocal(:)
	  integer                         :: localPet
          integer                         :: petCount
        end type

	integer :: SHYFEM_numOfImportFields
        integer :: SHYFEM_numOfExportFields

	type(SHYFEM_Metadata), allocatable, &
     &     save :: SHYFEM_FieldMetadata(:)

	type(SHYFEM_Mesh), save :: ShyfemToEsmf_Mesh

	public SetServices

	!-----------------------------------------------------------------------------
	contains
	!-----------------------------------------------------------------------------

        !! We have already seen how to create a component with the aid of the NUOPC
        !! layer.
	subroutine SetServices(model, rc)
	  type(ESMF_GridComp)  :: model
	  integer, intent(out) :: rc
        
	  rc = ESMF_SUCCESS

	  !! We do not comment again the derive/specialize commands, being
	  !! exactly equal to what we have seen for the earth component.

	  ! derive from NUOPC_Model
	  call NUOPC_CompDerive(model, modelSS, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &	    line=__LINE__,  &
     &	    file=__FILE__))  &
     &	    return  ! bail out

	  !! Instead we note that different methods are registered each one with a 
	  !! specific \textsf{specLabel} that denote a specific operation. The labels
	  !! can be regrouped into three phases: 
	  !! \begin{itemize}
  	  !! \item initialization: advertise and realize
	  !! \item run: advance
	  !! \item finalization: finalize
	  !! \end{itemize}

	  ! specialize model
	  call NUOPC_CompSpecialize(model, specLabel=label_Advertise, &
     &	    specRoutine=Advertise, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &      line=__LINE__, &
     &      file=__FILE__)) &
     &      return  ! bail out
	  call NUOPC_CompSpecialize(model,specLabel=label_RealizeProvided, &
     &      specRoutine=Realize, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &	    line=__LINE__, &
     &	    file=__FILE__)) &
     &	    return  ! bail out
	  call NUOPC_CompSpecialize(model, specLabel=label_Advance, &
     &	    specRoutine=Advance, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &	    line=__LINE__, &
     &	    file=__FILE__)) &
     &	    return  ! bail out
	  call NUOPC_CompSpecialize(model, specLabel=label_Finalize, &
     &	    specRoutine=Finalize, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &	    line=__LINE__, &
     &	    file=__FILE__)) &
     &	    return  ! bail out

	end subroutine

	!-----------------------------------------------------------------------------

	!! As done for the earth system component, the following subroutines register
	!! the members and methods for the the ocean model object that enters as first
	!! argument.
        !!
	!! \subsection{Advertise}
        !!
	!! The purpose of \textsf{Advertise()} is for your model to advertise its 
	!! import and export fields. This means that your model announces which model 
	!! variables it is capable of exporting (e.g., an ocean might export water 
	!! state at sea level) and which model variables it requires (e.g., an ocean 
	!! might require mass, momentum and heat fluxes as a boundary condition). The 
	!! reason there is an explicit advertise phase is because NUOPC dynamically 
	!! matches fields among all the models participating in a coupled simulation 
	!! during runtime. So, we need to collect the list of possible input and output 
	!! fields from all the models during their initialization. 
	!!  
	!! Advertising a Field does NOT allocate memory. Note that NUOPC does not 
	!! allocate memory for fields during the advertise phase or when 
	!! \textsf{NUOPC\_Advertise} is called. Instead, this is simply a way for 
	!! models to communicate the standard names of fields. During a later phase, 
	!! only those fields that are connected (e.g., a field exported from one model 
	!! that is imported by another) need to have memory allocated. Also, since ESMF
	!! will accept pointers to pre-allocated memory, it is usually not necessary to 
	!! change how memory is allocated for your model's variables. 
	subroutine Advertise(model, rc)

	  type(ESMF_GridComp)  :: model
	  integer, intent(out) :: rc

	  ! local variables
	  type(ESMF_State)        :: importState, exportState
	  integer                 :: var, num

	  rc = ESMF_SUCCESS

	  !! Since we are in the initialization phase, we use this subroutine 
	  !! to call also the SHYFEM initialization. Here the grid and data are read,
	  !! the initial conditions are set and data structures are
	  !! initialized. The \textsf{.false.} argument forces SHYFEM to not intialize MPI,
	  !! since this is left to the the coupler.
	  call shyfem_initialize(.false.)
	  call ESMF_LogWrite("Initialized OCN", ESMF_LOGMSG_INFO, rc=rc)

	  ! query for importState and exportState
	  call NUOPC_ModelGet(model, importState=importState, &
     &	    exportState=exportState, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &	    line=__LINE__, &
     &	    file=__FILE__)) &
     &	    return  ! bail out

	  !! The fields are organized with a precise order in a \textsf{SHYFEM\_FieldMetadata}
	  !! array. This ordering is beneficial to factorize the code and perform the 
	  !! many operations of the fields with do loops. You have to keep in mind that, 
	  !! first we have import fields, and only after we have export fields. 
	  allocate( SHYFEM_FieldMetadata(11) )

          !! The field list is defined here. You can modify this part adding new fields.
	  !! The metadata of each field is filled with a simple constructor statement.
	  !! The constructor is feeded with the field names which should correspond to the
	  !! the standard ones, for compatibility with the other components, see section 2.2.2 in:
	  !!
	  !! \texttt{https://earthsystemmodeling.org/docs/nightly/develop/NUOPC\_refdoc/node3.html}
	  !!
	  !! In the following we describe each field rapidly.
          num = 0

#ifdef WITHFIELDS_MOMENTUMFLUX
          !! \begin{itemize}
          !! \item atmospheric pressure at sea-level $p_a$ that
          !! act as a forcing term on the momentum equation.
          num = num + 1	  
          SHYFEM_FieldMetadata(1) = SHYFEM_Metadata(fieldName="pmsl", &
     &      longName="air_pressure_at_sea_level", &
     &      meshloc=ESMF_MESHLOC_NODE)

          !! \item ocean-atmosphere flux $F_{oa}(U^{n}_o,U^{n}_a)$ 
          !! that describes the flux across the surface for momentum (wind stress) 
          !! and temperature (heat flux). These fluxes are typically computed by the
          !! the atmospheric component with the aid of bulk formulae.
          !! For now we have added only momentum flux that has two components.
          num = num + 1	  
          SHYFEM_FieldMetadata(2) = SHYFEM_Metadata(fieldName="smes", &
     &      longName="surface_downward_eastward_stress", &
     &      meshloc=ESMF_MESHLOC_NODE)

          num = num + 1	  
          SHYFEM_FieldMetadata(3) = SHYFEM_Metadata(fieldName="smns", &
     &      longName="surface_downward_northward_stress", &
     &      meshloc=ESMF_MESHLOC_NODE)

#endif
#ifdef WITHFIELDS_HEATFLUX
	  !! \item downward-short-wave-radiation flux $R_{sw}$.
          num = num + 1	  
          SHYFEM_FieldMetadata(num) = SHYFEM_Metadata( &
     &      fieldName="rsns", &
     &      longName="surface_net_downward_shortwave_flux", &
     &      meshloc=ESMF_MESHLOC_NODE)

	  !! \item downward-long-wave-radiation flux $R_{lw}$.
          num = num + 1
          SHYFEM_FieldMetadata(num) = SHYFEM_Metadata( &
     &      fieldName="rlns", &
     &      longName="surface_net_downward_longwave_flux", &
     &      meshloc=ESMF_MESHLOC_NODE)

          !! \item sensible heat flux $Q_{sens}$.
          num = num + 1
          SHYFEM_FieldMetadata(num) = SHYFEM_Metadata( &
     &      fieldName="stsh", &
     &      longName="surface_downward_sensible_heat_flux_in_air", &
     &      meshloc=ESMF_MESHLOC_NODE)

          !! \item sensible heat flux $Q_{lat}$.
          num = num + 1
          SHYFEM_FieldMetadata(num) = SHYFEM_Metadata( &
     &      fieldName="stlh", &
     &      longName="surface_downward_latent_heat_flux_in_air", &
     &      meshloc=ESMF_MESHLOC_NODE)
#endif
#ifdef WITHFIELDS_MASSFLUX
          !! \item evaportaion/precipitation rate $P-E$.
          num = num + 1
	  SHYFEM_FieldMetadata(num) = SHYFEM_Metadata( &
     &      fieldName="prec", &
     &      longName="precipitation_flux", &
     &      meshloc=ESMF_MESHLOC_NODE)
#endif
          SHYFEM_numOfImportFields = num
#ifdef WITHFIELDS_SST
          num = num + 1
	  SHYFEM_FieldMetadata(num) = SHYFEM_Metadata( &
     &      fieldName="umap", &
     &      longName="atmospheric_unmapped_points", &
     &      meshloc=ESMF_MESHLOC_NODE)

	  !! \item sea surface temperature.
          num = num + 1
	  SHYFEM_FieldMetadata(num) = SHYFEM_Metadata( &
     &      fieldName="sst", &
     &      longName="sea_surface_temperature", &
     &      meshloc=ESMF_MESHLOC_NODE)
#endif
#ifdef WITHFIELDS_SEAROUGHNESS
	  !! \item sea roughness length
          num = num + 1
          SHYFEM_FieldMetadata(num) = SHYFEM_Metadata( &
     &      fieldName="z0", &
     &      longName="roughness_length", &
     &      meshloc=ESMF_MESHLOC_NODE)
#endif
          SHYFEM_numOfExportFields = num - SHYFEM_numOfImportFields

          !! \end{itemize}

	  !! With a do loop we advertise the import state in the 
	  !! ocean component.
	  do var=1,SHYFEM_numOfImportFields

	    call NUOPC_Advertise(importState, &
     &	      StandardName=SHYFEM_FieldMetadata(var)%longName, &
     &        name=SHYFEM_FieldMetadata(var)%fieldName, rc=rc) 
	    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &	      line=__LINE__, &
     &	      file=__FILE__)) &
     &	      return  ! bail out

	  enddo

	  !! As output the ocean state at the surface layer must be 
	  !! exorted. For now we have exported only the sea surface temperature.
          do var=SHYFEM_numOfImportFields+1, &
     &         SHYFEM_numOfImportFields+SHYFEM_numOfExportFields

	    call NUOPC_Advertise(exportState, &
     &        StandardName=SHYFEM_FieldMetadata(var)%longName, &
     &        name=SHYFEM_FieldMetadata(var)%fieldName, rc=rc)
	    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &	      line=__LINE__, &
     &	      file=__FILE__)) &
     &	      return  ! bail out

	  enddo

	end subroutine

	!-----------------------------------------------------------------------------

        !!
        !! \subsection{Realize}
        !!
	!! The following code fragment shows the \textsf{Realize} subroutine. During 
	!! this phase, fields that were previously advertised should now be realized. 
	!! Realizing a field means that an \textsf{ESMF\_Field} object is created and 
	!! it is added to the appropriate \textsf{ESMF\_State}, either import or export.
	!! To create a field it is necessary to import the mesh from SHYFEM.
	subroutine Realize(model, rc)
	  type(ESMF_GridComp)  :: model
	  integer, intent(out) :: rc

	  ! local variables
	  type(ESMF_State)        :: importState, exportState
	  type(ESMF_TimeInterval) :: stabilityTimeStep
	  type(ESMF_Field)        :: field
          type(ESMF_Mesh)         :: meshIn
	  double precision, pointer :: fieldPtr(:)
	  integer                 :: var
    	  type(ESMF_VM) :: vm


	  rc = ESMF_SUCCESS

	  ! query for importState and exportState
	  call NUOPC_ModelGet(model, importState=importState, &
     &      exportState=exportState, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &	    line=__LINE__, &
     &	    file=__FILE__)) &
     &	    return  ! bail out

	  !! The next call returns information about the ocean model,
	  !! depending on the argument. We ask for the local PET and for
	  !! the total number of PETs. These are stored in global variables 
	  !! since they stay constant along the simulations. For now
	  !! these info are used only to append to the name of the output file 
	  !! info about the PET.
	  
!          call ESMF_GridCompGet(model,
!     +      localPet=ShyfemToEsmf_Mesh%localPet,
!     +      petCount=ShyfemToEsmf_Mesh%petCount, rc=rc)
!          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
!     +      line=__LINE__,
!     +      file=__FILE__))
!     +      return  ! bail out

	  call ESMF_VMGetCurrent(vm, rc=rc)

	  call ESMF_VMGet(vm, localPet=ShyfemToEsmf_Mesh%localPet, &
     &      petCount=ShyfemToEsmf_Mesh%petCount, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &      line=__LINE__, &
     &      file=__FILE__)) &
     &      return  ! bail out

	  !! In order to create an \textsf{ESMF\_Field}, you'll first need to create one 
	  !! of the ESMF geometric types, \textsf{ESMF\_Grid}, \textsf{ESMF\_Mesh}, or 
	  !! \textsf{ESMF\_LocStream}. For 2D and 3D logically rectangular grids (such as
	  !! a lat-lon grid), the typical choice is \textsf{ESMF\_Grid}. For unstructured
	  !! grids, use an \textsf{ESMF\_Mesh}. 
	  !! SHYFEM has already read the mesh during the advertise phase. With the command
	  !! \textsf{SHYFEM\_MeshGet} the mesh is imported into a \textsf{ESMF\_Mesh} 
	  !! object and later we also write into vtk file to check that evreything is ok.
	  !! Vtk files can be opened with Paraview.
	  call SHYFEM_MeshGet(meshIn, rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &      line=__LINE__, &
     &      file=__FILE__)) &
     &      return
          call ESMF_LogWrite("  Get Mesh OCN", ESMF_LOGMSG_INFO, rc=rc)
          call ESMF_MeshWrite(meshIn, "shyfem_mesh", rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &      line=__LINE__, &
     &      file=__FILE__)) &
     &      return  ! bail out

          !! An ESMF Field represents a physical field, such as temperature. 
          !! The motivation for including Fields in ESMF is that bundles of Fields 
          !! are the entities that are normally exchanged when coupling Components.
          !! The ESMF Field class contains distributed and discretized field data, 
          !! a reference to its associated grid, and metadata. The Field class stores 
          !! the grid staggering for that physical field. This is the relationship 
          !! of how the data array of a field maps onto a grid (e.g. one item per cell 
          !! located at the cell center, one item per cell located at the NW corner, 
          !! one item per cell vertex, etc.). This means that different Fields which 
          !! are on the same underlying ESMF Grid but have different staggerings can 
          !! share the same Grid object without needing to replicate it multiple times.
          !! Fields can be added to States for use in inter-Component data 
          !! communications. Fields can also be added to FieldBundles, which are 
          !! groups of Fields on the same underlying Grid. One motivation for packing 
          !! Fields into FieldBundles is convenience; another is the ability to perform 
          !! optimized collective data transfers.
          !! Field communication capabilities include: data redistribution, regridding, 
          !! scatter, gather, sparse-matrix multiplication, and halo update. These are 
          !! discussed in more detail in the documentation for the specific method calls. 
          !! ESMF does not currently support vector fields, so the components of a 
          !! vector field must be stored as separate Field objects.
	  !!
	  !! With the following commands we create a field for the pressure. Fields
	  !! in ESMF are of type \textsf{ESMF\_field}.
	  !! The subroutine \textsf{ESMF\_FieldCreate} simply associates the data 
	  !! with the Grid. The keywords are very important. The keyword \textsf{staggerloc}
	  !! specifies how the data is attached to the grid. \textsf{typekind} tells about 
          !! the data type, in this case double precision. The name that have been advertised must be 
	  !! also appended.
	  !! 
          !! An \textsf{ESMF\_Field} is created by passing the field name 
          !! (should be the same as advertised), the grid, the staggering, and the data type of the 
          !! field to \textsf{ESMF\_FieldCreate}.
	  do var=1,SHYFEM_numOfImportFields

	    field = ESMF_FieldCreate(meshIn, &
     &	      ESMF_TYPEKIND_R8, &
     &        meshloc=SHYFEM_FieldMetadata(var)%meshloc,  &
     &        name=SHYFEM_FieldMetadata(var)%fieldName, &
     &        rc=rc)
	    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &	      line=__LINE__, &
     &	      file=__FILE__)) &
     &	      return  ! bail out

            !! Fields are put into import or export States by calling 
            !! \textsf{NUOPC\_Realize}.
	    call NUOPC_Realize(importState, field=field, rc=rc)
	    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &	      line=__LINE__, &
     &	      file=__FILE__)) &
     &	      return  ! bail out

	  enddo

          do var=SHYFEM_numOfImportFields+1, &
     &         SHYFEM_numOfImportFields+SHYFEM_numOfExportFields

	    field = ESMF_FieldCreate(meshIn, &
     &        ESMF_TYPEKIND_R8, &
     &        meshloc=SHYFEM_FieldMetadata(var)%meshloc, &
     &        name=SHYFEM_FieldMetadata(var)%fieldName, &
     &        rc=rc)
	    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &	      line=__LINE__, &
     &	      file=__FILE__)) &
     &	      return  ! bail out

	    call NUOPC_Realize(exportState, field=field, rc=rc)
	    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &	      line=__LINE__, &
     &        file=__FILE__)) &
     &        return  ! bail out

            !! The model has not been initialized. This is dangerous
	    !! because at the first timestep the atmospheric model would run with an
	    !! uninitialized ocean state.
            !! With the following commands we initialize the ocean component.
	    !! We retrieve the field pointer with \textsf{ESMF\_FieldGet}, that 
	    !! get a DE-local information (included the array pointer) 
	    !! from the field. Local DE is requested with \textsf{localDE} argument. In
	    !! our MPI parallelization, every process corresponds to one DE, which means
	    !! \textsf{localDeCount=1}. We have thus commented the loop
	    !! over the DE, being trivial. Then \textsf{localDe} argument may be also omitted, 
	    !! in which case it will default to localDe=0. Once we have
	    !! the field pointer we fill it in \textsf{SHYFEM\_FieldGet}.
	    !do j=0,localDECount-1   !uncomment in case of multiple DE
            call ESMF_FieldGet(field, localDe=0, farrayPtr=fieldPtr, &
     &        rc=rc)
	    if (ESMF_LogFoundError(rcToCheck=rc, &
     &        msg=ESMF_LOGERR_PASSTHRU, &
     &        line=__LINE__, &
     &        file=__FILE__)) &
     &        return  ! bail out

            call SHYFEM_FieldGet(fieldPtr, &
     &        SHYFEM_FieldMetadata(var)%fieldName, rc)
!     +        ShyfemToEsmf_Mesh, rc)
	    !enddo

	  enddo

	end subroutine

	!-----------------------------------------------------------------------------

        !!
        !! \subsection{Advance SHYFEM}
        !!
	!! The SHYFEM ocean model advances the shallow water multilayer equations for
	!! stratified flows of one time step:
	!! \[
	!! U^{n+1}_o = U^{n+1}_o + \Delta t_{ao}\left( L_o(U^{n+1}_o)
	!!           + F_{oa}(U^{n}_o,U^{n}_a) \right)
	!! \]
	!! where $U_o(x,t)$ is the ocean state, $\Delta t_{ao}$ is the ocean-atmosphere
	!! timestep. In multilayer models $U_o$ pile-up the ocean state of each layer
	!! as:
	!! \[
	!! U_o = \{\, \zeta,\,U_{o,1},\,...\,U_{o,\alpha},\,...U_{o,N} \,\}
	!! \]	
	!! with $\zeta(x,t)$  the free-surface, $U_{o,\alpha}$ the momentum, 
	!! temperature and salinity of layer $\alpha$. $N$ the number of layers.
	!! The operator $L_o$ is a finite element discretization of the shallow water 
	!! multilayer equations. $F_{oa}$ are the atmosphere-ocean fluxes.
	!! In this subrotine we do two things. We take the set the fluxes $F_{oa}$
	!! inside SHYFEM, through a quite complex chain of commands. Then we advance 
	!! SHYFEM.
	subroutine Advance(model, rc)
	  type(ESMF_GridComp)  :: model
	  integer, intent(out) :: rc

  	  !! In the ESMF vocabulary the import state are atmosphere-ocean
	  !! fluxes $F_{oa}$, the export state is the oceanic state at the
	  !! first layer $U_{o,1}$.
	  ! local variables
	  type(ESMF_Clock)            :: clock
    	  type(ESMF_State)            :: importState, exportState
	  type(ESMF_Time)             :: currTime
	  type(ESMF_TimeInterval)     :: timeStep
	  character(len=160)          :: msgString
          character(len=160)          :: dateString, petString
          double precision            :: timeStepSec

	  type(ESMF_Field)            :: field
	  double precision, pointer   :: fieldPtr(:)
!	  integer                     :: localPet, petCount
	  integer                     :: var
	  integer                     :: currYear, currMonth, currDay, &
     &                                   currHour, currMinute, &
     &                                   currSecond, currTimeSec
	  real, external              :: getpar	  
	  integer                     :: idtatm

	  integer, save               :: iuinfo = 0

#define NUOPC_TRACE__OFF
#ifdef NUOPC_TRACE
	  call ESMF_TraceRegionEnter("OCN:Advance")
#endif

	  rc = ESMF_SUCCESS

          !! The various objects are nested. We access the fields with three
          !! consecutive \textsf{Get} command. We start accessing the most
          !! exterior component, namely the model. We query for three
          !! objects: the \textsf{clock}, the \textsf{importState} 
          !! and \textsf{exportState}.
          call NUOPC_ModelGet(model, modelClock=clock, &
     &      importState=importState, exportState=exportState, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &      line=__LINE__, &
     &      file=__FILE__)) &
     &      return  ! bail out

          !! Before cascading from the model to the fields, we
	  !! access to the clock object. We do it in different ways: here we
          !! output the current time to the log file:
          call ESMF_ClockPrint(clock, options="currTime", &
     &	    preString="------>Advancing OCN from: ",unit=msgString,rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &      line=__LINE__, &
     &      file=__FILE__)) &
     &      return  ! bail out
          call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &      line=__LINE__, &
     &      file=__FILE__)) &
     &      return  ! bail out

          !! We retrieve two important members of the clock object,
          !! the \textsf{currTime} object and the \textsf{timeStep} object
          call ESMF_ClockGet(clock, currTime=currTime, &
     &      timeStep=timeStep, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &      line=__LINE__, &
     &      file=__FILE__)) &
     &      return  ! bail out

          !! SHYFEM advance of one time step $\Delta t_{ao}$, we need
          !! to access this information. From the \textsf{timeStep}
          !! object we query for the time step. The argument
          !! \textsf{s\_r8} tells to return it in second with double
          !! precision, which is the correct format for time variables
          !! in SHYFEM.
          call ESMF_TimeIntervalGet(timeStep, s_r8=timeStepSec, rc=rc)

	  !! To write the output fields we also write retrieve the current date in two
	  !! differen unit. One is used to name the output file in an
	  !! ordered fashion. The second one is to check if the current
	  !! date in seconds correspond to an output tick. In that case
	  !! we write the fields to file. 
	  call ESMF_TimeGet(currTime, yy=currYear, mm=currMonth, &
     &       dd=currDay, h=currHour, m=currMinute, s=currSecond, rc=rc) 
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &      line=__LINE__, &
     &      file=__FILE__)) &
     &      return  ! bail out
          call ESMF_TimeGet(currTime, yy=currYear, s=currTimeSec, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &      line=__LINE__, &
     &      file=__FILE__)) &
     &      return  ! bail out
	  write(dateString, "(I4,'-',I2.2,'-',I2.2, &
     &      '_',I2.2,':',I2.2,':',I2.2)") currYear, &
     &      currMonth, currDay, currHour, currMinute, currSecond
	  idtatm = nint( max(getpar('idtatm'), timeStepSec ))
          write(petString, "('.',I0,'.',I0)") &
     &     ShyfemToEsmf_Mesh%petCount,ShyfemToEsmf_Mesh%localPet

          idtatm = nint( max(getpar('idtatm'), timeStepSec ))

	  !! With a little bit of paranoia we print the import state object
	  !! (this needs to be removed or added only in a debug mode):
	  !call ESMF_StatePrint(importState, rc=rc)

	  !! Finally from the component, we retrieve the state and from the state
	  !! we retreive the field. The field data is typically a large array 
	  !! stored in memory, we use a pointer to it. Field can be written to a file.
	  !! with a  given frequency \textsf{idtatm}. Finally, after all these 
	  !! \textsf{Get} commands, a \textsf{Set} commands
	  !! copy the field pointer to the SHYFEM flux variables.
	  do var=1,SHYFEM_numOfImportFields

	    call ESMF_StateGet(importState, field=field, & 
     &        itemName=SHYFEM_FieldMetadata(var)%fieldName, rc=rc)
	    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &	      line=__LINE__, &
     &	      file=__FILE__)) &
     &	      return  ! bail out
          
	    call ESMF_FieldGet(field, localDe=0, farrayPtr=fieldPtr, &
     &	      rc=rc)
	    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &	      line=__LINE__, &
     &	      file=__FILE__)) &
     &	      return  ! bail out

	    if ( currTimeSec>0. .and. mod(currTimeSec,idtatm).lt.1 ) then
	      call SHYFEM_FieldWrite(field, &
     &		trim(SHYFEM_FieldMetadata(var)%fieldName)//trim(dateString) &
     &          //trim(petString)//trim(".vtk"), rc)
	      if(ESMF_LogFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU, &
     &          line=__LINE__, &
     &          file=__FILE__)) &
     &          return  ! bail out
	      !print *, fieldPtr !lrp: debug
	    endif

	    call SHYFEM_FieldSet(fieldPtr, &
     &        SHYFEM_FieldMetadata(var)%fieldName, rc)
            if(ESMF_LogFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU, &
     &        line=__LINE__, &
     &        file=__FILE__)) &
     &        return  ! bail out

	  enddo

	  !! The first task is completed. We move to the second one:
	  !! advancing the model of one timestep.
          !! This is the call to the SHYFEM subroutine that timesteps
          !! the ocean variables for one ocean-atmosphere timestep.
	  !! If everything goes well we print to the PET a message
	  !! saying that the ocean timestep has been completed
          if( shympi_is_master() ) then
            call getinfo(iuinfo)  !unit number of info file
          end if
	  if( iuinfo > 0 ) then	  
	    print * 
	    print *,'shyfem_component_run: start nuopc timestep at date ', &
     &        TRIM(dateString)
            print *
	  endif

	  call shyfem_run(timeStepSec)

	  call ESMF_TimePrint(currTime + timeStep, &
     &	    preString="---------------------> to: ",unit=msgString,rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &	    line=__LINE__, &
     &      file=__FILE__)) &
     &	    return  ! bail out
	  call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
 	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &	    line=__LINE__, &
     &	    file=__FILE__)) &
     &	    return  ! bail out

          !! The export state must be updated. A loop on the export fields is performed. 
          !! We can retrieve what is needed from the export state.
          !! To say that we want to retrieve a field item (of \textsf{ESMF\_Field} type) we 
          !! add the keyworld \textsf{field}. With \textsf{itemName} we choose
          !! the specific field among the ones that we have created. The next step is to run
          !! the toy model, that is we simply re-evaluate at the new time the arrays of fluxes. The
          !! fields are then assigned to the updated arrays:
          do var=SHYFEM_numOfImportFields+1, &
     &      SHYFEM_numOfImportFields+SHYFEM_numOfExportFields

            call ESMF_StateGet(exportState, field=field, &
     &        itemName=SHYFEM_FieldMetadata(var)%fieldName, rc=rc)
	    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &        line=__LINE__, &
     &        file=__FILE__)) &
     &        return  ! bail out

            call ESMF_FieldGet(field, localDe=0, farrayPtr=fieldPtr, &
     &        rc=rc)
	    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &        line=__LINE__, &
     &        file=__FILE__)) &
     &        return  ! bail out

            call SHYFEM_FieldGet(fieldPtr, &
     &        SHYFEM_FieldMetadata(var)%fieldName, rc)
!     &        ShyfemToEsmf_Mesh, rc)
 
            !print *, "OCN: ", fieldPtr !lrp: debug
 
          enddo


#ifdef NUOPC_TRACE
	  call ESMF_TraceRegionExit("OCN:Advance")
#endif
	end subroutine

	!-----------------------------------------------------------------------------

        !!
        !! \subsection{Finalize SHYFEM}
        !!
  	!! We register a method for the finalization of the ocean code. Here files are 
  	!! closed and the processes is killed with the memory cleaned through the 
	!! deallocations of the data structures. 
	subroutine Finalize(model, rc)
	  type(ESMF_GridComp)  :: model
	  integer, intent(out) :: rc

#define NUOPC_TRACE__OFF
#ifdef NUOPC_TRACE
	  call ESMF_TraceRegionEnter("OCN:Finalize")
#endif

	  rc = ESMF_SUCCESS

	  ! HERE THE MODEL IS FINALIZED: 
	  call shyfem_finalize

	  !deallocate(SHYFEM_FieldMetadata)
	  !deallocate(ShyfemToEsmf_Mesh%id_node_ghost)
	  !deallocate(ShyfemToEsmf_Mesh%ipv_ghost)
	  !deallocate(ShyfemToEsmf_Mesh%xgv_ghost)
          !deallocate(ShyfemToEsmf_Mesh%ygv_ghost)
	  !deallocate(ShyfemToEsmf_Mesh%nen3v_ghost)
	  !deallocate(ShyfemToEsmf_Mesh%table_ghostToLocal)

	  call ESMF_LogWrite("Finalized OCN", ESMF_LOGMSG_INFO, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &      line=__LINE__, &
     &      file=__FILE__)) &
     &	    return  ! bail out

#ifdef NUOPC_TRACE
	  call ESMF_TraceRegionExit("OCN:Finalize")
#endif
	end subroutine

	!-----------------------------------------------------------------------------

        !!
        !! \subsection{Get SHYFEM mesh}
        !!
	!! This section describes how to get the mesh from SHYFEM and set into a ESMF 
	!! Mesh class. It starts with an explanation of creating a Mesh and then goes 
	!! through other Mesh methods. To create a Mesh we need to set some properties 
	!! of the Mesh as a whole, some properties of each node in the mesh and then 
	!! some properties of each element which connects the nodes.
	!! For the Mesh as a whole we set its parametric dimension (\textsf{parametricDim}) 
	!! and spatial dimension (\textsf{spatialDim}). A Mesh's parametric dimension can 
	!! be thought of as the dimension of the elements which make up the Mesh.
	!! A Mesh's spatial dimension, on the other hand, is the is the number of coordinate 
	!! dimensions needed to describe the location of the nodes making up the Mesh. For
	!! SHYFEM both are always two, which means that the mesh is defined in parametric
	!! two-dimensional manifold. Even is spherical coordinates, SHYFEM use the spherical
	!! transformation to transform the differential operator on a two-dimensional. 
	!! Another important properties of SHYFEM mesh is that they are made of only 
	!! triangular elements. Without loss of generality we can assume that a mesh looks
	!! like the one in figure:
	!!
	!!\begin{minipage}{\linewidth} 
	!!\begin{verbatim}
	!!
	!! 
	!!  2.0   7 ------- 8 ------- 9
	!!        |  \   6  |  \   4  |
	!!        |    \    |    \    |
	!!        |  7    \ |  5    \ |
	!!  1.0   4 ------- 5 ------- 6
	!!        |  \   8  |  \   3  |
	!!        |    \    |    \    |
	!!        |  1    \ |  2   \  |
	!!  0.0   1 ------- 2 ------- 3
	!!
	!!       0.0       1.0        2.0 
	!! 
	!!        Node Id labels at corners
	!!       Element Id labels in centers
	!!       (Everything owned by PET 0) 
	!!
	!!\end{verbatim}
	!!\end{minipage}
	!!
	!! With this is mind let's have a look to the subroutine that get SHYFEM mesh.
	subroutine SHYFEM_MeshGet(SHYFEM_mesh, rc)

          type(ESMF_Mesh), intent(out)    :: SHYFEM_mesh
          integer, intent(out) 		  :: rc

	  !! First note the type \textsf{ESMF\_Mesh} which defines an unstructured grid
	  !! object. Other global and local properties of the mesh follows: 
	  type(ESMF_CoordSys_Flag)        :: SHYFEM_coordSys
	  integer                         :: SHYFEM_spatialDim = 2
	  integer                         :: SHYFEM_numberOfNodes
	  integer                         :: SHYFEM_numberOfElements
	  integer, allocatable            :: SHYFEM_nodeIds(:)
	  integer, allocatable            :: SHYFEM_nodeOwners(:)
          integer, allocatable            :: SHYFEM_nodeMask(:)
	  real(ESMF_KIND_R8), allocatable :: SHYFEM_nodeCoords(:)
          integer, allocatable            :: SHYFEM_elementIds(:)
	  integer, allocatable            :: SHYFEM_elementTypes(:)
	  integer, allocatable            :: SHYFEM_elementConn(:)
          real(ESMF_KIND_R8), allocatable :: SHYFEM_elementCoords(:)

	  integer                         :: i, ie
          real, external                  :: getpar

	  rc = ESMF_SUCCESS

!          call SHYFEM_ConvertMeshToEsmf(ShyfemToEsmf_Mesh)

	  !! We save the number of nodes and the number of elements
	  !! We got this from the well known global variables \textsf{nkn}
	  !! and \textsf{nel} of SHYFEM. These are available into the module
	  !! \textsf{grid} which is nested into \textsf{mod\_shyfem}
	  SHYFEM_numberOfNodes    = nkn_local
!     &     ShyfemToEsmf_Mesh%nkn_ghost
          SHYFEM_numberOfElements = nel_local
!     &     nel_unique

          allocate( SHYFEM_nodeIds(SHYFEM_numberOfNodes) )
          allocate( SHYFEM_nodeMask(SHYFEM_numberOfNodes) )
	  allocate( SHYFEM_nodeCoords(SHYFEM_numberOfNodes &
     &      *SHYFEM_spatialDim) )
          allocate( SHYFEM_nodeOwners(SHYFEM_numberOfNodes) )
	  allocate( SHYFEM_elementIds(SHYFEM_numberOfElements) )
          allocate( SHYFEM_elementTypes(SHYFEM_numberOfElements) )
          allocate( SHYFEM_elementConn(SHYFEM_numberOfElements*3) )
	  allocate( SHYFEM_elementCoords(SHYFEM_numberOfElements &
     &      *SHYFEM_spatialDim) )

	  !! The structure of the per node and element information used to 
	  !! create a Mesh is influenced by the Mesh distribution strategy. 
	  !! The Mesh class is distributed by elements. This means that a node 
	  !! must be present on any PET that contains an element associated with 
	  !! that node, but not on any other PET (a node can't be on a PET
	  !! without an element "home"). Since a node may be used by two or more 
	  !! elements located on different PETs, a node may be duplicated on 
	  !! multiple PETs. When a node is duplicated in this manner, 
	  !! one and only one of the PETs that contain the node must "own" the node. 
	  !! The user sets this ownership when they define the nodes during Mesh 
	  !! creation. For each node in the Mesh we set three properties:
	  !! the global id of the node (\textsf{nodeIds}), node coordinates 
	  !! (\textsf{nodeCoords}), and which PET owns the node ({\textsf{nodeOwners}).
	  !! The node id is a unique (across all PETs) integer attached to 
	  !! the particular node. It is used to indicate which nodes are the 
	  !! same when connecting together pieces of the Mesh on different 
	  !! processors. The node coordinates indicate the location of a node 
	  !! in space and are used in the \textsf{ESMF\_FieldRegrid()} functionality 
	  !! when interpolating. The node owner indicates which PET is in charge 
	  !! of the node. This is used when creating a Field on the Mesh to 
	  !! indicate which PET should contain a Field location for the data. 
	  !! The \textsf{nodeCoords} is an array containing the physical coordinates 
	  !! of the nodes to be created on this PET. This input consists of a 1D array 
	  !! the size of the number of nodes on this PET times the Mesh's spatial 
	  !! dimension (\textsf{spatialDim}). The coordinates in this array are ordered 
	  !! so that the coordinates for a node lie in sequence in memory. For the 
	  !! example shown above coordinates are in the following sequence
	  !! $\{x_1,\,y_1,\,...,\,x_9,\,y_9\}$. 
	  !SHYFEM_nodeOwners = 0
          do i=1,SHYFEM_numberOfNodes
            SHYFEM_nodeOwners(i) = id_node(i)
!     &        ShyfemToEsmf_Mesh%id_node_ghost(i)
            SHYFEM_nodeIds(i) = ipv(i)
!     &        ShyfemToEsmf_Mesh%ipv_ghost(i)
            SHYFEM_nodeMask(i) = 0
            SHYFEM_nodeCoords(i*SHYFEM_spatialDim-1) = xgv(i)
!     &        ShyfemToEsmf_Mesh%xgv_ghost(i)
            SHYFEM_nodeCoords(i*SHYFEM_spatialDim) = ygv(i)
!     &        ShyfemToEsmf_Mesh%ygv_ghost(i)
          enddo

	  !! For each element in the Mesh we set three properties: the global id 
	  !! of the element (\textsf{elementIds}), the topology type of the element 
	  !! (\textsf{elementTypes}), and which nodes are connected together to form 
	  !! the element (\textsf{elementConn}). The element id is a unique (across all PETs) 
	  !! integer attached to the particular element. The element type describes 
	  !! the topology of the element (e.g. a triangle vs. a quadrilateral). 
	  !! The range of choices for the topology of the elements in a Mesh are 
	  !! The element connectivity indicates which nodes are to be connected 
	  !! together to form the element. The number of nodes connected together for 
	  !! each element is implied by the elements topology type ({\textsf{elementTypes}). 
	  !! It is IMPORTANT to note, that 
	  !! the entries in this list are NOT the global ids of the nodes, but are indices 
	  !! into the PET local lists of node info used in the Mesh Create. In other words, 
	  !! the element connectivity isn't specified in terms of the global list of nodes, 
	  !! but instead is specified in terms of the locally described node info. One 
	  !! other important point about connectivities is that the order of the nodes in the 
	  !! connectivity list of an element is important. In general, when specifying an 
	  !! element with parametric dimension 2, the nodes should be given in 
	  !! counterclockwise order around the element. 
	  SHYFEM_elementTypes = ESMF_MESHELEMTYPE_TRI
	  SHYFEM_elementCoords = 0.
          do ie=1,SHYFEM_numberOfElements
           SHYFEM_elementIds(ie) = ipev(ie)
	    do i=1,3
              SHYFEM_elementConn((ie-1)*3+i) = &
!     &          ShyfemToEsmf_Mesh%nen3v_ghost(i,ie)
     &          nen3v(i,ie)
	      SHYFEM_elementCoords(ie*SHYFEM_spatialDim-1) = &
     &		SHYFEM_elementCoords(ie*SHYFEM_spatialDim-1) + xgv(nen3v(i,ie))
              SHYFEM_elementCoords(ie*SHYFEM_spatialDim)   = &
     &		SHYFEM_elementCoords(ie*SHYFEM_spatialDim)   + ygv(nen3v(i,ie))
	    enddo
          enddo
	  SHYFEM_elementCoords = 0.333333D0 * SHYFEM_elementCoords
	  !! The coordinate system is asked to the SHYFEM configuration file. Please
          !! note that the user must assure that a single coordinate system
          !! among all the the components. Two choice are possible in SHYFEM :
          !! \texttt{isphe=1} stands for spherical (in degrees) or 
	  !! \texttt{isphe=0} stands for cartesian. Please check the SHYFEM manual
	  !! for the details.
	  if ( nint(getpar('isphe')) == 1) then
	    SHYFEM_coordSys = ESMF_COORDSYS_SPH_DEG
	  else
            SHYFEM_coordSys = ESMF_COORDSYS_CART
	  endif

	  !! Once we have collected the mesh properties in the correct form, 
	  !! a final call create the mesh object at once. We also print to
	  !! file in vtk format, that can be visualize with Paraview.
	  SHYFEM_mesh = ESMF_MeshCreate(parametricDim=SHYFEM_spatialDim, &
     &      spatialDim=SHYFEM_spatialDim, &
     &      coordSys=SHYFEM_coordSys, &
     &      nodeIds=SHYFEM_nodeIds, nodeCoords=SHYFEM_nodeCoords, &
     &      nodeOwners=SHYFEM_nodeOwners, &
     &      elementIds=SHYFEM_elementIds, &
     &      elementTypes=SHYFEM_elementTypes, & 
     &	    elementConn=SHYFEM_elementConn, &
     &      rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &      line=__LINE__, &
     &      file=__FILE__)) &
     &      return  ! bail out
          call SHYFEM_MaskCreate(SHYFEM_mesh, SHYFEM_nodeMask, rc)

	  call ESMF_MeshFreeMemory(SHYFEM_mesh, rc)
          !! Once we have collected the mesh properties in the correct form,
          !! a final call create the mesh object at once. We also print to
          !! file in vtk format, that can be visualize with Paraview.
          SHYFEM_mesh = ESMF_MeshCreate(parametricDim=SHYFEM_spatialDim, &
     &      spatialDim=SHYFEM_spatialDim, &
     &      coordSys=SHYFEM_coordSys, &
     &      nodeIds=SHYFEM_nodeIds, nodeCoords=SHYFEM_nodeCoords, &
     &      nodeOwners=SHYFEM_nodeOwners, nodeMask=SHYFEM_nodeMask, &
     &      elementIds=SHYFEM_elementIds, &
     &      elementTypes=SHYFEM_elementTypes, &
     &      elementConn=SHYFEM_elementConn, &
     &      elementCoords=SHYFEM_elementCoords, &
     &      rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &      line=__LINE__, &
     &      file=__FILE__)) &
     &      return  ! bail out

          deallocate( SHYFEM_nodeIds, SHYFEM_nodeOwners, &
     &      SHYFEM_nodeMask, SHYFEM_nodeCoords )
	  deallocate( SHYFEM_elementIds, SHYFEM_elementTypes, & 
     &      SHYFEM_elementConn, SHYFEM_elementCoords )

	end subroutine

        !-----------------------------------------------------------------------------

        !! In parallel simulation we have to give to ESMF a distributed mesh
	!! in a different form from the many available in SHYFEM (unique, inner or local
	!! distribution).
	!! We start from the previous example with global nodes/elements numbering
	!! as follows:
        !!
	!!\begin{minipage}{\linewidth} 
        !!\begin{verbatim}
        !! 
        !!  2.0   7 ------- 8 ------- 9
        !!        |  \   6  |  \   4  |
        !!        |    \    |    \    |
        !!        |  7    \ |  5    \ |
        !!  1.0   4 ------- 5 ------- 6
        !!        |  \   8  |  \   3  |
        !!        |    \    |    \    |
        !!        |  1    \ |  2   \  |
        !!  0.0   1 ------- 2 ------- 3
        !!
        !!       0.0       1.0        2.0
        !!        Node Id labels at corners
        !!       Element Id labels in centers
	!!
        !!\end{verbatim}
        !!\end{minipage}
	!!
	!! The correct distributed mesh for ESMF on two PETs looks like:
        !!
        !!\begin{minipage}{\linewidth} 
        !!\begin{verbatim}
	!!
        !!  2.0   0 ------- 1 ------- 1
        !!        |  \   0  |  \   1  |
        !!        |    \    |    \    |
        !!        |  0    \ |  1    \ |
        !!  1.0   0 ------- 0 ------- 1
        !!        |  \   0  |  \   1  |
        !!        |    \    |    \    |
        !!        |  0    \ |  1   \  |
        !!  0.0   0 ------- 0 ------- 1
        !!
        !!       0.0       1.0        2.0
        !!        Node owners labels at corners
        !!       Element owners labels in centers
	!!
        !!\end{verbatim}
        !!\end{minipage}
	!!
	!! Please note that elements are unique to a single PET, while
	!! nodes can be shared among PETs. We consider the distributed
	!! mesh for PET0. The PET1 case is similar and it is not shown.
	!! SHYFEM knows the following decompositions. The
	!! so-called "inner nodes" appear at the beginning of the list
	!! and they corresponds to the node with \textsf{nodeOwner=0}
	!! from 1 to 5. Then in the local node list the other nodes
	!! appears without a predictable order:
	!!
	!!
        !!\begin{minipage}{\linewidth} 
        !!\begin{verbatim}
	!!	
        !!        5 ------- 8
        !!        |  \   4  | \
        !!        |    \    |    \
        !!        |  3    \ |  7    \
        !!        4 ------- 3 ------- 7
        !!        |  \   2  |  \   6  |
        !!        |    \    |    \    |
        !!        |  1    \ |  5   \  |
        !!        1 ------- 2 ------- 6
	!!
	!!        Local Node Id labels at corners
        !!       Local Element Id labels in centers
	!!		SHYFEM mesh on PET0
	!!
        !!\end{verbatim}
        !!\end{minipage}	
	!!
	!! ESMF, according to the tests, requests the following mesh on PET0
	!! which include only the so-called "inner elements", that is
	!! the elements that has \textsf{elementOwner=0}, and their
	!! vertices. Retrieve the elements from the previous mesh is
	!! trivial. For the nodes we need a consecutive list of node and
	!! we need to replace the local node index 8 with a local node
	!! index of 6.
	!!
        !!\begin{minipage}{\linewidth} 
        !!\begin{verbatim}
	!!
        !!        5 ------- 6
        !!        |  \   4  |
        !!        |    \    |
        !!        |  3    \ |
        !!        4 ------- 3
        !!        |  \   2  |
        !!        |    \    |
        !!        |  1    \ |
        !!        1 ------- 2
        !!
        !!        Local Node Id labels at corners
        !!       Local Element Id labels in centers	
        !!		ESMF mesh on PET0
	!!
        !!\end{verbatim}
        !!\end{minipage}
	!!	
	!! The reconstruction of the correct node list, and of the
	!! corresponding connetivity is what is done here:
        subroutine SHYFEM_ConvertMeshToEsmf(mesh)

	  type(SHYFEM_Mesh), intent(inout) :: mesh
          integer                          :: i, j, ie
	  integer                          :: nkn_ghost
	  integer                          :: id_ghost
          integer, allocatable             :: ipv_inner(:)
          logical, allocatable             :: nodeFlag(:)
!          integer, allocatable             :: table_ghostToLocal(:)
          integer, allocatable             :: table_localToGhost(:)
          integer, allocatable             :: tmp(:)
!
	  !! We have to define the set of ghost node.
	  !! From what seen before, inner and ghost nodes corresponds, i.e. they have the
	  !! same local numbering.
          allocate( ipv_inner(nkn_inner) )
          allocate( nodeFlag(nkn_local) )
          allocate( mesh%table_ghostToLocal(nkn_local) )
          allocate( table_localToGhost(nkn_local) )

          table_localToGhost = -1
          do i=1,nkn_inner
            ipv_inner(i)=ipv(i)
            mesh%table_ghostToLocal(i) = i
            table_localToGhost(i) = i
          enddo

	  !! The nodes after the inner set do not correspond. We need to
	  !! find the correct map between them. We loop over all the
	  !! unique element set and, through the global node numbering
	  !! we check if the node has been already added to the set of
	  !! ghost node or not. If not we assign a new local numbering
	  !! in an increasing fashion. The tables define the maps
	  !! between the inner set and the ghost set. The map between
	  !! that goes from the ghost set to the local set is saved for
	  !! future purposes.
          nkn_ghost=nkn_inner
          nodeFlag = .true.
          do ie=1,nel_unique
            do i=1,3

              id_ghost = ipv(nen3v(i,ie))
              if ( ALL(id_ghost.ne.ipv_inner) .and. &
     &             nodeFlag(nen3v(i,ie)) ) then
                nkn_ghost=nkn_ghost+1
                mesh%table_ghostToLocal(nkn_ghost) = nen3v(i,ie)
                table_localToGhost(nen3v(i,ie)) = nkn_ghost
                nodeFlag(nen3v(i,ie)) = .false.
              endif

            enddo
!	    print *, ShyfemToEsmf_Mesh%localPet,
!     &               ipev(ie)

          enddo

	  mesh%nkn_ghost = nkn_ghost

          allocate(tmp(nkn_ghost))
          tmp = mesh%table_ghostToLocal(1:nkn_ghost)
          call move_alloc(tmp, mesh%table_ghostToLocal)

          allocate( mesh%id_node_ghost(nkn_ghost) )
          allocate( mesh%ipv_ghost(nkn_ghost) )
          allocate( mesh%xgv_ghost(nkn_ghost) )
          allocate( mesh%ygv_ghost(nkn_ghost) )
          allocate( mesh%nen3v_ghost(3, nel_unique) )

          do i=1,nkn_inner
            mesh%id_node_ghost(i) = id_node(i)
            mesh%ipv_ghost(i) = ipv(i)
            mesh%xgv_ghost(i) = xgv(i)
            mesh%ygv_ghost(i) = ygv(i)
!	    print *, ShyfemToEsmf_Mesh%localPet,
!     &               id_node(i), ipv(i), i
          enddo

          do i=nkn_inner+1,nkn_ghost
            mesh%id_node_ghost(i) = id_node(mesh%table_ghostToLocal(i))
            mesh%ipv_ghost(i) = ipv(mesh%table_ghostToLocal(i))
            mesh%xgv_ghost(i) = xgv(mesh%table_ghostToLocal(i))
            mesh%ygv_ghost(i) = ygv(mesh%table_ghostToLocal(i))
!	    print *, ShyfemToEsmf_Mesh%localPet,
!     &               mesh%id_node_ghost(i),
!     &               mesh%ipv_ghost(i), i
          enddo

          do ie=1,nel_unique
            do i=1,3
              id_ghost = table_localToGhost(nen3v(i,ie))
              if (id_ghost < 1) then
                print *, "Bad conversion from"
                print *, "SHYFEM distributed mesh to the ESMF one"
                stop
              endif
              mesh%nen3v_ghost(i,ie) = id_ghost
            enddo
!	    print *, " con : ", ShyfemToEsmf_Mesh%localPet,
!     &                         ipev(ie),
!     &                          mesh%nen3v_ghost(1,ie),
!     &                          mesh%nen3v_ghost(2,ie),
!     &                          mesh%nen3v_ghost(3,ie)	    
          enddo
!	  print *, "Conversion to ESMF distributed mesh: done"
!          print *, "         nkn_inner  nkn_local  nkn_ghost"
!          print *, "     ",  nkn_inner, nkn_local, nkn_ghost, nel_unique

	  deallocate( ipv_inner, nodeFlag )
	  deallocate( table_localToGhost )

	end subroutine

        !-----------------------------------------------------------------------------

	!! We create the mask file. SHYFEM mesh need to be masked for
	!! points near land-sea boundaries. This is a known problem of
	!! the interpolation algorithm in atmosphere-ocean models. Here
	!! we use a trick to compute rapidly such mask.
	subroutine SHYFEM_MaskCreate(SHYFEM_mesh, nodeMask, rc)

          type(ESMF_Mesh), intent(in)           :: SHYFEM_mesh
	  integer, dimension(:), intent(inout)  :: nodeMask
	  integer, intent(out)                  :: rc

	  type(ESMF_Grid)                       :: gridForMask
	  type(ESMF_Field)                      :: onesField, maskField
	  double precision, pointer             :: fieldPtr2d(:,:)
	  double precision, pointer             :: fieldPtr1d(:)
	  type(ESMF_RouteHandle)                :: rh
          character(len=160)                    :: petString

	  rc = ESMF_SUCCESS

	  !! We need to flag land point of atmospheric grid. These are
	  !! computed by WRF but here I do not have access to this info
	  !! and I need to use a proxy.
	  !! We create a structured mesh which should be equal to the
	  !! atmospheric one. Then a real field is created on it.
	  !! A second constant field (one everywhere) is created attached to
	  !! the SHYFEM mesh.
	  gridForMask = ESMF_GridCreateNoPeriDimUfrm( &
     &      maxIndex=(/101, 101/), & !maxIndex=(/20, 20/),
     &      minCornerCoord=(/0.0D0,0.0D0/), & !minCornerCoord=(/9.0D0,38.0D0/),
     &      maxCornerCoord=(/3000000.0D0, 3000000.0D0/), & !maxCornerCoord=(/21.0D0, 46.0D0/),
     &      coordSys=ESMF_COORDSYS_CART, & !coordSys=ESMF_COORDSYS_SPH_DEG,
     &      staggerLocList=(/ESMF_STAGGERLOC_CENTER/), &
     &      rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &      line=__LINE__, &
     &      file=__FILE__)) &
     &      return  ! bail out

	  maskField = ESMF_FieldCreate(gridForMask, ESMF_TYPEKIND_R8, &
     &      staggerloc=ESMF_STAGGERLOC_CENTER, name="ones", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &      line=__LINE__, &
     &      file=__FILE__)) &
     &      return  ! bail out

	  onesField = ESMF_FieldCreate(SHYFEM_mesh,  ESMF_TYPEKIND_R8, &
     &      meshloc=ESMF_MESHLOC_NODE, name="mask", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &      line=__LINE__, &
     &      file=__FILE__)) &
     &      return  ! bail out
	  call ESMF_FieldGet(onesField, localDe=0, &
     &      farrayPtr=fieldPtr1d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &      line=__LINE__, &
     &      file=__FILE__)) &
     &      return  ! bail out
	  fieldPtr1d=1.

	  !! We reamp the constant field onto the structured grid.
	  !! Unmapped points (outise the boundaries of the SHYFEM mesh)
	  !! will be set to zero. These are considered the land points
	  !! of the atmospheric model. Of course points that are on
	  !! the sea but outside an open boundary are also flagged as
	  !! land points, which is clearly false but I don't know how to
	  !! do better.
	  call ESMF_FieldRegridStore(onesField, maskField, &
     &      unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
     &      routehandle=rh, &
     &      regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
     &      rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &      line=__LINE__, &
     &      file=__FILE__)) &
     &      return  ! bail out
	  call ESMF_FieldRegrid(onesField, maskField, rh, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &      line=__LINE__, &
     &      file=__FILE__)) &
     &      return  ! bail out

	  !! A second reamp from the strucured grid to the unstructured
	  !! one will leave a constant field only on such SHYFEM grid points
	  !! that does not fall in any atmospheric quad with even one
	  !! land point. We identify the constant field with a trivial
	  !! \textsf{where} statement. For these points we turn off the node mask.
	  call ESMF_FieldRegridStore(maskField, onesField, &
     &      unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
     &      routehandle=rh, &
     &      regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
     &      rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &      line=__LINE__, &
     &      file=__FILE__)) &
     &      return  ! bail out
	  call ESMF_FieldRegrid(maskField, onesField, rh, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &      line=__LINE__, &
     &      file=__FILE__)) &
     &      return  ! bail out

	  where( fieldPtr1d<1.D0) nodeMask=1

	  write(petString, "('.',I0,'.',I0)") ShyfemToEsmf_Mesh%petCount, &
     &                                        ShyfemToEsmf_Mesh%localPet
	  call SHYFEM_FieldWrite(onesField, &
     &      trim("shyfem_mask")//trim(petString)//trim(".vtk"), rc)

	end subroutine	

        !-----------------------------------------------------------------------------

        !! Eventually we write the regridded fields to files.
	!! This can be helpful for debugging and checking the interpolations.
	!! We can write such a file with ESMF subroutine \textsf{SHYFEM\_FieldWrite}
	!! but this works only with the third party library PARALLELIO (PIO).
	!! Moreover the only format allowed when this manual was written was netcdf
	!! (ugrid). We have preferred to use vtk format to visualize the data,
	!! as done with the mesh. This lead to only one type of file outputted.
	!! Vtk files can be visualize nicely with Paraview.
	!! A brief comment on parallel file writing. Each PET plots its
	!! own mesh and field. Since ESMF interpolates only the inner nodes, or
	!! the nodes with \textsf{nodeOwner=PET}, we have outputted only
	!! such nodes to file. In this way the locally unmapped nodes (computed
	!! correctly on other PETs) are not visualized. Collecting all
	!! the local file in paraview will allows to visualize all nodes
	!! with the correct remapped field but not all elements. A tiny band 
	!! of elements in between the inner nodes belongin to different PET 
	!! is skipped. This is not important because all interpolation fields are nodal
	!! and no useful information is attached to this elements.
        subroutine SHYFEM_FieldWrite(field, filename, rc)

          type(ESMF_Field), intent(in)  :: field
	  character(len=*), intent(in)  :: filename
          integer, intent(out)          :: rc

	  ! local variables
          double precision, pointer   :: ptr(:)
          character(20) :: str
          integer :: i, ie, ii

          rc = ESMF_SUCCESS

	  !! We recovet the field with the usual call
          call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, &
     &      rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     &      line=__LINE__, &
     &      file=__FILE__)) &
     &      return  ! bail out

	  !! We spend some word about the vtk format. We use
	  !! three digits to represents point coordinates, that means
	  !! that you can reach mm of resolution in you mesh but not beyond.
	  !! Paraview does not like zeros as float 0.00 and I had to manually
	  !! write zero as an integer. For the variable precision we assume
	  !! that three digits are also sufficient.
	  open(1, file = filename, status = 'unknown')
	  write(1,'(a)') "# vtk DataFile Version 3.0"
	  write(1,'(a)') "This file was generated by NUOPC"
	  write(1,'(a)') "ASCII"
	  write(1,'(a)') "DATASET UNSTRUCTURED_GRID"
	  write(1,'(a)', advance='no') "POINTS ";
!	  write(1, *) ShyfemToEsmf_Mesh%nkn_ghost, " double"
	  write(1, *) nkn_inner, " double"
	  do i=1,nkn_inner      !ShyfemToEsmf_Mesh%nkn_ghost
            write(str,'(f11.3)') xgv(i) !ShyfemToEsmf_Mesh%xgv_ghost(i)
	    str = trim(adjustl(str))
	    do ii = len_trim(str),1,-1
	      if (str(ii:ii)/="0") exit
	    enddo
            if (str(ii:ii)==".") ii=ii-1
	    write(1,'(a,a)', advance="no") str(1:ii), " "
            write(str,'(f11.3)') ygv(i) !ShyfemToEsmf_Mesh%ygv_ghost(i)
            str = trim(adjustl(str))
            do ii = len_trim(str),1,-1
              if (str(ii:ii)/="0") exit
            enddo
	    if (str(ii:ii)==".") ii=ii-1
            write(1,'(a,a,i1)') str(1:ii), " 0"
	  end do
	  write(1,'(a)', advance='no') "CELLS "
	  write(1, *) nel_inner, nel_inner*4  !nel_unique, nel_unique*4
          do ie=1,nel_inner     !nel_unique
            write(1,"(i6X,i6X,i6X,i6)") 3, &
!     &	      ShyfemToEsmf_Mesh%nen3v_ghost(1,ie)-1,
!     &       ShyfemToEsmf_Mesh%nen3v_ghost(2,ie)-1, 
!     &       ShyfemToEsmf_Mesh%nen3v_ghost(3,ie)-1
     &        nen3v(1,ie)-1, &
     &        nen3v(2,ie)-1, &
     &        nen3v(3,ie)-1
          end do
          write(1,'(a)', advance="no") "CELL_TYPES "
	  write(1, *) nel_inner !nel_unique
          do ie=1,nel_inner     !nel_unique
            write(1,*) 5
          end do
	  write(1,'(a)', advance="no") "POINT_DATA "
	  write(1, *) nkn_inner !ShyfemToEsmf_Mesh%nkn_ghost
	  write(1,'(a)') "SCALARS _NODE_NUM double 1"
	  write(1,'(a)') "LOOKUP_TABLE default"
          do i=1,nkn_inner      !ShyfemToEsmf_Mesh%nkn_ghost
            write(str,'(f11.4)') ptr(i)
            str = trim(adjustl(str))
            do ii = len_trim(str),1,-1
              if (str(ii:ii)/="0") exit
            enddo
            if (str(ii:ii)==".") ii=ii-1
            write(1,'(a,a)', advance="no") str(1:ii), " "
          end do
	  close(1)

        end subroutine

        !-----------------------------------------------------------------------------

	!! At the end of the cap layer, we come to the most important
	!! part. Where we update SHYFEM global (flux) variables copying them
	!! to the pointer value. From now on, the (flux) variables are
	!! syncronyzed with the atmospheric model. The right-hand-side
	!! of the expression is the SHYFEM variable name. SHYFEM in fact
	!! store as global variable the following variables:
	!! \begin{itemize}
	!! \item \textsf{ppv} the atmospheric pressure, in $[Pa]$
	!! \item \textsf{wxv} and \textsf{wyv} are the 10m wind
	!! components, in $[m/s]$
	!! \item \textsf{mettair} is the air temperature, in $[C]$
	!! \item \textsf{methum} is the air humidity $[0-1]$
	!! \item \textsf{metrad} is incoming short-wave-radiation
	!! $[W/m^2]$
	!! \item \textsf{metcc} is the cloud cover $[0-1]$
	!! \end{itemize}
	!! The name of the variable can be confusing. In fact depending on the 
	!! parametrizations SHYFEM can take as input also
	!! the wind stress and fluxes. In this
	!! case, these are stored in the same variables which must be
	!! intended as:
        !! \begin{itemize}      
        !! \item \textsf{ppv} the atmospheric pressure, in $[Pa]$
        !! \item \textsf{wxv} and \textsf{wyv} are the 10m wind
        !! stress components, in $[N/m^2]$ 
        !! \item \textsf{mettair} is the upward sensible heat flux $Q_{sens} $in $[W/m^2]$
        !! \item \textsf{methum} is the upward latent heat flux $Q_{lat} in $ $[W/m^2]$
	!! \item \textsf{metrad} is downward short-wave-radiation $R_{sw}$
	!! in $[W/m^2]$
        !! \item \textsf{metcc} is outgoing (upward) long-wave-radiation $[W/m^2]$
        !! \end{itemize}
	!! Also the sign should not be surprising: SHYFEM reverts the
	!! signs of all upward flux internally (positive if downward,
	!! negative if upward).
	!!
	!! The rest of the algorithm is unsurprising. We copy the inner
        !! nodes only with a loop since these are the ones owned by the
	!! local PET and remapped correctly. The remaining part of the
	!! arrays are copied from other PETs thanks to an MPI
	!! syncronization, which is encapsulated in the SHYFEM call
	!! \textsf{shympi\_exchange\_2d\_node}.
        subroutine SHYFEM_FieldSet(fieldPtr, fieldName, rc)

          double precision, pointer, intent(in) :: fieldPtr(:)
          character(len=*), intent(in)          :: fieldName
          integer, intent(out)                  :: rc

          rc = ESMF_SUCCESS

	  select case (fieldName)
	    case ("pmsl")
	      ppv(1:nkn_inner) = fieldPtr(1:nkn_inner)
	      call shympi_exchange_2d_node(ppv)
	    case ("smes")
	      wxv(1:nkn_inner) = fieldPtr(1:nkn_inner)
	      call shympi_exchange_2d_node(wxv)
            case ("smns")
	      wyv(1:nkn_inner) = fieldPtr(1:nkn_inner)
	      call shympi_exchange_2d_node(wyv)
            case ("stsh")
              mettair(1:nkn_inner) = fieldPtr(1:nkn_inner)
	      call shympi_exchange_2d_node(mettair)
            case ("stlh")
              methum(1:nkn_inner) = fieldPtr(1:nkn_inner)
	      call shympi_exchange_2d_node(methum)
	    case ("rsns")
	      metrad(1:nkn_inner) = fieldPtr(1:nkn_inner)
	      call shympi_exchange_2d_node(metrad)
            case ("rlns")
              metcc(1:nkn_inner) = fieldPtr(1:nkn_inner)
	      call shympi_exchange_2d_node(metcc)
	    case ("prec")
	      metrain(1:nkn_inner) = fieldPtr(1:nkn_inner)
	      call shympi_exchange_2d_node(metrain)
	    case default
              call ESMF_LogWrite("  OCN unknown field name", &
     &          ESMF_LOGMSG_INFO, rc=rc)
	     rc = 1
	  end select

	end subroutine

	!-----------------------------------------------------------------------------

        !! Next we syncronise the variable to be exported. At the
	!! opposite of the set routine, here we "get" the SHYFEM variable
        !! to be exported in order to pass them to the export state. 
	!! Here a list of such variables that are global
	!! in SHYFEM:
        !! \begin{itemize}
        !! \item \textsf{tempv} the sea temperature, in $[C]$
        !! \item \textsf{charn} the sea roughness expressed
        !!       with the Charnock parameter, dimensionless
        !! \end{itemize}
	!! We pass to the coupler only inner nodes. To realize
	!! interpolation I expect to pass also ghost nodes but
	!! a segmentation fault occurs. By the way passing inner
	!! nodes seems enough for ESMF: remapped sst field are correct.
        subroutine SHYFEM_FieldGet(fieldPtr, fieldName, rc)

          double precision, pointer, intent(inout) :: fieldPtr(:)
          character(len=*), intent(in)             :: fieldName
!          type(SHYFEM_Mesh), intent(inout)         :: mesh
	  integer, intent(out)                     :: rc

	  integer				   :: i

          rc = ESMF_SUCCESS

          select case (fieldName)
            case ("umap")
	      fieldPtr(1:nkn_inner) = 1.D0
!              do i=1,mesh%nkn_ghost
!                fieldPtr(i) = 1.0
!              enddo	  
            case ("sst")
	      fieldPtr(1:nkn_inner) = tempv(1,1:nkn_inner)
!              do i=1,mesh%nkn_ghost
!	        fieldPtr(i) = !tempv(1,mesh%table_ghostToLocal(i))
!	      enddo !lrp
            case ("z0")
              fieldPtr(1:nkn_inner) = charn(1:nkn_inner)
!              do i=1,mesh%nkn_ghost
!               fieldPtr(i) = charn(mesh%table_ghostToLocal(i))
!             enddo !lrp
            case default
              call ESMF_LogWrite("  OCN unknown field name", &
     &          ESMF_LOGMSG_INFO, rc=rc)
              rc = 1
          end select

        end subroutine


	end module
