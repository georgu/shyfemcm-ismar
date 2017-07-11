  module mod_enkf

  use mod_states
  use mod_observations

  implicit none

  character (len=80), save :: basfile, obsfile, rstfile

  ! parameters for the initial ensemble of states
  !
  integer, save :: bnew_ens
  integer nx_in,ny_in		!number of x and y grid points
  integer fmult_in		!mult factor to determine the supersampling
  real theta_in			!rotation of the random fields (0 East, anticlockwise)
  real sigma_in			!standard deviation of the fields (level)
  ! parameters for the computation of the model error
  !
  integer, save :: bmod_err
  integer nx_er,ny_er		!number of x and y grid points
  integer fmult_er		!mult factor to determine the supersampling 
  real theta_er			!rotation of the random fields (0 East, anticlockwise)
  real sigma_er			!standard deviation of the fields (level)
  double precision dt_er	!time between 2 analysis steps
  double precision tau_er	!time decorrelation (>=dt) of the old error: dq/dt = -(1/tau)*q

  integer, save :: nrens, nanal
  ! Observations
  !
  type(files), allocatable, dimension(:), private ::  ofile
  integer, private :: n_0dlev,n_2dvel		! number of obs files of different type
  double precision, save :: tobs
  integer, save :: nobs_tot			! total number of obs with status = 0 (good)
  type(levels), save, allocatable :: o0dlev(:)
  type(currentf), save, allocatable :: o2dvel(:)

  type(states), allocatable, save  :: A(:) 	! ensemble states
  type(states), save  :: Am			! average state
  type(states), allocatable :: qA(:)		! model error
  type(dstates), allocatable, save  :: Aaug(:) 	! double state with model error

  real, save, allocatable :: D(:,:)		! matrix holding perturbed measurments
  real, save, allocatable :: D1(:,:)		! Innovation vectors
  real, save, allocatable :: E(:,:)		! matrix holding perturbations (mode=?3)
  real, save, allocatable :: R(:,:)		! Obs error cov matrix

  real, save, allocatable :: S(:,:)		! matrix holding HA' mod perturbations
  real, save, allocatable :: innov(:)		! innovation vector holding d-H*mean(A)
  real, save, allocatable :: HA(:,:)		! matrix HA, ens model on obs space
  ! min-max values for the observations (by topaz)
  !
  real, parameter, private :: TEM_MIN = 0.0d0
  real, parameter, private :: TEM_MAX = 50.0d0
  real, parameter, private :: SAL_MIN = 0.0d0
  real, parameter, private :: SAL_MAX = 40.0d0
  real, parameter, private :: SSH_MIN = -3.0d0
  real, parameter, private :: SSH_MAX = 3.0d0
  real, parameter, private :: VEL_MIN = 0.0d0
  real, parameter, private :: VEL_MAX = 4.0d0
  ! multiplication factor for the ens std, to resize
  ! the obs std. Set <= 0 to disable it.
  !
  real, parameter, private :: KSTD = 2.
  ! decay time for the red noise of the observations
  !
  double precision, parameter, private :: TTAU = 2. * 86400.
  ! time interval to select the observations
  !
  double precision, parameter, private :: TEPS = 300.
  ! standard flag for bad data
  !
  real, parameter, private :: OFLAG = -999.

  contains


!********************************************************
!********************************************************
!********************************************************

  subroutine read_info

  implicit none

  integer n

  open(20, file='analysis.info', status='old')

  read(20,*) nrens	! number of ens members
  read(20,*) nanal		! analysis step
  read(20,*) basfile	! name of bas file (no extension)
  read(20,*) tobs	! current time of the observations
  read(20,*) obsfile	! name of obs file list
  read(20,*) bnew_ens	! 1 to create a new initial ens of states
  read(20,*) bmod_err	! 1 to use an augmented state with mod err

  close(20)

  if( mod(nrens,2).eq.0 ) error stop 'read_info: n of ens members must be odd, with the control as first.'

  if( bnew_ens.eq.1 ) then
    open(21, file='init_ens.info', status='old')
    read(21,*) nx_in,ny_in,fmult_in,theta_in,sigma_in	
    close(21)
  end if
  if( bmod_err.eq.1 ) then
    open(22, file='mod_err.info', status='old')
    read(22,*) nx_er,ny_er,fmult_er,theta_er,sigma_er,dt_er,tau_er
    close(22)
  end if

  write(*,*) 'time: ',tobs
  write(*,*) 'n. of ens members: ',nrens
  write(*,*) 'bnew_ens: ',bnew_ens
  write(*,*) 'bmod_err: ',bmod_err
  
  end subroutine read_info

!********************************************************

  subroutine set_shyfem_vars

  use basin
  implicit none

  open(21, file=basfile, status='old', form='unformatted')
  call basin_read_by_unit(21)
  close(21)

  if (( nkn.ne.nnkn ).or.( nel.ne.nnel)) error stop "read_basin: dim error"

  ! init modules
  !
  call init_shyfem_vars(nnkn,nnel,nnlv)

  ! addpar for restart
  !
  call add_rst_pars

  end subroutine set_shyfem_vars


!********************************************************

  subroutine read_obs
  implicit none

  integer ios
  character(len=80) :: line

  integer n,nfile
  integer kinit,kend
  logical linit
  integer nobs

     write(*,*) 'Observation file list: ',trim(obsfile)

  !-------------------------------
  ! Read a list of obs files
  !-------------------------------
     n = 1
     open(25,file = obsfile, status = 'old', form = 'formatted', iostat = ios)
     if( ios.ne.0 ) error stop 'read_obs: error opening file'
 88   read(25,*,end=98) line
      n = n + 1
     goto 88
 98  continue
     nfile = n - 1
     rewind(unit=25)

     allocate(ofile(nfile))
     do n = 1,nfile
        read(25,*,err=95) ofile(n)
     end do
     close(25)

  !-------------------------------
  ! Read every file one time to find the number and the type of each obs
  !-------------------------------
  n_0dlev = 0
  n_2dvel = 0
  do n = 1,nfile

    select case (trim(ofile(n)%ty))

        case default
             error stop 'read_obs: Unknown file type'

        ! Level timeseries
        !
        case ('0DLEV')
             linit = .true.
             kinit = n_0dlev
             call read_level(linit,trim(ofile(n)%name),tobs,TEPS,&
                             ofile(n)%id,kinit,kend,nobs)
             if (kend > kinit) n_0dlev = n_0dlev + 1

        ! 2d current fem files
        !
        case ('2DVEL')
             n_2dvel = n_2dvel + 1

    end select

  end do


  !-------------------------------
  ! allocate obs vars
  !-------------------------------
  if (n_0dlev > 0) allocate(o0dlev(n_0dlev))
  if (n_2dvel > 0) allocate(o2dvel(n_2dvel))

  !-------------------------------
  ! read a second time and store
  !-------------------------------
  n_2dvel = 0
  kinit = 0
  nobs_tot = 0
  do n = 1,nfile

    select case (trim(ofile(n)%ty))
        case default
             error stop 'read_obs: Unknown file type'
        case ('0DLEV')
             linit = .false.
             call read_level(linit,trim(ofile(n)%name),tobs,TEPS,&
                             ofile(n)%id,kinit,kend,nobs)
             kinit = kend
             nobs_tot = nobs_tot + nobs
        case ('2DVEL')
             n_2dvel = n_2dvel + 1
             call read_2dvel(trim(ofile(n)%name),ofile(n)%id,tobs,n_2dvel,TEPS,nobs)
             nobs_tot = nobs_tot + 2 * nobs ! u and v components
    end select

  end do

  return

 95 write(*,*) 'read_obs error reading file: ',trim(ofile(n)%name)
    stop

  end subroutine read_obs

!******************************************************

  subroutine read_level(linit,filin,tobs,eps,id,kinit,kend,nobs)

  implicit none
  logical,intent(in)           :: linit
  character(len=*),intent(in)  :: filin
  double precision, intent(in) :: tobs
  double precision, intent(in) :: eps
  integer, intent(in)	       :: id
  integer, intent(in)          :: kinit
  integer, intent(out)         :: kend
  integer, intent(out)         :: nobs
  integer ios
  double precision tt
  real x, y, v, stdv
  integer k
  integer ostatus

  nobs = 0

  ! read info file for the obs file with coords and std
  !
  open(27,file=trim(filin)//'.info', status = 'old', form = 'formatted', iostat = ios)
  if (ios /= 0) error stop 'read_level: error opening info file'
  read(27,*) x,y,stdv
  close(27)

  ! read obs file
  !
  open(26,file=trim(filin), status = 'old', form = 'formatted', iostat = ios)
  if (ios /= 0) error stop 'read_level: error opening file'

  select case(linit)
  
     ! just look if there is a valid obs inside
     !
     case (.true.)
          k = kinit
 90       read(26,*,end=100) tt, v
          ! Take only records with times near tobs
          !
          if (abs(tt - tobs) < eps) then
             ostatus = 0
             ! check if the obs value is out of range
             !
             call check_obs('0DLEV',v,v,OFLAG,ostatus)
             k = k + 1
          end if
          goto 90
 100      close(26)
          kend = k

     ! store the obs
     !
     case (.false.)
          k = kinit
 91       read(26,*,end=101) tt, v
          ostatus = .false.
          ! Take only records with times near tobs
          !
          if (abs(tt - tobs) < eps) then
             ostatus = 0
             ! check if the obs value is out of range
             !
             call check_obs('0DLEV',v,v,OFLAG,ostatus)
             k = k + 1
             o0dlev(k)%t = tt
             o0dlev(k)%x = x
             o0dlev(k)%y = y
             o0dlev(k)%val = v
             o0dlev(k)%std = stdv
             o0dlev(k)%status = ostatus
             o0dlev(k)%id = id

             nobs = nobs + 1

             goto 101	!take just 1 lev obs
          end if
          goto 91
 101      close(26)
          kend = k

  end select

  end subroutine read_level

!********************************************************

  subroutine read_2dvel(filin,fid,tobs,nrec,eps,nobs)
  implicit none

  character(len=*),intent(in)  :: filin
  integer,intent(in)           :: fid
  double precision,intent(in)  :: tobs
  integer, intent(in)          :: nrec
  double precision,intent(in)  :: eps
  integer, intent(out)          :: nobs

  integer ios
  integer np,iformat,iunit
  integer irec,i,ii
  double precision tt
  integer nvers           !version of file format
  integer lmax            !vertical values
  integer nvar            !number of variables to write
  integer ntype           !type of information contained
  integer datetime(2)     !date and time information
  integer ierr            !return error code
  integer nlvddi
  real*4,allocatable :: hhlv(:)   !vertical structure
  real*4 regpar(7)               !regular array params
  integer nx,ny
  real flag,dx,dy,x0,y0
  integer,allocatable :: ilhkv(:)
  real*4,allocatable :: hd(:)
  real*4,allocatable :: data(:,:),dataens(:,:)
  character(len=50) :: string
  integer ostatus
  logical bdata
  real x,y,uu,vv,ostd
  integer ix,iy

  nobs = 0

  ! read info file for the obs file with std
  !
  open(27,file=trim(filin)//'.info', status = 'old', form = 'formatted', iostat = ios)
  if (ios /= 0) error stop 'read_2dvel: error opening file'
  read(27,*) ostd
  close(27)

  ! open fem file
  !
  np = 0
  call fem_file_read_open(trim(filin),np,iformat,iunit)

  irec = 0
  bdata = .false.
  do
    irec = irec + 1

    ! Reads headers
    !
    call fem_file_read_params(iformat,iunit,tt,nvers,np,lmax,nvar,ntype,datetime,ierr)
    if( ierr .lt. 0 ) exit

    allocate(hhlv(lmax))
    nlvddi = lmax
    call fem_file_read_2header(iformat,iunit,ntype,lmax,hhlv,regpar,ierr)

    nx = nint(regpar(1))
    ny = nint(regpar(2))
    x0 = regpar(3)
    y0 = regpar(4)
    dx = regpar(5)
    dy = regpar(6)
    flag = regpar(7)

    if (flag /= OFLAG) error stop 'read_2dvel: bad flag'

    if (abs(tt - tobs) > eps) then
       do i=1,nvar
          call fem_file_skip_data(iformat,iunit,nvers,np,lmax,string,ierr)
          if (ierr /= 0) error stop 'read_2dvel: error reading file'
       end do
       cycle
    end if

    allocate(ilhkv(np),hd(np),data(nx,ny))
    allocate(o2dvel(nrec)%x(nx),o2dvel(nrec)%y(ny),&
             o2dvel(nrec)%u(nx,ny),o2dvel(nrec)%v(nx,ny),&
             o2dvel(nrec)%std(nx,ny),o2dvel(nrec)%status(nx,ny))

    do i = 1,nvar
       call fem_file_read_data(iformat,iunit,nvers,np,lmax,string,ilhkv,hd,nlvddi,&
                               data,ierr)

       select case (i)
              case default
                   error stop 'read_2dvel: too many variables'
              case (1)
                   o2dvel(nrec)%u = data
              case (2)
                   o2dvel(nrec)%v = data
       end select
    end do
    deallocate(ilhkv,hd,data)

    ! find and assign the coords
    !
    do ii = 1,nx
       o2dvel(nrec)%x(ii) = x0 + dx * (ii-1)	! x coords
    end do
    do ii = 1,ny
       o2dvel(nrec)%y(ii) = y0 + dy * (ii-1)	! y coords
    end do
    o2dvel(nrec)%z = 0.		! depth

    ! check observation values
    !
    do ix = 1,nx
    do iy = 1,ny
       uu = o2dvel(nrec)%u(ix,iy)
       vv = o2dvel(nrec)%v(ix,iy)

       call check_obs('2DVEL',uu,vv,flag,ostatus)
       o2dvel(nrec)%status = ostatus

       if (ostatus == 0) then
          bdata = .true. !at least 1 rec
          nobs = nobs + 1
       end if

    end do
    end do

    ! last variables to assign
    !
    o2dvel(nrec)%nx = nx	! n of x points
    o2dvel(nrec)%ny = ny	! n od y points
    o2dvel(nrec)%std = ostd	! standard dev
    o2dvel(nrec)%id = fid	! id of the original file
    o2dvel(nrec)%status = ostatus ! status of the obs

    ! exit after finding 1 field
    !
    exit

  end do
  close(iunit)

  if (.not.bdata) error stop 'read_2dvel: 2dvel file without valid data'

  end subroutine read_2dvel

!********************************************************

  subroutine average_mat(date,time,rstw)

  use mod_hydro
  use mod_ts
  implicit none
  integer,intent(in) :: date,time
  integer,intent(in) :: rstw
  integer ne
  character(len=16) :: rstname
  character(len=3) :: nal

  type(states4) :: A4

  ! make the average
  !
  Am = 0.
  do ne = 1,nrens
    Am = Am + A(ne)
  end do
  Am = states_real_mult(Am,1./float(nrens))

  ! write in a file
  !
  A4 = Am
  call pull_state(A4)
  call num2str(nanal,nal)
  if (rstw == -1) then
        write(*,*) 'Writing average background state...'
        rstname = 'an'//nal//'_enavrb.rst'
        call rst_write(rstname,tobs,date,time)
  elseif (rstw == -2) then
        write(*,*) 'Writing average analysis state...'
        rstname = 'an'//nal//'_enavra.rst'
        call rst_write(rstname,tobs,date,time)
  end if

  end subroutine average_mat


!********************************************************

  subroutine make_matrices
! R (only used if mode=?1 or ?2) (not with low rank R: (N-1) R = EE')
!
  implicit none

  integer nook

  allocate(D(nobs_tot,nrens),E(nobs_tot,nrens),&
           R(nobs_tot,nobs_tot))
  allocate(S(nobs_tot,nrens),innov(nobs_tot))
  allocate(HA(nobs_tot,nrens),D1(nobs_tot,nrens))

  R(:,:) = 0.
  nook = 0
  ! This is for the levels
  !
  if (n_0dlev > 0) then
     call fill_levels(n_0dlev,nook)
  end if

  ! This is for the currents
  !
  if (n_2dvel > 0) then
     call fill_scurrents(n_2dvel,nook)
  end if

  end subroutine make_matrices

!********************************************************

  subroutine fill_levels(nfile,nook)
  implicit none
  integer, intent(in) :: nfile
  integer, intent(inout) :: nook
  integer nf,ne
  real x,y
  integer iemin,kmin
  real oval,stdv
  integer ostatus
  real inn1,inn2
  real pvec(nrens)
  character(len=3) :: nal,nfil

  do nf = 1,nfile 

     ! create a red noise random vector with mean 0 and std 1
     !
     call make_0Dpert('z',nrens,nanal,o0dlev(nf)%id,pvec,tobs,TTAU)

     ! next if the observation is not good
     !
     if (o0dlev(nf)%status /= 0) cycle

     nook = nook + 1

     x = o0dlev(nf)%x
     y = o0dlev(nf)%y
     call find_el_node(x,y,iemin,kmin)

     ! compute the observation errors R
     !
     !R(nook,nook) = o0dlev(nf)%std**2

     ! compute the model perturbed values, S = HA' and HA
     ! Remember for enKF: Aa = Af + A' [HA']^t [ U L^-1 U^t ] D' and D' = D-HA
     !
     do ne = 1,nrens
        !S(nook,ne) = sum( A(ne)%ze(:,iemin) - Am%ze(:,iemin) )/3. !average of the three vertexes
        !HA(nook,ne) = sum( A(ne)%ze(:,iemin) )/3. !average of the three vertexes
        S(nook,ne) = A(ne)%z(kmin) - Am%z(kmin)
        HA(nook,ne) = A(ne)%z(kmin)
     end do

     ! compute the innovation vector
     !
     oval = o0dlev(nf)%val
     ostatus = o0dlev(nf)%status
     stdv = o0dlev(nf)%std
     call check_obs_inn('0DLEV',x,y,0.,oval,oval,stdv,inn1,inn2,ostatus)
     o0dlev(nf)%std = stdv
     innov(nook) = inn1
     write(*,'(a,i5,3f8.4)') ' nobs, vobs, vmod, innov: ',nf,o0dlev(nf)%val,Am%z(kmin),inn1
 
     ! compute the perturbations E, the perturbed observations D
     ! and the innovation vectors D1
     !
     E(nook,:) = o0dlev(nf)%std * pvec
     D(nook,:) = o0dlev(nf)%val + (o0dlev(nf)%std * pvec)
     D1(nook,:) = D(nook,:) - HA(nook,:)
 
  end do
  end subroutine fill_levels

!********************************************************

  subroutine fill_scurrents(nfile,nook)
  implicit none
  integer, intent(in) :: nfile
  integer, intent(inout) :: nook
  integer nf,ix,iy,ne
  real pvec1(nrens),pvec2(nrens)
  real x,y
  integer iemin,kmin
  real uu,vv,stdv
  integer ostatus
  real inn1,inn2

 
  do nf = 1,nfile 
     
    ! create a red noise random vector with mean 0 and std 1
    !
    call make_0Dpert('u',nrens,nanal,o2dvel(nf)%id,pvec1,tobs,TTAU)
    call make_0Dpert('v',nrens,nanal,o2dvel(nf)%id,pvec2,tobs,TTAU)

    do iy = 1,o2dvel(nf)%ny
    do ix = 1,o2dvel(nf)%nx

     ! next if the observation is not good
     !
     if (o2dvel(nf)%status(ix,iy) /= 0) cycle

     nook = nook + 2	! 2 obs, u and v components

     x = o2dvel(nf)%x(ix)
     y = o2dvel(nf)%y(iy)
     call find_el_node(x,y,iemin,kmin)

     ! compute the observation errors R
     !
     !R(nook-1,nook-1) = o2dvel(nf)%std(ix,iy)**2
     !R(nook,nook) = o2dvel(nf)%std(ix,iy)**2

     ! compute the model perturbed values, S = HA' and HA
     ! Remember for enKF: Aa = Af + A' [HA']^t [ U L^-1 U^t ] D' and D' = D-HA
     !
     do ne = 1,nrens
        S(nook-1,ne) = A(ne)%u(1,iemin) - Am%u(1,iemin)
        S(nook,ne) = A(ne)%v(1,iemin) - Am%v(1,iemin)
        HA(nook-1,ne) = A(ne)%u(1,iemin)
        HA(nook,ne) = A(ne)%v(1,iemin)
     end do

     ! compute the innovation vector
     !
     ostatus = o2dvel(nf)%status(ix,iy)
     uu = o2dvel(nf)%u(ix,iy)
     vv = o2dvel(nf)%v(ix,iy)
     stdv = o2dvel(nf)%std(ix,iy)
     call check_obs_inn('2DVEL',x,y,0.,uu,vv,stdv,inn1,inn2,ostatus)
     o2dvel(nf)%std(ix,iy) = stdv

     innov(nook-1) = inn1
     innov(nook) = inn2
 
     ! compute the perturbations E, the perturbed observations D
     ! and the innovation vectors D1
     !
     ! u component
     !
     E(nook-1,:) = o2dvel(nf)%std(ix,iy) * pvec1
     D(nook-1,:) = o2dvel(nf)%u(ix,iy) + (o2dvel(nf)%std(ix,iy) * pvec1)
     D1(nook-1,:) = D(nook-1,:) - HA(nook-1,:)
     ! v component
     !
     E(nook,:) = o2dvel(nf)%std(ix,iy) * pvec2
     D(nook,:) = o2dvel(nf)%v(ix,iy) + (o2dvel(nf)%std(ix,iy) * pvec2)
     D1(nook,:) = D(nook,:) - HA(nook,:)
 
    end do
    end do

  end do

  end subroutine fill_scurrents

!********************************************************

  subroutine read_ensemble(date,time)

   implicit none

   integer,intent(out) :: date,time

   type(states4) :: A4
   type(states) :: Ap(nrens-1)
   character(len=3) :: nrel,nal
   character(len=16) rstname
   integer ne

   ! Allocates the state A to store the ens states
   allocate(A(nrens))

   call num2str(nanal,nal)

   if ((bnew_ens.eq.0) .or. (nanal.gt.1)) then

     write(*,*) 'Loading an ensemble of initial states'
     do ne = 1,nrens
        call num2str(ne-1,nrel)
        rstname='an'//nal//'_'//'en'//nrel//'b.rst'
        call rst_read(rstname,tobs,date,time)
        call push_state(A4)
        A(ne)=A4
     end do

   else if ((bnew_ens.eq.1) .and. (nanal.eq.1)) then

     write(*,*) 'Creating a new ensemble of initial states'
     !read an input restart file
     call num2str(0,nrel)
     rstname = 'an'//nal//'_'//'en'//nrel//'b.rst'
     call rst_read(rstname,tobs,date,time)

     !push the vars into the state and makes the ens
     call push_state(A4)
     call make_init_ens(A4)
     
     !save the initial ens in new restart files
     call num2str(nanal,nal)
     do ne = 1,nrens
        call num2str(ne-1,nrel)
        A4 = A(ne)
        call pull_state(A4)
        rstname='an'//nal//'_'//'en'//nrel//'b.rst'
        call rst_write(rstname,tobs,date,time)
     end do

   else

     write(*,*) 'Not a valid option for bnew_ens'
     error stop

   end if

   return
  end subroutine read_ensemble

!********************************************************

  subroutine write_ensemble(date,time)
   implicit none
   integer,intent(in) :: date,time

   type(states4) :: A4
   character(len=3) :: nrel,nal
   character(len=16) rstname
   integer ne
   integer dt,tm

   call num2str(nanal,nal)
   do ne = 1,nrens
      call num2str(ne-1,nrel)
      rstname='an'//nal//'_'//'en'//nrel//'b.rst'
      call rst_read(rstname,tobs,dt,tm) !This is to load var not present 
                                        ! in the ens state. It should be removed.
      A4 = A(ne)
      call pull_state(A4)
      rstname='an'//nal//'_'//'en'//nrel//'a.rst'
      call rst_write(rstname,tobs,date,time)
   end do
  end subroutine write_ensemble

!********************************************************

   subroutine push_state(AA)
    use mod_hydro
    use mod_hydro_vel
    use mod_ts
    use mod_conz
    implicit none
 
    type(states4) :: AA

    AA%u = utlnv
    AA%v = vtlnv
    AA%ze = zenv
    AA%z = znv
    AA%t = tempv
    AA%s = saltv
   
   end subroutine push_state

!********************************************************

   subroutine pull_state(AA)
    use mod_hydro
    use mod_hydro_vel
    use mod_ts
    use mod_conz
    implicit none
 
    type(states4) :: AA

    utlnv = AA%u
    vtlnv = AA%v
    zenv = AA%ze
    znv = AA%z
    tempv = AA%t
    saltv = AA%s

   end subroutine pull_state


!********************************************************

  subroutine make_init_ens(A4)
   use basin
   implicit none
   type(states4),intent(in) :: A4

   type(states), allocatable, save :: Apert(:)
   type(states) :: Aaux
   real kvec(nnkn,nrens-1),evec(nnel,nrens-1)

   integer ne,n,ie,k

   Aaux = A4

   allocate(Apert(nrens-1))

   ! perturbation for ze
   call make_pert(kvec,nnkn,nrens-1,fmult_in,theta_in,nx_in,ny_in)

   do ne = 1,nrens-1
     call assign_states(Apert(ne),0.)

     do ie = 1,nnel
        do n = 1,3
           k = nen3v(n,ie)
           Apert(ne)%ze(n,ie) = kvec(k,ne) * sigma_in
        end do
     end do

   end do

   ! The first state is unperturbed
   A(1) = Aaux
   do ne = 2,nrens
      A(ne) = add_states(Aaux,Apert(ne-1))
   end do

  end subroutine make_init_ens

!********************************************************

  subroutine push_aug
   use basin
   implicit none
   integer nst 		!= nanal-1	number of time steps from the begin of the assimilation
   double precision alpha,rho

   type(states) :: A1
   real kvec(nnkn,nrens)

   ! Parameters for the sample fields
   integer nx,ny 	!grid dimension
   integer fmult	!the start ensemble is fmult*nrens
   real theta		!rotation of the fields (0 East, anticlockwise)

   real mfact
   integer ne,ie,n,k

   write(*,*) '*********************************************'
   write(*,*) 'Creating an augmented state with model errors'
   write(*,*) '*********************************************'

   !---------------------------------------
   ! defines parameters for the model error
   !---------------------------------------
   nst = nanal 
 
   if( tau_er.lt.dt_er ) error stop 'make_aug: parameter error'
 
   alpha = 1. - (dt_er/tau_er)
   rho=sqrt( (1.0-alpha)**2 /&
     (dt_er*(float(nst) - 2.0*alpha - float(nst)*alpha**2 + 2.0*alpha**(nst+1))) )
 
   !---------------------------------------
   !makes the new white noise field
   !---------------------------------------
   call make_pert(kvec,nnkn,nrens,fmult_er,theta_er,nx_er,ny_er)
 
   allocate(qA(nrens))
   do ne = 1,nrens
      call assign_states(qA(ne),0.)
 
      do ie = 1,nnel
        do n = 1,3
          k = nen3v(n,ie)
          qA(ne)%ze(n,ie) = kvec(k,ne)
        end do
      end do

   end do
 
   !---------------------------------------
   !if exist old error q0 load it and add to qA (q1 = alpha*q0 + sqrt(1-alpha**2)*w)
   !---------------------------------------
   call load_error(alpha,tobs-dt_er)
 
   !---------------------------------------
   !compute the new state qA = A + sqrt(dt)*sigma*rho*q1
   !---------------------------------------
   ! change this if you want errors not only in zeta
   mfact = sqrt(dt_er) * sigma_er * rho
   do ne = 1,nrens
      A1 = states_real_mult(qA(ne),mfact)
      A(ne) = A(ne) + A1
   end do
    
   !---------------------------------------
   !make the augmented state Aaug = (A,qA)
   !---------------------------------------
   allocate(Aaug(nrens))
   do ne = 1,nrens
      call push_dstate(A(ne),qA(ne),Aaug(ne))
   end do
   deallocate(A,qA)
 
  end subroutine push_aug

!********************************************************

  subroutine load_error(alpha,tt)

   implicit none
   real, intent(in) :: alpha
   double precision, intent(in) :: tt

   logical :: bfile
   integer ne
   character(len=3) :: nrel,nal
   character(len=19) rstname
   type(states4) :: A4
   type(states) :: A1,A2
   integer date,time

   real mfact

   ! Old analysis step
   call num2str(nanal-1,nal)

   ! Check if error files exist
   do ne = 1,nrens
      call num2str(ne-1,nrel)
      rstname='an'//nal//'_'//'en'//nrel//'_err.rst'
      inquire(file=rstname, exist=bfile)
      if(.not. bfile) goto 777
   end do

   ! Add the old error to the new one
   write(*,*) '********'
   write(*,*) 'Loading model error from files'
   do ne = 1,nrens
      call num2str(ne-1,nrel)
      rstname='an'//nal//'_'//'en'//nrel//'_err.rst'
      call rst_read(rstname,tt,date,time)
      call push_state(A4)

      ! q1 = alpha*q0 + sqrt(1-alpha**2)*w
      A1 = A4
      A1 = states_real_mult(A1,alpha)

      mfact = sqrt(1 - alpha**2)
      A2 = states_real_mult(qA(ne),mfact)

      qA(ne) = A1 + A2

   end do

   return

 777 continue

   write(*,*) '********'
   write(*,*) 'Model error files not found'

  end subroutine load_error

!********************************************************

  subroutine pull_aug(date,time)
  implicit none
  integer,intent(in) :: date,time
  character(len=3) :: nrel,nal
  character(len=19) rstname
  type(states4) :: A4
  integer ne

  write(*,*) '********'
  write(*,*) 'Saving model errors'

  allocate(A(nrens),qA(nrens))
  do ne = 1,nrens
     call pull_dstate(A(ne),qA(ne),Aaug(ne))
  end do
  deallocate(Aaug)

  call num2str(nanal,nal)

  do ne = 1,nrens
     A4 = qA(ne)
     call pull_state(A4)

     call num2str(ne-1,nrel)
     rstname='an'//nal//'_'//'en'//nrel//'_err.rst'
     call rst_write(rstname,tobs,date,time)
  end do

  end subroutine pull_aug

!********************************************************

  subroutine check_obs(ty,v1,v2,flag,status)

  implicit none
  character(len=*), intent(in) :: ty
  real, intent(in) :: v1,v2
  real, intent(in) :: flag
  integer, intent(out) :: status

  real vmin,vmax,v

  status = 0

  if (v1 == flag .or. v2 == flag) then
     status = 2
     return
  end if

  if (trim(ty) == '0DLEV') then
     vmin = SSH_MIN
     vmax = SSH_MAX
     v = v1
  else if (trim(ty) == '2DVEL') then
     vmin = VEL_MIN
     vmax = VEL_MAX
     v = sqrt(v1**2 + v2**2)
  else
     write(*,*) 'Observation not implemented yet'
     error stop
  end if

  if (v < vmin .or. v > vmax) then
     status = 1
  end if

  end subroutine check_obs

!********************************************************

  ! Compute the innovation vectors and change the std if
  ! it is too low
  !
  subroutine check_obs_inn(ty,x,y,z,v1,v2,stdv,d1,d2,status)

  use levels
  implicit none
  character(len=*), intent(in) :: ty
  real, intent(in) :: x,y,z
  integer, intent(in) :: status
  real, intent(in) :: v1,v2
  real, intent(inout) :: stdv
  real, intent(out) :: d1,d2

  real vmod,vmod_ens(nrens)
  real v
  integer ie,ik,ne
  real h_1st_layer
  real ens_std,stdv_new
  real inn

  ! Exit if it is a bad obs, not to be assimilated
  !
  if (status /= 0) return

  ! Find the nearest node/element
  !
  call find_el_node(x,y,ie,ik)

  select case (trim(ty))

    case default
     write(*,*) 'Observation not implemented yet'
     error stop

    case ('0DLEV')

     vmod = Am%z(ik)
     vmod_ens = A(:)%z(ik)
     v = v1
     d1 = v - vmod
     d2 = v - vmod
     inn = d1

    case ('2DVEL')

     if (size(hlv) <= 1)&
     error stop 'check_obs_inn: a 3D sim is necessary to assimilate surface currents'

     vmod = sqrt(Am%u(1,ie)**2 + Am%v(1,ie)**2)
     vmod_ens = sqrt(A(:)%u(1,ie)**2 + A(:)%v(1,ie)**2)
     h_1st_layer = hlv(1) + Am%z(ik)

     v = sqrt(v1**2 + v2**2)
     d1 = v1*h_1st_layer - Am%u(1,ie)
     d2 = v2*h_1st_layer - Am%v(1,ie)
     inn = sqrt(d1**2 + d2**2)

  end select

  ! Compute the ensemble spread (std)
  !
  ens_std = 0.0d0
  do ne = 1,nrens
     ens_std = ens_std + (vmod_ens(ne) - vmod)**2
  end do
  ens_std = sqrt(ens_std/float(nrens-1))

  ! Formula from: Sakov, 2012 (topaz)
  !
  if ((KSTD > 0) .and. (abs(inn) > KSTD * ens_std)) then
     stdv_new = sqrt( sqrt( (ens_std**2 + stdv**2)**2 +&
                 (1/KSTD * ens_std * inn)**2 ) - ens_std**2 )

     write(*,*) 'innovation too high. Changing std:'
     write(*,'(a5,2x,4f8.4)') trim(ty),inn,ens_std,stdv,stdv_new

     stdv = stdv_new
  end if

  end subroutine check_obs_inn

  end module mod_enkf
