!
! Copyright (C) 2020, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
!------------------------------------------------------------------------------
! Creates an ensemble of files with perturbed BC from an initial BC file.
! Can handle both lateral and surface BC.
!------------------------------------------------------------------------------
program perturbeBC

  use m_sample2D
  implicit none

  character(len=6) :: arg1
  character(len=3) :: arg2
  character(len=3) :: arg3
  character(len=80) :: arg4
  character(len=3) :: arg5
  character(len=12) :: arg6
  character(len=12) :: arg7
  character(len=80) :: arg8

  integer nrens, var_dim
  character(len=80) :: filein
  character(len=3) :: filety
  character(len=80) :: variable
  real var_std
  double precision mem_time
  integer pert_type

  ! FEM files
  integer iunit,iformat
  integer nvers           !version of file format
  integer np              !size of data (horizontal, nodes or elements)
  integer nvar            !number of variables to write
  integer ntype           !type of information contained
  integer datetime(2)     !date and time information
  integer ierr            !return error code
  real*4 regpar(7)        !regular array params
  integer lmax            !vertical values
  real*4,allocatable :: hlv(:)     !vertical structure
  character(len=50),allocatable :: vstring(:)
  integer,allocatable :: ilhkv(:)
  real*4,allocatable :: hd(:)
  real x0,y0,dx,dy,flag
  integer nx,ny
  integer nlvddi
  !---

  integer           :: nrec_file

  real              :: var
  double precision  :: tnew,told,dtime
  character(len=80) :: dstring

  real, allocatable :: pvec0(:),pvec0_old(:)
  real, allocatable :: pvec1(:,:),pvec1_old(:,:)
  real, allocatable :: pmat(:,:,:,:),pmat_old(:,:,:,:)
  real, allocatable :: amat(:,:,:)

  ! 2D perturbation parameters
  real rx,ry,theta
  integer fmult
  logical samp_fix,verbose

  ! geostrophic vars
  real               :: flat
  logical, parameter :: bpress = .false.
  real, parameter    :: sigmaP = 100.	!press std in Pa

  ! SHYFEM single precision variables
  real*4, allocatable :: var2d(:,:,:),var2d_ens(:,:,:)

  integer n,i,ne
  integer fid

  ! read input
  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  call get_command_argument(4, arg4)
  call get_command_argument(5, arg5)
  call get_command_argument(6, arg6)
  call get_command_argument(7, arg7)
  call get_command_argument(8, arg8)

  if (trim(arg8) .eq. '') then
          write(*,*) ''
          write(*,*) 'Usage:'
          write(*,*) ''
          write(*,*) 'perturbeBC [nrens] [file_type] [pert_type] [variable] [var_dim] [var_std] [mem_time] [input]'
          write(*,*) ''
          write(*,*) '[nrens] is the n. of ens members, control included.'
          write(*,*) '[file_type] can be fem or ts (timeseries).'
	  write(*,*) '[pert_type] type of 2D perturbations. See below.'
          write(*,*) '[variable] name of the dynamic variable.'
          write(*,*) '[var_dim] spatial dimension of the variable (0->3).'
          write(*,*) '[var_std] standard deviation of the ensemble distribution (same units of the var).'
          write(*,*) '[mem_time] time correlation of the perturbations (red noise) in seconds.'
          write(*,*) '[input] is the name of the input unperturbed BC file.'
          write(*,*) ''
          write(*,*) 'pert_type can be:'
          write(*,*) ''
          write(*,*) '1- Spatially constant perturbation of each variable (no pressure).'
	  write(*,*) '2- 2D pseudo Gaussian perturbation of each variable (no pressure).'
          write(*,*) '3- 2D pseudo Gaussian geostrophic perturbation of wind and pressure.'
          write(*,*) '4- Same as 3, but using two perturbation fields. Pressure not perturbed.'
          write(*,*) '5- Only the wind speed is perturbed (same directions). Pressure not perturbed.'
	  write(*,*) '6- Honestly I do not remember.. but it is not working.'
          write(*,*) ''
          stop
  end if

  read(arg1,*) nrens
  filety = arg2
  read(arg3,*) pert_type
  variable = arg4
  read(arg5,*) var_dim
  read(arg6,*) var_std
  read(arg7,*) mem_time
  filein = arg8

  ! some checks on inputs
  if (( nrens < 3 ) .or. (mod(nrens,2) == 0) .or. (nrens > 1000)) error stop 'perturbeBC: bad nrens'
  if (( var_dim > 3 ) .or. ( var_dim < 0 )) error stop 'perturbeBC: bad dimension'
  if (( trim(filety) /= 'ts' ) .and. ( trim(filety) /= 'fem' )) error stop 'perturbeBC: bad file_type'
  if (mem_time <= 0) write(*,*) 'Warning: zero or negative mem_time. Setting to 0 (white noise).'
  if ( pert_type > 6 ) error stop 'perturbeBC: bad pert_type'

  ! open and close to read informations
  select case(trim(filety))
  case('ts')
     call read_ts(.true.,filein,nrec_file,var,tnew,dstring)
  case('fem')
     np = 0
     y0 = 0.
     x0 = 0.
     dx = 0.
     dy = 0.
     nx = 0
     ny = 0
     call fem_file_read_open(filein,np,iformat,iunit)
  end select

  iunit = 20
  call open_files(iunit,nrens,filein,filety)

  told = -999.
  ! time loop
  n = 0
  do

    n = n + 1

    ! read record
    select case(trim(filety))
    case('ts')
	if (n > nrec_file) exit
        call read_ts(.false.,filein,nrec_file,var,tnew,dstring)
    case('fem')

	! read 1st header
	call fem_file_read_params(iformat,iunit,dtime &
               ,nvers,np,lmax,nvar,ntype,datetime,ierr)
        if( ierr .lt. 0 ) exit
	call dts_convert_to_atime(datetime,dtime,tnew)

	! read 2nd header
	if (.not.allocated(hlv)) allocate(hlv(lmax))
        if (.not.allocated(ilhkv)) allocate(ilhkv(np))
	if (.not.allocated(hd)) allocate(hd(np))
	if (.not.allocated(vstring)) allocate(vstring(nvar))

        nlvddi = lmax
        call fem_file_read_2header(iformat,iunit,ntype,lmax &
             ,hlv,regpar,ierr)
        nx = nint(regpar(1))
        ny = nint(regpar(2))
	x0 = regpar(3)
	y0 = regpar(4)
        dx = regpar(5)
        dy = regpar(6)
        flag = regpar(7)

	! read variables
	do i = 1,nvar
	   if (.not.allocated(var2d)) allocate(var2d(nvar,nx,ny))
           call fem_file_read_data(iformat,iunit,nvers,np,lmax, &
		   vstring(i),ilhkv,hd,nlvddi,var2d(i,:,:),ierr)
	   if (n==1) write(*,*) 'Reading: ',vstring(i)
        end do
    end select

    ! generate perturbations and write the files
    select case(var_dim)
    case(0)	! 0D variable

	if (.not.allocated(pvec0)) allocate(pvec0(nrens-1))
	if (.not.allocated(pvec0_old)) allocate(pvec0_old(nrens-1))
        call perturbe_0d(nrens-1,pvec0)
        call red_noise_0d(told,pvec0_old,tnew,pvec0,nrens-1,mem_time)
        call write_record_0d(iunit,told,nrens,dstring,var,var_std,pvec0)

    case(1)	! 1D variable

        print*, 'todo'
        !call perturbe_1d

    case(2)	! 2D variable

	if (trim(filety) /= 'fem') error stop 'Bad file format'

	select case(pert_type)
	case(1)

            if (n==1) write(*,*) 'Case 1: spatially constant perturbations'

	    allocate(pvec0(nrens-1),pvec0_old(nrens-1))
	    do i=1,nvar
	       if (.not.allocated(pvec1)) allocate(pvec1(nvar,nrens-1))
	       if (.not.allocated(pvec1_old)) allocate(pvec1_old(nvar,nrens-1))
	       pvec0 = pvec1(i,:)
	       pvec0_old = pvec1_old(i,:)
               call perturbe_0d(nrens-1,pvec0)
               call red_noise_0d(told,pvec0_old,tnew,pvec0,nrens-1,mem_time)
	       pvec1(i,:) = pvec0
	       pvec1_old(i,:) = pvec0_old
	    end do
	    deallocate(pvec0,pvec0_old)

	    do ne=1,nrens
	      ! create member
	      if (.not.allocated(var2d_ens)) allocate(var2d_ens(nvar,nx,ny))
              call make_member(nvar,nrens,nx,ny,variable,ne,pvec1,var2d,var2d_ens,var_std)

	      ! write file
              fid = iunit + 10 + ne
              call fem_file_write_header(iformat,fid,dtime,nvers,np,lmax &
                     ,nvar,ntype,nlvddi,hlv,datetime,regpar)
              do i = 1,nvar
                 call fem_file_write_data(iformat,fid,nvers,np,lmax &
                       ,vstring(i),ilhkv,hd,nlvddi,var2d_ens(i,:,:))
              end do
	    end do

	case(2,3,4,5,6)

	    if (n==1) write(*,*) 'Case 2: 2D perturbations'
	    if (n==1) write(*,*) 'Warning: check pressure parameters inside the code'

	    ! set parameters
	    verbose = .true.
	    samp_fix = .true.
	    theta = 0.
	    fmult = 8
	    call set_decorrelation(nx,dx,rx)
	    call set_decorrelation(ny,dy,ry)
	    flat = y0 + (ny/2 * dy)

	    if (n==1) write(*,*) 'dx,dy,rx,ry,theta,samp_fix: ',dx,dy,rx,ry,theta,samp_fix

            ! Make the random fields
            allocate(amat(nx,ny,nrens-1))
	    do i=1,nvar
               if(.not.allocated(pmat)) allocate(pmat(nvar,nx,ny,nrens-1))
               call sample2D(amat,nx,ny,nrens-1,fmult,dx,dy,rx,ry,theta &
                    ,samp_fix,verbose)
               pmat(i,:,:,:) = amat
	    end do
            deallocate(amat)

            if(.not.allocated(pmat_old)) allocate(pmat_old(nvar,nx,ny,nrens-1))
            call red_noise_2d(nx,ny,nvar,told,pmat_old,tnew,pmat,nrens-1,mem_time)

	    do ne=1,nrens

	      if (.not.allocated(var2d_ens)) allocate(var2d_ens(nvar,nx,ny))

	      if (pert_type == 2) then
                 call make_member_2D(nvar,nrens,nx,ny,variable,ne,pmat,var2d,var2d_ens,var_std)

	      else if (pert_type == 3) then
	         if (trim(variable) /= 'wind') error stop 'Bad varable. Use just wind.'
	         call make_geo_field(nvar,ne,nrens,nx,ny,dx,dy,pmat,var2d,var2d_ens,flag,var_std, &
                                sigmaP,flat,bpress)

	      else if (pert_type == 4) then
	         if (trim(variable) /= 'wind') error stop 'Bad varable. Use just wind.'
	         call make_2geo_field(nvar,ne,nrens,nx,ny,dx,dy,pmat,var2d,var2d_ens,flag,var_std, &
                                sigmaP,flat)

	      else if (pert_type == 5) then
	         if (trim(variable) /= 'wind') error stop 'Bad varable. Use just wind.'
		 call make_ws_pert(nvar,ne,nrens,nx,ny,pmat,var2d,var2d_ens,flag,var_std)

	      else if (pert_type == 6) then
	         if (trim(variable) /= 'wind') error stop 'Bad varable. Use just wind.'
		 call make_ind_field_ws(nvar,ne,nrens,nx,ny,pmat,var2d,var2d_ens,flag,var_std)

	      end if

	      !if (trim(variable) == 'wind') call moderate_max_wind

	      ! write file
              fid = iunit + 10 + ne
              call fem_file_write_header(iformat,fid,dtime,nvers,np,lmax &
                     ,nvar,ntype,nlvddi,hlv,datetime,regpar)
              do i = 1,nvar
                 call fem_file_write_data(iformat,fid,nvers,np,lmax &
                       ,vstring(i),ilhkv,hd,nlvddi,var2d_ens(i,:,:))
              end do
	    end do

	end select

    case(3)	! 3D variable

        print*, 'todo'
        !call perturbe_3d

    end select

    told = tnew
  ! end loop
  end do

  ! close files
  call close_files(iunit,nrens)

end program perturbeBC

!------------------------------------------------------------------
! subroutines
!------------------------------------------------------------------

!-----------------------------------------------
  subroutine read_ts(linit,filein,kdim,v,vtime,dstring)
!-----------------------------------------------

  use iso8601
  implicit none

  logical, intent(in)          :: linit
  character(len=*),intent(in)  :: filein
  integer, intent(inout)       :: kdim
  real, intent(out)            :: v
  double precision, intent(out):: vtime
  character*80, intent(out) :: dstring
  integer ios
  integer ierr
  integer date, time
  integer k

  v = -999.
  vtime = -999.

  select case(linit)
  
     ! just check values and see the lenght
     !
     case (.true.)
          open(26,file=trim(filein), status = 'old', form = 'formatted', iostat = ios)
          if (ios /= 0) error stop 'read_ts: error opening file'

          k = 0
 90       read(26,*,end=100) dstring,v
          if (isnan(v)) error stop 'perturbeBC: input file contains nans'
          call string2date(trim(dstring),date,time,ierr)
          if (ierr /= 0) error stop "read_ts: error reading string"
          k = k + 1
          goto 90

 100      continue
          kdim = k
          rewind(26, iostat = ios)
          if (ios /= 0) error stop 'read_ts: error in file'
          return

     ! read the values
     !
     case (.false.)

          read(26,*,end=101) dstring,v
          call string2date(trim(dstring),date,time,ierr)
          if (ierr /= 0) error stop "read_ts: error reading string"
          call dts_to_abs_time(date,time,vtime)

  end select

  return

  101  close(26)

  end subroutine read_ts


!-----------------------------------------------
  subroutine open_files(iunit,nrens,fin,ftype)
!-----------------------------------------------
  implicit none

  integer, intent(in) :: nrens,iunit
  character(len=80), intent(in) :: fin
  character(len=3), intent(in) :: ftype

  integer n
  integer fid
  integer :: ppos
  character(len=80) :: bname,fname
  character(len=3) :: nlab

  ! find basename
  ppos = scan(trim(fin),".", BACK= .true.)
  if ( ppos > 0 ) bname = fin(1:ppos-1)
  
  if (trim(ftype) == 'ts') then
     do n = 1,nrens
        fid = iunit + 10 + n
        write(nlab,'(i3.3)') n-1
        fname = trim(bname)//'_'//nlab//'.dat'
        open(fid,file=fname,status='unknown')
     end do
  !else
  !   error stop 'open file error'
  end if

  end subroutine open_files

!-----------------------------------------------
  subroutine close_files(iunit,nrens)
!-----------------------------------------------
  implicit none
  integer, intent(in) :: iunit
  integer, intent(in) :: nrens
  integer n

  do n = 1,nrens
     close(iunit + 10 + n)
  end do

  end subroutine close_files

!-----------------------------------------------
  subroutine perturbe_0d(nrensp,pvec0)
!-----------------------------------------------
  use m_random
  implicit none
  integer, intent(in) :: nrensp
  real, intent(out) :: pvec0(nrensp)

  integer n
  real aaux,ave

  call random(pvec0,nrensp)

  ! remove outlayers
  do n = 1,nrensp
     aaux = pvec0(n)
     if( abs(aaux).ge.3. ) then
        aaux = sign(1.,aaux) * (abs(aaux)-floor(abs(aaux)) + 1.)
     end if
     pvec0(n) = aaux
  end do
  ! set mean eq to zero
  ave = sum(pvec0)/float(nrensp)
  pvec0 = pvec0 - ave

  end subroutine perturbe_0d

!-----------------------------------------------
  subroutine red_noise_0d(told,pvec0_old,tnew,pvec0,nrensp,tau)
!-----------------------------------------------
  implicit none
  double precision, intent(in) :: tnew,tau
  double precision, intent(inout) :: told
  integer, intent(in) :: nrensp
  real, intent(inout) :: pvec0(nrensp),pvec0_old(nrensp)
  integer n

  double precision alpha

  if (told < 0) return

  if (tau > 0) then
     alpha = 1. -  (tnew - told)/tau
  else
     alpha = 0.
  end if

  do n=1,nrensp
     pvec0(n) = alpha * pvec0_old(n) + sqrt(1 - alpha**2) * pvec0(n)
  end do

  told = tnew
  pvec0_old = pvec0
  
  end subroutine red_noise_0d

!-----------------------------------------------
  subroutine red_noise_2d(nx,ny,nvar,told,pmat_old,tnew,pmat,nrensp,tau)
!-----------------------------------------------
  implicit none
  integer, intent(in) :: nx,ny,nvar
  double precision, intent(in) :: tnew,tau
  double precision, intent(inout) :: told
  integer, intent(in) :: nrensp
  real, intent(inout) :: pmat(nvar,nx,ny,nrensp),pmat_old(nvar,nx,ny,nrensp)
  integer n

  double precision alpha

  if (told < 0) return

  if (tau > 0) then
     alpha = 1. -  (tnew - told)/tau
  else
     alpha = 0.
  end if

  pmat = alpha * pmat_old + sqrt(1 - alpha**2) * pmat

  told = tnew
  pmat_old = pmat
  
  end subroutine red_noise_2d


!-----------------------------------------------
  subroutine write_record_0d(iunit,told,nrens,dstring,var0,var_std,pvec0)
!-----------------------------------------------
  implicit none
  double precision, intent(in) :: told
  integer, intent(in) :: nrens,iunit
  character(len=80), intent(in) :: dstring
  real, intent(in) :: var0,var_std
  real, intent(in) :: pvec0(nrens-1)

  integer n,fid
  real var


  do n = 1,nrens
     fid = iunit + 10 + n

     if (n > 1) then
        var = (pvec0(n-1) * var_std) + var0
     else
        var = var0
     end if

     if (told < 0) var = var0

     write(fid,*) trim(dstring),var
  end do

  end subroutine write_record_0d

!-----------------------------------------------
  subroutine make_member(nvar,nrens,nx,ny,variable,ne,vec1,var2d,var2d_ens,var_std)
!-----------------------------------------------
  implicit none
  integer, intent(in) :: nvar,nrens,nx,ny,ne
  character(len=80), intent(in) :: variable
  real, intent(in) :: vec1(nvar,nrens-1),var_std
  real*4, intent(in) :: var2d(nvar,nx,ny)
  real*4, intent(out) :: var2d_ens(nvar,nx,ny)

  integer i
  real vec11(nvar,nrens)

  vec11(:,1) = 0.
  vec11(:,2:-1) = vec1

  ! does not perturbe the pressure (3rd variable)
  if (trim(variable) == 'wind') vec11(3,:) = 0.

  do i=1,nvar
     var2d_ens(i,:,:) = var2d(i,:,:) + var_std * vec11(i,ne)
  end do

  end subroutine make_member

!-----------------------------------------------
  subroutine make_member_2D(nvar,nrens,nx,ny,variable,ne,pmat,var2d,var2d_ens,var_std)
!-----------------------------------------------
  implicit none
  integer, intent(in) :: nvar,nrens,nx,ny,ne
  character(len=80), intent(in) :: variable
  real, intent(in) :: pmat(nvar,nx,ny,nrens-1),var_std
  real*4, intent(in) :: var2d(nvar,nx,ny)
  real*4, intent(out) :: var2d_ens(nvar,nx,ny)

  integer i
  real pmat1(nvar,nx,ny,nrens)

  pmat1(:,:,:,1) = 0.
  pmat1(:,:,:,2:-1) = pmat

  ! does not perturbe the pressure (3rd variable)
  if (trim(variable) == 'wind') pmat1(3,:,:,:) = 0.

  do i=1,nvar
     var2d_ens(i,:,:) = var2d(i,:,:) + var_std * pmat1(i,:,:,ne)
  end do

  end subroutine make_member_2D

!-----------------------------------------------
  subroutine set_decorrelation(ntot,delta,lrange)
!-----------------------------------------------
! For pseudo 2D fields, set rx,ry, the horizontal decorrelation length. Suppose to use geo degree
  implicit none
  real delta,lrange,lrange_low
  integer ntot
  real ltot
  integer ires
  integer, parameter :: nmin = 12	!minimum delta to resolve the range

  ! total length
  ltot = delta * (ntot-1)

  ! minimum range allowed
  lrange_low = delta * nmin

  if (lrange_low < 4. ) then
	  lrange = 4.
  else 
	  lrange = lrange_low
  end if

  end subroutine set_decorrelation


!--------------------------------------------------
	subroutine make_geo_field(nvar,iens,nrens,nx,ny,dx,dy, &
     		pmat,datain,dataout,flag,err,sigmaP,flat,bpress)
!--------------------------------------------------
	implicit none

	integer,intent(in) :: nvar,iens
	integer,intent(in) :: nrens,nx,ny
	real,intent(in) :: dx,dy
	real,intent(in) :: err,sigmaP,flag
	real,intent(in) :: flat
        logical,intent(in) :: bpress
	real,intent(in) :: pmat(nvar,nx,ny,nrens-1)
	real*4,intent(in) :: datain(nvar,nx,ny)
	real*4,intent(out) :: dataout(nvar,nx,ny)

	real :: mat(nvar,nx,ny,nrens)

	integer ix,iy,ivar
	real Fp1,Fp2,Up,Vp
	real sigmaU, sigmaV
	real dxm,dym

	real, parameter :: pi = acos(-1.)
	real, parameter :: rhoa = 1.2041
	real, parameter :: er1 = 63781370. !max earth radius
	real, parameter :: er2 = 63567523. !min earth radius
	real fcor,er,theta

        mat(:,:,:,1) = 0.
        mat(:,:,:,2:-1) = pmat
  
	theta = flat * pi/180.
	fcor = 2. * sin(theta) * (2.* pi / 86400.)
	! earth radius with latitude
	er = sqrt( ( (er1**2 * cos(theta))**2 + &
     			(er2**2 * sin(theta))**2 ) / &
     			( (er1 * cos(theta))**2 + &
     			(er2 * sin(theta))**2 ) ) 

	dxm = dx * pi/180. * er
	dym = dy * pi/180. * er

	do ivar = 1,nvar

          select case(ivar)
	      case(1)	!u-wind

		do ix = 1,nx
		do iy = 2,ny

		  Fp2 = mat(1,ix,iy,iens)
		  Fp1 = mat(1,ix,iy-1,iens)
		  ! if err is relative use this
		  !sigmaU = err * abs(datain(ivar,ix,iy))
		  sigmaU = err

		  !Up = - (((Fp2 - Fp1)/dym) * sigmaP) / (rhoa * fcor)
		  Up = - (Fp2 - Fp1) * sigmaU

		  !dataout(ivar,ix,iy) = Up
		  dataout(ivar,ix,iy) = datain(ivar,ix,iy) + Up

		  if (datain(ivar,ix,iy) == flag) dataout(ivar,ix,iy) = flag

		end do
		end do

		! First row unchanged
		dataout(ivar,:,1) = datain(ivar,:,1)
		
	    case(2)	!v-wind

		do iy = 1,ny
		do ix = 2,nx

		  Fp2 = mat(1,ix,iy,iens)
		  Fp1 = mat(1,ix-1,iy,iens)
		  ! if err is relative use this
		  !sigmaV = err * abs(datain(ivar,ix,iy))
		  sigmaV = err

		  !Vp = (((Fp2 - Fp1)/dxm) * sigmaP) / (rhoa * fcor)
		  Vp = (Fp2 - Fp1) * sigmaV

		  !dataout(ivar,ix,iy) =  Vp
		  dataout(ivar,ix,iy) = datain(ivar,ix,iy) + Vp

		  if (datain(ivar,ix,iy) == flag) &
     			dataout(ivar,ix,iy) = flag

		end do
		end do

		! First row unchanged
		dataout(ivar,1,:) = datain(ivar,1,:)

	    case(3)	!pressure

	      if (bpress) then
		do iy = 1,ny
		do ix = 1,nx
		     dataout(ivar,ix,iy) = datain(ivar,ix,iy) + sigmaP &
     					* mat(1,ix,iy,iens)
!		  dataout(ivar,ix,iy) = sigmaP * mat(1,ix,iy,iens)
		  if (datain(ivar,ix,iy) == flag) dataout(ivar,ix,iy) = flag
		end do
		end do
              else
		write(*,*) 'pressure not perturbed'
		dataout = datain ! no perturbation
              end if

	  end select
	end do
	
	end subroutine make_geo_field


!--------------------------------------------------
	subroutine make_2geo_field(nvar,iens,nrens,nx,ny,dx,dy, &
     		pmat,datain,dataout,flag,err,sigmaP,flat)
!--------------------------------------------------
	implicit none

	integer,intent(in) :: nvar,iens
	integer,intent(in) :: nrens,nx,ny
	real,intent(in) :: dx,dy
	real,intent(in) :: err,sigmaP,flag
	real,intent(in) :: flat
	real,intent(in) :: pmat(nvar,nx,ny,nrens-1)
	real*4,intent(in) :: datain(nvar,nx,ny)
	real*4,intent(out) :: dataout(nvar,nx,ny)

	real :: mat(nvar,nx,ny,nrens)

	integer ix,iy,ivar
	real F1p,F2p,Up,Vp
	real sigmaU, sigmaV
	real dxm,dym

	real, parameter :: pi = acos(-1.)
	real, parameter :: rhoa = 1.2041
	real, parameter :: er1 = 63781370. !max earth radius
	real, parameter :: er2 = 63567523. !min earth radius
	real fcor,er,theta

        mat(:,:,:,1) = 0.
        mat(:,:,:,2:-1) = pmat

	theta = flat * pi/180.
	fcor = 2. * sin(theta) * (2.* pi / 86400.)
	! earth radius with latitude
	er = sqrt( ( (er1**2 * cos(theta))**2 + &
     			(er2**2 * sin(theta))**2 ) / &
     			( (er1 * cos(theta))**2 + &
     			(er2 * sin(theta))**2 ) ) 

	dxm = dx * pi/180. * er
	dym = dy * pi/180. * er

	do ivar = 1,nvar

          select case(ivar)
	      case(1)	!u-wind

		do ix = 1,nx
		do iy = 2,ny

		  F1p = - (mat(1,ix,iy,iens) - mat(1,ix,iy-1,iens))
		  F2p = (mat(2,ix,iy,iens) - mat(2,ix,iy-1,iens))

		  ! if err is relative use this
		  !sigmaU = err * abs(datain(ivar,ix,iy))
		  sigmaU = err

		  !Up = - (((F1p2 - F1p1)/dym) * sigmaP) / (rhoa * fcor)
		  Up = (F1p + F2p) * sigmaU

		  !dataout(ivar,ix,iy) = Up
		  dataout(ivar,ix,iy) = datain(ivar,ix,iy) + Up

		  if( datain(ivar,ix,iy).eq.flag ) dataout(ivar,ix,iy) = flag

		end do
		end do

		! First row unchanged
		dataout(ivar,:,1) = datain(ivar,:,1)
		
	    case(2)	!v-wind

		do iy = 1,ny
		do ix = 2,nx

		  F1p = mat(1,ix,iy,iens) - mat(1,ix-1,iy,iens)
		  F2p = - (mat(2,ix,iy,iens) - mat(2,ix-1,iy,iens))
		  ! if err is relative use this
		  !sigmaV = err * abs(datain(ivar,ix,iy))
		  sigmaV = err

		  !Vp = (((F1p2 - F1p1)/dxm) * sigmaP) / (rhoa * fcor)
		  Vp = (F1p + F2p) * sigmaV

		  !dataout(ivar,ix,iy) =  Vp
		  dataout(ivar,ix,iy) = datain(ivar,ix,iy) + Vp

		  if( datain(ivar,ix,iy).eq.flag ) dataout(ivar,ix,iy) = flag

		end do
		end do

		! First row unchanged
		dataout(ivar,1,:) = datain(ivar,1,:)

	    case(3)	!pressure

		dataout(ivar,:,:) = datain(ivar,:,:) ! no perturbation

	  end select
	end do
	
	end subroutine make_2geo_field


!--------------------------------------------------
	subroutine make_ws_pert(nvar,iens,nrens,nx,ny,pmat, &
		datain,dataout,flag,err)
!--------------------------------------------------
	implicit none
        integer,intent(in) :: nvar,iens
        integer,intent(in) :: nrens,nx,ny
        real,intent(in) :: pmat(nvar,nx,ny,nrens-1)
        real*4,intent(in) :: datain(nvar,nx,ny)
        real*4,intent(out) :: dataout(nvar,nx,ny)
        real,intent(in) :: err,flag

	real :: mat(nvar,nx,ny,nrens)

	real wso(nx,ny),wse(nx,ny)
	integer ivar,ix,iy

        mat(:,:,:,1) = 0.
        mat(:,:,:,2:-1) = pmat

	wso = sqrt(datain(1,:,:)**2 + datain(2,:,:)**2)
	wse = wso + err * mat(1,:,:,iens)

	do ivar = 1,nvar
	   select case (ivar)
	     case default	!wind
		do ix = 1,nx
		do iy = 1,ny
	     	  dataout(ivar,ix,iy) = datain(ivar,ix,iy) + &
                                    wse(ix,iy)/wso(ix,iy)
		  if (datain(ivar,ix,iy) == flag) dataout(ivar,ix,iy) = flag
		end do
		end do
	     case (3)	!pressure
	       dataout(ivar,:,:) = datain(ivar,:,:)
	   end select
	end do


	end subroutine make_ws_pert


!--------------------------------------------------
	subroutine make_ind_field_ws(nvar,iens,nrens,nx,ny,pmat, &
     				datain,dataout,flag,err1)
!--------------------------------------------------
! The idea was to perturb u and to use a coefficient k around 1
! k = wsp/ws
! With these conditions the perturbation for v is:
! dv=(-b +/- sqrt(b**2 - 4.*c))/2.
! with b=2.*v and c=-(k**2 - 1.)*ws**2 + du**2 + 2.*u*du
! but the sqrt must be positive so:
! k>=k0==sqrt((du**2 + 2.*u*du -v**2)/ws**2 + 1.)
! Anyway there is something wrong... NaNs and other..
	implicit none

	integer,intent(in) :: nvar,iens
	integer,intent(in) :: nrens,nx,ny
	real,intent(in) :: err1,flag
	real,intent(in) :: pmat(nvar,nx,ny,nrens-1)
	real*4,intent(in) :: datain(nvar,nx,ny)
	real*4,intent(out) :: dataout(nvar,nx,ny)
        real, allocatable :: u(:,:),du(:,:),v(:,:),dv1(:,:),k(:,:), &
                            ws(:,:),uu(:,:),vv(:,:),b(:,:),c(:,:), &
                            dv2(:,:)
	real :: mat(nvar,nx,ny,nrens)
        real k0

        mat(:,:,:,1) = 0.
        mat(:,:,:,2:-1) = pmat

        if (nvar /= 3) error stop 'dimension error'

        allocate(u(nx,ny),du(nx,ny),v(nx,ny),dv1(nx,ny),k(nx,ny), &
                ws(nx,ny),uu(nx,ny),vv(nx,ny),b(nx,ny),c(nx,ny), &
                dv2(nx,ny))

        u = datain(1,:,:)
        v = datain(2,:,:)
        ws = sqrt(u**2 + v**2)
        du = mat(1,:,:,iens) * (err1 * abs(u))

        k0 = maxval(sqrt((du**2 + 2.*u*du -v**2)/ws**2 + 1.))
        k = k0 + mat(2,:,:,iens)/3. !Gaussian with av 1 and 0 and 2 as limits

        ! u component
        uu = flag
        where (u /= flag)
          uu = u + du
        end where
        dataout(1,:,:) = uu

	! v component
        b = 2. * v
        c = -(k**2 - 1.)*ws**2 + du**2 + 2.*u*du
        dv1 = (-b + sqrt(b**2 - 4.*c))/2.
        dv2 = (-b - sqrt(b**2 - 4.*c))/2.
        vv = flag
        where (v /= flag)
          vv = v + dv2
        end where
        dataout(2,:,:) = vv

	! pressure
	dataout(3,:,:) = datain(3,:,:)

	end subroutine make_ind_field_ws


!--------------------------------------------------
	subroutine check_wind(nx,ny,wind,ierr,flag)
!--------------------------------------------------
        implicit none
        integer,intent(in) :: nx,ny
        real,intent(in) :: wind(2,nx,ny)
        integer,intent(out) :: ierr
	real,intent(in) :: flag
        real,parameter :: wsmax = 45.
        real ws
        integer ix,iy

        ierr = 0
        do ix = 1,nx
        do iy = 1,ny
	   if ((wind(1,ix,iy) == flag) .or. &
              (wind(2,ix,iy) == flag)) cycle

           ws = sqrt(wind(1,ix,iy)**2 + wind(2,ix,iy)**2)
           if (ws > wsmax) then
              write(*,*) 'wind too high (ix,iy,ws): ',ix,iy,ws
              ierr = 1
           end if
        end do
        end do
              
	end subroutine check_wind


!--------------------------------------------------
	subroutine correct_wind(bcorr,ne,nvar,nx,ny,err,omet,emet,flag)
!--------------------------------------------------
        implicit none
	logical, intent(in) :: bcorr
	integer,intent(in) :: ne
        integer,intent(in) :: nx,ny,nvar
        real,intent(in) :: err
        real,intent(in) :: omet(nvar,nx,ny)
        real,intent(inout) :: emet(nvar,nx,ny)
	real,intent(in) :: flag
        real ws,u,v
        real ws0,u0,v0,wsmax
        integer ix,iy
	real wmax

	if (.not.bcorr) return

        wmax = 30.
        do ix = 1,nx
        do iy = 1,ny
           if ((emet(1,ix,iy) == flag) .or.(emet(2,ix,iy) == flag)) cycle

           ! Reduce ws of 5%
           emet(1,ix,iy) = 0.95 * emet(1,ix,iy)
           emet(2,ix,iy) = 0.95 * emet(2,ix,iy)

           u = emet(1,ix,iy)
           v = emet(2,ix,iy)
           ws = sqrt(u**2 + v**2)
           u0 = omet(1,ix,iy)
           v0 = omet(2,ix,iy)
           ws0 = sqrt(u0**2 + v0**2)
           wsmax = ws0 + err * ws0
	   if (wsmax > wmax) wsmax = wmax
           if (ws > wsmax) then
	      write(*,*) 'Limiting max wind speed: ',ws,wsmax
              emet(1,ix,iy) = u/ws * wsmax
              emet(2,ix,iy) = v/ws * wsmax
              !write(*,*) 'moderating wind. Ens member:',ne
              !write(*,*) 'u->um v->vm: ',u,emet(1,ix,iy),v,emet(2,ix,iy)
              !write(*,*) 'ws->wsm: ',ws,ws0,wsmax
           end if
        end do
        end do
              
	end subroutine correct_wind
