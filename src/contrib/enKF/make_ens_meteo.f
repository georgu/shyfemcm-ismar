	program make_ens_meteo
	use m_sample2D
	implicit none

	integer iunit, iformat
	integer ounit
	character(len=80) :: filein
	double precision tt  !time stamp
	double precision tt_old
	integer nvers           !version of file format
	integer np              !size of data (horizontal, nodes or elements)
	integer nvar            !number of variables to write
        integer ntype           !type of information contained
        integer datetime(2)     !date and time information
        integer ierr		!return error code
	real*4 regpar(7)          !regular array params
	integer lmax		!vertical values
	real*4,allocatable :: hlv(:)          !vertical structure
	character*50 string
	real flag
	integer,allocatable :: ilhkv(:)
	real*4,allocatable :: hd(:)
	integer nlvddi
	real*4,allocatable :: data(:,:),dataens(:,:)
        real,allocatable :: metaux(:,:,:,:)
	real sigmaUV,sigmaP	!u-v relative error, std of pressure
	integer nrens
	integer irec
	integer i,nr

	integer pert_type	! type of perturbation:	1 = 0D spatially constant
				!			2 = 2D u and v indipendently
				!			3 = Pressure pert and geostr wind pert
	logical bpress
	real tau_er	! e-folding time for the error memory
	integer nx,ny
	real dx,dy
	integer fmult	!Mult factor for the super-sampling
	real rx,ry	!x and y length scales for the perturbations
	real theta	!Rotation of the random fields (theta=0 is east, rotation anticlocwise)
	logical verbose,samp_fix
	real,allocatable :: pmat(:,:,:,:)
	real,allocatable :: pmat_old(:,:,:,:)
	real,allocatable :: pvec(:,:)
	real,allocatable :: pvec_old(:,:)
	real,allocatable :: amat(:,:,:)
	integer ix,iy
	real dval
	real flat	!Latitude for the Coriolis factor
        real wsmax	!maximum wind speed

	! n. of ens members
	!
	nrens = 30

	! type of perturbed field
	!
	pert_type = 2

	! relative error
	!
	sigmaUV = 0.6
	!sigmaUV = 3.

	! maximum wind speed. If more resize it
	!
	wsmax = 40.

	! false to remove pressure perturbation. Only if pert_type = 3
	!
	bpress = .true.

	! pressure standard deviation
	!(1mbar->100Pa) Used only for the pressure, see the routine
	!
	sigmaP = 100. * 2.

	! decorrelation e-folding time
	!
	tau_er = 86400*2

	! Average latitude for the Coriolis factor. Used only with pert_type = 3
	!
	flat = 40.

	! 2d pseudo random fields params
	!
	fmult = 8
	theta = 0.
	!theta = 135	!prevalent wind direction in the Med
	rx = 4.
	ry = 4.
	verbose = .false.
	samp_fix = .true.     !keep true

	! input file
	!
	filein = 'wind.fem'
	!filein = 'hfr.fem'


	np = 0
	call fem_file_read_open(filein,np,iformat,iunit)

        !--------------------------------------------------
	! Reads and writes records
        !--------------------------------------------------
	irec = 0
	do

	  irec = irec + 1
          !--------------------------------------------------
	  ! Reads headers
          !--------------------------------------------------
 	  call fem_file_read_params(iformat,iunit,tt
     +                          ,nvers,np,lmax,nvar,ntype,datetime,ierr)
	  if( ierr .lt. 0 ) exit
	  write(*,*) 'time ',tt
	  allocate(hlv(lmax))
	  nlvddi = lmax
	  call fem_file_read_2header(iformat,iunit,ntype,lmax
     +                  ,hlv,regpar,ierr)

	  nx = nint(regpar(1))
	  ny = nint(regpar(2))
	  flag = regpar(7)
	  dx = regpar(5)
	  dy = regpar(6)

	  allocate(ilhkv(np),hd(np),data(nx,ny))
	  allocate(dataens(nx,ny))

          !--------------------------------------------------
	  ! Makes 2D perturbed normalised fields and update them
          !--------------------------------------------------

	  select case(pert_type)

	    case(1)		!constant field

		write(*,*) 'Case 1: spatially constant perturbations'

		if(.not.allocated(pvec)) allocate(pvec(2,nrens))
		if(.not.allocated(pvec_old)) 
     +			allocate(pvec_old(2,nrens))

        	call make_random_0D(pvec(1,:),nrens)
	        call make_random_0D(pvec(2,:),nrens)

		call merge_old_vec(irec,tt,tt_old,tau_er,
     +                  nrens,2,pvec,pvec_old)

	    case(2)		! u and v indipendent

		write(*,*) 'Case 2: indipendent wind perturbations'

	        if(.not.allocated(pmat)) allocate(pmat(2,nx,ny,nrens))
	        if(.not.allocated(pmat_old))
     +			allocate(pmat_old(2,nx,ny,nrens))
	        if(.not.allocated(amat)) allocate(amat(nx,ny,nrens))

		! dx must resolve rx, but you need enough nx.
		! The same for y.
		!xlength = (nx-1) * dx
		!ylength = (ny-1) * dy
		!rx = ((nx-1)/4) * dx
		!ry = ((ny-1)/4) * dy

		! Make the random fields
		call sample2D(amat,nx,ny,nrens,fmult,dx,dy,rx,ry,theta
     +               ,samp_fix,verbose)
		pmat(1,:,:,:) = amat
		call sample2D(amat,nx,ny,nrens,fmult,dx,dy,rx,ry,theta
     +               ,samp_fix,verbose)
		pmat(2,:,:,:) = amat

		! Merge with the old fields in order to have time correlation
		call merge_old_field(irec,tt,tt_old,tau_er,nx,ny,
     +			2,nrens,pmat,pmat_old)

	    case(3)		! perturbed p, u and v with geostroph. pert.
		
		write(*,*) 'Case 3: geostrophic wind perturbations'

	        if(.not.allocated(pmat)) allocate(pmat(1,nx,ny,nrens))
	        if(.not.allocated(pmat_old))
     +			allocate(pmat_old(1,nx,ny,nrens))
	        if(.not.allocated(amat)) allocate(amat(nx,ny,nrens))

		call sample2D(amat,nx,ny,nrens,fmult,dx,dy,rx,ry,theta
     +               ,samp_fix,verbose)
		pmat(1,:,:,:) = amat

		call merge_old_field(irec,tt,tt_old,tau_er,nx,ny,
     +                  1,nrens,pmat,pmat_old)

	  end select 

          !--------------------------------------------------
	  ! Writes headers
          !--------------------------------------------------
	  do nr = 1,nrens

	     ounit = iunit + 10 + nr
	     call fem_file_write_header(iformat,ounit,tt
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype
     +                          ,nlvddi,hlv,datetime,regpar)

	  end do


	  !--------------------------------------------------
	  ! Loop on variables
	  !--------------------------------------------------
          allocate(metaux(3,nx,ny,nrens))
 	  do i = 1,nvar

	     ! Reads variable
	     call fem_file_read_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvddi,data
     +                          ,ierr)
	     write(*,*) 'Reading: ',string
	  
             ! loop on ens members
	     do nr = 1,nrens

		! Select perturbation type
                select case (pert_type)
		   case (1)
			call make_const_field(i,nr,2,nrens,nx,ny,pvec,
     +					data,dataens,flag,sigmaUV)
		   case (2)
			call make_ind_field(i,nr,2,nrens,nx,ny,pmat,
     +					data,dataens,flag,sigmaUV)
		   case (3)
			call make_geo_field(i,nr,1,nrens,nx,ny,dx,dy,
     +			pmat,data,dataens,flag,sigmaUV,sigmaP,flat,
     +                  bpress)
		end select

                metaux(i,:,:,nr) = dataens

	     end do !------loop on ens members----
	  end do !------loop on variables----

          ! Write record
	  do nr = 1,nrens

             !call check_wind(nx,ny,nrens,metaux,ierr)
             !if (ierr /= 0) error stop 'wind too high'
             call correct_wind(nx,ny,metaux(:,:,:,nr),wsmax)

	     ounit = iunit + 10 + nr
	     do i = 1,nvar
	        call fem_file_write_data(iformat,ounit
     +                    ,nvers,np,lmax
     +                    ,string
     +                    ,ilhkv,hd
     +                    ,nlvddi,metaux(i,:,:,nr))
             end do

          end do

	  deallocate(ilhkv)
	  deallocate(hd)
	  deallocate(data)
	  deallocate(dataens)
	  deallocate(hlv)
          deallocate(metaux)

	end do 	!-------loop on records-----

	! closes files
	close(iunit)
	do nr = 1,nrens
	   ounit = iunit + 10 + nr
	   close(ounit)
        end do

	end program make_ens_meteo

!***********************************************************
!***********************************************************
!***********************************************************
	
!--------------------------------------------------
	subroutine merge_old_field(irec,tt,tt_old,tau_er,nx,ny,
     +                  nvar,nrens,pmat,pmat_old)
!--------------------------------------------------
	implicit none
	integer irec
	double precision tt,tt_old,tau_er
	integer nx,ny,nrens,nvar
	real pmat(nvar,nx,ny,nrens),pmat_old(nvar,nx,ny,nrens)

	double precision dt_er,alpha

        if ( irec.gt.1 ) then
	  dt_er = tt - tt_old
	  alpha = 1. - (dt_er/tau_er)
	  pmat = alpha * pmat_old + sqrt(1 - alpha**2) * pmat
        end if

	! Update old fields
	pmat_old = pmat
	tt_old = tt

	end subroutine merge_old_field

!--------------------------------------------------
	subroutine merge_old_vec(irec,tt,tt_old,tau_er,
     +                  nrens,nvar,pvec,pvec_old)
!--------------------------------------------------
	implicit none
	integer irec
	double precision tt,tt_old,tau_er
	integer nx,ny,nrens,nvar
	real pvec(nvar,nrens),pvec_old(nvar,nrens)

	double precision dt_er,alpha

        if ( irec.gt.1 ) then
	  dt_er = tt - tt_old
	  alpha = 1. - (dt_er/tau_er)
	  pvec = alpha * pvec_old + sqrt(1 - alpha**2) * pvec
        end if

	! Update old fields
	pvec_old = pvec
	tt_old = tt

	end subroutine merge_old_vec


!--------------------------------------------------
	subroutine make_random_0D(vec,nvec)
!--------------------------------------------------
        use m_random
	implicit none
	integer, intent(in) :: nvec
	real, intent(out) :: vec(nvec)
	integer n
	real aaux,ave
	
        call random(vec,nvec)

        ! remove outlayers
        do n = 1,nvec
           aaux = vec(n)
           if( abs(aaux).ge.3. ) then
             aaux = aaux/abs(aaux) * (abs(aaux)-floor(abs(aaux)) + 1.)
           end if
           vec(n) = aaux
        end do
        ! set mean eq to zero
        ave = sum(vec)/float(nvec)
        vec = vec - ave

	end subroutine make_random_0D


!--------------------------------------------------
	subroutine make_const_field(ivar,iens,vardim,nrens,nx,ny,vec,
     +					data,dataens,flag,err)
!--------------------------------------------------
	implicit none

	integer,intent(in) :: ivar,iens
	integer,intent(in) :: vardim,nrens,nx,ny
	real,intent(in) :: err,flag
	real,intent(in) :: vec(vardim,nrens)
	real*4,intent(in) :: data(nx,ny)
	real*4,intent(out) :: dataens(nx,ny)

	integer ix,iy

	select case (ivar)

	  case default

	    do iy = 1,ny
	    do ix = 1,nx
! if err is relative use this
		dataens(ix,iy) = data(ix,iy) + (err * abs(data(ix,iy))) 
     +					* vec(ivar,iens)
!		dataens(ix,iy) = data(ix,iy) + (err * vec(ivar,iens))
		if( data(ix,iy).eq.flag ) dataens(ix,iy) = flag
	    end do
	    end do

	  case (3)

	    dataens = data

	end select

	end subroutine make_const_field


!--------------------------------------------------
	subroutine make_ind_field(ivar,iens,vardim,nrens,nx,ny,mat,
     +				data,dataens,flag,err)
!--------------------------------------------------
	implicit none

	integer,intent(in) :: ivar,iens
	integer,intent(in) :: vardim,nrens,nx,ny
	real,intent(in) :: err,flag
	real,intent(in) :: mat(vardim,nx,ny,nrens)
	real*4,intent(in) :: data(nx,ny)
	real*4,intent(out) :: dataens(nx,ny)

	integer ix,iy

	select case (ivar)

	  case default

	    do iy = 1,ny
	    do ix = 1,nx
! if err is relative use this
		dataens(ix,iy) = data(ix,iy) + (err*abs(data(ix,iy)))
     +					* mat(ivar,ix,iy,iens)
!		dataens(ix,iy) = data(ix,iy) + 
!     +				(err * mat(ivar,ix,iy,iens))
		if( data(ix,iy).eq.flag ) dataens(ix,iy) = flag
	    end do
	    end do

	  case (3)

	    dataens = data

	end select

	end subroutine make_ind_field



!--------------------------------------------------
	subroutine make_geo_field(ivar,iens,vardim,nrens,nx,ny,dx,dy,
     +			mat,data,dataens,flag,err,sigmaP,flat,bpress)
!--------------------------------------------------
	implicit none

	integer,intent(in) :: ivar,iens
	integer,intent(in) :: vardim,nrens,nx,ny
	real,intent(in) :: dx,dy
	real,intent(in) :: err,sigmaP,flag
	real,intent(in) :: flat
        logical,intent(in) :: bpress
	real,intent(in) :: mat(vardim,nx,ny,nrens)
	real*4,intent(in) :: data(nx,ny)
	real*4,intent(out) :: dataens(nx,ny)

	integer ix,iy
	real Fp1,Fp2,Up,Vp
	real sigmaU, sigmaV
	real dxm,dym

	real, parameter :: pi = acos(-1.)
	real, parameter :: rhoa = 1.2041
	real, parameter :: er1 = 63781370. !max earth radius
	real, parameter :: er2 = 63567523. !min earth radius
	real fcor,er,theta

	theta = flat * pi/180.
	fcor = 2. * sin(theta) * (2.* pi / 86400.)
	! earth radius with latitude
	er = sqrt( ( (er1**2 * cos(theta))**2 + 
     +			(er2**2 * sin(theta))**2 ) /
     +			( (er1 * cos(theta))**2 + 
     +			(er2 * sin(theta))**2 ) ) 

	dxm = dx * pi/180. * er
	dym = dy * pi/180. * er

        select case(ivar)
	  case(1)	!u-wind

		do ix = 1,nx
		do iy = 2,ny

		  Fp2 = mat(1,ix,iy,iens)
		  Fp1 = mat(1,ix,iy-1,iens)
		  ! if err is relative use this
		  sigmaU = err * abs(data(ix,iy))
		  !sigmaU = err

		  !Up = - (((Fp2 - Fp1)/dym) * sigmaP) / (rhoa * fcor)
		  Up = - (Fp2 - Fp1) * sigmaU

		  !dataens(ix,iy) = Up
		  dataens(ix,iy) = data(ix,iy) + Up

		  if( data(ix,iy).eq.flag ) dataens(ix,iy) = flag

		end do
		end do

		! First row unchanged
		dataens(:,1) = data(:,1)
		
	  case(2)	!v-wind

		do iy = 1,ny
		do ix = 2,nx

		  Fp2 = mat(1,ix,iy,iens)
		  Fp1 = mat(1,ix-1,iy,iens)
		  ! if err is relative use this
		  sigmaV = err * abs(data(ix,iy))
		  !sigmaV = err

		  !Vp = (((Fp2 - Fp1)/dxm) * sigmaP) / (rhoa * fcor)
		  Vp = (Fp2 - Fp1) * sigmaV

		  !dataens(ix,iy) =  Vp
		  dataens(ix,iy) = data(ix,iy) + Vp

		  if( data(ix,iy).eq.flag ) dataens(ix,iy) = flag

		end do
		end do

		! First row unchanged
		dataens(1,:) = data(1,:)

	  case(3)	!pressure

	      if (bpress) then
		do iy = 1,ny
		do ix = 1,nx
		     dataens(ix,iy) = data(ix,iy) + sigmaP *
     +					mat(1,ix,iy,iens)
!		  dataens(ix,iy) = sigmaP * mat(1,ix,iy,iens)
		  if (data(ix,iy) == flag) dataens(ix,iy) = flag
		end do
		end do
              else
		write(*,*) 'pressure not perturbed'
		dataens = data ! no perturbation
              end if

	end select
	
	end subroutine make_geo_field


!--------------------------------------------------
	subroutine check_wind(nx,ny,wind,ierr)
!--------------------------------------------------
        implicit none
        integer,intent(in) :: nx,ny
        real,intent(in) :: wind(2,nx,ny)
        integer,intent(out) :: ierr
        real,parameter :: wsmax = 45.
        real ws
        integer ix,iy

        ierr = 0
        do ix = 1,nx
        do iy = 1,ny
           ws = sqrt(wind(1,ix,iy)**2 + wind(2,ix,iy)**2)
           if (ws > wsmax) then
              write(*,*) 'wind too high (ix,iy,ws): ',ix,iy,ws
              ierr = 1
           end if
        end do
        end do
              
	end subroutine check_wind


!--------------------------------------------------
	subroutine correct_wind(nx,ny,met,wsmax)
!--------------------------------------------------
        implicit none
        integer,intent(in) :: nx,ny
        real,intent(in) :: wsmax
        real,intent(inout) :: met(3,nx,ny)
        real ws,u,v
        integer ix,iy

        do ix = 1,nx
        do iy = 1,ny
           u = met(1,ix,iy)
           v = met(2,ix,iy)
           ws = sqrt(u**2 + v**2)
           if (ws > wsmax) then
              met(1,ix,iy) = u/ws * wsmax
              met(2,ix,iy) = v/ws * wsmax
              write(*,*) 'moderating wind u: ',u,met(1,ix,iy)
              write(*,*) 'moderating wind v: ',v,met(2,ix,iy)
           end if
        end do
        end do
              
	end subroutine correct_wind
