!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
	program make_ens_meteo
	use m_sample2D
	implicit none

	integer iunit, iformat
	integer ounit
	character(len=80) :: filein
	double precision atime,dtime,tt  !time stamp
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
	character*50,allocatable :: string(:)
	real flag
	integer,allocatable :: ilhkv(:)
	real*4,allocatable :: hd(:)
	integer nlvddi
	real*4,allocatable :: data(:,:),dataens(:,:)
        real,allocatable :: ens_met(:,:,:),ctrl_met(:,:,:)
	real sigmaUV,sigmaP,sigmaWS	!u-v relative error, std of pressure
	integer nrens
	integer irec
	integer i,nr

	integer pert_type	! type of perturbation:	1 = 0D spatially constant, no press
				!			2 = 2D u and v indipendently, no press
				!			3 = Pressure pert and geostr wind pert
				!			4 = Wind speed perturbation, no press
				!			5 = geostr pert in both senses
				!			6 = u perturbed and v perturbed according to ws
	logical bpress
	double precision tau_er	! e-folding time for the error memory
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
	logical bcorr

	! n. of ens members
	!
	nrens = 50

	! type of perturbed field
	!
	pert_type = 2

	! correction for extreme wind. Use it if pert_type = 2
	!
	bcorr = .true.
	!bcorr = .false.

	! relative error
	!
	!sigmaUV = 0.4	!40%
	sigmaUV = 0.5
	!sigmaWS = 0.3

	! false to remove pressure perturbation. Only if pert_type = 3
	!
	bpress = .false.

	! pressure standard deviation
	!(1mbar->100Pa) Used only for the pressure, see the routine
	!
	sigmaP = 100. * 2.

	! decorrelation e-folding time
	!
	tau_er = 2*86400.

	! Average latitude for the Coriolis factor. Used only with pert_type = 3
	!
	flat = 40.

	! 2d pseudo random fields params
	!
	fmult = 8
	theta = 0.
	!theta = 135	!prevalent wind direction in the Med
	rx = 5
	ry = 5
	verbose = .false.
	samp_fix = .true.     !keep true

	! input file
	!
	filein = 'zbound.fem'
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
 	  call fem_file_read_params(iformat,iunit,dtime &
               ,nvers,np,lmax,nvar,ntype,datetime,ierr)
	  if( ierr .lt. 0 ) exit

          call dts_convert_to_atime(datetime,dtime,atime)

          tt = atime
	  write(*,*) 'time ',tt

	  allocate(hlv(lmax))
	  nlvddi = lmax
	  call fem_file_read_2header(iformat,iunit,ntype,lmax &
	       ,hlv,regpar,ierr)

	  nx = nint(regpar(1))
	  ny = nint(regpar(2))
	  flag = regpar(7)
	  dx = regpar(5)
	  dy = regpar(6)

	  allocate(ilhkv(np),hd(np),data(nx,ny))

          !--------------------------------------------------
	  ! Makes 2D perturbed normalised fields and update them
          !--------------------------------------------------

	  select case(pert_type)

	    case(1)		!constant field

		write(*,*) 'Case 1: spatially constant perturbations'

		if(.not.allocated(pvec)) allocate(pvec(2,nrens))
		if(.not.allocated(pvec_old)) &
			allocate(pvec_old(2,nrens))

        	call make_random_0D(pvec(1,:),nrens)
	        call make_random_0D(pvec(2,:),nrens)

		call merge_old_vec(irec,tt,tt_old,tau_er, &
                       nrens,2,pvec,pvec_old)

	    case(2,5,6)		! u and v indipendent

		write(*,*) 'Case 2: two perturbations'

	        if(.not.allocated(pmat)) allocate(pmat(2,nx,ny,nrens))
	        if(.not.allocated(pmat_old)) &
     			allocate(pmat_old(2,nx,ny,nrens))
		allocate(amat(nx,ny,nrens))

		! dx must resolve rx, but you need enough nx.
		! The same for y.
		!xlength = (nx-1) * dx
		!ylength = (ny-1) * dy
		!rx = ((nx-1)/4) * dx
		!ry = ((ny-1)/4) * dy

		! Make the random fields
		call sample2D(amat,nx,ny,nrens,fmult,dx,dy,rx,ry,theta &
                    ,samp_fix,verbose)
		pmat(1,:,:,:) = amat
		call sample2D(amat,nx,ny,nrens,fmult,dx,dy,rx,ry,theta &
                    ,samp_fix,verbose)
		pmat(2,:,:,:) = amat

	        deallocate(amat)

		! Merge with the old fields in order to have time correlation
		call merge_old_field(irec,tt,tt_old,tau_er,nx,ny, &
     			2,nrens,pmat,pmat_old)

	    case(3,4)		! perturbed p, u and v with geostroph. pert.
		
		write(*,*) 'Case 3: one perturbation'

	        if(.not.allocated(pmat)) allocate(pmat(1,nx,ny,nrens))
	        if(.not.allocated(pmat_old)) &
     			allocate(pmat_old(1,nx,ny,nrens))
		allocate(amat(nx,ny,nrens))

		call sample2D(amat,nx,ny,nrens,fmult,dx,dy,rx,ry,theta &
                    ,samp_fix,verbose)
		pmat(1,:,:,:) = amat

		call merge_old_field(irec,tt,tt_old,tau_er,nx,ny, &
                       1,nrens,pmat,pmat_old)

	        deallocate(amat)

	  end select 


	  !--------------------------------------------------
	  ! Loop on variables
	  !--------------------------------------------------
          allocate(ctrl_met(nvar,nx,ny))
	  allocate(string(nvar))
 	  do i = 1,nvar

	     ! Reads variable
	     call fem_file_read_data(iformat,iunit &
                               ,nvers,np,lmax &
                               ,string(i) &
                               ,ilhkv,hd &
                               ,nlvddi,data &
                               ,ierr)
	     write(*,*) 'Reading: ',string(i)

	     ! store control meteo
	     ctrl_met(i,:,:) = data
	  
	  end do !------loop on variables----

          ! loop on ens members
	  do nr = 1,nrens

             allocate(ens_met(nvar,nx,ny))
	     ! Select perturbation type
             select case (pert_type)
	        case (1)
	  	  call make_const_field(nr,nvar,nrens,nx,ny,pvec, &
     				ctrl_met,ens_met,flag,sigmaUV)
		case (2)
		  call make_ind_field(nvar,nr,nrens,nx,ny,pmat, &
     		        ctrl_met,ens_met,flag,sigmaUV,sigmaP)

		case (3)
		  call make_geo_field(nvar,nr,nrens,nx,ny,dx,dy, &
     				pmat,ctrl_met,ens_met,flag,sigmaUV, &
     				sigmaP,flat,bpress)
		case (4)
		  call make_ws_pert(nvar,nr,nrens,nx,ny,pmat, &
                      ctrl_met,ens_met,flag,sigmaUV)

		case (5)
		  call make_2geo_field(nvar,nr,nrens,nx,ny,dx,dy, &
     				pmat,ctrl_met,ens_met,flag,sigmaUV, &
     				sigmaP,flat)
		case (6)
		  call make_ind_field_ws(nvar,nr,nrens,nx,ny,pmat, &
                               ctrl_met,ens_met,flag,sigmaUV,sigmaWS)
	     end select

             !call check_wind(nx,ny,nrens,ens_met,ierr,flag)
             !if (ierr /= 0) error stop 'wind too high'
             call correct_wind(bcorr,nr,nvar,nx,ny,sigmauv,ctrl_met, &
                              ens_met,flag)

	     ounit = iunit + 10 + nr
	     call fem_file_write_header(iformat,ounit,dtime &
                               ,nvers,np,lmax &
                               ,nvar,ntype &
                               ,nlvddi,hlv,datetime,regpar)


	     allocate(dataens(nx,ny))
	     do i = 1,nvar
	        dataens = ens_met(i,:,:) !write in real4
	        call fem_file_write_data(iformat,ounit &
                         ,nvers,np,lmax &
                         ,string(i) &
                         ,ilhkv,hd &
                         ,nlvddi,dataens)
             end do

             deallocate(ens_met,dataens)

          end do

	  deallocate(ilhkv)
	  deallocate(hd)
	  deallocate(data)
	  deallocate(hlv)
          deallocate(string)
          deallocate(ctrl_met)

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
	subroutine merge_old_field(irec,tt,tt_old,tau_er,nx,ny, &
                       nvar,nrens,pmat,pmat_old)
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
	subroutine merge_old_vec(irec,tt,tt_old,tau_er, &
                       nrens,nvar,pvec,pvec_old)
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
	subroutine make_const_field(iens,nvar,nrens,nx,ny,vec, &
     					datain,dataout,flag,err)
!--------------------------------------------------
	implicit none

	integer,intent(in) :: nvar,iens
	integer,intent(in) :: nrens,nx,ny
	real,intent(in) :: err,flag
	real,intent(in) :: vec(2,nrens)
	real,intent(in) :: datain(nvar,nx,ny)
	real,intent(out) :: dataout(nvar,nx,ny)

	integer ix,iy,ivar

        do ivar = 1,nvar

 	   select case (ivar)

	     case default

	       do iy = 1,ny
	       do ix = 1,nx
		! if err is relative use this
!		dataout(ivar,ix,iy) = datain(ivar,ix,iy) + &
!                      (err * abs(datain(ivar,ix,iy))) * vec(ivar,iens)
		dataout(ivar,ix,iy) = datain(ivar,ix,iy) + &
                      err * vec(ivar,iens)
		if (datain(ivar,ix,iy) == flag) &
                      dataout(ivar,ix,iy) = flag
	       end do
	       end do

	     case (3)

	       dataout(ivar,:,:) = datain(ivar,:,:)

	   end select

	end do

	end subroutine make_const_field


!--------------------------------------------------
	subroutine make_ind_field(nvar,iens,nrens,nx,ny,mat, &
     		        datain,dataout,flag,err,errp)
!--------------------------------------------------
	implicit none

	integer,intent(in) :: nvar,iens
	integer,intent(in) :: nrens,nx,ny
	real,intent(in) :: err,errp,flag
	real,intent(in) :: mat(2,nx,ny,nrens)
	real,intent(in) :: datain(nvar,nx,ny)
	real,intent(out) :: dataout(nvar,nx,ny)
        real ppert

	integer ix,iy,ivar

        do ivar = 1,nvar
	  select case (ivar)

	    case default

	      do iy = 1,ny
	      do ix = 1,nx
	        ! if err is relative use this
		dataout(ivar,ix,iy) = datain(ivar,ix,iy) + &
                 (err*abs(datain(ivar,ix,iy))) * mat(ivar,ix,iy,iens)
		! dataout(ivar,ix,iy) = datain(ivar,ix,iy) + &
!     				(err * mat(ivar,ix,iy,iens))
		if (datain(ivar,ix,iy) == flag) &
     			dataout(ivar,ix,iy) = flag
	      end do
	      end do

	    case (3)

	      dataout(ivar,:,:) = datain(ivar,:,:)
!	      ! pressure pert is a linear combination
!	      do iy = 1,ny
!	      do ix = 1,nx
!
!                 ppert = 1/2. * mat(1,ix,iy,iens) + &
!                        1/2. * mat(2,ix,iy,iens)
!
!	         dataout(ivar,ix,iy) = datain(ivar,ix,iy) + &
!                   ppert * errp
!
!		 if (datain(ivar,ix,iy) == flag) &
!     			dataout(ivar,ix,iy) = flag
!
!	      end do
!	      end do

	  end select

	end do
	!write(*,*) iens,datain(:,5,5),dataout(:,5,5)

	end subroutine make_ind_field



!--------------------------------------------------
	subroutine make_geo_field(nvar,iens,nrens,nx,ny,dx,dy, &
     		mat,datain,dataout,flag,err,sigmaP,flat,bpress)
!--------------------------------------------------
	implicit none

	integer,intent(in) :: nvar,iens
	integer,intent(in) :: nrens,nx,ny
	real,intent(in) :: dx,dy
	real,intent(in) :: err,sigmaP,flag
	real,intent(in) :: flat
        logical,intent(in) :: bpress
	real,intent(in) :: mat(1,nx,ny,nrens)
	real,intent(in) :: datain(nvar,nx,ny)
	real,intent(out) :: dataout(nvar,nx,ny)

	integer ix,iy,ivar
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
		  sigmaU = err * abs(datain(ivar,ix,iy))
		  !sigmaU = err

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
		  sigmaV = err * abs(datain(ivar,ix,iy))
		  !sigmaV = err

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
     		mat,datain,dataout,flag,err,sigmaP,flat)
!--------------------------------------------------
	implicit none

	integer,intent(in) :: nvar,iens
	integer,intent(in) :: nrens,nx,ny
	real,intent(in) :: dx,dy
	real,intent(in) :: err,sigmaP,flag
	real,intent(in) :: flat
	real,intent(in) :: mat(2,nx,ny,nrens)
	real,intent(in) :: datain(nvar,nx,ny)
	real,intent(out) :: dataout(nvar,nx,ny)

	integer ix,iy,ivar
	real F1p,F2p,Up,Vp
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
		  sigmaU = err * abs(datain(ivar,ix,iy))
		  !sigmaU = err

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
		  sigmaV = err * abs(datain(ivar,ix,iy))
		  !sigmaV = err

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
	subroutine make_ws_pert(nvar,iens,nrens,nx,ny,mat, &
		datain,dataout,flag,err)
!--------------------------------------------------
	implicit none
        integer,intent(in) :: nvar,iens
        integer,intent(in) :: nrens,nx,ny
        real,intent(in) :: mat(1,nx,ny,nrens)
        real,intent(in) :: datain(nvar,nx,ny)
        real,intent(out) :: dataout(nvar,nx,ny)
        real,intent(in) :: err,flag

	real wso(nx,ny),wse(nx,ny)
	integer ivar,ix,iy

	wso = sqrt(datain(1,:,:)**2 + datain(2,:,:)**2)
	wse = wso + err * wso * mat(1,:,:,iens)

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
	subroutine make_ind_field_ws(nvar,iens,nrens,nx,ny,mat, &
     				datain,dataout,flag,err1,err2)
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
	real,intent(in) :: err1,err2,flag
	real,intent(in) :: mat(2,nx,ny,nrens)
	real,intent(in) :: datain(nvar,nx,ny)
	real,intent(out) :: dataout(nvar,nx,ny)
        real, allocatable :: u(:,:),du(:,:),v(:,:),dv1(:,:),k(:,:), &
                            ws(:,:),uu(:,:),vv(:,:),b(:,:),c(:,:), &
                            dv2(:,:)
        real k0

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



