!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
program make_eof_ens_meteo

 implicit none

 integer i,irec,nrec

 integer iunit, iformat
 integer ounit
 character(len=80) :: filein
 double precision, allocatable :: tt(:)  !time stamp
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
 character*80,allocatable :: string(:)
 real flag
 integer,allocatable :: ilhkv(:)
 real*4,allocatable :: hd(:)
 integer nlvddi
 real*4,allocatable :: idata(:,:)

 integer nx,ny
 real dx,dy

 real, allocatable :: ddata(:,:)
 real, allocatable :: A(:,:)
 real, allocatable :: press(:,:)
 real, allocatable :: av(:)

 integer nrens,nr
 integer ix,iy,k

 character jobu,jobvt
 integer m,n
 integer mnmin,mnmax
 integer lda,ldu,ldvt
 integer lwork
 real, allocatable :: work(:), L(:)
 real, allocatable :: U(:,:),VT(:,:)
 real, allocatable :: PC(:,:),COV(:,:),UL(:,:)
 real, allocatable :: matA(:,:),matB(:,:)
 real Ltot,Lnorm
 integer nmode
 real, allocatable :: pert(:,:),vec(:)
 real, allocatable :: PCe(:,:),VTe(:,:)

 real perc,shfact

 !----------------------------------
 ! number of records, set this!
 !
 nrec = 120

 ! n. of ens members
 !
 nrens = 30
 nx = 0
 ny = 0

 ! Minimum percentage of variance to keep a PC component
 ! in the creation of the ensemble members
 !
 !perc = 1.

 ! shrinking factor for the distribution of random numbers
 ! use 3 at least, in order to have positive numbers.
 !
 shfact = 3.
 !----------------------------------

 ! input file
 !
 filein = 'wind.fem'
 !filein = 'hfr.fem'

 np = 0
 call fem_file_read_open(filein,np,iformat,iunit)

 !----------------------------------
 ! read wind fields
 !
 allocate(tt(nrec))
 do irec = 1,nrec

   call fem_file_read_params(iformat,iunit,tt(irec),nvers,np,&
                             lmax,nvar,ntype,datetime,ierr)
 
   if( ierr .lt. 0 ) exit
   write(*,*) 'time ',tt(irec) 
   if (.not. allocated(hlv)) allocate(hlv(lmax))
   nlvddi = lmax
 
   call fem_file_read_2header(iformat,iunit,ntype,lmax,&
                               hlv,regpar,ierr)
   nx = nint(regpar(1))
   ny = nint(regpar(2))
   flag = regpar(7)
   dx = regpar(5)
   dy = regpar(6)
   
   ! Allocate
   !
   if (.not. allocated(ilhkv)) allocate(ilhkv(np))
   if (.not. allocated(hd)) allocate(hd(np))
   if (.not. allocated(string)) allocate(string(nvar))
   if (.not. allocated(A)) allocate(A(nrec,2*nx*ny))
   if (.not. allocated(press)) allocate(press(nrec,nx*ny))
   allocate(idata(nx,ny))
   allocate(ddata(nx,ny))
 
   ! Read variables
   do i = 1,nvar
 
      call fem_file_read_data(iformat,iunit,nvers,np,lmax,&
                      string(i),ilhkv,hd,nlvddi,idata,ierr)
      ddata = idata
 
      write(*,*) 'Reading: ',string(i)
      
      select case (i)
        case (1)
	   k = 0
           do ix = 1,nx
           do iy = 1,ny
              k = k + 1
              A(irec,k) = ddata(ix,iy)
           end do
           end do
        case (2)
	   k = 0
           do ix = 1,nx
           do iy = 1,ny
              k = k + 1
              A(irec,nx*ny + k) = ddata(ix,iy)
           end do
           end do
        case (3)
	   k = 0
           do ix = 1,nx
           do iy = 1,ny
              k = k + 1
              press(irec,k) = ddata(ix,iy)    
           end do
           end do
      end select
 
   end do

   deallocate(idata)
   deallocate(ddata)

 end do
 close(iunit)

 !----------------------------------
 ! SVD
 ! some math: A = ULV^t -> COV = AA^t/(m-1) = 
 ! = UL^2U^t/(m-1)
 ! Y = U1^tA = U1^t(ULV^t) = (U1^tU)L V^t = L V^t
 !----------------------------------

 !----------------------------------
 ! detrend data
 !
 n = 2*nx*ny
 allocate(av(n))
 do k = 1,n
    av(k) = sum(A(:,k))/float(nrec)
    A(:,k) = A(:,k) - av(k)
 end do

 !----------------------------------
 ! settings
 !
 jobu = 'S'
 jobvt = 'S'
 m = nrec
 n = 2*nx*ny
 mnmin = min(m,n)
 mnmax = max(m,n)
 lda = m
 ldu = m
 
 allocate(L(mnmin))
 L = 0.
 if (jobu == 'A') then
    allocate(U(ldu,m)) 
    ldvt = n
 else if (jobu == 'S') then
    allocate(U(ldu,mnmin))
    ldvt = mnmin
 else
    error stop
 end if
 allocate(VT(ldvt,n)) 
 lwork = 2*max(3*mnmin+mnmax,5*mnmin)
 allocate(work(lwork))
 
 !----------------------------------
 ! SVD
 !
 call dgesvd(jobu, jobvt, m, n, A, lda, L, U,&
             ldu, VT, ldvt, work, lwork, ierr)

 if (ierr > 0) error&
     stop 'The algorithm computing SVD failed to converge.'

 !----------------------------------
 ! Principal components and covariance
 !
 allocate(PC(ldu,mnmin))
 allocate(UL(ldu,mnmin))
 allocate(COV(mnmin,mnmin))
 do ix = 1,ldu
 do iy = 1,mnmin
    PC(ix,iy) = U(ix,iy) * L(iy)
    UL(ix,iy) = U(ix,iy) * L(iy)**2
 end do
 end do

 !----------------------------------
 ! Compute the variance of the modes
 !
 Ltot = sum(L(:)**2)
 do k = 1,m
    Lnorm = L(k)**2/Ltot * 100.
    write(*,'(a,i3,2x,f6.3)') 'mode n, var(%): ',k,Lnorm
    ! set nmode if variance is less than 1%
    !if (Lnorm < perc) nmode = k
 end do
 nmode = m 

 ! covariance: UL U
 !
 call dgemm('n','t',ldu,mnmin,ldu,1.0,UL,ldu,U,ldu,0.0,COV,mnmin)
 COV = COV/(mnmin-1)


 !----------------------------------
 !----------------------------------
 
 !----------------------------------
 ! Make random numbers
 !
 write(*,*) 'Using ',nmode,' PCs'
 allocate(pert(nrens,nmode),vec(nrens))
 do k = 1,nmode
    call make_random_0D(vec,nrens) !gaussian, 0 mean, 1 std
    pert(:,k) = vec
 end do

 allocate(VTe(nmode,n))
 VTe = VT(1:nmode,:)
 !----------------------------------
 ! Loop on ensemble members
 !
 do nr = 1,nrens

   !----------------------------------
   ! make PCs
   if (.not. allocated(PCe)) allocate(PCe(ldu,nmode))
   do k = 1,nmode
      PCe(:,k) = PC(:,k) * ((pert(nr,k)/shfact) + 1)
   end do

   !----------------------------------
   ! rebuild wind A = PC VT
   !
   call dgemm('n','n',ldu,n,nmode,1.0,PCe,ldu,VTe,nmode,0.0,A,lda)
 
   !----------------------------------
   ! save just one mode
   !
   !number of mode
   !nmode = 1
   !allocate(matA(m,1),matB(1,n))
   !matA(:,1) = PC(:,nmode)
   !matB(1,:) = VT(nmode,:)
   !call dgemm('n','n',ldu,n,1,1.0,matA,ldu,matB,1,0.0,A,lda)
   !write(*,*) 'PCs ',matA
 

   ! add the average value
   do k = 1,n
      A(:,k) = A(:,k) + av(k)
   end do


   !----------------------------------
   ! write output
   !
   ounit = iunit + 10 + nr
   do irec = 1,nrec

      call fem_file_write_header(iformat,ounit,tt(irec),nvers,np,&
                  lmax,nvar,ntype,nlvddi,hlv,datetime,regpar)
  
      ! U
      k = 0
      allocate(idata(nx,ny))
      do ix = 1,nx
      do iy = 1,ny
         k = k + 1
         idata(ix,iy) = A(irec,k)
      end do
      end do
      call fem_file_write_data(iformat,ounit,nvers,np,lmax,&
                        string(1),ilhkv,hd,nlvddi,idata)
      deallocate(idata)

      ! V
      k = 0
      allocate(idata(nx,ny))
      do ix = 1,nx
      do iy = 1,ny
         k = k + 1
         idata(ix,iy) = A(irec,nx*ny + k)
      end do
      end do
      call fem_file_write_data(iformat,ounit,nvers,np,lmax,&
                        string(2),ilhkv,hd,nlvddi,idata)
      deallocate(idata)

      ! pressure
      k = 0
      allocate(idata(nx,ny))
      do ix = 1,nx
      do iy = 1,ny
         k = k + 1
         idata(ix,iy) = press(irec,k)
      end do
      end do
      call fem_file_write_data(iformat,ounit,nvers,np,lmax,&
                        string(3),ilhkv,hd,nlvddi,idata)
      deallocate(idata)
   end do

   close(ounit)

 end do

end program make_eof_ens_meteo


!***********************************************************
!***********************************************************
!***********************************************************


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
