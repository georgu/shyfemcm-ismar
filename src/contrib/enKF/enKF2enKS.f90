! ======================================================================
!  MODULE: mod_enks_data  (used internally)
! ======================================================================
module mod_enks_data
   use iso_fortran_env, only: dp=>real64
   implicit none
   save

   type t_x5_record
      real(dp)              :: tt
      real(dp), allocatable :: mat(:,:)
   end type t_x5_record

   type(t_x5_record), allocatable :: all_x5(:)
   integer :: total_x5_records = 0
   real(dp) :: dt_x5 = -1.0_dp

contains

! ======================================================================
! Sort X5 timestamps (selection sort, deterministic)

subroutine check_sort_x5()
   implicit none
   integer :: i, j, kmin
   type(t_x5_record) :: tmp

   if (total_x5_records <= 1) return

   ! Quick check
   do i=1,total_x5_records-1
      if (all_x5(i+1)%tt < all_x5(i)%tt) exit
   end do
   if (i==total_x5_records) return   ! already sorted

   ! Selection sort
   do i=1,total_x5_records-1
      kmin = i
      do j=i+1,total_x5_records
         if (all_x5(j)%tt < all_x5(kmin)%tt) kmin=j
      end do
      if (kmin /= i) then
         tmp = all_x5(i)
         all_x5(i) = all_x5(kmin)
         all_x5(kmin) = tmp
      end if
   end do
end subroutine check_sort_x5

! ======================================================================
! Estimate dt_x5

real(dp) function estimate_dt_x5()
   implicit none
   integer :: i
   if (total_x5_records < 2) then
      estimate_dt_x5 = 0.0_dp
      return
   end if

   do i=2,total_x5_records
      if (abs(all_x5(i)%tt - all_x5(i-1)%tt) > 1.0e-8_dp) then
         estimate_dt_x5 = all_x5(i)%tt - all_x5(i-1)%tt
         return
      end if
   end do

   estimate_dt_x5 = 0.0_dp
end function estimate_dt_x5

! ======================================================================
! float comparison

logical function equal_time(t1,t2)
   real(dp), intent(in) :: t1,t2
   equal_time = abs(t1-t2) < 1.0e-6_dp
end function equal_time

! ======================================================================
subroutine allocate_states(A, Amean, Astd, n, nrens)
   use mod_hydro
   use mod_hydro_vel
   use mod_ts
   use mod_restart, only : ibarcl_rst
   implicit none

   real(dp), allocatable, intent(out) :: A(:,:), Amean(:), Astd(:)
   integer, intent(out)               :: n
   integer, intent(in)                :: nrens

   n = size(utlnv) + size(vtlnv) + size(znv)
   if (ibarcl_rst /= 0) n = n + size(tempv) + size(saltv)

   allocate(A(n,nrens))
   allocate(Amean(n), Astd(n))
end subroutine allocate_states

end module mod_enks_data

! ======================================================================
! SUBROUTINES
! ======================================================================

! ======================================================================
subroutine load_all_x5(nrens)
   use iso_fortran_env, only : dp=>real64
   use mod_enks_data
   implicit none

   integer, intent(in) :: nrens

   integer :: ios, k, nren_tmp, nrobs_tmp
   real(dp) :: tt_tmp
   character(len=6) :: label_tmp
   character(len=2) :: tag_tmp

   open(15, file="X5_tot.uf", status="old", form="unformatted", iostat=ios)
   if (ios /= 0) error stop "Cannot open X5_tot.uf"

   ! First pass: count X5
   total_x5_records = 0
   do
      read(15, iostat=ios) tt_tmp, label_tmp, tag_tmp
      if (ios /= 0) exit
      if (tag_tmp == 'X5') then
         read(15) nren_tmp
         read(15)
         total_x5_records = total_x5_records + 1
      else
         if (tag_tmp == 'X3') then
            read(15) nren_tmp, nrobs_tmp
            read(15)
         else
            read(15) nren_tmp
            read(15)
         end if
      end if
   end do

   allocate(all_x5(total_x5_records))
   rewind(15)

   ! Second pass: load actual matrices
   k = 0
   do
      read(15, iostat=ios) tt_tmp, label_tmp, tag_tmp
      if (ios /= 0) exit
      if (tag_tmp == 'X5') then
         read(15) nren_tmp
         if (nren_tmp /= nrens) error stop "X5 size mismatch"
         k = k + 1
         all_x5(k)%tt = tt_tmp
         allocate(all_x5(k)%mat(nrens, nrens))
         read(15) all_x5(k)%mat
      else
         if (tag_tmp == 'X3') then
            read(15) nren_tmp, nrobs_tmp
            read(15)
         else
            read(15) nren_tmp
            read(15)
         end if
      end if
   end do
   close(15)

   call check_sort_x5()
   dt_x5 = estimate_dt_x5()
   if (dt_x5 <= 0.0_dp) then
      write (*,*) 'Analysis timestep negative. Setting to 3600s.'
      dt_x5 = 3600.0_dp
   end if

   write(*,*) "Loaded X5 records:", total_x5_records
   write(*,*) "Estimated dt_x5 =", dt_x5

end subroutine load_all_x5

! ======================================================================
subroutine init_shyfem(basinf, nnlv)
   use basin
   use levels
   use mod_geom_dynamic
   use mod_hydro
   use mod_hydro_vel
   use mod_ts
   use mod_restart
   use shympi
   implicit none
      
   character(len=*), intent(in) :: basinf
   integer, intent(in) :: nnlv
   integer :: ios
      
   open(21,file=basinf,status="old",form="unformatted",iostat=ios)
   if (ios/=0) error stop "Cannot open basin"
   call basin_read_by_unit(21)
   close(21)
      
   nlv = nnlv
   nlvdi = nnlv

   call mod_geom_dynamic_init(nkn,nel)
   call mod_hydro_init(nkn,nel,nlv)
   call mod_hydro_vel_init(nkn,nel,nlv)
   call mod_ts_init(nkn,nlv)
   call levels_init(nkn,nel,nlv)

   nkn_global=nkn
   nel_global=nel
   nlv_global=nlv
end subroutine init_shyfem

! ======================================================================
subroutine num2str(num, str)
   implicit none
   integer, intent(in)        :: num
   character(len=5), intent(out) :: str

   if (num < 0 .or. num > 99999) then 
      error stop "num2str: argument out of range"
   end if   
            
   write(str, '(I5.5)') num
end subroutine num2str

! ======================================================================
subroutine rst_read(rstname, atimea)
   use iso_fortran_env, only : dp=>real64
   use mod_restart
   use levels, only : nlvdi, nlv, hlv, ilhv, ilhkv
   use shympi
   use mod_enks_data, only : equal_time
   implicit none

   character(len=*), intent(in) :: rstname
   real(dp), intent(in) :: atimea
      
   integer :: ierr, iflag, ios
   real(dp) :: atimef
   integer, save :: icall = 0
   real*4 :: ibarcl4, iconz4, imerc4, iturb4, iwvert_rst4, ieco_rst4, zero4
   zero4 = 0.0
      
   open(24,file=trim(rstname),status='old',form='unformatted',action='read',iostat=ios)
   if (ios /= 0) error stop "rst_read: cannot open"
      
   do 
      call rst_read_record(24, atimef, iflag, ierr)
      if (ierr /= 0) then
         close(24)
         error stop "rst_read: time not found"
      end if
      if (equal_time(atimef, atimea)) exit
   end do
   close(24)

   if (icall==0) then
      hlv        = hlvrst
      hlv_global = hlvrst
      ilhv       = ilhrst
      ilhkv      = ilhkrst

      ibarcl4     = ibarcl_rst
      iwvert_rst4 = iwvert_rst
      ieco_rst4   = ieco_rst
      iconz4      = iconz_rst
      imerc4      = imerc_rst
      iturb4      = iturb_rst

      call addpar("ibarcl",ibarcl4)
      call addpar("iconz", iconz4)
      call addpar("iwvert",iwvert_rst4)
      call addpar("ieco"  ,ieco_rst4)
      call addpar("ibio"  ,zero4)
      call addpar("ibfm"  ,zero4)
      call addpar("imerc" ,imerc4)
      call addpar("iturb" ,iturb4)

      call daddpar("date",0.0_dp)
      call daddpar("time",0.0_dp)
      
      write(*,*) 'SHYFEM flags from restart:'
      write(*,*) 'nvers = ', nvers_rst
      write(*,*) 'nvmax = ', nvmax
      write(*,*) 'ibarcl = ', ibarcl_rst
      write(*,*) 'iconz  = ', iconz_rst
      write(*,*) 'iwvert = ', iwvert_rst
      write(*,*) 'ieco   = ', ieco_rst
      write(*,*) 'imerc  = ', imerc_rst
      write(*,*) 'iturb  = ', iturb_rst
      write(*,*) 'nlv   = ', nlv
      !write(*,*) 'hlvrst = ', hlvrst(1:nlv)
      write(*,*) 'hlv    = ', hlv(1:nlv)

   end if

   icall = icall + 1
end subroutine rst_read

! ======================================================================
subroutine push_matrix(sdim, nrens, nre, Amat)
   use iso_fortran_env, only: dp=>real64
   use mod_hydro
   use mod_hydro_vel
   use mod_ts
   use mod_restart, only : ibarcl_rst
   implicit none
      
   integer, intent(in) :: sdim, nrens, nre
   real(dp), intent(inout) :: Amat(sdim,nrens)
   integer :: dimuv, dimz, dimts
      
   dimz  = size(znv)
   dimuv = size(utlnv)
   dimts = size(znv)*size(utlnv,1)
      
   Amat(1:dimuv, nre) = reshape(real(utlnv,dp), [dimuv])
   Amat(dimuv+1:2*dimuv, nre) = reshape(real(vtlnv,dp), [dimuv])
   Amat(2*dimuv+1:2*dimuv+dimz, nre) = real(znv,dp)

   if (ibarcl_rst /= 0) then
      Amat(2*dimuv+dimz+1 : 2*dimuv+dimz+dimts, nre) = reshape(real(tempv,dp), [dimts])
      Amat(2*dimuv+dimz+dimts+1 : 2*dimuv+dimz+2*dimts, nre) = reshape(real(saltv,dp), [dimts])
   end if
end subroutine push_matrix

! ======================================================================
! pull_matrix 
subroutine pull_matrix(sdim, nrens, nre, Amat)
   use iso_fortran_env, only: dp=>real64
   use mod_hydro
   use mod_hydro_vel
   use mod_ts
   use mod_restart, only : ibarcl_rst
   implicit none
      
   integer, intent(in) :: sdim, nrens, nre
   real(dp), intent(in) :: Amat(sdim,nrens)
   integer :: dimuv, dimz, dimts
   integer :: nnkn, nnel, nnlv
      
   dimz  = size(znv)
   dimuv = size(utlnv)
   dimts = size(znv)*size(utlnv,1)
   nnkn  = size(znv)
   nnel  = size(utlnv,2)
   nnlv  = size(utlnv,1)

   utlnv = reshape(real(Amat(1:dimuv, nre)), [nnlv, nnel])
   vtlnv = reshape(real(Amat(dimuv+1:2*dimuv, nre)), [nnlv, nnel])
   znv   = real(Amat(2*dimuv+1:2*dimuv+dimz, nre))
   
   if (ibarcl_rst /= 0) then
      tempv = reshape(real(Amat(2*dimuv+dimz+1:2*dimuv+dimz+dimts,nre)), [nnlv, nnkn])
      saltv = reshape(real(Amat(2*dimuv+dimz+dimts+1:2*dimuv+dimz+2*dimts,nre)), [nnlv, nnkn])
   end if
end subroutine pull_matrix

! ======================================================================
! Lagged smoother
subroutine make_analysis(atime, sdim, nrens, Amat, nlag)
   use iso_fortran_env, only : dp=>real64
   use mod_enks_data
   implicit none

   real(dp), intent(in)    :: atime
   integer, intent(in)     :: sdim, nrens, nlag
   real(dp), intent(inout) :: Amat(sdim,nrens)

   real(dp) :: t_start, t_end
   integer :: i, count
   real(dp), allocatable :: Xacc(:,:), Xtmp(:,:)
   real(dp), allocatable :: Xmat_tmp(:,:)

   t_start = atime
   if (nlag == -1) then
      t_end = huge(1.0_dp)
   else
      t_end = atime + nlag * dt_x5
   end if

   allocate(Xacc(nrens,nrens), Xtmp(nrens,nrens))
   allocate(Xmat_tmp(sdim,nrens))

   Xacc = 0.0_dp
   do i=1,nrens
      Xacc(i,i) = 1.0_dp
   end do

   count = 0
   do i=1,total_x5_records
      if (all_x5(i)%tt >= t_start .and. all_x5(i)%tt <= t_end) then
         call dgemm("N", "N", nrens, nrens, nrens, 1.0_dp, Xacc, nrens, &
                    all_x5(i)%mat, nrens, 0.0_dp, Xtmp, nrens)
         Xacc = Xtmp
         count = count + 1
      end if
   end do

   if (count>0) then
      call dgemm("N", "N", sdim, nrens, nrens, 1.0_dp, Amat, sdim, &
                 Xacc, nrens, 0.0_dp, Xmat_tmp, sdim)
      Amat = Xmat_tmp
   end if

   deallocate(Xacc, Xtmp, Xmat_tmp)
end subroutine make_analysis

! ======================================================================
! Stats
subroutine make_mn_std(n, nens, Amat, mean_vec, std_vec)
    use iso_fortran_env, only: dp => real64
    implicit none

    ! Input/Output
    integer, intent(in) :: n, nens
    real(dp), intent(in) :: Amat(n, nens)
    real(dp), intent(out) :: mean_vec(n)
    real(dp), intent(out) :: std_vec(n)

    integer :: i

    do i = 1, n
        mean_vec(i) = sum(Amat(i, :)) / real(nens, dp)
    end do

    do i = 1, n
        std_vec(i) = sqrt(sum((Amat(i, :) - mean_vec(i))**2) / real(nens, dp))
    end do
end subroutine make_mn_std

! ======================================================================
! PROGRAM enKF2enKS
!
! Compute the ensemble mean and std of the KS from the KF and the X5 matrices.
! All the members are necessary.
! It is easy to implement the output for each member.
! ======================================================================
program enKF2enKS

    use iso_fortran_env, only: dp=>real64
    use mod_enks_data
    implicit none
 
    !---------------------------- CLI arguments ----------------------------
    character(len=80) :: basinf           ! basin filename (with extension)
    character(len=3)  :: lnnlv            ! number of vertical levels (as text)
    character(len=6)  :: lnrens           ! number of ensemble members (as text)
    character(len=6)  :: lnlag            ! smoother lag steps, -1 => all

    integer :: nnlv                       ! parsed nlv
    integer :: nrens                      ! parsed ensemble size
    integer :: nlag                       ! smoother lag steps, -1 => all

    !---------------------------- program variables ----------------------------
    integer           :: rrec, nre
    character(len=5)  :: nrel
    real(dp)          :: atime
    integer           :: sdim
    real(dp), allocatable :: Astate(:,:), AmeanKS(:), AstdKS(:)

! ======================================================================
! Input data
    call get_command_argument(1,basinf)
    call get_command_argument(2,lnnlv)
    call get_command_argument(3,lnrens)
    call get_command_argument(4,lnlag)

    if (lnlag == '') then
       write(*,*) ''
       write(*,*) 'Usage: enKF2enKS [basinf] [nlv] [nrens] [nlag]'
       write(*,*) ''
       write(*,*) '[basinf] bas file with the extension.'
       write(*,*) '[nlv]     number of vertical levels.'
       write(*,*) '[nrens]   number of ensemble members, control included.'
       write(*,*) '[nlag]    number of forward analysis steps to consider (-1 = all)'
       write(*,*) ''
       stop
    end if

    read(lnnlv, *) nnlv
    if (nnlv < 1) error stop "Invalid nlv"

    read(lnrens, *) nrens
    if (nrens < 2) error stop "nrens must be >= 2"

    read(lnlag, *) nlag
    if ((nlag < 2) .and. (nlag /= -1)) error stop "Invalid smoother lag"

    write(*,*) "EnKS with lag =", nlag

! ======================================================================
! WARNING
  write(*,*) 'WARNING!!! SOMETHING WRONG, RST FILES ARE TOO SMALL AND RESULTS BAD'

! ======================================================================
! Load, sort and store X5 (nrens x nrens) matricies
    call load_all_x5(nrens)

! ======================================================================
! Initialise SHYFEM data
    call init_shyfem(basinf,nnlv)

! ======================================================================
! Opening outputs
   open(18,file="analKS_mean.rst",status="unknown",form="unformatted")
   open(19,file="analKS_std.rst" ,status="unknown",form="unformatted")

! ======================================================================
! Time loop
   do rrec = 1, total_x5_records

      atime = all_x5(rrec)%tt

      ! -----------------------------------
      ! Read all members
      ! -----------------------------------
      do nre = 1, 1!nrens

         call num2str(nre-1,nrel)
         call rst_read("analKF_en"//nrel//".rst", atime)

         ! allocate the matrices
         if (.not. allocated(Astate)) call allocate_states(Astate, AmeanKS, AstdKS, sdim, nrens)

         ! From restart variables to matrix
         call push_matrix(sdim,nrens,nre,Astate)
      end do

      ! -----------------------------------
      ! Smoother
      ! -----------------------------------
      call make_analysis(atime, sdim, nrens, Astate, nlag)

      ! -----------------------------------
      ! Stats
      ! -----------------------------------
      call make_mn_std(sdim, nrens, Astate, AmeanKS, AstdKS)

      ! From matrix to restart variables
      call pull_matrix(sdim, nrens, 1, AmeanKS)
      call rst_write_record(atime, 18)

      call pull_matrix(sdim, nrens, 1, AstdKS)
      call rst_write_record(atime, 19)
  
      write(*,*) "Record, time:", rrec, atime

   end do
  
! ======================================================================
! Closing output
   close(18)
   close(19)

end program enKF2enKS
