!======================================================================
!  File: subenkf.f90
!  Purpose:
!    Helper subroutines for EnKF workflows:
!      - Restart read/write wrappers
!      - Grid search utilities (nodes/elements/weights)
!      - Random vector generation and spatial perturbations (0D/2D)
!      - Super-observation construction
!      - Archival of X5 matrices for EnKS
!
!  Copyright:
!    (C) 2017, Marco Bajo, CNR-ISMAR Venice. All rights reserved.
!    Updated comments and corrections (2026-02-13).
!
!======================================================================

!----------------------------------------------------------------------
! Read the SHYFEM restart at the target absolute time and, on first call,
! propagate restart flags/levels into global parameters.
!----------------------------------------------------------------------
subroutine rst_read(rstname, atimea)
  use mod_restart
  use levels, only : nlvdi, nlv, hlv, ilhv, ilhkv
  use shympi
  implicit none
  character(len=*), intent(in) :: rstname
  double precision, intent(in) :: atimea

  integer :: ierr, iflag, ios
  double precision :: atimef
  integer, save :: icall = 0

  ! Single-precision parameters expected by addpar/daddpar
  real*4 :: ibarcl4, iconz4, imerc4, iturb4, iwvert_rst4, ieco_rst4, zero4
  zero4 = 0.0

  open(24, file=trim(rstname), status='old', form='unformatted', action='read', iostat=ios)
  if (ios /= 0) error stop 'rst_read: error opening restart file'

  !atimef = atimea
  do
    call rst_read_record(24, atimef, iflag, ierr)
    !call rst_read_dummy(24, atimef, iflag, ierr)
    if (ierr /= 0) then
      close(24)
      write(*,*) 'Error in the restart file. Is the analysis time present among restart records?'
      error stop
    end if
    if (atimef == atimea) exit    ! exact match as in original logic
  end do

  close(24)

  ! On first call, publish restart meta/flags
  if (icall == 0) then
    hlv  = hlvrst
    hlv_global = hlvrst
    ilhv = ilhrst
    ilhkv = ilhkrst

    ibarcl4     = ibarcl_rst
    iwvert_rst4 = iwvert_rst
    ieco_rst4   = ieco_rst
    iconz4      = iconz_rst
    imerc4      = imerc_rst
    iturb4      = iturb_rst

    call addpar('ibarcl', ibarcl4)
    call addpar('iconz',  iconz4)
    call addpar('iwvert', iwvert_rst4)
    call addpar('ieco',   ieco_rst4)
    call addpar('ibio',   zero4)
    call addpar('ibfm',   zero4)
    call addpar('imerc',  imerc4)
    call addpar('iturb',  iturb4)
    call daddpar('date',  0.0d0)
    call daddpar('time',  0.0d0)

    write(*,*) 'SHYFEM flags from restart:'
    write(*,*) 'nvers	= ',nvers_rst
    write(*,*) 'nvmax	= ',nvmax
    write(*,*) 'ibarcl = ', ibarcl_rst
    write(*,*) 'iconz  = ', iconz_rst
    write(*,*) 'iwvert = ', iwvert_rst
    write(*,*) 'ieco   = ', ieco_rst
    write(*,*) 'imerc  = ', imerc_rst
    write(*,*) 'iturb  = ', iturb_rst
    write(*,*) 'nlv	= ',nlv
    !write(*,*) 'hlvrst	= ',hlvrst(1:nlv)
    write(*,*) 'hlv	= ',hlv(1:nlv)
  end if

  icall = icall + 1
end subroutine rst_read

!----------------------------------------------------------------------
! Write a minimalist restart record at time atimea.
!----------------------------------------------------------------------
subroutine rst_write(rstname, atimea)
  implicit none
  character(len=*), intent(in) :: rstname
  double precision, intent(in) :: atimea
  integer :: ios

  open(34, file=trim(rstname), form='unformatted', status='replace', action='write', iostat=ios)
  if (ios /= 0) error stop 'rst_write: error opening restart for write'

  call rst_write_record(atimea, 34)
  close(34)
end subroutine rst_write

!----------------------------------------------------------------------
! Find the nearest node (ik) of triangle element ie to point (x,y).
!----------------------------------------------------------------------
subroutine find_node(x, y, ie, ik)
  use shyfile
  use basin
  implicit none
  real,    intent(in)  :: x, y
  integer, intent(in)  :: ie
  integer, intent(out) :: ik

  integer :: i, ikmin
  real    :: d, dmin

  ikmin = -999
  dmin  = huge(dmin)

  do i = 1, 3
    ik = nen3v(i, ie)
    d  = (x - xgv(ik))**2 + (y - ygv(ik))**2
    if (d <= dmin) then
      dmin = d
      ikmin = ik
    end if
  end do

  ik = ikmin
end subroutine find_node

!----------------------------------------------------------------------
! Set layer thickness proxy zenv from znv, ensuring minimum water depth.
!----------------------------------------------------------------------
subroutine layer_thick(nelem)
  use mod_hydro
  use basin
  implicit none
  integer, intent(in) :: nelem

  integer :: ie, ii, k
  real*4  :: z, h
  real*4  :: zmin

  zmin = 0.01

  do ie = 1, nelem
    do ii = 1, 3
      k = nen3v(ii, ie)
      zenv(ii, ie) = znv(k)
      z = zenv(ii, ie)
      h = hm3v(ii, ie)
      if (z + h < zmin) then
        write(*,*) 'Warning: small layer thickness. Correcting: ', z, zmin
        z = zmin
        zenv(ii, ie) = z
      end if
    end do
  end do
end subroutine layer_thick

!----------------------------------------------------------------------
! Convert INTEGER to zero-padded 5-char string: "00000" .. "99999".
!----------------------------------------------------------------------
subroutine num2str(num, str)
  implicit none
  integer,          intent(in)  :: num
  character(len=5), intent(out) :: str
  if (num < 0 .or. num > 99999) error stop 'num2str: num out of range'
  write(str, '(I5.5)') num
end subroutine num2str

!----------------------------------------------------------------------
! Generate a random vector with zero mean. v(1)=0; v(2:) ~ centered RNG.
! Outliers |a|>=3 are compressed to keep finite variance.
!----------------------------------------------------------------------
subroutine random_vec(v, vdim)
  use m_random
  implicit none
  integer, intent(in)  :: vdim
  real,    intent(out) :: v(vdim)

  real :: vaux(max(1, vdim-1))
  real :: aaux, ave
  integer :: n

  if (vdim < 1) then
    return
  else if (vdim == 1) then
    v(1) = 0.0
    return
  end if

  call random(vaux, vdim-1)

  ! Compress outliers to avoid extremely large magnitudes
  do n = 1, vdim-1
    aaux = vaux(n)
    if (abs(aaux) >= 3.0) then
      aaux = sign(1.0, aaux) * (abs(aaux) - floor(abs(aaux)) + 1.0)
    end if
    vaux(n) = aaux
  end do

  ! Zero mean and insert v(1)=0
  ave = sum(vaux) / real(vdim-1)
  vaux = vaux - ave
  v(1) = 0.0
  v(2:vdim) = vaux
end subroutine random_vec

!----------------------------------------------------------------------
! Temporal red-noise perturbation for 0D observations (z/t/s).
! If tau <= 0, returns white noise. Otherwise blends with previous state.
!   vec <- alpha * vec_old + sqrt(1 - alpha^2) * vec,  alpha=1 - dt/tau
!----------------------------------------------------------------------
subroutine make_0Dpert(vflag, n, na, id, vec, t, tau)
  implicit none
  character(len=5), intent(in)  :: vflag
  integer,         intent(in)   :: n, na, id
  real,            intent(out)  :: vec(n)
  double precision,intent(in)   :: t, tau

  character(len=5)  :: nal, idl
  character(len=21) :: pfile
  logical :: bfile
  integer :: nf, ios
  real, allocatable :: vec_old(:)
  double precision :: t_old, alpha, dt, s1
  character(len=1)  :: vfl

  ! New random perturbation
  call random_vec(vec, n)

  ! White noise if no temporal correlation requested
  if (tau <= 0.0d0) return

  select case (vflag)
  case ('0DLEV'); vfl = 'z'
  case ('0DTEM'); vfl = 't'
  case ('0DSAL'); vfl = 's'
  case default
    error stop 'make_0Dpert: unknown vflag'
  end select

  ! Load and blend previous perturbation, if available
  call num2str(na-1, nal)
  call num2str(id,   idl)
  pfile = vfl // 'pert_' // nal // '_' // idl // '.bin'

  inquire(file=pfile, exist=bfile)
  if (bfile) then
    open(22, file=pfile, status='old', form='unformatted', action='read', iostat=ios)
    if (ios /= 0) error stop 'make_0Dpert: cannot open previous perturbation'
    read(22, iostat=ios) nf;     if (ios /= 0) error stop 'make_0Dpert: read nf failed'
    if (nf /= n) error stop 'make_0Dpert: dimension mismatch'
    read(22, iostat=ios) t_old;  if (ios /= 0) error stop 'make_0Dpert: read time failed'
    allocate(vec_old(nf))
    read(22, iostat=ios) vec_old; if (ios /= 0) error stop 'make_0Dpert: read vector failed'
    close(22)

    dt    = t - t_old
    alpha = 1.0d0 - (dt / tau)
    if (alpha < 0.0d0) alpha = 0.0d0
    if (alpha > 1.0d0) alpha = 1.0d0
    s1    = sqrt(max(0.0d0, 1.0d0 - alpha*alpha))

    vec = real(alpha) * vec_old + real(s1) * vec
    deallocate(vec_old)
  end if

  ! Save the latest perturbation
  call num2str(na, nal)
  pfile = vfl // 'pert_' // nal // '_' // idl // '.bin'
  open(32, file=pfile, form='unformatted', status='replace', action='write', iostat=ios)
  if (ios /= 0) error stop 'make_0Dpert: cannot open output perturbation'
  write(32) n
  write(32) t
  write(32) vec
  close(32)
end subroutine make_0Dpert

!----------------------------------------------------------------------
! Generate 2D spatial perturbations on a regular grid, then interpolate
! to the FEM grid. Sampling density depends on domain size.
!----------------------------------------------------------------------
subroutine make_2Dpert(vec, n, nens)
  use basin
  use m_sample2D
  use mod_para
  implicit none
  integer, intent(in)  :: n, nens
  real,    intent(out) :: vec(n, nens)

  integer :: nx, ny
  real :: x1, x2, y1, y2, xlength, ylength
  real :: dx, dy, rx, ry
  real*4 :: sdx, sdy, sx0, sy0
  real*4, parameter :: sflag = -999.0
  real, allocatable :: mat(:,:,:)
  real*4, allocatable :: mat4(:,:), vec4fem(:)
  integer :: ne

  ! Domain extents
  x1 = floor(minval(xgv)); x2 = ceiling(maxval(xgv))
  y1 = floor(minval(ygv)); y2 = ceiling(maxval(ygv))
  xlength = x2 - x1; ylength = y2 - y1

  if ((xlength > 180.0) .or. (ylength > 90.0)) error stop 'make_2Dpert: coordinates not geographical'

  if (xlength < 4.0) then
    dx = 0.05; dy = 0.05
    rx = 2.0;  ry = 2.0
  else
    dx = 0.10; dy = 0.10
    rx = 4.0;  ry = 4.0
  end if
  nx = int(xlength/dx) + 1
  ny = int(ylength/dy) + 1

  if (verbose) then
    write(*,'(a20,2f8.4,i5,f8.4)') 'x1,xlength,nx,dx: ', x1, xlength, nx, dx
    write(*,'(a20,2f8.4,i5,f8.4)') 'y1,ylength,ny,dy: ', y1, ylength, ny, dy
    write(*,'(a14,2f10.4,1x,i3)') 'rx,ry,fmult: ', rx, ry, fmult_init
  end if

  ! Create samples on regular grid
  allocate(mat(nx, ny, nens))
  call sample2D(mat, nx, ny, nens, fmult_init, dx, dy, rx, ry, theta_init, sample_fix_init, verbose)

  ! Interpolate to FEM grid
  write(*,*) 'Interpolating 2D field over the FEM grid...'
  sdx = dx; sdy = dy; sx0 = x1; sy0 = y1
  call setgeo(sx0, sy0, sdx, sdy, sflag)

  allocate(mat4(nx, ny), vec4fem(n))
  do ne = 1, nens
    mat4 = mat(:,:,ne)
    call am2av(mat4, vec4fem, nx, ny)
    vec(:, ne) = vec4fem
  end do

  deallocate(mat, mat4, vec4fem)
end subroutine make_2Dpert

!----------------------------------------------------------------------
! Find the element iie and its nearest node ik to point (x,y).
! Works with single-precision geometry (FEM utilities).
!----------------------------------------------------------------------
subroutine find_el_node(x, y, iie, ik)
  use basin
  implicit none
  real,    intent(in)  :: x, y
  integer, intent(out) :: iie, ik

  real*4 :: x4, y4
  real :: dst, dstmax
  integer :: iik, ii, k

  x4 = x; y4 = y
  call find_element(x4, y4, iie)

  dstmax = 1.0e15

  if (iie /= 0) then
    do ii = 1, 3
      iik = nen3v(ii, iie)
      dst = sqrt( (xgv(iik)-x4)**2 + (ygv(iik)-y4)**2 )
      if (dst < dstmax) then
        dstmax = dst
        ik = iik
      end if
    end do
  else
    ! Outside the grid: choose nearest global node, then element around it
    do k = 1, nkn
      dst = sqrt( (xgv(k)-x4)**2 + (ygv(k)-y4)**2 )
      if (dst < dstmax) then
        dstmax = dst
        ik = k
      end if
    end do
    call find_element(xgv(ik), ygv(ik), iie)
  end if
end subroutine find_el_node

!----------------------------------------------------------------------
! Gaspari–Cohn taper weight w(rho, dst) with compact support [0,2*rho].
! s = dst / rho.
!----------------------------------------------------------------------
subroutine find_weight_GC(rho, dst, w)
  implicit none
  real, intent(in)  :: rho, dst
  real, intent(out) :: w
  real :: s

  if (rho <= 0.0) then
    w = 0.0
    return
  end if

  s = dst / rho

  if ((s >= 0.0) .and. (s < 1.0)) then
    w = 1.0 - (5.0/3.0)*s**2 + (5.0/8.0)*s**3 + 0.5*s**4 - 0.25*s**5
  else if ((s >= 1.0) .and. (s < 2.0)) then
    w = - (2.0/3.0)*s**(-1) + 4.0 - 5.0*s + (5.0/3.0)*s**2 + (5.0/8.0)*s**3  &
        - 0.5*s**4 + (1.0/12.0)*s**5
  else
    w = 0.0
  end if
end subroutine find_weight_GC

!----------------------------------------------------------------------
! Build super-observations by averaging observations that fall into the
! same FEM element. First one becomes the super-ob, others are marked merged.
!----------------------------------------------------------------------
subroutine superobs_horiz_el(no, x, y, ostatus, val1, val2)
  use basin
  implicit none
  integer, intent(in)    :: no
  real,    intent(in)    :: x(no), y(no)
  integer, intent(inout) :: ostatus(no)
  real,    intent(inout) :: val1(no), val2(no)

  integer :: nobs(nel)
  integer, allocatable :: ieobs(:)
  integer, allocatable :: oindex(:,:)
  integer :: omax
  real*4 :: x4, y4
  integer :: n, ie, nn
  real :: avval

  allocate(ieobs(no))
  ieobs = -999
  nobs  = 0

  ! For each observation, find containing element (if any)
  do n = 1, no
    if (ostatus(n) > 0) cycle
    x4 = x(n)
    y4 = y(n)
    call find_element(x4, y4, ie)
    ieobs(n) = ie
    if (ie >= 1 .and. ie <= nel) nobs(ie) = nobs(ie) + 1
  end do

  ! Maximum number of obs per element
  omax = 0
  do ie = 1, nel
    omax = max(omax, nobs(ie))
  end do

  ! Index table: oindex(0,ie)=count; oindex(1..count,ie)=obs indices
  allocate(oindex(0:omax, nel))
  oindex = 0

  do n = 1, no
    if (ostatus(n) > 0) cycle
    ie = ieobs(n)
    if (ie < 1 .or. ie > nel) cycle
    nn = oindex(0, ie) + 1
    oindex(nn, ie) = n
    oindex(0,  ie) = nn
  end do

  ! Average values within each element; mark statuses
  do ie = 1, nel
    nn = oindex(0, ie)
    if (nn == 0) cycle
    write(*,*) 'making super-observation...'
    avval = sum(val1(oindex(1:nn, ie))) / real(nn)
    val1(oindex(1, ie)) = avval
    avval = sum(val2(oindex(1:nn, ie))) / real(nn)
    val2(oindex(1, ie)) = avval
    ostatus(oindex(1, ie)) = 1
    if (nn > 1) ostatus(oindex(2:nn, ie)) = 2
  end do

  deallocate(ieobs, oindex)
end subroutine superobs_horiz_el

!----------------------------------------------------------------------
! Append the current (global/local) analysis transform matrices to file
! 'X5_tot.uf' for EnKS use. Supports 'X3' (+S) or 'X5' encodings.
!----------------------------------------------------------------------
subroutine save_X5(alabel, tt)
  implicit none
  character(len=*), intent(in) :: alabel   ! e.g., 'GLOBAL'/'LOCAL' (<=6 chars)
  double precision, intent(in) :: tt

  integer :: nrens, nrobs, ios
  real, allocatable :: X5(:,:), X3(:,:), S(:,:)
  character(len=2) :: tag2

  ! Read X5.uf (either 'X3' or 'X5' format)
  open(unit=25, file='X5.uf', form='unformatted', status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'save_X5: cannot open X5.uf'

  read(25, iostat=ios) tag2
  if (ios /= 0) then
    close(25)
    error stop 'save_X5: failed reading tag'
  end if
  rewind(25)

  if (tag2 == 'X3') then
    read(25, iostat=ios) tag2, nrens, nrobs
    if (ios /= 0) then
      close(25); error stop 'save_X5: failed reading header (X3)'
    end if
    allocate(X3(nrens, nrens), S(nrobs, nrens))
    rewind(25)
    read(25, iostat=ios) tag2, nrens, nrobs, X3, S
    if (ios /= 0) then
      close(25); error stop 'save_X5: failed reading matrices (X3)'
    end if
  else
    read(25, iostat=ios) tag2, nrens
    if (ios /= 0) then
      close(25); error stop 'save_X5: failed reading header (X5)'
    end if
    allocate(X5(nrens, nrens))
    rewind(25)
    read(25, iostat=ios) tag2, nrens, X5
    if (ios /= 0) then
      close(25); error stop 'save_X5: failed reading matrix (X5)'
    end if
  end if
  close(25)

  ! Append to cumulative file
  write(*,*) 'saving X5 matrix for EnKS'
  open(unit=35, file='X5_tot.uf', form='unformatted', status='unknown', access='append', &
       action='write', iostat=ios)
  if (ios /= 0) error stop 'save_X5: cannot open X5_tot.uf'

  write(35) tt, alabel, tag2
  if (tag2 == 'X3') then
    write(35) nrens, nrobs
    write(35) X3, S
  else
    write(35) nrens
    write(35) X5
  end if
  close(35)

  if (allocated(X3)) deallocate(X3)
  if (allocated(S))  deallocate(S)
  if (allocated(X5)) deallocate(X5)
end subroutine save_X5
