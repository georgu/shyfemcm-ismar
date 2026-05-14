!======================================================================
! File: subenkf.F90
! Purpose:
!   Helper subroutines for EnKF workflows:
!   - Restart read/write wrappers
!   - Grid search utilities (nodes/elements/weights)
!   - Random vector generation and spatial perturbations (0D/2D)
!   - Super-observation construction
!   - Archival of X5 matrices for EnKS
!
! Copyright:
!   (C) 2017, Marco Bajo, CNR-ISMAR Venice. All rights reserved.
!   Updated comments and corrections (2026-02-13).
!======================================================================

!----------------------------------------------------------------------
! Read the SHYFEM restart at the target absolute time and, on first call,
! propagate restart flags/levels into global parameters.
! NOTE:
!  - Keep restart flags/fields in single precision (file format).
!  - Internal time variables are in double precision via dp kind.
!----------------------------------------------------------------------
subroutine rst_read(rstname, atimea)
  use iso_fortran_env, only : dp => real64
  use mod_restart
  use levels, only : nlvdi, nlv, hlv, ilhv, ilhkv
  use shympi
  implicit none
  character(len=*), intent(in) :: rstname
  real(dp),        intent(in) :: atimea
  integer :: ierr, iflag, ios
  real(dp) :: atimef
  integer, save :: icall = 0
  ! Single-precision parameters expected by addpar/daddpar (restart flags)
  real*4 :: ibarcl4, iconz4, imerc4, iturb4, iwvert_rst4, ieco_rst4, zero4

  zero4 = 0.0

  open(24, file=trim(rstname), status='old', form='unformatted', action='read', iostat=ios)
  if (ios /= 0) error stop 'rst_read: error opening restart file'

  do
     call rst_read_record(24, atimef, iflag, ierr)
     if (ierr /= 0) then
        close(24)
        write(*,*) 'Error in the restart file. Is the analysis time present among restart records?'
        error stop
     end if
     if (abs(atimef - atimea) < epsilon(atimea)) exit   ! exact match
  end do
  close(24)

  ! On first call, publish restart meta/flags
  if (icall == 0) then
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

     call addpar('ibarcl', ibarcl4)
     call addpar('iconz' , iconz4 )
     call addpar('iwvert', iwvert_rst4)
     call addpar('ieco'  , ieco_rst4)
     call addpar('ibio'  , zero4)
     call addpar('ibfm'  , zero4)
     call addpar('imerc' , imerc4)
     call addpar('iturb' , iturb4)

     call daddpar('date', 0.0_dp)
     call daddpar('time', 0.0_dp)

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

!----------------------------------------------------------------------
! Write a restart record at time atimea.
!----------------------------------------------------------------------
subroutine rst_write(rstname, atimea)
  use iso_fortran_env, only : dp => real64
  use basin, only : nel, nkn
  use levels, only : nlv
  use mod_hydro
  use mod_hydro_vel
  use mod_layer_thickness
  implicit none
  character(len=*), intent(in) :: rstname
  real(dp),        intent(in) :: atimea
  integer :: ios

  call mod_layer_thickness_init(nkn,nel,nlv)

  call setzev_enkf

  zov = znv
  utlov = utlnv
  vtlov = vtlnv
  zeov = zenv

  open(34, file=trim(rstname), form='unformatted', status='replace', action='write', iostat=ios)
  if (ios /= 0) error stop 'rst_write: error opening restart for write'
  call rst_write_record(atimea, 34)
  close(34)

end subroutine rst_write

!----------------------------------------------------------------------
! Duplicate of setzev with some modifications
!----------------------------------------------------------------------
subroutine setzev_enkf

  use mod_geom_dynamic
  use mod_hydro
  use basin

  implicit none

! local
  integer ie,ii

  iwetv = 0

  do ie=1,nel
      do ii=1,3
        zenv(ii,ie)=znv(nen3v(ii,ie))
        if ((zenv(ii,ie) + hm3v(ii,ie)) < 0.05) then
	  iwegv(ie) = iwegv(ie) + 1
	  zenv(ii,ie) = 0.02 - hm3v(ii,ie)
          znv(nen3v(ii,ie)) = zenv(ii,ie)
	end if
      end do
  end do

end subroutine setzev_enkf

!----------------------------------------------------------------------
! Find the nearest node (ik) of triangle element ie to point (x,y).
! Single-precision geometry (mesh).
!----------------------------------------------------------------------
subroutine find_node(x, y, ie, ik)
  use shyfile
  use basin
  implicit none
  real*4,    intent(in)  :: x, y
  integer, intent(in)  :: ie
  integer, intent(out) :: ik
  integer :: i, ikmin
  real*4 :: d, dmin

  ikmin = -999
  dmin  = huge(dmin)

  do i = 1, 3
     ik = nen3v(i, ie)
     d  = (x - xgv(ik))**2 + (y - ygv(ik))**2
     if (d <= dmin) then
        dmin  = d
        ikmin = ik
     end if
  end do
  ik = ikmin
end subroutine find_node

!----------------------------------------------------------------------
! Convert INTEGER to zero-padded 5-char string: "00000" .. "99999".
!----------------------------------------------------------------------
subroutine num2str(num, str)
  implicit none
  integer,        intent(in)  :: num
  character(len=5), intent(out) :: str

  if (num < 0 .or. num > 99999) error stop 'num2str: num out of range'
  write(str, '(I5.5)') num
end subroutine num2str

!----------------------------------------------------------------------
! Generate a random vector with zero mean. v(1)=0; v(2:) ~ centered RNG.
! Outliers |a|>=3 are compressed to keep finite variance.
!----------------------------------------------------------------------
subroutine random_vec(v, vdim)
  use iso_fortran_env, only : dp => real64
  use m_random
  implicit none
  integer, intent(in)  :: vdim
  real(dp),    intent(out) :: v(vdim)
  real(dp) :: vaux(max(1, vdim-1))
  real(dp) :: aaux, ave
  integer :: n

  if (vdim < 1) then
     return
  else if (vdim == 1) then
     v(1) = 0.0_dp
     return
  end if

  call random(vaux, vdim-1)

  ! Compress outliers to avoid extremely large magnitudes
  do n = 1, vdim-1
     aaux = vaux(n)
     if (abs(aaux) >= 3.0_dp) then
        ! Use 1.0_dp to match the kind of aaux
        aaux = sign(1.0_dp, aaux) * (abs(aaux) - floor(abs(aaux)) + 1.0_dp)
     end if
     vaux(n) = aaux
  end do

  ! Zero mean and insert v(1)=0
  ave = sum(vaux) / real(vdim-1, dp)
  vaux = vaux - ave
  v(1) = 0.0_dp
  v(2:vdim) = vaux
end subroutine random_vec

!----------------------------------------------------------------------
! Temporal red-noise perturbation for 0D observations (z/t/s).
! If tau <= 0, returns white noise. Otherwise blends with previous state.
! vec <- alpha * vec_old + sqrt(1 - alpha^2) * vec, alpha = 1 - dt/tau
!----------------------------------------------------------------------
subroutine make_0Dpert(vflag, n, na, id, vec, t, tau)
  use iso_fortran_env, only : dp => real64
  implicit none
  character(len=5), intent(in) :: vflag
  integer,         intent(in)  :: n, na, id
  real(dp),            intent(out) :: vec(n)
  real(dp),        intent(in)  :: t, tau

  character(len=5)  :: nal, idl
  character(len=21) :: pfile
  logical :: bfile
  integer :: nf, ios
  real(dp), allocatable :: vec_old(:)
  real(dp) :: t_old, alpha, dt, s1
  character(len=1) :: vfl

  ! New random perturbation
  call random_vec(vec, n)

  ! White noise if no temporal correlation requested
  if (tau <= 0.0_dp) return

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
     read(22, iostat=ios) nf;    if (ios /= 0) error stop 'make_0Dpert: read nf failed'
     if (nf /= n) error stop 'make_0Dpert: dimension mismatch'
     read(22, iostat=ios) t_old; if (ios /= 0) error stop 'make_0Dpert: read time failed'
     allocate(vec_old(nf))
     read(22, iostat=ios) vec_old; if (ios /= 0) error stop 'make_0Dpert: read vector failed'
     close(22)

     dt    = t - t_old
     alpha = 1.0_dp - (dt / tau)
     if (alpha < 0.0_dp) alpha = 0.0_dp
     if (alpha > 1.0_dp) alpha = 1.0_dp
     s1 = sqrt( max(0.0_dp, 1.0_dp - alpha*alpha) )

     vec = real(alpha, dp) * vec_old + real(s1, dp) * vec
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
  use iso_fortran_env, only : dp => real64
  use basin
  use m_sample2D
  use mod_para
  implicit none

  ! Arguments
  integer, intent(in)      :: n, nens
  real(dp), intent(out)    :: vec(n, nens) ! Final ensemble perturbations

  ! Local variables
  integer                  :: nx, ny, ne
  real(dp)                 :: x1, x2, y1, y2, xlength, ylength
  real(dp)                 :: dx, dy, rx, ry
  real(dp)                 :: dx_m, dy_m, rx_m, ry_m
  real*4                   :: sdx, sdy, sx0, sy0
  real*4, parameter        :: sflag = -999.0
  real(dp), allocatable    :: mat(:,:,:)
  real*4,   allocatable    :: mat4(:,:), vec4fem(:)

  ! 1. Determine domain extents from FEM grid (xgv, ygv)
  x1 = floor(minval(xgv)); x2 = ceiling(maxval(xgv))
  y1 = floor(minval(ygv)); y2 = ceiling(maxval(ygv))
  xlength = x2 - x1; ylength = y2 - y1

  ! Safety check for geographical coordinates
  if ((xlength > 180.0) .or. (ylength > 90.0)) then
     error stop 'make_2Dpert: coordinates do not appear to be geographical (Lat/Lon)'
  end if

  ! 2. Set grid resolution and decorrelation scales based on domain size
  if (xlength < 4.0) then
     dx = 0.05; dy = 0.05
     rx = 2.0;  ry = 2.0
  else
     dx = 0.10; dy = 0.10
     rx = 4.0;  ry = 4.0
  end if

  ! Use nint to avoid truncation errors in grid size calculation
  nx = nint(xlength/dx) + 1
  ny = nint(ylength/dy) + 1

  if (verbose) then
     write(*,'(a20,2f8.4,i5,f8.4)') 'x1,xlength,nx,dx: ', x1, xlength, nx, dx
     write(*,'(a20,2f8.4,i5,f8.4)') 'y1,ylength,ny,dy: ', y1, ylength, ny, dy
     write(*,'(a14,2f10.4,1x,i3)') 'rx,ry,fmult: ', rx, ry, fmult_init
  end if

  ! 3. Generate random samples on a regular grid (requires dp)
  allocate(mat(nx, ny, nens))

  ! Convert degrees to meters at the reference latitude (y1)
  call deg2meters(x1, y1, dx, .true.,  dx_m)
  call deg2meters(x1, y1, dy, .false., dy_m)
  call deg2meters(x1, y1, rx, .true., rx_m)
  call deg2meters(x1, y1, ry, .false., ry_m)

  ! Call the external library (seed is pre-set externally)
  call sample2D(mat, nx, ny, nens, fmult_init, dx_m, dy_m, rx_m, ry_m, &
                theta_init, sample_fix_init, verbose)

  ! 4. Interpolate regular grid to SHYFEM (unstructured) grid
  write(*,*) 'Interpolating 2D field over the FEM grid...'

  ! Casting to real*4 as required by SHYFEM interpolation routines
  sdx = real(dx, 4); sdy = real(dy, 4)
  sx0 = real(x1, 4); sy0 = real(y1, 4)

  ! Initialize the interpolation engine
  call setgeo(sx0, sy0, sdx, sdy, sflag)

  allocate(mat4(nx, ny), vec4fem(n))

  do ne = 1, nens
     ! Copy dp slice to real*4 temporary array
     mat4 = real(mat(:,:,ne), 4)

     ! am2av: Interpolation from regular grid (mat4) to FEM nodes (vec4fem)
     call am2av(mat4, vec4fem, nx, ny)

     ! Store result back into the output array (auto-cast to dp)
     vec(:, ne) = vec4fem
  end do

  ! 5. Clean up memory
  deallocate(mat, mat4, vec4fem)

end subroutine make_2Dpert


!-----------------------------------------------------------------------
! Convert a distance expressed in degrees (lon/lat direction) into meters
! at a given geographic position.
!-----------------------------------------------------------------------
subroutine deg2meters(lon, lat, dist_deg, is_lon, dist_m)
  use iso_fortran_env, only : dp => real64
  implicit none
  real(dp), intent(in)  :: lon, lat, dist_deg
  logical, intent(in)   :: is_lon
  real(dp), intent(out) :: dist_m
  real(dp), parameter :: deg2rad = 3.141592653589793_dp / 180.0_dp
  real(dp), parameter :: R_earth = 6371000.0_dp  ! mean Earth radius [m]

  if (is_lon) then
     ! distance along longitude depends on latitude
     dist_m = dist_deg * deg2rad * R_earth * cos(lat * deg2rad)
  else
     ! distance along latitude is constant
     dist_m = dist_deg * deg2rad * R_earth
  end if
end subroutine deg2meters

!----------------------------------------------------------------------
! Find the element iie and its nearest node ik to point (x,y).
! Works with single-precision geometry (FEM utilities).
!----------------------------------------------------------------------
subroutine find_el_node(x, y, iie, ik)
  use iso_fortran_env, only : dp => real64
  use basin
  implicit none
  real(dp),    intent(in)  :: x, y
  integer, intent(out) :: iie, ik
  real*4 :: x4, y4
  real(dp) :: dst, dstmax
  integer :: iik, ii, k

  x4 = x; y4 = y
  call find_element(x4, y4, iie)

  dstmax = huge(0.0_dp)
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
! Gaspari–Cohn taper weight w(rho, dst) with compact support [0, 2*rho].
! s = dst / rho.
!----------------------------------------------------------------------
subroutine find_weight_GC(rho, dst, w)
  use iso_fortran_env, only : dp => real64
  implicit none
  real(dp), intent(in)  :: rho, dst
  real(dp), intent(out) :: w
  real(dp) :: s

  if (rho <= 0.0) then
     w = 0.0
     return
  end if

  s = dst / rho
  if ((s >= 0.0) .and. (s < 1.0)) then
     w = 1.0 - (5.0/3.0)*s**2 + (5.0/8.0)*s**3 + 0.5*s**4 - 0.25*s**5
  else if ((s >= 1.0) .and. (s < 2.0)) then
     w = - (2.0/3.0)*s*(-1) + 4.0 - 5.0*s + (5.0/3.0)*s**2 + (5.0/8.0)*s**3 &
         - 0.5*s**4 + (1.0/12.0)*s**5
  else
     w = 0.0
  end if
end subroutine find_weight_GC

!----------------------------------------------------------------------
! Append the current (global/local) analysis transform matrices to file
! 'X5_tot.uf' for EnKS use. Supports 'X3' (+S) or 'X5' encodings.
!----------------------------------------------------------------------
subroutine save_X5(tt)
  use iso_fortran_env, only : dp => real64
  implicit none

  real(dp),        intent(in) :: tt       ! Current time

  integer :: nrens, nrobs, ios, u_in, u_out
  real(dp), allocatable :: X5(:,:), X3(:,:), S(:,:)
  character(len=2) :: tag2

  ! 1. Open the temporary file generated by dumpX3/dumpX5
  open(newunit=u_in, file='X5.uf', form='unformatted', status='old', action='read', iostat=ios)
  if (ios /= 0) return ! Exit quietly if file doesn't exist (e.g., local analysis skipped)

  ! Identify the format (X3 or X5)
  read(u_in, iostat=ios) tag2
  rewind(u_in)

  if (tag2 == 'X3') then
     ! CRITICAL FIX: X3 dimensions must be (nrobs, nrens) as defined in analysis/dumpX3
     read(u_in) tag2, nrens, nrobs
     allocate(X3(nrobs, nrens), S(nrobs, nrens))
     rewind(u_in)
     read(u_in) tag2, nrens, nrobs, X3, S
  else
     read(u_in) tag2, nrens
     allocate(X5(nrens, nrens))
     rewind(u_in)
     read(u_in) tag2, nrens
     read(u_in) X5
  end if
  close(u_in)


  ! 2. Append to the cumulative Smoother file (X5_tot.uf)
  ! Using 'access=stream' or standard sequential append
  write(*,*) 'save_X5: Saving weights for EnKS at time ', tt

  open(newunit=u_out, file='X5_tot.uf', form='unformatted', status='unknown', &
       position='append', action='write', iostat=ios)
  if (ios /= 0) error stop 'save_X5: cannot open X5_tot.uf'

  ! Header for this specific time step
  write(u_out) tt, tag2

  if (tag2 == 'X3') then
     write(u_out) nrens, nrobs
     write(u_out) X3
     write(u_out) S
  else
     write(u_out) nrens
     write(u_out) X5
  end if
  close(u_out)

  ! 3. Cleanup
  if (allocated(X3)) deallocate(X3)
  if (allocated(S )) deallocate(S)
  if (allocated(X5)) deallocate(X5)
end subroutine save_X5
