module mod_pestimation

  use iso_fortran_env, only: dp => real64
  use mod_init_enkf,   only: nanal, nrens ! Contains 'nanal' and 'nrens'

  implicit none
  private
  
  ! Public variables and routines available to other modules/main
  public :: pe_mat, pe_meta, read_pe_files, write_pe_files, check_pe_bounds

  ! Structure to hold metadata for each parameter (from member 1)
  type :: pe_metadata
     character(len=32) :: label
     real(dp)          :: vstd
     real(dp)          :: vmin
     real(dp)          :: vmax
  end type pe_metadata

  ! Module Arrays
  real(dp), allocatable         :: pe_mat(:,:)   ! Dimensions: (ndim_pe, nrens)
  type(pe_metadata), allocatable :: pe_meta(:)    ! Dimensions: (ndim_pe)

contains

 !=======================================================================
 ! READ PE FILES (*b.dat) AND POPULATE MODULE DATA
 !=======================================================================
 subroutine read_pe_files(mode_an, ndim_pe)
    integer, intent(inout) :: mode_an
    integer, intent(out)   :: ndim_pe

    character(len=33) :: filename
    character(len=5)  :: nrel, nal
    integer           :: ne, k, ios, iu
    logical           :: bexist

    ! Temporary reading variables
    character(len=32) :: tmp_label
    real(dp)          :: tmp_val, tmp_vstd, tmp_vmin, tmp_vmax

    ndim_pe = 0
    k = 0

    call num2str(nanal, nal)

    ! --- PHASE 1: Verify files and determine ndim_pe (using member 1) ---
    do ne = 1, nrens
      call num2str(ne-1, nrel)
      filename = 'pe_parameters_an'//trim(nal)//'_en'//trim(nrel)//'.dat'
      !write(*,*) 'Reading file: ', filename

      inquire(file=filename, exist=bexist)
      if (.not. bexist) then
        write(6,*) 'PE files missing. Disabling PE.'
        mode_an = 0
        ndim_pe = 0
        if (allocated(pe_mat))  deallocate(pe_mat)
        if (allocated(pe_meta)) deallocate(pe_meta)
        return
      end if

      if (ne == 1) then
        open(newunit=iu, file=filename, status='old', form='formatted')
        do
          read(iu, *, iostat=ios) tmp_label, tmp_val, tmp_vstd, tmp_vmin, tmp_vmax
          if (ios /= 0) exit
          k = k + 1
        end do
        ndim_pe = k
        close(iu)
      end if
    end do

    ! --- PHASE 2: Allocate module arrays ---
    if (allocated(pe_mat))  deallocate(pe_mat)
    if (allocated(pe_meta)) deallocate(pe_meta)
    
    allocate(pe_mat(ndim_pe, nrens))
    allocate(pe_meta(ndim_pe))

    ! --- PHASE 3: Read data and store values & metadata ---
    do ne = 1, nrens
      call num2str(ne-1, nrel)
      filename = 'pe_parameters_an'//trim(nal)//'_en'//trim(nrel)//'.dat'

      open(newunit=iu, file=filename, status='old', form='formatted')
      do k = 1, ndim_pe
        read(iu, *) tmp_label, pe_mat(k, ne), tmp_vstd, tmp_vmin, tmp_vmax
        
        ! Store metadata only once from the first ensemble member
        if (ne == 1) then
           pe_meta(k)%label = tmp_label
           pe_meta(k)%vstd  = tmp_vstd
           pe_meta(k)%vmin  = tmp_vmin
           pe_meta(k)%vmax  = tmp_vmax
        end if
      end do
      close(iu)
    end do

    mode_an = 2

    write(6,*) 'PE matrix and metadata allocated. Dimensions: ', ndim_pe, 'x', nrens
 end subroutine read_pe_files


 !=======================================================================
 ! WRITE PE ANALYSIS FILES (*a.dat) USING STORED METADATA
 !=======================================================================
 subroutine write_pe_files(ndim_pe)
    integer, intent(in) :: ndim_pe

    character(len=33) :: filename
    character(len=5)  :: nrel, nal
    integer           :: nanal_new
    integer           :: ne, k, iu

    ! Safety check
    if (ndim_pe <= 0 .or. .not. allocated(pe_mat) .or. .not. allocated(pe_meta)) then
       write(6,*) 'PE not active or data unallocated. Skipping write.'
       return
    end if

    nanal_new = nanal + 1
    call num2str(nanal_new, nal)

    ! --- Write *a.dat files for all ensemble members ---
    do ne = 1, nrens
      call num2str(ne-1, nrel)
      filename = 'pe_parameters_an'//trim(nal)//'_en'//trim(nrel)//'.dat'
      !write(*,*) 'Writing file: ', filename

      open(newunit=iu, file=filename, status='replace', form='formatted')
      do k = 1, ndim_pe
         ! Write back using metadata from pe_meta and updated value from pe_mat
         write(iu, '(A10, 4(1X, F14.6))') &
               trim(pe_meta(k)%label), &
               real(pe_mat(k, ne)), &
               real(pe_meta(k)%vstd), &
               real(pe_meta(k)%vmin), &
               real(pe_meta(k)%vmax)

      end do
      close(iu)
    end do

    write(6,*) 'PE analysis files successfully written.' 
 end subroutine write_pe_files

 !=======================================================================
 ! BOUNDS CHECK: CLIP VALUES TO STORED MIN/MAX LIMITS
 !=======================================================================
 subroutine check_pe_bounds(ndim_pe)
    integer, intent(in) :: ndim_pe

    integer :: ne, k
    integer :: count_min, count_max

    ! Safety check to ensure arrays are allocated and active
    if (ndim_pe <= 0 .or. .not. allocated(pe_mat) .or. .not. allocated(pe_meta)) then
       write(6,*) 'PE bounds check skipped: matrix or metadata not allocated.'
       return
    end if

    count_min = 0
    count_max = 0

    ! Loop through all parameters and all ensemble members
    do ne = 1, nrens
       do k = 1, ndim_pe
          
          ! Check lower bound
          if (pe_mat(k, ne) < pe_meta(k)%vmin) then
             pe_mat(k, ne) = pe_meta(k)%vmin
             count_min = count_min + 1
          end if

          ! Check upper bound
          if (pe_mat(k, ne) > pe_meta(k)%vmax) then
             pe_mat(k, ne) = pe_meta(k)%vmax
             count_max = count_max + 1
          end if

       end do
    end do

    ! Print summary info to standard output
    if (count_min > 0 .or. count_max > 0) then
       write(6,'(A,I6,A,I6,A)') 'PE Bounds Correction: clipped ', &
            count_min, ' values to MIN and ', count_max, ' values to MAX.'
    else
       write(6,*) 'PE Bounds Correction: all values are within limits.'
    end if

 end subroutine check_pe_bounds

 !=======================================================================
 ! AUXILIARY ROUTINE: INTEGER TO STRING
 !=======================================================================
 subroutine num2str(num, str)
   integer, intent(in)           :: num
   character(len=5), intent(out) :: str

   if (num < 0 .or. num > 99999) then
      error stop "num2str: argument out of range"
   end if
   write(str, '(I5.5)') num
 end subroutine num2str

end module mod_pestimation
