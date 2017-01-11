!********************************************************

  subroutine random2(work1,n)
! Returns a vector of random values N(variance=1,mean=0)
! From Evensen's code 
   implicit none
   integer, intent(in) :: n
   real,   intent(out) :: work1(n)
   real,   allocatable :: work2(:)
   real, parameter :: pi=3.14159253589

   allocate (work2(n))

   call random_number(work1)
   call random_number(work2)
   work1= sqrt(-2.0*log(work1))*cos(2.0*pi*work2)

   deallocate(work2)

  end subroutine random2

!********************************************************

  subroutine check_innov_val(inn,typeo)
  implicit none
  double precision, intent(inout) :: inn
  character(len=*), intent(in) :: typeo

  if( trim(typeo).eq.'level' ) then
    if ( abs(inn) .gt. 1. ) then
       write(*,*) 'Warning: level innovation too large: ',inn
       write(*,*) '         Setting it to zero...'
       inn = 0.
    end if
  else
    write(*,*) 'check_innov_val: warning observation not implemented'
  end if

  end subroutine check_innov_val

!********************************************************

  subroutine num2str(num,str)
  implicit none
  integer, intent(in) :: num
  character(len=3), intent(out) :: str

  if( (num.ge.0).and.(num.lt.10) ) then
    write(str,'(a2,i1)') '00',num
  elseif( (num.ge.10).and.(num.lt.100) ) then
    write(str,'(a1,i2)') '0',num
  elseif( (num.ge.100).and.(num.lt.1000) ) then
    write(str,'(i3)') num
  else
    stop 'num2str: num out of range'
  end if

  end subroutine num2str
