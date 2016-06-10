module m_sub1
contains
subroutine sub1()
   use m_sub2
   implicit none
#ifdef DEBUG
   print *,'Sub1'
#endif
   call sub2
end subroutine sub1
end module m_sub1
