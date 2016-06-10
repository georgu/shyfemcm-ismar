












module m_newton2D
contains
subroutine newton2D(r1,r2,n1,n2,dx,dy,rx,ry,lconv,verbose)
   use m_newtonfunc2D
   implicit none
   integer                    n1,n2
   logical, intent(in)     :: verbose
   real,    intent(inout)  :: r1,r2
   real,    intent(in)     :: dx,dy
   real,    intent(in)     :: rx,ry

   logical lconv
   real, parameter :: eps=1.0E-05
   real inc1,inc2,err1,err2
   integer i,j
   real f,g,f1,g1,f2,g2

   real r1ini,r2ini,gamma


   r1ini=r1
   r2ini=r2
   gamma=1.25

   lconv=.false.

   open(10,file='newton.cnv')
   do j=1,10
      r1=real(j)*r1ini
      r2=real(j)*r2ini
      gamma=gamma-0.25/real(j)

      do i=1,100

         call newtonfunc2D(f,g,f1,g1,f2,g2,r1,r2,n1,n2,dx,dy,rx,ry)

         inc1 =  (f*g2-f2*g)/(f1*g2-f2*g1)
         inc2 =  (f1*g-f*g1)/(f1*g2-f2*g1)

         r1 = r1 - gamma*inc1
         r2 = r2 - gamma*inc2

         r1=max(r1,0.000001)
         r2=max(r2,0.000001)

         r1=min(r1,1.0)
         r2=min(r2,1.0)

         err1 = inc1/(abs(r1)+eps)
         err2 = inc2/(abs(r2)+eps)

         write(10,'(i5,10g13.5)')i,r1,r2,f,g,f1,g1,f2,g2,err1,err2

         if (abs(err1)+abs(err2) < eps) then
            close(10)
            if (verbose) print *,' Newton converged in iteration ',j,i
            lconv=.true.
            exit
         endif

      enddo
      if (lconv) exit

   enddo
   close(10)

   if (.not. lconv) then
      print *,'Newton did not converge.'
      print *,'Probably an error in input parameters.'
      print *,'Is the dx resolving the rx scale?'
   endif

end subroutine newton2D
end module m_newton2D




