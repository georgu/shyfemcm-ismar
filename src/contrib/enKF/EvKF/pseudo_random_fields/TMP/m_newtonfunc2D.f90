












module m_newtonfunc2D
contains
subroutine newtonfunc2D(f,g,f1,g1,f2,g2,r1,r2,n1,n2,dx,dy,rx,ry)
   implicit none
   integer                 n1,n2
   real,    intent(out) :: f,g,f1,g1,f2,g2
   real,    intent(in)  :: r1,r2
   real,    intent(in)  :: dx,dy
   real,    intent(in)  :: rx,ry
   real, parameter   :: pi=3.141592653589

   integer l,p
   real pi2,e,kappa,kappa2,lambda,lambda2


   pi2=2.0*pi

   kappa=pi2/(real(n1)*dx)
   kappa2=kappa**2

   lambda=pi2/(real(n2)*dy)
   lambda2=lambda**2

   f=0.0; g=0.0; f1=0.0; g1=0.0; f2=0.0; g2=0.0

   do p=-n2/2+1,n2/2
   do l=-n1/2+1,n1/2
      if (p == 0 .and. l == 0) cycle

      e=exp( -2.0*( kappa2*real(l*l)/r1**2 + lambda2*float(p*p)/r2**2 ) )

      f=f + e * ( cos(kappa *real(l)*rx) - exp(-1.0) )
      g=g + e * ( cos(lambda*real(p)*ry) - exp(-1.0) )

      f1=f1 + e * (4.0*kappa2 *real(l*l)/r1**3) * ( cos(kappa *float(l)*rx) - exp(-1.0) )
      g1=g1 + e * (4.0*kappa2 *real(l*l)/r1**3) * ( cos(lambda*float(p)*ry) - exp(-1.0) )
      f2=f2 + e * (4.0*lambda2*real(p*p)/r2**3) * ( cos(kappa *float(l)*rx) - exp(-1.0) )
      g2=g2 + e * (4.0*lambda2*real(p*p)/r2**3) * ( cos(lambda*float(p)*ry) - exp(-1.0) )
      
   enddo
   enddo

end subroutine newtonfunc2D
end module m_newtonfunc2D




