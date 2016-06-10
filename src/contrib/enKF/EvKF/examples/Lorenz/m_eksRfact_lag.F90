module m_eksRfact_lag
contains
subroutine eksRfact_lag(S,RR,bb,nsing,obsvar)
   use mod_dimensions
   implicit none
   real, intent(in)        :: S(nrobs,nrsamp)
   real, intent(out)       :: RR(nrobs,nrobs)
   real, intent(inout)     :: bb(nrobs,0:nrsamp)
   integer, intent(out)    :: nsing
   real, intent(inout)     :: obsvar

   integer i,j,m,mm,info,i1,i2
   real rcond,summ,summ2

   integer, allocatable :: ipvt(:)
   real, allocatable :: ss(:)
   real, allocatable :: H(:,:)
   real, allocatable :: RT(:)
   real, allocatable :: ZZ(:,:)
   real, allocatable :: w(:)

   real, allocatable :: TTT(:,:)

   allocate (ipvt(nrobs))
   allocate (ss(nrobs))
   allocate (H(nrobs,3))
   allocate (RT(nrobs*nrobs))
   allocate (ZZ(nrobs,nrobs))
   allocate (w(8*nrobs))

   print *,'start eksRfact'

!  Representer matrix */
   
   allocate (TTT(nrsamp,nrobs))
   TTT=transpose(S)
   RR=matmul(S,TTT)
   deallocate(TTT)

   print *,'RR'

   RR=RR/float(nrsamp)
   print *,'RR2'
!     print '(8f10.4)',RR

!  Measurement error covariance matrix RR */
   H=0.0
   do m=1,3
      H(m,m)=1.0
   enddo

   do i=1,nrobs
      RR(i,i)=RR(i,i)+obsvar
   enddo

   print '(3f10.3)',RR
   print *,'calling rsm'
!  Calculating the Eigendecomposition of RR using eispack */
   call rsm(nrobs,nrobs,RR,ss,nrobs,ZZ,w,ipvt,info)
   print *,'rsm ok'
   RR=ZZ

   print *,'ss'
   print '(10g11.3)',ss

   summ=sum(ss)
   summ2=0.0
   do i=nrobs,1,-1
      summ2=summ2+ss(i)
      if (summ2/summ > 0.99999)  then
         m=nrobs-i+1
         print *,'Number of singular values= ',m
         exit
      endif
   enddo
   nsing=min(m,nrobs)
   print *,'Number of singular values= ',nsing

   open(10,file='eig.dat')
   do i=1,nrobs
      write(10,'(i4,g13.5)')i,ss(i)/summ
   enddo
   close(10)


!     print *,'b1'
!     print '(8g10.2)',b

   ZZ=transpose(RR)
   bb=matmul(ZZ,bb)
      print *,'b2'
!      print '(8g10.2)',b

   do i=nrobs,nrobs-nsing+1,-1
      bb(i,:)=bb(i,:)/ss(i)
   enddo
      print *,'b3'
!      print '(8g10.2)',b


   deallocate (ipvt)
   deallocate (ss)
   deallocate (H)
   deallocate (RT)
   deallocate (ZZ)
   deallocate (w)
  
end subroutine eksRfact_lag
end module m_eksRfact_lag
