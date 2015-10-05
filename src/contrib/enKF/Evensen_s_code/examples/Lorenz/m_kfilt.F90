module m_kfilt
contains
subroutine kfilt(A,z,p,vars)
   use mod_dimensions
   use m_ranmean2
   use m_random2
   implicit none
   type(variances), intent(in) :: vars
   real, intent(inout) :: A(3,nrsamp)
   real, intent(in)    :: z(nrobs)
   real, intent(inout) :: p(3)

   integer i,j,m,mm,info
   real rcond,summ,summ2
   logical zlz,svd,lud
   real dd(3)          ! workspace

   integer, allocatable :: ipvt(:)
   real, allocatable :: K(:,:)
   real, allocatable :: H(:,:)
   real, allocatable :: RR(:,:)
   real, allocatable :: RT(:)
   real, allocatable :: w(:)
   real, allocatable :: RS(:,:)
   real, allocatable :: ZZ(:,:)
   real, allocatable :: tmp(:,:)
   real, allocatable :: d(:)
   real, allocatable :: ss(:)
   real, allocatable :: beta(:)

   zlz=.true.
   svd=.false.
   lud=.false.


   allocate (ipvt(3))
   allocate (K(nrobs,3), H(nrobs,3), RR(nrobs,nrobs))
   allocate (w(3), RS(nrobs,nrobs), ZZ(nrobs,nrobs), tmp(nrobs,nrsamp))
   allocate (d(nrobs), ss(nrobs), beta(nrobs))
   allocate (RT(nrobs*nrobs))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   RR=0.0
   H=0.0
   do m=1,nrobs
      RR(m,m)=vars%mes   ! 3x3 diagonal measurment err covar mat.
      H(m,m)=1.0        ! 3x3 identity matrix
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculating $K=PH' = AA'H'/nrsamp$, i.e. the representers
! where $A=[ens(1)-m,...,ens(nrsamp)-m]$ and m is the mean 


   call ranmean2(A,w)
   do i=1,nrsamp
      A(:,i)=A(:,i)-w(:)
   enddo

   tmp=matmul(H,A)
   print *,'kfilt tmp'
!      print '(4g12.4)',tmp(1:nrobs,1:20)


!      K=matmul(tmp,transpose(A))
   call sgemm('n','t',nrobs,3,nrsamp,1.0,tmp,nrobs,A,3,0.0,K,nrobs)

   print *,'kfilt K'

   K=(1.0/float(nrsamp))*K

   do i=1,nrsamp
      A(:,i)=A(:,i)+w(:)
   enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   RR=RR+(1.0/float(nrsamp))*matmul(tmp,transpose(tmp))
   print *,'kfilt RR'
!      print '(4g12.4)',RR(1:nrobs,1:nrobs)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  if(zlz) then
! Calculating the Eigendecomposition of RR using eispack */
!      call rsm(nrobs,nrobs,RR,ss,nrobs,ZZ,w,ipvt,info)
!!!!   m=0
!!!!   do i=1,nrobs
!!!!   do j=i,nrobs
!!!!      m=m+1
!!!!      RT(m)=RR(j,i)
!!!!   enddo
!!!!   enddo
!!!!   call sspev(1,RT,ss,ZZ,nrobs,nrobs,w,3)
!!!!   print *,'sigma'
!!!!   print '(10g11.3)',ss
!!!!   K=matmul(transpose(ZZ),K)
!!!!
!!!!
!!!!   summ=sum(ss)
!!!!   summ2=0.0
!!!!   do i=nrobs,1,-1
!!!!      summ2=summ2+ss(i)
!!!!      if (summ2/summ > 0.95)  then
!!!!         m=nrobs-i+1
!!!!         print *,'Number of singular values= ',m
!!!!         exit
!!!!      endif
!!!!   enddo
!!!!
!!!!   open(10,file='eig.dat')
!!!!   do i=1,nrobs
!!!!      write(10,'(i4,g13.5)')i,ss(i)/summ
!!!!   enddo
!!!!   close(10)
!!!!
!!!!
!!!!   d=z-matmul(H,p)
!!!!   b=matmul(d,ZZ)
!!!!   do i=nrobs,nrobs-m+1,-1
!!!!      beta(i)=beta(i)/ss(i)
!!!!      p=p+beta(i)*K(i,:)
!!!!   enddo
!!!!
!!!!   do j=1,nrsamp
!!!!      d=z-matmul(H,A(:,j))
!!!!      beta=matmul(d,ZZ)
!!!!      do i=nrobs,nrobs-m+1,-1
!!!!         beta(i)=beta(i)/ss(i)
!!!!         A(:,j)=A(:,j)+beta(i)*K(i,:)
!!!!      enddo
!!!!   enddo
!!!!
!!!!   elseif (svd) then
!!!!! Calculating the SVD of RR using ssvdc from Linpack */
!!!!   mm=nrobs
!!!!   call ssvdc(RR,mm,mm,mm,ss,beta,ZZ,mm,RS,mm,w,11,info)
!!!!   if (info != 0) then
!!!!      print *,'kfilt: info=',info
!!!!   endif
!!!!
!!!!
!!!!   K=matmul(transpose(RS),K)
!!!!
!!!!   j=0
!!!!   do i=1,nrobs
!!!!      w(:)=K(i,:)
!!!!      write(tag3,'(i2.2)')i
!!!!      call prt(pdx,w,'datak'\/tag)
!!!!   enddo
!!!!
!!!!
!!!!   summ=sum(ss)
!!!!   summ2=0.0
!!!!   do i=1,nrobs
!!!!      summ2=summ2+ss(i)
!!!!      if (summ2/summ > 0.95)  then
!!!!         m=i
!!!!         print *,'Number of singular values= ',m
!!!!         exit
!!!!      endif
!!!!   enddo
!!!!
!!!!   open(10,file='eig.dat')
!!!!   do i=1,nrobs
!!!!      write(10,'(i4,g13.5)')i,ss(i)/summ
!!!!   enddo
!!!!   close(10)
!!!!
!!!!   
!!!!   d=z-matmul(H,p)
!!!!   beta=matmul(d,ZZ)
!!!!
!!!!   do i=1,m
!!!!      beta(i)=beta(i)/ss(i)
!!!!      p=p+beta(i)*K(i,:)
!!!!   enddo
!!!!
!!!!   do j=1,nrsamp
!!!!      d=z-matmul(H,A(:,j))
!!!!      beta=matmul(d,ZZ)
!!!!      do i=1,m
!!!!         beta(i)=beta(i)/ss(i)
!!!!         A(:,j)=A(:,j)+beta(i)*K(i,:)
!!!!      enddo
!!!!   enddo
!!!!
!!!!
!!!!
!!!!  elseif (lud) then
! Standard solver for Coefficients beta */
   call spoco(RR,nrobs,nrobs,rcond,w,info)
   write(*,*)'rcond=',rcond

   d=z-matmul(H,p)
   call sposl(RR,nrobs,nrobs,d)
   p=p+matmul(d,K) 


   do j=1,nrsamp
      call random2(d,3)
      d=z+sqrt(vars%mes)*d-matmul(H,A(:,j))
      call sposl(RR,nrobs,nrobs,d)
      A(:,j)=A(:,j)+matmul(d,K) 
   enddo
!!!!   endif


   deallocate (RT)
   deallocate (ipvt)
   deallocate (K, H, RR)
   deallocate (w, RS, ZZ, tmp)
   deallocate (d, ss, beta)

end subroutine kfilt
end module m_kfilt


