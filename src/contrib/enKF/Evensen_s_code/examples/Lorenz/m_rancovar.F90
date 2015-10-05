module m_rancovar
contains
subroutine rancovar(ens,nn,ave,var,cov,corr)
   use mod_dimensions
   real, intent(in)  :: ens(3,ndim,nrsamp)
   real, intent(in)  :: ave(3,ndim)
   real, intent(in)  :: var(3,ndim)
   real, intent(out) :: cov(3,ndim)
   real, intent(out) :: corr(3,ndim)
   integer, intent(in)  :: nn

   integer i,j
   cov=0.0
   do j=1,nrsamp
      do i=1,ndim
         cov(:,i)=cov(:,i)+(ens(:,i,j)-ave(:,i))*(ens(:,nn,j)-ave(:,nn))
      enddo
   enddo
   cov=(1.0/real(nrsamp-1))*cov
   do i=1,ndim
      corr(:,i)=cov(:,i)/(sqrt(var(:,i)*var(:,nn))+0.0000001)
   enddo

end subroutine rancovar
end module
