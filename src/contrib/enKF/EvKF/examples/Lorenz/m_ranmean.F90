module m_ranmean
contains
subroutine ranmean(ens,ave)
   use mod_dimensions
   real, intent(in)  :: ens(3*ndim,nrsamp)
   real, intent(out) :: ave(3*ndim)
   integer j

   ave(:)=ens(:,1)
   do j=2,nrsamp
      ave(:)=ave(:)+ens(:,j)
   enddo
   ave=(1.0/real(nrsamp))*ave
end subroutine ranmean
end module
