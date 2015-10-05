module m_ranmean2
contains
subroutine ranmean2(ens,ave)
   use mod_dimensions
   real, intent(in)  :: ens(3,nrsamp)
   real, intent(out) :: ave(3)
   integer j

   ave(:)=ens(:,1)
   do j=2,nrsamp
      ave(:)=ave(:)+ens(:,j)
   enddo
   ave=(1.0/real(nrsamp))*ave
end subroutine ranmean2
end module
