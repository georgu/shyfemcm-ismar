module m_ranvar
contains
subroutine ranvar(ens,ave,var)
   use mod_dimensions
   real, intent(in)  :: ens(3*ndim,nrsamp)
   real, intent(in)  :: ave(3*ndim)
   real, intent(out) :: var(3*ndim)
   integer j

   var=0.0
   do j=1,nrsamp
       var(:)=var(:)+(ens(:,j)-ave(:))*(ens(:,j)-ave(:))
   enddo
   var=(1.0/real(nrsamp-1))*var
end subroutine ranvar
end module
