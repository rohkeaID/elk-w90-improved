
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmdedn
! !INTERFACE:
subroutine rdmdedn(dedn)
! !USES:
use modmain
use modrdm
! !INPUT/OUTPUT PARAMETERS:
!   dedn : free energy derivative (out,real(nstsv,nkpt))
! !DESCRIPTION:
!   Calculates the negative of the derivative of total free energy w.r.t.
!   occupation numbers.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! arguments
real(8), intent(out) :: dedn(nstsv,nkpt)
! allocatable
complex(8), allocatable :: evecsv(:,:),c(:,:)
integer ik,ist
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evecsv,c,ist)
!$OMP DO
do ik=1,nkpt
  allocate(evecsv(nstsv,nstsv),c(nstsv,nstsv))
! get eigenvectors from file
  call getevecsv(filext,vkl(:,ik),evecsv)
! kinetic and Coulomb potential contribution
  call zgemm('C','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,dkdc(:,:,ik),nstsv, &
   zzero,c,nstsv)
  do ist=1,nstsv
    dedn(ist,ik)=-(dble(c(ist,ist))+dble(vclmat(ist,ist,ik)))
  end do
  deallocate(evecsv,c)
end do
!$OMP END DO
!$OMP END PARALLEL
! add exchange correlation contribution
call rdmdexcdn(dedn)
! add entropic contribution if needed
if (rdmtemp.gt.0.d0) call rdmdtsdn(dedn)
return
end subroutine
!EOC
