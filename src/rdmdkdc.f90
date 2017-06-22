
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmdkdc
! !INTERFACE:
subroutine rdmdkdc
! !USES:
use modmain
use modrdm
! !DESCRIPTION:
!   Calculates the derivative of kinetic energy w.r.t. the second-variational
!   coefficients {\tt evecsv}.
!
! !REVISION HISTORY:
!   Created October 2008 (Sharma)
!EOP
!BOC
implicit none
! allocatable arrays
complex(8), allocatable :: evecsv(:,:),kmat(:,:)
integer ik
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evecsv,kmat)
!$OMP DO
do ik=1,nkpt
  allocate(evecsv(nstsv,nstsv),kmat(nstsv,nstsv))
  call getevecsv(filext,vkl(:,ik),evecsv)
  call getkmat(ik,kmat)
  call zgemm('N','N',nstsv,nstsv,nstsv,zone,kmat,nstsv,evecsv,nstsv,zzero, &
   dkdc(:,:,ik),nstsv)
  deallocate(evecsv,kmat)
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine
!EOC
