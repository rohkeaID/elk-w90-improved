
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writevclijji
! !INTERFACE:
subroutine writevclijji
! !USES:
use modmain
use modmpi
! !DESCRIPTION:
!   Generates Coulomb matrix elements of the type $(i-jj-i)$ and outputs them to
!   the file {\tt VCLIJJI.OUT}.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! allocatable arrays
real(8), allocatable :: vclijji(:,:,:)
integer recl,ik
! determine record length for vclijji and open file
allocate(vclijji(nstsv,nstsv,nkpt))
inquire(iolength=recl) vclijji
deallocate(vclijji)
open(100,file='VCLIJJI.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
 recl=recl)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(vclijji)
!$OMP DO
do ik=1,nkptnr
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(vclijji(nstsv,nstsv,nkpt))
!$OMP CRITICAL
  write(*,'("Info(writevclijji): ",I6," of ",I6," k-points")') ik,nkptnr
!$OMP END CRITICAL
! calculate Coulomb matrix elements of the type (i-jj-i)
  call genvclijji(ik,vclijji)
!$OMP CRITICAL
  write(100,rec=ik) vclijji
!$OMP END CRITICAL
  deallocate(vclijji)
end do
!$OMP END DO
!$OMP END PARALLEL
close(100)
return
end subroutine
!EOC
