
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
use modomp
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
integer recl,ik,nthd
! determine record length for vclijji and open file
allocate(vclijji(nstsv,nstsv,nkpt))
inquire(iolength=recl) vclijji
deallocate(vclijji)
open(260,file='VCLIJJI.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
call omp_hold(nkptnr/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(vclijji) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ik=1,nkptnr
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(vclijji(nstsv,nstsv,nkpt))
!$OMP CRITICAL(writevclijji_)
  write(*,'("Info(writevclijji): ",I6," of ",I6," k-points")') ik,nkptnr
!$OMP END CRITICAL(writevclijji_)
! calculate Coulomb matrix elements of the type (i-jj-i)
  call genvclijji(ik,vclijji)
!$OMP CRITICAL(u260)
  write(260,rec=ik) vclijji
!$OMP END CRITICAL(u260)
  deallocate(vclijji)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
close(260)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
return
end subroutine
!EOC

