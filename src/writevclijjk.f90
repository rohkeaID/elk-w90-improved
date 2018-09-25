
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writevclijjk
! !INTERFACE:
subroutine writevclijjk
! !USES:
use modmain
use modmpi
use modomp
! !DESCRIPTION:
!   Generates Coulomb matrix elements of the type $(i-jj-k)$ and outputs them to
!   the file {\tt VCLIJJK.OUT}. Also writes the real diagonal of this matrix,
!   $(i-jj-i)$, to {\tt VCLIJJI.OUT}.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! local variables
integer ik,ist,recl,nthd
! allocatable arrays
real(8), allocatable :: vclijji(:,:,:)
complex(8), allocatable :: vclijjk(:,:,:,:)
! determine record length for vclijji and open file
allocate(vclijji(nstsv,nstsv,nkpt))
inquire(iolength=recl) vclijji
deallocate(vclijji)
open(260,file='VCLIJJI.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
! determine record length for vclijjk and open file
allocate(vclijjk(nstsv,nstsv,nstsv,nkpt))
inquire(iolength=recl) vclijjk
deallocate(vclijjk)
open(262,file='VCLIJJK.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
call omp_hold(nkptnr/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(vclijji,vclijjk,ist) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ik=1,nkptnr
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(vclijji(nstsv,nstsv,nkpt))
  allocate(vclijjk(nstsv,nstsv,nstsv,nkpt))
!$OMP CRITICAL(writevclijjk_)
  write(*,'("Info(writevclijjk): ",I6," of ",I6," k-points")') ik,nkptnr
!$OMP END CRITICAL(writevclijjk_)
! calculate Coulomb matrix elements of the type (i-jj-k)
  call genvclijjk(ik,vclijjk)
! make a copy of the diagonal elements (i-jj-i)
  do ist=1,nstsv
    vclijji(ist,:,:)=dble(vclijjk(ist,ist,:,:))
  end do
!$OMP CRITICAL(u260)
  write(260,rec=ik) vclijji
!$OMP END CRITICAL(u260)
!$OMP CRITICAL(u262)
  write(262,rec=ik) vclijjk
!$OMP END CRITICAL(u262)
  deallocate(vclijji,vclijjk)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
close(260)
close(262)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
return
end subroutine
!EOC

