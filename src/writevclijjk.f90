
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
integer ik,ist,recl
! allocatable arrays
real(8), allocatable :: vclijji(:,:,:)
complex(8), allocatable :: vclijjk(:,:,:,:)
! determine record length for vclijji and open file
allocate(vclijji(nstsv,nstsv,nkpt))
inquire(iolength=recl) vclijji
deallocate(vclijji)
open(100,file='VCLIJJI.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
 recl=recl)
! determine record length for vclijjk and open file
allocate(vclijjk(nstsv,nstsv,nstsv,nkpt))
inquire(iolength=recl) vclijjk
deallocate(vclijjk)
open(101,file='VCLIJJK.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
 recl=recl)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(vclijji,vclijjk,ist)
!$OMP DO
do ik=1,nkptnr
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(vclijji(nstsv,nstsv,nkpt))
  allocate(vclijjk(nstsv,nstsv,nstsv,nkpt))
!$OMP CRITICAL
  write(*,'("Info(writevclijjk): ",I6," of ",I6," k-points")') ik,nkptnr
!$OMP END CRITICAL
! calculate Coulomb matrix elements of the type (i-jj-k)
  call genvclijjk(ik,vclijjk)
! make a copy of the diagonal elements (i-jj-i)
  do ist=1,nstsv
    vclijji(ist,:,:)=dble(vclijjk(ist,ist,:,:))
  end do
!$OMP CRITICAL
  write(100,rec=ik) vclijji
  write(101,rec=ik) vclijjk
!$OMP END CRITICAL
  deallocate(vclijji,vclijjk)
end do
!$OMP END DO
!$OMP END PARALLEL
close(100)
close(101)
return
end subroutine
!EOC
