
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmldbse
use modmain
use modmpi
implicit none
! local variables
integer ik2
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
do ik2=1,nkptnr
! distribute among MPI processes
  if (mod(ik2-1,np_mpi).ne.lp_mpi) cycle
!$OMP CRITICAL
  write(*,'("Info(hmldbse): ",I6," of ",I6," k-points")') ik2,nkptnr
!$OMP END CRITICAL
  call hmldbsek(ik2)
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine

