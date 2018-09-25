
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmlxbse
use modmain
use modmpi
use modomp
implicit none
! local variables
integer ik2,nthd
call omp_hold(nkptnr/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ik2=1,nkptnr
! distribute among MPI processes
  if (mod(ik2-1,np_mpi).ne.lp_mpi) cycle
!$OMP CRITICAL(hmlxbse_)
  write(*,'("Info(hmlxbse): ",I6," of ",I6," k-points")') ik2,nkptnr
!$OMP END CRITICAL(hmlxbse_)
  call hmlxbsek(ik2)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
return
end subroutine

