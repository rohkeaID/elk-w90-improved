
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genpmat(tfv,tsv)
use modmain
use modmpi
use modomp
implicit none
! arguments
logical, intent(in) :: tfv,tsv
! local variables
integer ik,nthd
if (mp_mpi) write(*,*)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
call omp_hold(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
!$OMP CRITICAL(genpmat_)
  write(*,'("Info(genpmat): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL(genpmat_)
  call putpmat(tfv,tsv,ik)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
return
end subroutine

