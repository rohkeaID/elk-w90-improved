
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genevfsv
use modmain
use modmpi
implicit none
! local variables
integer ik,lp
! allocatable arrays
real(8), allocatable :: evalfv(:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
! begin parallel loop over k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evalfv,evecfv,evecsv)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(evalfv(nstfv,nspnfv))
  allocate(evecfv(nmatmax,nstfv,nspnfv))
  allocate(evecsv(nstsv,nstsv))
! solve the first- and second-variational eigenvalue equations
  call eveqn(ik,evalfv,evecfv,evecsv)
! write the eigenvalues/vectors to file
  call putevalfv(filext,ik,evalfv)
  call putevalsv(filext,ik,evalsv(:,ik))
  call putevecfv(filext,ik,evecfv)
  call putevecsv(filext,ik,evecsv)
  deallocate(evalfv,evecfv,evecsv)
end do
!$OMP END DO
!$OMP END PARALLEL
! synchronise MPI processes
call mpi_barrier(mpi_comm_kpt,ierror)
! broadcast eigenvalue array to every process
do ik=1,nkpt
  lp=mod(ik-1,np_mpi)
  call mpi_bcast(evalsv(:,ik),nstsv,mpi_double_precision,lp,mpi_comm_kpt,ierror)
end do
return
end subroutine

