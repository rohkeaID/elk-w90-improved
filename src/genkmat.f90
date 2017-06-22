
! Copyright (C) 2007-2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genkmat
! !INTERFACE:
subroutine genkmat(tfv,tvclcr)
! !USES:
use modmain
use modmpi
! !INPUT/OUTPUT PARAMETERS:
!   tfv    : .true. if the matrix elements are to be expressed in the
!            first-variational basis; second-variational otherwise (in,logical)
!   tvclvr : .true. if the non-local Coulomb potential from the core states is
!            to be included in the kinetic matrix elements (in,logical)
! !DESCRIPTION:
!   Computes the kinetic matrix elements in the first- or second-variational
!   basis and stores them in the file {\tt KMAT.OUT}. See routine {\tt putkmat}.
!
! !REVISION HISTORY:
!   Created January 2007 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tfv,tvclcr
! local variables
integer ik,is,ias
! allocatable arrays
real(8), allocatable :: vmt(:,:,:),vir(:)
allocate(vmt(lmmaxvr,nrcmtmax,natmtot),vir(ngtot))
! convert muffin-tin Kohn-Sham potential to spherical coordinates
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(is)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  call rbsht(nrcmt(is),nrcmtinr(is),lradstp,vsmt(:,:,ias),1,vmt(:,:,ias))
end do
!$OMP END DO
!$OMP END PARALLEL
! multiply Kohn-Sham interstitial potential by characteristic function
vir(:)=vsir(:)*cfunir(:)
if (mp_mpi) write(*,*)
! loop over k-points
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
!$OMP CRITICAL
  write(*,'("Info(genkmat): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL
  call putkmat(tfv,tvclcr,ik,vmt,vir)
end do
!$OMP END DO
!$OMP END PARALLEL
deallocate(vmt,vir)
! synchronise MPI processes
call mpi_barrier(mpi_comm_kpt,ierror)
return
end subroutine
!EOC

