
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine current
use modmain
use modtddft
use modmpi
implicit none
! local variables
integer ik,ist,i
real(8) t1
! allocatable arrays
complex(8), allocatable :: evecsvt(:,:),pmat(:,:,:)
complex(8), allocatable :: a(:,:),b(:,:)
! zero the total current
curtot(:)=0.d0
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evecsvt,pmat,a,b,i,ist)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(evecsvt(nstsv,nstsv))
  allocate(pmat(nstsv,nstsv,3))
  allocate(a(nstsv,nstsv),b(nstsv,nstsv))
! get the time-dependent Kohn-Sham eigenvectors from file
  call getevecsv(filext,vkl(:,ik),evecsvt)
! get the first-variational basis momentum matrix elements from file
  call getpmat(.true.,vkl(:,ik),pmat)
  do i=1,3
    call zgemm('N','N',nstsv,nstsv,nstsv,zone,pmat(:,:,i),nstsv,evecsvt,nstsv, &
     zzero,a,nstsv)
    call zgemm('C','N',nstsv,nstsv,nstsv,zone,evecsvt,nstsv,a,nstsv,zzero,b, &
     nstsv)
!$OMP CRITICAL
    do ist=1,nstsv
      curtot(i)=curtot(i)+wkpt(ik)*occsv(ist,ik)*dble(b(ist,ist))
    end do
!$OMP END CRITICAL
  end do
  deallocate(evecsvt,pmat,a,b)
end do
!$OMP END DO
!$OMP END PARALLEL
curtot(:)=curtot(:)/omega
! add currents from all mpi process and redistribute
if (np_mpi.gt.1) then
  call mpi_allreduce(mpi_in_place,curtot,3,mpi_double_precision,mpi_sum, &
   mpi_comm_kpt,ierror)
end if
! add vector potential contribution to make current gauge invariant
t1=chgval/(omega*solsc)
curtot(:)=curtot(:)-t1*afieldt(:,itimes)
return
end subroutine

