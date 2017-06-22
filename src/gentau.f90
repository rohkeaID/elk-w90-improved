
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gentau(taumt,tauir)
use modmain
use modmpi
implicit none
! arguments
real(8), intent(out) :: taumt(lmmaxvr,nrmtmax,natmtot,nspinor)
real(8), intent(out) :: tauir(ngtot,nspinor)
! local variables
integer ik,ispn,is,ias
integer ir,itp,n
! allocatable arrays
real(8), allocatable :: rfmt(:,:,:),rfir(:)
real(8), allocatable :: rvfmt(:,:,:,:),rvfir(:,:)
! set the kinetic energy density to zero
taumt(:,:,:,:)=0.d0
tauir(:,:)=0.d0
! if wavefunctions do not exist tau cannot be computed
if (iscl.le.1) return
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  call gentauk(ik,taumt,tauir)
end do
!$OMP END DO
!$OMP END PARALLEL
allocate(rfmt(lmmaxvr,nrmtmax,natmtot))
! convert taumt to spherical harmonics
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(is,ispn,ir)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  do ispn=1,nspinor
    do ir=1,nrmt(is),lradstp
      rfmt(:,ir,ias)=taumt(:,ir,ias,ispn)
    end do
    call rfsht(nrcmt(is),nrcmtinr(is),lradstp,rfmt(:,:,ias),lradstp, &
     taumt(:,:,ias,ispn))
  end do
end do
!$OMP END DO
!$OMP END PARALLEL
! symmetrise tau
if (spinpol) then
! spin-polarised case: convert to scalar-vector form
  allocate(rfir(ngtot))
  allocate(rvfmt(lmmaxvr,nrmtmax,natmtot,ndmag))
  allocate(rvfir(ngtot,ndmag))
  rvfmt(:,:,:,1:ndmag-1)=0.d0
  do ias=1,natmtot
    is=idxis(ias)
    do ir=1,nrmt(is),lradstp
      rfmt(:,ir,ias)=taumt(:,ir,ias,1)+taumt(:,ir,ias,2)
      rvfmt(:,ir,ias,ndmag)=taumt(:,ir,ias,1)-taumt(:,ir,ias,2)
    end do
  end do
  rfir(:)=tauir(:,1)+tauir(:,2)
  rvfir(:,1:ndmag-1)=0.d0
  rvfir(:,ndmag)=tauir(:,1)-tauir(:,2)
  call symrf(lradstp,rfmt,rfir)
  call symrvf(lradstp,rvfmt,rvfir)
  do ias=1,natmtot
    is=idxis(ias)
    do ir=1,nrmt(is),lradstp
      taumt(:,ir,ias,1)=0.5d0*(rfmt(:,ir,ias)+rvfmt(:,ir,ias,ndmag))
      taumt(:,ir,ias,2)=0.5d0*(rfmt(:,ir,ias)-rvfmt(:,ir,ias,ndmag))
    end do
  end do
  tauir(:,1)=0.5d0*(rfir(:)+rvfir(:,ndmag))
  tauir(:,2)=0.5d0*(rfir(:)-rvfir(:,ndmag))
  deallocate(rfir,rvfmt,rvfir)
else
! spin-unpolarised case
  call symrf(lradstp,taumt,tauir)
end if
! convert taumt to spherical coordinates
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(is,ispn,ir)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  do ispn=1,nspinor
    do ir=1,nrmt(is),lradstp
      rfmt(:,ir,ias)=taumt(:,ir,ias,ispn)
    end do
    call rbsht(nrcmt(is),nrcmtinr(is),lradstp,rfmt(:,:,ias),lradstp, &
     taumt(:,:,ias,ispn))
  end do
end do
!$OMP END DO
!$OMP END PARALLEL
! convert taumt from a coarse to a fine radial mesh
do ispn=1,nspinor
  call rfmtctof(taumt(:,:,:,ispn))
end do
! add tau from each process and redistribute
if (np_mpi.gt.1) then
  n=lmmaxvr*nrmtmax*natmtot*nspinor
  call mpi_allreduce(mpi_in_place,taumt,n,mpi_double_precision,mpi_sum, &
   mpi_comm_kpt,ierror)
  n=ngtot*nspinor
  call mpi_allreduce(mpi_in_place,tauir,n,mpi_double_precision,mpi_sum, &
   mpi_comm_kpt,ierror)
end if
! add the core contribution
call gentaucr(taumt)
! make sure tau is positive everywhere
do ispn=1,nspinor
  do ias=1,natmtot
    is=idxis(ias)
    do ir=1,nrmt(is)
      do itp=1,lmmaxvr
        if (taumt(itp,ir,ias,ispn).lt.0.d0) taumt(itp,ir,ias,ispn)=0.d0
      end do
    end do
  end do
  do ir=1,ngtot
    if (tauir(ir,ispn).lt.0.d0) tauir(ir,ispn)=0.d0
  end do
end do
deallocate(rfmt)
return
end subroutine

