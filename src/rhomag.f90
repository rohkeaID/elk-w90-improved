
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rhomag
use modmain
use modmpi
implicit none
! local variables
integer ik,ispn,idm,n
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
! set the charge density and magnetisation to zero
rhomt(:,:,:)=0.d0
rhoir(:)=0.d0
if (spinpol) then
  magmt(:,:,:,:)=0.d0
  magir(:,:)=0.d0
end if
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(apwalm,evecfv,evecsv,ispn)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
  allocate(evecfv(nmatmax,nstfv,nspnfv))
  allocate(evecsv(nstsv,nstsv))
! get the eigenvectors from file
  call getevecfv(filext,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(filext,vkl(:,ik),evecsv)
! find the matching coefficients
  do ispn=1,nspnfv
    call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
     sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
  end do
! add to the density and magnetisation
  call rhomagk(ngk(:,ik),igkig(:,:,ik),wkpt(ik),occsv(:,ik),apwalm,evecfv, &
   evecsv)
  deallocate(apwalm,evecfv,evecsv)
end do
!$OMP END DO
!$OMP END PARALLEL
! convert muffin-tin density/magnetisation to spherical harmonics
call rhomagsh
! symmetrise the density
call symrf(lradstp,rhomt,rhoir)
! symmetrise the magnetisation
if (spinpol) call symrvf(lradstp,magmt,magir)
! convert the density from a coarse to a fine radial mesh
call rfmtctof(rhomt)
! convert the magnetisation from a coarse to a fine radial mesh
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
do idm=1,ndmag
  call rfmtctof(magmt(:,:,:,idm))
end do
!$OMP END DO
!$OMP END PARALLEL
! add densities from each process and redistribute
if (np_mpi.gt.1) then
  n=lmmaxvr*nrmtmax*natmtot
  call mpi_allreduce(mpi_in_place,rhomt,n,mpi_double_precision,mpi_sum, &
   mpi_comm_kpt,ierror)
  call mpi_allreduce(mpi_in_place,rhoir,ngtot,mpi_double_precision,mpi_sum, &
   mpi_comm_kpt,ierror)
  if (spinpol) then
    n=n*ndmag
    call mpi_allreduce(mpi_in_place,magmt,n,mpi_double_precision,mpi_sum, &
     mpi_comm_kpt,ierror)
    n=ngtot*ndmag
    call mpi_allreduce(mpi_in_place,magir,n,mpi_double_precision,mpi_sum, &
     mpi_comm_kpt,ierror)
  end if
end if
! synchronise MPI processes
call mpi_barrier(mpi_comm_kpt,ierror)
! add the core density to the total density
call rhocore
! calculate the charges
call charge
! calculate the moments
if (spinpol) call moment
! apply smoothing operation to the density and magnetisation if required
call rhomagsm
! normalise the density
call rhonorm
return
end subroutine

