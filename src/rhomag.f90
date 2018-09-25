
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rhomag
use modmain
use modmpi
use modomp
implicit none
! local variables
integer ik,ispn,idm
integer is,ias,n,nthd
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:),evecfv(:,:,:),evecsv(:,:)
! set the charge density and magnetisation to zero
do ias=1,natmtot
  is=idxis(ias)
  rhomt(1:npcmt(is),ias)=0.d0
end do
rhoir(:)=0.d0
do idm=1,ndmag
  do ias=1,natmtot
    is=idxis(ias)
    magmt(1:npcmt(is),ias,idm)=0.d0
  end do
  magir(:,idm)=0.d0
end do
call omp_hold(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(apwalm,evecfv,evecsv,ispn) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
  allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
! get the eigenvectors from file
  call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(filext,ik,vkl(:,ik),evecsv)
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
call omp_free(nthd)
! convert muffin-tin density/magnetisation to spherical harmonics
call rhomagsh
! symmetrise the density
call symrf(nrcmt,nrcmti,npcmt,npmtmax,rhomt,rhoir)
! symmetrise the magnetisation
if (spinpol) call symrvf(.true.,ncmag,nrcmt,nrcmti,npcmt,npmtmax,magmt,magir)
! convert the density from a coarse to a fine radial mesh
call rfmtctof(rhomt)
! convert the magnetisation from a coarse to a fine radial mesh
call omp_hold(ndmag,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do idm=1,ndmag
  call rfmtctof(magmt(:,:,idm))
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
! add densities from each process and redistribute
if (np_mpi.gt.1) then
  n=npmtmax*natmtot
  call mpi_allreduce(mpi_in_place,rhomt,n,mpi_double_precision,mpi_sum,mpicom, &
   ierror)
  call mpi_allreduce(mpi_in_place,rhoir,ngtot,mpi_double_precision,mpi_sum, &
   mpicom,ierror)
  if (spinpol) then
    n=n*ndmag
    call mpi_allreduce(mpi_in_place,magmt,n,mpi_double_precision,mpi_sum, &
     mpicom,ierror)
    n=ngtot*ndmag
    call mpi_allreduce(mpi_in_place,magir,n,mpi_double_precision,mpi_sum, &
     mpicom,ierror)
  end if
end if
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! add the core density to the total density
call rhocore
! calculate the charges
call charge
! calculate the moments
if (spinpol) call moment
! normalise the density
call rhonorm
return
end subroutine

