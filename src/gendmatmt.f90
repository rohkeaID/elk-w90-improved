
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gendmatmt
use modmain
use moddftu
use modmpi
implicit none
! local variables
integer ik,ispn,ist,ias,n
real(8) t1
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
complex(8), allocatable :: dmat(:,:,:,:,:)
! zero the density matrix
dmatmt(:,:,:,:,:)=0.d0
! begin parallel loop over k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(apwalm,evecfv,evecsv,dmat) &
!$OMP PRIVATE(ispn,ias,ist,t1)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
  allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
  allocate(dmat(lmmaxdm,nspinor,lmmaxdm,nspinor,nstsv))
! find the matching coefficients
  do ispn=1,nspnfv
    call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
     sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
  end do
! get the eigenvectors and occupancies from file
  call getevecfv(filext,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(filext,vkl(:,ik),evecsv)
  call getoccsv(filext,vkl(:,ik),occsv(:,ik))
! begin loop over atoms and species
  do ias=1,natmtot
    call gendmat(.false.,.false.,0,lmaxdm,ias,ngk(:,ik),apwalm,evecfv,evecsv, &
     lmmaxdm,dmat)
    do ist=1,nstsv
      t1=wkpt(ik)*occsv(ist,ik)
!$OMP CRITICAL
      dmatmt(:,:,:,:,ias)=dmatmt(:,:,:,:,ias)+t1*dmat(:,:,:,:,ist)
!$OMP END CRITICAL
    end do
  end do
  deallocate(apwalm,evecfv,evecsv,dmat)
end do
!$OMP END DO
!$OMP END PARALLEL
! add density matrices from each process and redistribute
if (np_mpi.gt.1) then
  n=((lmmaxdm*nspinor)**2)*natmtot
  call mpi_allreduce(mpi_in_place,dmatmt,n,mpi_double_complex,mpi_sum, &
   mpi_comm_kpt,ierror)
end if
! initialise with symmetry-breaking tensor moments
if (ftmtype.lt.0) then
  dmftm=dmftm*reducebf
  dmatmt=dmatmt+dmftm
endif
! symmetrise the density matrix
call symdmat(lmaxdm,lmmaxdm,dmatmt)
return
end subroutine

