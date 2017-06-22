
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvmat(vmt,vir,vmat)
! generates potential matrix elements for all states and k-points
use modmain
use modmpi
implicit none
! arguments
real(8), intent(in) :: vmt(lmmaxvr,nrmtmax,natmtot),vir(ngtot)
complex(8), intent(out) :: vmat(nstsv,nstsv,nkpt)
! local variables
integer ik,ist,ispn
integer is,ias,n,lp
! automatic arrays
integer idx(nstsv)
! allocatable arrays
real(8), allocatable :: vmt1(:,:,:),vir1(:)
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:,:),wfir(:,:,:)
! allocate local arrays
allocate(vmt1(lmmaxvr,nrcmtmax,natmtot),vir1(ngtot))
! convert muffin-tin potential to spherical coordinates
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(is)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  call rbsht(nrcmt(is),nrcmtinr(is),lradstp,vmt(:,:,ias),1,vmt1(:,:,ias))
end do
!$OMP END DO
!$OMP END PARALLEL
! multiply interstitial potential by characteristic function
vir1(:)=vir(:)*cfunir(:)
! index to all states
do ist=1,nstsv
  idx(ist)=ist
end do
! loop over k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(apwalm,evecfv,evecsv) &
!$OMP PRIVATE(wfmt,wfir,ispn)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
  allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
  allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
  allocate(wfir(ngkmax,nspinor,nstsv))
! find the matching coefficients
  do ispn=1,nspnfv
    call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
     sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
  end do
! get the eigenvectors from file
  call getevecfv(filext,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(filext,vkl(:,ik),evecsv)
! calculate the wavefunctions for all states of the input k-point
  call genwfsv(.false.,.true.,nstsv,idx,ngk(:,ik),igkig(:,:,ik),apwalm, &
   evecfv,evecsv,wfmt,ngkmax,wfir)
  call genvmatk(vmt1,vir1,ngk(:,ik),igkig(:,:,ik),wfmt,ngkmax,wfir,vmat(:,:,ik))
  deallocate(apwalm,evecfv,evecsv,wfmt,wfir)
end do
!$OMP END DO
!$OMP END PARALLEL
! broadcast matrix elements to every process
if (np_mpi.gt.1) then
  n=nstsv*nstsv
  do ik=1,nkpt
    lp=mod(ik-1,np_mpi)
    call mpi_bcast(vmat(:,:,ik),n,mpi_double_complex,lp,mpi_comm_kpt,ierror)
  end do
end if
deallocate(vmt1,vir1)
return
end subroutine

