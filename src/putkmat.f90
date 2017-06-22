
! Copyright (C) 2007-2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putkmat(tfv,tvclcr,ik,vmt,vir)
use modmain
use modmpi
implicit none
! arguments
logical, intent(in) :: tfv,tvclcr
integer, intent(in) :: ik
real(8), intent(in) :: vmt(lmmaxvr,nrcmtmax,natmtot),vir(ngtot)
! local variables
integer ist,ispn,recl
! automatic arrays
integer idx(nstsv)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: kmat(:,:),a(:,:)
! get the eigenvalues/vectors from file for input k-point
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
call getevalsv(filext,vkl(:,ik),evalsv(:,ik))
call getevecfv(filext,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
call getevecsv(filext,vkl(:,ik),evecsv)
! find the matching coefficients
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
do ispn=1,nspnfv
  call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
   sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
end do
! index to all states
do ist=1,nstsv
  idx(ist)=ist
end do
! calculate the wavefunctions for all states of the input k-point
allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir(ngkmax,nspinor,nstsv))
call genwfsv(.false.,.true.,nstsv,idx,ngk(:,ik),igkig(:,:,ik),apwalm,evecfv, &
 evecsv,wfmt,ngkmax,wfir)
deallocate(apwalm,evecfv)
! compute Kohn-Sham potential matrix elements
allocate(kmat(nstsv,nstsv))
if (spinpol) then
  call genvbmatk(vmt,vir,bsmt,bsir,ngk(:,ik),igkig(:,:,ik),wfmt,ngkmax,wfir, &
   kmat)
else
  call genvmatk(vmt,vir,ngk(:,ik),igkig(:,:,ik),wfmt,ngkmax,wfir,kmat)
end if
deallocate(wfir)
! add the DFT+U matrix elements if required
call genvmmtsv(wfmt,kmat)
! negate the potential matrix elements because we have to subtract them
kmat(:,:)=-kmat(:,:)
! add second-variational eigenvalues along the diagonal
do ist=1,nstsv
  kmat(ist,ist)=kmat(ist,ist)+evalsv(ist,ik)
end do
allocate(a(nstsv,nstsv))
! add the Coulomb core matrix elements if required
if (tvclcr) then
  call vclcore(wfmt,a)
  kmat(:,:)=kmat(:,:)+a(:,:)
end if
! rotate kinetic matrix elements to first-variational basis if required
if (tfv) then
  call zgemm('N','C',nstsv,nstsv,nstsv,zone,kmat,nstsv,evecsv,nstsv,zzero,a, &
   nstsv)
  call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,a,nstsv,zzero,kmat, &
   nstsv)
end if
deallocate(evecsv,wfmt,a)
! determine the record length
inquire(iolength=recl) vkl(:,1),nstsv,kmat
!$OMP CRITICAL
open(85,file='KMAT.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
 recl=recl)
write(85,rec=ik) vkl(:,ik),nstsv,kmat
close(85)
!$OMP END CRITICAL
deallocate(kmat)
return
end subroutine

