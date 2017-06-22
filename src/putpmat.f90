
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putpmat(tfv,tsv,ik)
use modmain
use modmpi
implicit none
! arguments
logical, intent(in) :: tfv,tsv
integer, intent(in) :: ik
! local variables
integer ist,ispn,i,recl
! automatic arrays
integer idx(nstsv)
! allocatable arrays
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: wfmt(:,:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: pmat(:,:,:),a(:,:)
! index to all states
do ist=1,nstsv
  idx(ist)=ist
end do
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
! get the eigenvectors from file
call getevecfv(filext,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
call getevecsv(filext,vkl(:,ik),evecsv)
! find the matching coefficients
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
do ispn=1,nspnfv
  call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
   sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
end do
! calculate the wavefunctions for all states
allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir(ngkmax,nspinor,nstsv))
call genwfsv(.true.,.true.,nstsv,idx,ngk(:,ik),igkig(:,:,ik),apwalm,evecfv, &
 evecsv,wfmt,ngkmax,wfir)
deallocate(evecfv,apwalm)
! calculate the momentum matrix elements
allocate(pmat(nstsv,nstsv,3))
call genpmatk(ngk(:,ik),igkig(:,:,ik),vgkc(:,:,:,ik),wfmt,wfir,pmat)
deallocate(wfmt,wfir)
! determine the record length
inquire(iolength=recl) vkl(:,1),nstsv,pmat
! write the matrix elements in the second-variational basis if required
if (tsv) then
!$OMP CRITICAL
  open(85,file='PMAT.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
   recl=recl)
  write(85,rec=ik) vkl(:,ik),nstsv,pmat
  close(85)
!$OMP END CRITICAL
end if
! write matrix elements in first-variational basis if required
if (tfv) then
  allocate(a(nstsv,nstsv))
  do i=1,3
    call zgemm('N','C',nstsv,nstsv,nstsv,zone,pmat(:,:,i),nstsv,evecsv,nstsv, &
     zzero,a,nstsv)
    call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,a,nstsv,zzero, &
     pmat(:,:,i),nstsv)
  end do
  deallocate(a)
!$OMP CRITICAL
  open(85,file='PMATFV.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
   recl=recl)
  write(85,rec=ik) vkl(:,ik),nstsv,pmat
  close(85)
!$OMP END CRITICAL
end if
deallocate(evecsv,pmat)
return
end subroutine

