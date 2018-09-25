
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
real(8), intent(in) :: vmt(npcmtmax,natmtot),vir(ngtot)
! local variables
integer ist,ispn,recl,i
! automatic arrays
integer idx(nstsv)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: kmat(:,:),a(:,:)
! get the eigenvalues/vectors from file for input reduced k-point
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
call getevecsv(filext,ik,vkl(:,ik),evecsv)
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
allocate(wfmt(npcmtmax,natmtot,nspinor,nstsv),wfir(ngkmax,nspinor,nstsv))
call genwfsv(.false.,.true.,nstsv,idx,ngridg,igfft,ngk(:,ik),igkig(:,:,ik), &
 apwalm,evecfv,evecsv,wfmt,ngkmax,wfir)
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
! add scissor correction if required
if (scissor.ne.0.d0) then
  do ist=1,nstsv
    if (evalsv(ist,ik).gt.efermi) then
      kmat(ist,ist)=kmat(ist,ist)+scissor
    end if
  end do
end if
! add the Coulomb core matrix elements if required
if (tvclcr) call vclcore(wfmt,kmat)
! rotate kinetic matrix elements to first-variational basis if required
if (tfv) then
  allocate(a(nstsv,nstsv))
  call zgemm('N','C',nstsv,nstsv,nstsv,zone,kmat,nstsv,evecsv,nstsv,zzero,a, &
   nstsv)
  call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,a,nstsv,zzero,kmat, &
   nstsv)
  deallocate(a)
end if
deallocate(evecsv,wfmt)
! determine the record length
inquire(iolength=recl) vkl(:,1),nstsv,kmat
!$OMP CRITICAL(u140)
do i=1,2
  open(140,file='KMAT.OUT',form='UNFORMATTED',access='DIRECT',recl=recl,err=10)
  write(140,rec=ik,err=10) vkl(:,ik),nstsv,kmat
  close(140)
  exit
10 continue
  if (i.eq.2) then
    write(*,*)
    write(*,'("Error(putkmat): unable to write to KMAT.OUT")')
    write(*,*)
    stop
  end if
  close(140)
end do
!$OMP END CRITICAL(u140)
deallocate(kmat)
return
end subroutine

