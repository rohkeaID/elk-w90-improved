
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine curdenk(ik)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
! local variables
integer ispn,jspn,nst,ist,jst
integer is,ia,ias,nrc,nrci,npc
integer igk,ifg,i
real(8) wo
complex(8) z1
! automatic arrays
integer idx(nstsv)
! allocatable arrays
real(8), allocatable :: rfmt(:)
complex(8), allocatable :: apwalm(:,:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfgp(:,:,:)
complex(8), allocatable :: gzfmt(:,:),zfmt1(:),zfmt2(:)
complex(8), allocatable :: zfft1(:),zfft2(:)
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
! find the matching coefficients
do ispn=1,nspnfv
  call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
   sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
end do
! get the eigenvectors from file
call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
call getevecsv(filext,ik,vkl(:,ik),evecsv)
! count and index the occupied states
nst=0
do ist=1,nstsv
  if (abs(occsv(ist,ik)).lt.epsocc) cycle
  nst=nst+1
  idx(nst)=ist
end do
! calculate the second-variational wavefunctions for occupied states
allocate(wfmt(npcmtmax,natmtot,nspinor,nst),wfgp(ngkmax,nspinor,nst))
call genwfsv(.true.,.true.,nst,idx,ngridg,igfft,ngk(:,ik),igkig(:,:,ik), &
 apwalm,evecfv,evecsv,wfmt,ngkmax,wfgp)
deallocate(apwalm,evecfv,evecsv)
!------------------------------------!
!     muffin-tin current density     !
!------------------------------------!
allocate(rfmt(npcmtmax))
allocate(gzfmt(npcmtmax,3),zfmt1(npcmtmax),zfmt2(npcmtmax))
do ist=1,nst
  jst=idx(ist)
  wo=wkpt(ik)*occsv(jst,ik)
  do ispn=1,nspinor
    do is=1,nspecies
      nrc=nrcmt(is)
      nrci=nrcmti(is)
      npc=npcmt(is)
      do ia=1,natoms(is)
        ias=idxas(ia,is)
! compute the gradient of the wavefunction
        call gradzfmt(nrc,nrci,rcmt(:,is),wfmt(:,ias,ispn,ist),npcmtmax,gzfmt)
! convert wavefunction to spherical coordinates and conjugate
        call zbsht(nrc,nrci,wfmt(:,ias,ispn,ist),zfmt1)
        zfmt1(1:npc)=conjg(zfmt1(1:npc))
        do i=1,3
! convert wavefunction gradient to spherical coordinates
          call zbsht(nrc,nrci,gzfmt(:,i),zfmt2)
! compute the partial current density
          rfmt(1:npc)=aimag(zfmt1(1:npc)*zfmt2(1:npc))
!$OMP CRITICAL(currentk_1)
          call daxpy(npc,wo,rfmt,1,cdmt(:,ias,i),1)
!$OMP END CRITICAL(currentk_1)
        end do
      end do
    end do
  end do
end do
deallocate(wfmt,rfmt,gzfmt,zfmt1,zfmt2)
!--------------------------------------!
!     interstitial current density     !
!--------------------------------------!
allocate(zfft1(ngtot),zfft2(ngtot))
do ist=1,nst
  jst=idx(ist)
  wo=wkpt(ik)*occsv(jst,ik)/omega
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
! Fourier transform to real-space and conjugate
    zfft1(:)=0.d0
    do igk=1,ngk(jspn,ik)
      ifg=igfft(igkig(igk,jspn,ik))
      zfft1(ifg)=wfgp(igk,ispn,ist)
    end do
    call zfftifc(3,ngridg,1,zfft1)
    zfft1(:)=conjg(zfft1(:))
    do i=1,3
! compute the gradient of the wavefunction
      zfft2(:)=0.d0
      do igk=1,ngk(jspn,ik)
        ifg=igfft(igkig(igk,jspn,ik))
        z1=wfgp(igk,ispn,ist)
        zfft2(ifg)=vgkc(i,igk,jspn,ik)*cmplx(-aimag(z1),dble(z1),8)
      end do
      call zfftifc(3,ngridg,1,zfft2)
!$OMP CRITICAL(currentk_2)
      cdir(:,i)=cdir(:,i)+wo*aimag(zfft1(:)*zfft2(:))
!$OMP END CRITICAL(currentk_2)
    end do
  end do
end do
deallocate(wfgp,zfft1,zfft2)
return
end subroutine

