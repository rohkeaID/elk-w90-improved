
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gentauk(ik)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
! local variables
integer ispn,jspn,nst,ist,jst
integer is,ias,nrc,nrci
integer npc,igk,ifg,i
real(8) wo
complex(8) z1
! automatic arrays
integer idx(nstsv)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfgp(:,:,:)
complex(8), allocatable :: gzfmt(:,:),zfmt(:),zfft(:)
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
!-------------------------!
!     muffin-tin part     !
!-------------------------!
allocate(gzfmt(npcmtmax,3),zfmt(npcmtmax))
do ist=1,nst
  jst=idx(ist)
  wo=0.5d0*wkpt(ik)*occsv(jst,ik)
  do ispn=1,nspinor
    do ias=1,natmtot
      is=idxis(ias)
      nrc=nrcmt(is)
      nrci=nrcmti(is)
      npc=npcmt(is)
! compute the gradient of the wavefunction
      call gradzfmt(nrc,nrci,rcmt(:,is),wfmt(:,ias,ispn,ist),npcmtmax,gzfmt)
      do i=1,3
! convert gradient to spherical coordinates
        call zbsht(nrc,nrci,gzfmt(:,i),zfmt)
! add to total in muffin-tin
!$OMP CRITICAL(gentauk_1)
        taumt(1:npc,ias,ispn)=taumt(1:npc,ias,ispn) &
         +wo*(dble(zfmt(1:npc))**2+aimag(zfmt(1:npc))**2)
!$OMP END CRITICAL(gentauk_1)
      end do
    end do
  end do
end do
deallocate(wfmt,gzfmt,zfmt)
!---------------------------!
!     interstitial part     !
!---------------------------!
allocate(zfft(ngtot))
do ist=1,nst
  jst=idx(ist)
  wo=0.5d0*wkpt(ik)*occsv(jst,ik)/omega
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
    do i=1,3
      zfft(:)=0.d0
      do igk=1,ngk(jspn,ik)
        ifg=igfft(igkig(igk,jspn,ik))
        z1=wfgp(igk,ispn,ist)
        zfft(ifg)=vgkc(i,igk,jspn,ik)*cmplx(-aimag(z1),dble(z1),8)
      end do
      call zfftifc(3,ngridg,1,zfft)
!$OMP CRITICAL(gentauk_2)
      tauir(:,ispn)=tauir(:,ispn)+wo*(dble(zfft(:))**2+aimag(zfft(:))**2)
!$OMP END CRITICAL(gentauk_2)
    end do
  end do
end do
deallocate(wfgp,zfft)
return
end subroutine

