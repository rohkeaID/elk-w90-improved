
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvmatk(vmt,vir,ngp,igpig,wfmt,ld,wfir,vmat)
use modmain
implicit none
! arguments
real(8), intent(in) :: vmt(lmmaxvr,nrcmtmax,natmtot),vir(ngtot)
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
complex(8), intent(in) :: wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv)
integer, intent(in) :: ld
complex(8), intent(in) :: wfir(ld,nspinor,nstsv)
complex(8), intent(out) :: vmat(nstsv,nstsv)
! local variables
integer ist,jst,ispn,jspn
integer is,ias,nrc,nrci
integer igp,ifg
complex(8) z1
! allocatable arrays
complex(8), allocatable :: wfmt1(:,:),wfir1(:),z(:)
! external functions
complex(8) zfcmtinp,zdotc
external zfcmtinp,zdotc
! zero the matrix elements
vmat(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfmt1,ias,is,nrc,nrci) &
!$OMP PRIVATE(ispn,ist,z1)
!$OMP DO
do jst=1,nstsv
  allocate(wfmt1(lmmaxvr,nrcmtmax))
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmtinr(is)
    do ispn=1,nspinor
! apply potential to wavefunction
      call vmtapp(nrc,nrci,vmt(:,:,ias),wfmt(:,:,ias,ispn,jst),wfmt1)
      do ist=1,jst
! compute inner product (functions are in spherical coordinates)
        z1=zfcmtinp(nrc,nrci,rcmt(:,is),r2cmt(:,is),wfmt(:,:,ias,ispn,ist), &
         wfmt1)
        vmat(ist,jst)=vmat(ist,jst)+z1
      end do
    end do
  end do
  deallocate(wfmt1)
end do
!$OMP END DO
!$OMP END PARALLEL
!---------------------------!
!     interstitial part     !
!---------------------------!
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfir1,z,ispn,jspn) &
!$OMP PRIVATE(igp,ifg,ist,z1)
!$OMP DO
do jst=1,nstsv
  allocate(wfir1(ngtot),z(ngkmax))
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
! Fourier transform wavefunction to real-space
    wfir1(:)=0.d0
    do igp=1,ngp(jspn)
      ifg=igfft(igpig(igp,jspn))
      wfir1(ifg)=wfir(igp,ispn,jst)
    end do
    call zfftifc(3,ngridg,1,wfir1)
! apply potential to wavefunction
    wfir1(:)=vir(:)*wfir1(:)
! Fourier transform to G+p-space
    call zfftifc(3,ngridg,-1,wfir1)
    do igp=1,ngp(jspn)
      ifg=igfft(igpig(igp,jspn))
      z(igp)=wfir1(ifg)
    end do
    do ist=1,jst
! compute inner product
      z1=zdotc(ngp(jspn),wfir(:,ispn,ist),1,z,1)
      vmat(ist,jst)=vmat(ist,jst)+z1
    end do
  end do
  deallocate(wfir1,z)
end do
!$OMP END DO
!$OMP END PARALLEL
! lower triangular part
do ist=1,nstsv
  do jst=1,ist-1
    vmat(ist,jst)=conjg(vmat(jst,ist))
  end do
end do
return

contains

subroutine vmtapp(nr,nri,vmt,wfmt1,wfmt2)
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: vmt(lmmaxvr,nr)
complex(8), intent(in) :: wfmt1(lmmaxvr,nr)
complex(8), intent(out) :: wfmt2(lmmaxvr,nr)
! local variables
integer iro
wfmt2(1:lmmaxinr,1:nri)=vmt(1:lmmaxinr,1:nri)*wfmt1(1:lmmaxinr,1:nri)
iro=nri+1
wfmt2(:,iro:nr)=vmt(:,iro:nr)*wfmt1(:,iro:nr)
return
end subroutine

end subroutine

