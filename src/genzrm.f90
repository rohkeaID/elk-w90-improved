
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genzrm(wfmt11,wfmt12,wfir11,wfir12,wfmt21,wfmt22,wfir21,wfir22, &
 zrhomt,zrhoir,zmagmt,zmagir)
use modmain
implicit none
! arguments
complex(8), intent(in) ::  wfmt11(lmmaxvr,nrcmtmax,natmtot)
complex(8), intent(in) ::  wfmt12(lmmaxvr,nrcmtmax,natmtot)
complex(8), intent(in) ::  wfir11(ngtot),wfir12(ngtot)
complex(8), intent(in) ::  wfmt21(lmmaxvr,nrcmtmax,natmtot)
complex(8), intent(in) ::  wfmt22(lmmaxvr,nrcmtmax,natmtot)
complex(8), intent(in) ::  wfir21(ngtot),wfir22(ngtot)
complex(8), intent(out) :: zrhomt(lmmaxvr,nrcmtmax,natmtot),zrhoir(ngtot)
complex(8), intent(out) :: zmagmt(lmmaxvr,nrcmtmax,natmtot,ndmag)
complex(8), intent(out) :: zmagir(ngtot,ndmag)
! local variables
integer is,ias,ir
complex(8) z11,z12,z21,z22,z1,z2
!-------------------------!
!     muffin-tin part     !
!-------------------------!
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(is)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  call genzrmmt(nrcmt(is),nrcmtinr(is),wfmt11(:,:,ias),wfmt12(:,:,ias), &
   wfmt21(:,:,ias),wfmt22(:,:,ias),zrhomt(:,:,ias),natmtot,zmagmt(:,:,ias,1))
end do
!$OMP END DO
!$OMP END PARALLEL
!---------------------------!
!     interstitial part     !
!---------------------------!
if (ncmag) then
! non-collinear case
  do ir=1,ngtot
    z11=wfir11(ir)
    z12=wfir12(ir)
    z21=wfir21(ir)
    z22=wfir22(ir)
! up-dn spin density
    z1=conjg(z11)*z22
! dn-up spin density
    z2=conjg(z12)*z21
! x-component: up-dn + dn-up
    zmagir(ir,1)=z1+z2
! y-component: i*(dn-up - up-dn)
    z1=z2-z1
    zmagir(ir,2)=cmplx(-aimag(z1),dble(z1),8)
    z1=conjg(z11)*z21
    z2=conjg(z12)*z22
! z-component: up-up - dn-dn
    zmagir(ir,3)=z1-z2
! density: up-up + dn-dn
    zrhoir(ir)=z1+z2
  end do
else
! collinear case
  do ir=1,ngtot
    z1=conjg(wfir11(ir))*wfir21(ir)
    z2=conjg(wfir12(ir))*wfir22(ir)
    zmagir(ir,1)=z1-z2
    zrhoir(ir)=z1+z2
  end do
end if
return
end subroutine

