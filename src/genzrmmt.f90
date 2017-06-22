
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genzrmmt(nrc,nrci,wfmt11,wfmt12,wfmt21,wfmt22,zrhomt,ld,zmagmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nrc,nrci
complex(8), intent(in) :: wfmt11(lmmaxvr,nrcmtmax),wfmt12(lmmaxvr,nrcmtmax)
complex(8), intent(in) :: wfmt21(lmmaxvr,nrcmtmax),wfmt22(lmmaxvr,nrcmtmax)
complex(8), intent(out) :: zrhomt(lmmaxvr,nrcmtmax)
integer, intent(in) :: ld
complex(8), intent(out) :: zmagmt(lmmaxvr,nrcmtmax,ld,ndmag)
! local variables
integer irc,itp
complex(8) z11,z12,z21,z22,z1,z2
if (ncmag) then
! non-collinear case
  do irc=1,nrci
    do itp=1,lmmaxinr
      z11=wfmt11(itp,irc)
      z12=wfmt12(itp,irc)
      z21=wfmt21(itp,irc)
      z22=wfmt22(itp,irc)
! up-dn spin density
      z1=conjg(z11)*z22
! dn-up spin density
      z2=conjg(z12)*z21
! x-component: up-dn + dn-up
      zmagmt(itp,irc,1,1)=z1+z2
! y-component: i*(dn-up - up-dn)
      z1=z2-z1
      zmagmt(itp,irc,1,2)=cmplx(-aimag(z1),dble(z1),8)
      z1=conjg(z11)*z21
      z2=conjg(z12)*z22
! z-component: up-up - dn-dn
      zmagmt(itp,irc,1,3)=z1-z2
! density: up-up + dn-dn
      zrhomt(itp,irc)=z1+z2
    end do
  end do
  do irc=nrci+1,nrc
    do itp=1,lmmaxvr
      z11=wfmt11(itp,irc)
      z12=wfmt12(itp,irc)
      z21=wfmt21(itp,irc)
      z22=wfmt22(itp,irc)
      z1=conjg(z11)*z22
      z2=conjg(z12)*z21
      zmagmt(itp,irc,1,1)=z1+z2
      z1=z2-z1
      zmagmt(itp,irc,1,2)=cmplx(-aimag(z1),dble(z1),8)
      z1=conjg(z11)*z21
      z2=conjg(z12)*z22
      zmagmt(itp,irc,1,3)=z1-z2
      zrhomt(itp,irc)=z1+z2
    end do
  end do
else
! collinear case
  do irc=1,nrci
    do itp=1,lmmaxinr
      z1=conjg(wfmt11(itp,irc))*wfmt21(itp,irc)
      z2=conjg(wfmt12(itp,irc))*wfmt22(itp,irc)
      zmagmt(itp,irc,1,1)=z1-z2
      zrhomt(itp,irc)=z1+z2
    end do
  end do
  do irc=nrci+1,nrc
    do itp=1,lmmaxvr
      z1=conjg(wfmt11(itp,irc))*wfmt21(itp,irc)
      z2=conjg(wfmt12(itp,irc))*wfmt22(itp,irc)
      zmagmt(itp,irc,1,1)=z1-z2
      zrhomt(itp,irc)=z1+z2
    end do
  end do
end if
return
end subroutine

