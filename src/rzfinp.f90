
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

complex(8) function rzfinp(rfmt,rfir,zfmt,zfir)
use modmain
implicit none
! arguments
real(8), intent(in) :: rfmt(lmmaxvr,nrcmtmax,natmtot),rfir(ngtot)
complex(8), intent(in) :: zfmt(lmmaxvr,nrcmtmax,natmtot),zfir(ngtot)
! local variables
integer ias,is,ir
complex(8) zsum
! external functions
complex(8) rzfmtinp
external rzfmtinp
zsum=0.d0
! interstitial contribution
do ir=1,ngtot
  zsum=zsum+(cfunir(ir)*rfir(ir))*zfir(ir)
end do
zsum=zsum*(omega/dble(ngtot))
! muffin-tin contribution
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(is) REDUCTION(+:zsum)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  zsum=zsum+rzfmtinp(nrcmt(is),nrcmtinr(is),rcmt(:,is),r2cmt(:,is), &
   rfmt(:,:,ias),zfmt(:,:,ias))
end do
!$OMP END DO
!$OMP END PARALLEL
rzfinp=zsum
return
end function

