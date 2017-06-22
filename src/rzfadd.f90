
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rzfadd(za,zfmt,zfir,rfmt,rfir)
use modmain
implicit none
! arguments
complex(8), intent(in) :: za
complex(8), intent(in) :: zfmt(lmmaxvr,nrcmtmax,natmtot),zfir(ngtot)
real(8), intent(inout) :: rfmt(lmmaxvr,nrcmtmax,natmtot),rfir(ngtot)
! local variables
integer ias,is
! add in muffin-tin region
do ias=1,natmtot
  is=idxis(ias)
  call rzfmtadd(nrcmt(is),nrcmtinr(is),za,zfmt(:,:,ias),rfmt(:,:,ias))
end do
! add in interstitial region
rfir(:)=rfir(:)+dble(za*zfir(:))
return
end subroutine

