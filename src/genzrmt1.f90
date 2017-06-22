
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genzrmt1(nrc,nrci,wfmt1,wfmt2,zrhomt)
use modmain
implicit none
! arguments
integer, intent(in) :: nrc,nrci
complex(8), intent(in) :: wfmt1(lmmaxvr,nrcmtmax),wfmt2(lmmaxvr,nrcmtmax)
complex(8), intent(out) :: zrhomt(lmmaxvr,nrcmtmax)
! local variables
integer irc
do irc=1,nrci
  zrhomt(1:lmmaxinr,irc)=conjg(wfmt1(1:lmmaxinr,irc))*wfmt2(1:lmmaxinr,irc)
end do
do irc=nrci+1,nrc
  zrhomt(:,irc)=conjg(wfmt1(:,irc))*wfmt2(:,irc)
end do
return
end subroutine

