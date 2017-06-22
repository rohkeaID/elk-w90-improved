
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genzrmt2(nrc,nrci,wfmt11,wfmt12,wfmt21,wfmt22,zrhomt)
use modmain
implicit none
! arguments
integer, intent(in) :: nrc,nrci
complex(8), intent(in) :: wfmt11(lmmaxvr,nrcmtmax),wfmt12(lmmaxvr,nrcmtmax)
complex(8), intent(in) :: wfmt21(lmmaxvr,nrcmtmax),wfmt22(lmmaxvr,nrcmtmax)
complex(8), intent(out) :: zrhomt(lmmaxvr,nrcmtmax)
! local variables
integer irc
do irc=1,nrci
  zrhomt(1:lmmaxinr,irc)=conjg(wfmt11(1:lmmaxinr,irc))*wfmt21(1:lmmaxinr,irc) &
                        +conjg(wfmt12(1:lmmaxinr,irc))*wfmt22(1:lmmaxinr,irc)
end do
do irc=nrci+1,nrc
  zrhomt(:,irc)=conjg(wfmt11(:,irc))*wfmt21(:,irc) &
               +conjg(wfmt12(:,irc))*wfmt22(:,irc)
end do
return
end subroutine

