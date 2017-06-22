
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rzfmtadd(nr,nri,za,zfmt,rfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
complex(8), intent(in) :: za
complex(8), intent(in) :: zfmt(lmmaxvr,nr)
real(8), intent(inout) :: rfmt(lmmaxvr,nr)
! local variables
integer ir
! add on inner part of muffin-tin
do ir=1,nri
  rfmt(1:lmmaxinr,ir)=rfmt(1:lmmaxinr,ir)+dble(za*zfmt(1:lmmaxinr,ir))
end do
! add on outer part of muffin-tin
do ir=nri+1,nr
  rfmt(:,ir)=rfmt(:,ir)+dble(za*zfmt(:,ir))
end do
return
end subroutine

