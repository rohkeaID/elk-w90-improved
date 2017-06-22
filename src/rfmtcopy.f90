
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfmtcopy(nr,nri,lrstp,rfmt1,rfmt2)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri,lrstp
real(8), intent(in) :: rfmt1(lmmaxvr,nr)
real(8), intent(out) :: rfmt2(lmmaxvr,nr)
! local variables
integer ir
! inner part of muffin-tin
do ir=1,nri,lrstp
  call dcopy(lmmaxinr,rfmt1(:,ir),1,rfmt2(:,ir),1)
end do
! outer part of muffin-tin
do ir=nri+lrstp,nr,lrstp
  call dcopy(lmmaxvr,rfmt1(:,ir),1,rfmt2(:,ir),1)
end do
return
end subroutine

