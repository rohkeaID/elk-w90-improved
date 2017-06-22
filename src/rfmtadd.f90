
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfmtadd(nr,nri,lrstp,rfmt1,rfmt2)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri,lrstp
real(8), intent(in) :: rfmt1(lmmaxvr,nr)
real(8), intent(inout) :: rfmt2(lmmaxvr,nr)
! local variables
integer ir
! inner part of muffin-tin
do ir=1,nri,lrstp
  rfmt2(1:lmmaxinr,ir)=rfmt2(1:lmmaxinr,ir)+rfmt1(1:lmmaxinr,ir)
end do
! outer part of muffin-tin
do ir=nri+lrstp,nr,lrstp
  rfmt2(:,ir)=rfmt2(:,ir)+rfmt1(:,ir)
end do
return
end subroutine

