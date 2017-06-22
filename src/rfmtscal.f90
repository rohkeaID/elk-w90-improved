
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfmtscal(nr,nri,lrstp,c,rfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
integer, intent(in) :: lrstp
real(8), intent(in) :: c
real(8), intent(inout) :: rfmt(lmmaxvr,nr)
! local variables
integer ir
! inner part of muffin-tin
do ir=1,nri,lrstp
  call dscal(lmmaxinr,c,rfmt(:,ir),1)
end do
! outer part of muffin-tin
do ir=nri+lrstp,nr,lrstp
  call dscal(lmmaxvr,c,rfmt(:,ir),1)
end do
return
end subroutine

