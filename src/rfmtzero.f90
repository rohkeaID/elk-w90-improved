
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfmtzero(nr,nri,lrstp,rfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
integer, intent(in) :: lrstp
real(8), intent(out) :: rfmt(lmmaxvr,nr)
! local variables
integer ir
! inner part of muffin-tin
do ir=1,nri,lrstp
  rfmt(1:lmmaxinr,ir)=0.d0
end do
! outer part of muffin-tin
do ir=nri+lrstp,nr,lrstp
  rfmt(:,ir)=0.d0
end do
return
end subroutine

