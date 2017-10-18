
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rotrfmt(rot,nr,nri,rfmt1,rfmt2)
use modmain
implicit none
! arguments
real(8), intent(in) :: rot(3,3)
integer, intent(in) :: nr,nri
real(8), intent(in) :: rfmt1(*)
real(8), intent(out) :: rfmt2(*)
! local variables
integer nro,i
! inner part of muffin-tin
call rotrflm(rot,lmaxi,nri,lmmaxi,rfmt1,rfmt2)
! outer part of muffin-tin
nro=nr-nri
i=lmmaxi*nri+1
call rotrflm(rot,lmaxo,nro,lmmaxo,rfmt1(i),rfmt2(i))
return
end subroutine

