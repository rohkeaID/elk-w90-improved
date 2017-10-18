
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfsht(nr,nri,rfmt1,rfmt2)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: rfmt1(*)
real(8), intent(out) :: rfmt2(*)
! local variables
integer nro,i
! transform the inner part of the muffin-tin
call dgemm('N','N',lmmaxi,nri,lmmaxi,1.d0,rfshti,lmmaxi,rfmt1,lmmaxi,0.d0, &
 rfmt2,lmmaxi)
! transform the outer part of the muffin-tin
nro=nr-nri
if (nro.eq.0) return
i=lmmaxi*nri+1
call dgemm('N','N',lmmaxo,nro,lmmaxo,1.d0,rfshto,lmmaxo,rfmt1(i),lmmaxo,0.d0, &
 rfmt2(i),lmmaxo)
return
end subroutine

