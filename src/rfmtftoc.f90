
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfmtftoc(nr,nri,rfmt1,rfmt2)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: rfmt1(*)
real(8), intent(out) :: rfmt2(*)
! local variables
integer ir,i1,i2,j1
i1=1
i2=1
j1=lmmaxi*lradstp
do ir=1,nri,lradstp
  call dcopy(lmmaxi,rfmt1(i1),1,rfmt2(i2),1)
  i1=i1+j1
  i2=i2+lmmaxi
end do
i1=i1+(lradstp-1)*(lmmaxo-lmmaxi)
j1=lmmaxo*lradstp
do ir=nri+lradstp,nr,lradstp
  call dcopy(lmmaxo,rfmt1(i1),1,rfmt2(i2),1)
  i1=i1+j1
  i2=i2+lmmaxo
end do
return
end subroutine

