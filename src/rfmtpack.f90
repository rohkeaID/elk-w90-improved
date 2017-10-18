
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfmtpack(tpack,nr,nri,rfmt1,rfmt2)
use modmain
implicit none
! arguments
logical, intent(in) :: tpack
integer, intent(in) :: nr,nri
real(8), intent(in) :: rfmt1(*)
real(8), intent(out) :: rfmt2(*)
! local variables
integer ir,i,j,k
i=1
j=1
if (tpack) then
  do ir=1,nri
    call dcopy(lmmaxi,rfmt1(i),1,rfmt2(j),1)
    i=i+lmmaxo
    j=j+lmmaxi
  end do
else
  do ir=1,nri
    call dcopy(lmmaxi,rfmt1(i),1,rfmt2(j),1)
    i=i+lmmaxi
    k=j+lmmaxi
    j=j+lmmaxo
    rfmt2(k:j-1)=0.d0
  end do
end if
k=lmmaxo*(nr-nri)
call dcopy(k,rfmt1(i),1,rfmt2(j),1)
return
end subroutine

