
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfmtsm(m,lmmax,nr,ld,rfmt)
implicit none
! arguments
integer, intent(in) :: m,lmmax,nr,ld
real(8), intent(inout) :: rfmt(ld,*)
! local variables
integer lm
do lm=1,lmmax
  call fsmooth(m,nr,ld,rfmt(lm,1))
end do
return
end subroutine

