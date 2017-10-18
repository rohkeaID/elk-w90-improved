
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rtozfmt(nr,nri,rfmt,zfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: rfmt(*)
complex(8), intent(out) :: zfmt(*)
! local variables
integer ir,i
i=1
do ir=1,nri
  call rtozflm(lmaxi,rfmt(i),zfmt(i))
  i=i+lmmaxi
end do
do ir=nri+1,nr
  call rtozflm(lmaxo,rfmt(i),zfmt(i))
  i=i+lmmaxo
end do
return
end subroutine

