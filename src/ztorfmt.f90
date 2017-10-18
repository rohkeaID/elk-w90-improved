
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ztorfmt(nr,nri,zfmt,rfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
complex(8), intent(in) :: zfmt(*)
real(8), intent(out) :: rfmt(*)
! local variables
integer ir,i
i=1
do ir=1,nri
  call ztorflm(lmaxi,zfmt(i),rfmt(i))
  i=i+lmmaxi
end do
do ir=nri+1,nr
  call ztorflm(lmaxo,zfmt(i),rfmt(i))
  i=i+lmmaxo
end do
return
end subroutine

