
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rtozfmt(nr,nri,lrstp1,rfmt,lrstp2,zfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
integer, intent(in) :: lrstp1
real(8), intent(in) :: rfmt(lmmaxvr,lrstp1,nr)
integer, intent(in) :: lrstp2
complex(8), intent(out) :: zfmt(lmmaxvr,lrstp2,nr)
! local variables
integer ir
do ir=1,nri
  call rtozflm(lmaxinr,rfmt(:,1,ir),zfmt(:,1,ir))
end do
do ir=nri+1,nr
  call rtozflm(lmaxvr,rfmt(:,1,ir),zfmt(:,1,ir))
end do
return
end subroutine

