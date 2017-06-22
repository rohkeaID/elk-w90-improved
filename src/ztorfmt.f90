
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ztorfmt(nr,nri,lrstp1,zfmt,lrstp2,rfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
integer, intent(in) :: lrstp1
complex(8), intent(in) :: zfmt(lmmaxvr,lrstp1,nr)
integer, intent(in) :: lrstp2
real(8), intent(out) :: rfmt(lmmaxvr,lrstp2,nr)
! local variables
integer ir
do ir=1,nri
  call ztorflm(lmaxinr,zfmt(:,1,ir),rfmt(:,1,ir))
end do
do ir=nri+1,nr
  call ztorflm(lmaxvr,zfmt(:,1,ir),rfmt(:,1,ir))
end do
return
end subroutine

