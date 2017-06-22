
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dolpistl(ngp,ngpq,igpig,igpqig,ld,od)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: ngp,ngpq
integer, intent(in) :: igpig(ngkmax),igpqig(ngkmax)
integer, intent(in) :: ld
complex(8), intent(inout) :: od(ld,*)
! local variables
integer iv(3),jv(3),i,j
do j=1,ngp
  jv(:)=ivg(:,igpig(j))
  do i=1,ngpq
    iv(:)=ivg(:,igpqig(i))-jv(:)
    od(i,j)=od(i,j)+dcfunig(ivgig(iv(1),iv(2),iv(3)))
  end do
end do
return
end subroutine

