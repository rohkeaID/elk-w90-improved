
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dhmlistl(ngp,ngpq,igpig,igpqig,vgpc,vgpqc,ld,dh)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: ngp,ngpq
integer, intent(in) :: igpig(ngkmax),igpqig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax),vgpqc(3,ngkmax)
integer, intent(in) :: ld
complex(8), intent(inout) :: dh(ld,*)
! local variables
integer iv(3),jv(3),ig,i,j
real(8) vj(3),t1
do j=1,ngp
  jv(:)=ivg(:,igpig(j))
  vj(:)=vgpc(:,j)
  do i=1,ngpq
    iv(:)=ivg(:,igpqig(i))-jv(:)
    ig=ivgig(iv(1),iv(2),iv(3))
    t1=0.5d0*(vgpqc(1,i)*vj(1)+vgpqc(2,i)*vj(2)+vgpqc(3,i)*vj(3))
    dh(i,j)=dh(i,j)+dvsig(ig)+t1*dcfunig(ig)
  end do
end do
return
end subroutine

