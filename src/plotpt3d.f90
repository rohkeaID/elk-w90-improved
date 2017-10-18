
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine plotpt3d(vpl)
use modmain
implicit none
! arguments
real(8), intent(out) :: vpl(3,np3d(1)*np3d(2)*np3d(3))
! local variables
integer ip,i1,i2,i3
real(8) v1(3),v2(3),v3(3)
real(8) t1,t2,t3
! generate 3D grid from corner vectors
v1(:)=vclp3d(:,1)-vclp3d(:,0)
v2(:)=vclp3d(:,2)-vclp3d(:,0)
v3(:)=vclp3d(:,3)-vclp3d(:,0)
ip=0
do i3=0,np3d(3)-1
  t3=dble(i3)/dble(np3d(3))
  do i2=0,np3d(2)-1
    t2=dble(i2)/dble(np3d(2))
    do i1=0,np3d(1)-1
      t1=dble(i1)/dble(np3d(1))
      ip=ip+1
      vpl(:,ip)=t1*v1(:)+t2*v2(:)+t3*v3(:)+vclp3d(:,0)
    end do
  end do
end do
return
end subroutine

