
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine plotpt2d(cvec,cinv,vpnl,vpl,vppc)
use modmain
implicit none
! arguments
real(8), intent(in) :: cvec(3,3),cinv(3,3)
real(8), intent(out) :: vpnl(3)
real(8), intent(out) :: vpl(3,np2d(1)*np2d(2))
real(8), intent(out) :: vppc(2,np2d(1)*np2d(2))
! local variables
integer ip,i1,i2
real(8) vl1(3),vl2(3)
real(8) vc1(3),vc2(3),vc3(3)
real(8) d1,d2,d12,t1,t2
vl1(:)=vclp2d(:,1)-vclp2d(:,0)
vl2(:)=vclp2d(:,2)-vclp2d(:,0)
call r3mv(cvec,vl1,vc1)
call r3mv(cvec,vl2,vc2)
d1=sqrt(vc1(1)**2+vc1(2)**2+vc1(3)**2)
d2=sqrt(vc2(1)**2+vc2(2)**2+vc2(3)**2)
if ((d1.lt.epslat).or.(d2.lt.epslat)) then
  write(*,*)
  write(*,'("Error(plotpt2d): zero length plotting vectors")')
  write(*,*)
  stop
end if
d12=(vc1(1)*vc2(1)+vc1(2)*vc2(2)+vc1(3)*vc2(3))/(d1*d2)
! vector normal to plane
call r3cross(vc1,vc2,vc3)
t1=sqrt(vc3(1)**2+vc3(2)**2+vc3(3)**2)
if (t1.lt.epslat) then
  write(*,*)
  write(*,'("Error(plotpt2d): 2D plotting plane vectors are collinear")')
  write(*,*)
  stop
end if
vc3(:)=vc3(:)/t1
call r3mv(cinv,vc3,vpnl)
ip=0
do i2=0,np2d(2)-1
  do i1=0,np2d(1)-1
    ip=ip+1
    t1=dble(i1)/dble(np2d(1))
    t2=dble(i2)/dble(np2d(2))
! plot points in 3D space
    vpl(:,ip)=t1*vl1(:)+t2*vl2(:)+vclp2d(:,0)
! plot points on the plane
    vppc(1,ip)=t1*d1+t2*d2*d12
    vppc(2,ip)=t2*d2*sqrt(abs(1.d0-d12**2))
  end do
end do
return
end subroutine

