
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: axangrot
! !INTERFACE:
subroutine axangrot(v,th,rot)
! !INPUT/OUTPUT PARAMETERS:
!   v   : axis vector (in,real)
!   th  : rotation angle (in,real)
!   rot : rotation matrix (out,real(3,3))
! !DESCRIPTION:
!   Determines the $3\times 3$ rotation matrix of a rotation specified by an
!   axis-angle pair following the `right-hand rule'. The axis vector need not be
!   normalised. See {\tt rotaxang} for details.
!
! !REVISION HISTORY:
!   Created February 2014 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: v(3),th
real(8), intent(out) :: rot(3,3)
! local variables
real(8) x,y,z,x2,y2,z2
real(8) xy,xz,yz,cs,sn,t1
x=v(1)
y=v(2)
z=v(3)
t1=sqrt(x**2+y**2+z**2)
! if the axis has zero length then assume the identity
if (t1.lt.1.d-14) then
  rot(:,:)=0.d0
  rot(1,1)=1.d0
  rot(2,2)=1.d0
  rot(3,3)=1.d0
  return
end if
t1=1.d0/t1
x=x*t1
y=y*t1
z=z*t1
x2=x**2
y2=y**2
z2=z**2
xy=x*y
xz=x*z
yz=y*z
cs=cos(th)
sn=sin(th)
t1=1.d0-cs
rot(1,1)=cs+x2*t1
rot(2,1)=xy*t1+z*sn
rot(3,1)=xz*t1-y*sn
rot(1,2)=xy*t1-z*sn
rot(2,2)=cs+y2*t1
rot(3,2)=yz*t1+x*sn
rot(1,3)=xz*t1+y*sn
rot(2,3)=yz*t1-x*sn
rot(3,3)=cs+z2*t1
return
end subroutine
!EOC

