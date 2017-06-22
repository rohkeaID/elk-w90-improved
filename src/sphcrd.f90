
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: sphcrd
! !INTERFACE:
subroutine sphcrd(v,r,tp)
! !INPUT/OUTPUT PARAMETERS:
!   v  : input vector (in,real(3))
!   r  : length of v (out,real)
!   tp : (theta, phi) coordinates (out,real(2))
! !DESCRIPTION:
!   Returns the spherical coordinates $(r,\theta,\phi)$ of a vector
!   $$ {\bf v}=(r\sin(\theta)\cos(\phi), r\sin(\theta)\sin(\phi),
!    r\cos(\theta)). $$
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: v(3)
real(8), intent(out) :: r,tp(2)
! local variables
real(8), parameter :: pi=3.1415926535897932385d0
real(8), parameter :: eps=1.d-14
real(8) t1
r=sqrt(v(1)**2+v(2)**2+v(3)**2)
if (r.gt.eps) then
  t1=v(3)/r
  if (t1.ge.1.d0) then
    tp(1)=0.d0
  else if (t1.le.-1.d0) then
    tp(1)=pi
  else
    tp(1)=acos(t1)
  end if
  if ((abs(v(1)).gt.eps).or.(abs(v(2)).gt.eps)) then
    tp(2)=atan2(v(2),v(1))
  else
    tp(2)=0.d0
  end if
else
  tp(1)=0.d0
  tp(2)=0.d0
end if
return
end subroutine
!EOC
