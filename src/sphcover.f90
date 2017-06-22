
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: sphcover
! !INTERFACE:
subroutine sphcover(n,tp)
! !INPUT/OUTPUT PARAMETERS:
!   n  : number of required points (in,integer)
!   tp : (theta, phi) coordinates (out,real(2,n))
! !DESCRIPTION:
!   Produces a set of $N$ points which cover the unit sphere nearly optimally.
!   The points in spherical $(\theta,\phi)$ coordinates are generated using the
!   explicit `golden section' formula:
!   \begin{align*}
!    \theta_k&=\arccos\left[1-\left(k-\tfrac{1}{2}\right)\delta z\right] \\
!    \phi_k&=(k-1)\delta\phi,
!   \end{align*}
!   where $\delta z=2/n$ and $\delta\phi=\pi(1-\sqrt{5})$.
!
! !REVISION HISTORY:
!   Created April 2008 (JKD)
!   Improved covering, October 2009 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(out) :: tp(2,n)
! local variables
integer k
real(8), parameter :: pi=3.1415926535897932385d0
real(8) z,dz,p,dp
if (n.le.0) then
  write(*,*)
  write(*,'("Error(sphcover): n <= 0 : ",I8)') n
  write(*,*)
  stop
end if
dz=2.d0/dble(n)
z=1.d0-dz/2.d0
tp(1,1)=acos(z)
dp=pi*(1.d0-sqrt(5.d0))
p=0.d0
tp(2,1)=p
do k=2,n
  z=z-dz
  tp(1,k)=acos(z)
  p=p+dp
  tp(2,k)=mod(p,2.d0*pi)
end do
return
end subroutine
!EOC
