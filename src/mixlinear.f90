
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine mixlinear(iscl,beta,n,nu,mu,d)
implicit none
! arguments
integer, intent(in) :: iscl
real(8), intent(in) :: beta
integer, intent(in) :: n
real(8), intent(inout) :: nu(n),mu(n)
real(8), intent(out) :: d
! local variables
integer i
real(8) t0,t1
if (n.le.0) return
! initialise mixer
if (iscl.le.0) then
  mu(:)=nu(:)
  d=1.d0
  return
end if
t0=1.d0-beta
d=0.d0
do i=1,n
  t1=nu(i)-mu(i)
  nu(i)=beta*nu(i)+t0*mu(i)
  d=d+t1**2
  mu(i)=nu(i)
end do
d=sqrt(d/dble(n))
return
end subroutine

