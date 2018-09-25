
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine r3vo(x,y)
implicit none
! arguments
real(8), intent(in) :: x(3)
real(8), intent(inout) :: y(3)
! local variables
real(8) t1,t2
t1=x(1)**2+x(2)**2+x(3)**2
if (t1.lt.1.d-8) return
t2=(x(1)*y(1)+x(2)*y(2)+x(3)*y(3))/t1
y(:)=y(:)-t2*x(:)
return
end subroutine
