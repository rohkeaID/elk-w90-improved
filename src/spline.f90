
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: spline
! !INTERFACE:
subroutine spline(n,x,ld,f,cf)
! !INPUT/OUTPUT PARAMETERS:
!   n  : number of points (in,integer)
!   x  : abscissa array (in,real(n))
!   ld : leading dimension (in,integer)
!   f  : input data array (in,real(ld,n))
!   cf : cubic spline coefficients (out,real(3,n))
! !DESCRIPTION:
!   Calculates the coefficients of a cubic spline fitted to input data. In other
!   words, given a set of data points $f_i$ defined at $x_i$, where
!   $i=1\ldots n$, the coefficients $c_j^i$ are determined such that
!   $$ y_i(x)=f_i+c_1^i(x-x_i)+c_2^i(x-x_i)^2+c_3^i(x-x_i)^3, $$
!   is the interpolating function for $x\in[x_i,x_{i+1})$. The coefficients are
!   determined piecewise by fitting a cubic polynomial to adjacent points.
!
! !REVISION HISTORY:
!   Created November 2011 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: x(n)
integer, intent(in) :: ld
real(8), intent(in) :: f(ld,n)
real(8), intent(out) :: cf(3,n)
! local variables
integer i
real(8) x0,x1,x2,x3,y0,y1,y2,y3
real(8) c1,c2,c3,t0,t1,t2,t3,t4,t5,t6
if (n.le.0) then
  write(*,*)
  write(*,'("Error(spline): n <= 0 : ",I8)') n
  write(*,*)
  stop
end if
if (n.eq.1) then
  cf(:,1)=0.d0
  return
end if
if (n.eq.2) then
  cf(1,1)=(f(1,2)-f(1,1))/(x(2)-x(1))
  cf(2:3,1)=0.d0
  cf(1,2)=cf(1,1)
  cf(2:3,2)=0.d0
  return
end if
if (n.eq.3) then
  x0=x(1)
  x1=x(2)-x0
  x2=x(3)-x0
  y0=f(1,1)
  y1=f(1,2)-y0
  y2=f(1,3)-y0
  t0=1.d0/(x1*x2*(x2-x1))
  t1=x1*y2
  t2=x2*y1
  c1=t0*(x2*t2-x1*t1)
  c2=t0*(t1-t2)
  cf(1,1)=c1
  cf(2,1)=c2
  cf(3,1)=0.d0
  t3=2.d0*c2
  cf(1,2)=c1+t3*x1
  cf(2,2)=c2
  cf(3,2)=0.d0
  cf(1,3)=c1+t3*x2
  cf(2,3)=c2
  cf(3,3)=0.d0
  return
end if
x0=x(1)
x1=x(2)-x0
x2=x(3)-x0
x3=x(4)-x0
y0=f(1,1)
y1=f(1,2)-y0
y2=f(1,3)-y0
y3=f(1,4)-y0
t0=1.d0/(x1*x2*x3*(x1-x2)*(x1-x3)*(x2-x3))
t1=x1*x2*y3
t2=x2*x3*y1
t3=x3*x1*y2
t4=x1**2
t5=x2**2
t6=x3**2
y1=t3*t6-t1*t5
y3=t2*t5-t3*t4
y2=t1*t4-t2*t6
c1=t0*(x1*y1+x2*y2+x3*y3)
c2=-t0*(y1+y2+y3)
c3=t0*(t1*(x1-x2)+t2*(x2-x3)+t3*(x3-x1))
cf(1,1)=c1
cf(2,1)=c2
cf(3,1)=c3
cf(1,2)=c1+2.d0*c2*x1+3.d0*c3*t4
cf(2,2)=c2+3.d0*c3*x1
cf(3,2)=c3
if (n.eq.4) then
  cf(1,3)=c1+2.d0*c2*x2+3.d0*c3*t5
  cf(2,3)=c2+3.d0*c3*x2
  cf(3,3)=c3
  cf(1,4)=c1+2.d0*c2*x3+3.d0*c3*t6
  cf(2,4)=c2+3.d0*c3*x3
  cf(3,4)=c3
  return
end if
do i=3,n-2
  x0=x(i)
  x1=x(i-1)-x0
  x2=x(i+1)-x0
  x3=x(i+2)-x0
  y0=f(1,i)
  y1=f(1,i-1)-y0
  y2=f(1,i+1)-y0
  y3=f(1,i+2)-y0
  t1=x1*x2*y3
  t2=x2*x3*y1
  t3=x3*x1*y2
  t0=1.d0/(x1*x2*x3*(x1-x2)*(x1-x3)*(x2-x3))
  c3=t0*(t1*(x1-x2)+t2*(x2-x3)+t3*(x3-x1))
  t4=x1**2
  t5=x2**2
  t6=x3**2
  y1=t3*t6-t1*t5
  y2=t1*t4-t2*t6
  y3=t2*t5-t3*t4
  cf(1,i)=t0*(x1*y1+x2*y2+x3*y3)
  cf(2,i)=-t0*(y1+y2+y3)
  cf(3,i)=c3
end do
c1=cf(1,n-2)
c2=cf(2,n-2)
c3=cf(3,n-2)
cf(1,n-1)=c1+2.d0*c2*x2+3.d0*c3*t5
cf(2,n-1)=c2+3.d0*c3*x2
cf(3,n-1)=c3
cf(1,n)=c1+2.d0*c2*x3+3.d0*c3*t6
cf(2,n)=c2+3.d0*c3*x3
cf(3,n)=c3
return
end subroutine
!EOC

