
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: fderiv
! !INTERFACE:
subroutine fderiv(m,n,x,f,g)
! !INPUT/OUTPUT PARAMETERS:
!   m : order of derivative (in,integer)
!   n : number of points (in,integer)
!   x : abscissa array (in,real(n))
!   f : function array (in,real(n))
!   g : (anti-)derivative of f (out,real(n))
! !DESCRIPTION:
!   Given function $f$ defined on a set of points $x_i$ then if $m\ge 0$ this
!   routine computes the $m$th derivative of $f$ at each point. If $m<0$ the
!   anti-derivative of $f$ given by
!   $$ g(x_i)=\int_{x_1}^{x_i} f(x)\,dx $$
!   is calculated. If $m=-1$ then an accurate integral is computed by fitting
!   the function to a clamped cubic spline. When $m=-3$ the fast but low
!   accuracy trapezoidal integration method is used. Simpson's integration,
!   which is slower but more accurate than the trapezoidal method, is used if
!   $m=-2$.
!
! !REVISION HISTORY:
!   Created May 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: m,n
real(8), intent(in) :: x(n),f(n)
real(8), intent(out) :: g(n)
! local variables
integer i
real(8) x0,x1,x2,dx
! automatic arrays
real(8) cf(3,n)
if (n.le.0) then
  write(*,*)
  write(*,'("Error(fderiv): invalid number of points : ",I8)') n
  write(*,*)
  stop
end if
! low accuracy trapezoidal integration
if (m.eq.-3) then
  g(1)=0.d0
  do i=1,n-1
    g(i+1)=g(i)+0.5d0*(x(i+1)-x(i))*(f(i+1)+f(i))
  end do
  return
end if
! medium accuracy Simpson integration
if (m.eq.-2) then
  g(1)=0.d0
  do i=1,n-2
    x0=x(i)
    x1=x(i+1)
    x2=x(i+2)
    g(i+1)=g(i)+(x0-x1)*(f(i+2)*(x0-x1)**2+f(i+1)*(x2-x0)*(x0+2.d0*x1-3.d0*x2) &
     +f(i)*(x2-x1)*(2.d0*x0+x1-3.d0*x2))/(6.d0*(x0-x2)*(x1-x2))
  end do
  x0=x(n)
  x1=x(n-1)
  x2=x(n-2)
  g(n)=g(n-1)+(x1-x0)*(f(n-2)*(x1-x0)**2+f(n)*(x1-x2)*(3.d0*x2-x1-2.d0*x0) &
   +f(n-1)*(x0-x2)*(3.d0*x2-2.d0*x1-x0))/(6.d0*(x2-x1)*(x2-x0))
  return
end if
if (m.eq.0) then
  g(:)=f(:)
  return
end if
if (m.ge.4) then
  g(:)=0.d0
  return
end if
! high accuracy integration/differentiation from spline interpolation
call spline(n,x,1,f,cf)
select case(m)
case(:-1)
  g(1)=0.d0
  do i=1,n-1
    dx=x(i+1)-x(i)
    g(i+1)=g(i)+(((0.25d0*cf(3,i)*dx+0.3333333333333333333d0*cf(2,i))*dx &
     +0.5d0*cf(1,i))*dx+f(i))*dx
  end do
case(1)
  g(:)=cf(1,:)
case(2)
  g(:)=2.d0*cf(2,:)
case(3)
  g(:)=6.d0*cf(3,:)
end select
return
end subroutine
!EOC

