
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: polynom
! !INTERFACE:
real(8) function polynom(m,np,xa,ya,c,x)
! !INPUT/OUTPUT PARAMETERS:
!   m  : order of derivative (in,integer)
!   np : number of points to fit (in,integer)
!   xa : abscissa array (in,real(np))
!   ya : ordinate array (in,real(np))
!   c  : work array (out,real(np))
!   x  : evaluation abscissa (in,real)
! !DESCRIPTION:
!   Fits a polynomial of order $n_p-1$ to a set of $n_p$ points. If $m\ge 0$ the
!   function returns the $m$th derviative of the polynomial at $x$, while for
!   $m<0$ the integral of the polynomial from the first point in the array to
!   $x$ is returned.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! argmuments
integer, intent(in) :: m
integer, intent(in) :: np
real(8), intent(in) :: xa(np),ya(np)
real(8), intent(out) :: c(np)
real(8), intent(in) :: x
! local variables
integer i,j,k
real(8) x0,x1,x2,x3,y0,y1,y2,y3
real(8) t0,t1,t2,t3,t4,t5,t6
real(8) c1,c2,c3,sum
! fast evaluations for small np
select case(np)
case(1)
  select case(m)
  case(:-1)
    polynom=ya(1)*(x-xa(1))
  case(0)
    polynom=ya(1)
  case default
    polynom=0.d0
  end select
  return
case(2)
  c1=(ya(2)-ya(1))/(xa(2)-xa(1))
  t1=x-xa(1)
  select case(m)
  case(:-1)
    polynom=t1*(ya(1)+0.5d0*c1*t1)
  case(0)
    polynom=c1*t1+ya(1)
  case(1)
    polynom=c1
  case default
    polynom=0.d0
  end select
  return
case(3)
  x0=xa(1)
  x1=xa(2)-x0
  x2=xa(3)-x0
  y0=ya(1)
  y1=ya(2)-y0
  y2=ya(3)-y0
  t0=1.d0/(x1*x2*(x2-x1))
  t1=x1*y2
  t2=x2*y1
  c1=x2*t2-x1*t1
  c2=t1-t2
  t1=x-x0
  select case(m)
  case(:-1)
    polynom=t1*(y0+t0*t1*(0.5d0*c1+0.3333333333333333333d0*c2*t1))
  case(0)
    polynom=y0+t0*t1*(c1+c2*t1)
  case(1)
    polynom=t0*(2.d0*c2*t1+c1)
  case(2)
    polynom=t0*2.d0*c2
  case default
    polynom=0.d0
  end select
  return
case(4)
  x0=xa(1)
  x1=xa(2)-x0
  x2=xa(3)-x0
  x3=xa(4)-x0
  y0=ya(1)
  y1=ya(2)-y0
  y2=ya(3)-y0
  y3=ya(4)-y0
  t0=1.d0/(x1*x2*x3*(x1-x2)*(x1-x3)*(x2-x3))
  t1=x1*x2*y3
  t2=x2*x3*y1
  t3=x3*x1*y2
  c3=t1*(x1-x2)+t2*(x2-x3)+t3*(x3-x1)
  t6=x3**2
  t5=x2**2
  t4=x1**2
  c2=t1*(t5-t4)+t2*(t6-t5)+t3*(t4-t6)
  c1=t1*(x2*t4-x1*t5)+t2*(x3*t5-x2*t6)+t3*(x1*t6-x3*t4)
  t1=x-x0
  select case(m)
  case(:-1)
    polynom=t1*(y0+t0*t1*(0.5d0*c1+t1*(0.3333333333333333333d0*c2 &
     +0.25d0*c3*t1)))
  case(0)
    polynom=y0+t0*t1*(c1+t1*(c2+c3*t1))
  case(1)
    polynom=t0*(c1+t1*(2.d0*c2+3.d0*c3*t1))
  case(2)
    polynom=t0*(6.d0*c3*t1+2.d0*c2)
  case(3)
    polynom=t0*6.d0*c3
  case default
    polynom=0.d0
  end select
  return
end select
if (np.le.0) then
  write(*,*)
  write(*,'("Error(polynom): np <= 0 : ",I8)') np
  write(*,*)
  stop
end if
if (m.ge.np) then
  polynom=0.d0
  return
end if
! find the polynomial coefficients in divided differences form
c(:)=ya(:)
do i=2,np
  do j=np,i,-1
    c(j)=(c(j)-c(j-1))/(xa(j)-xa(j+1-i))
  end do
end do
! special case m=0
if (m.eq.0) then
  sum=c(1)
  t1=1.d0
  do i=2,np
    t1=t1*(x-xa(i-1))
    sum=sum+c(i)*t1
  end do
  polynom=sum
  return
end if
x0=xa(1)
! convert to standard form
do j=1,np-1
  do i=1,np-j
    k=np-i
    c(k)=c(k)+(x0-xa(k-j+1))*c(k+1)
  end do
end do
if (m.gt.0) then
! take the m th derivative
  do j=1,m
    do i=m+1,np
      c(i)=c(i)*dble(i-j)
    end do
  end do
  t1=c(np)
  t2=x-x0
  do i=np-1,m+1,-1
    t1=t1*t2+c(i)
  end do
  polynom=t1
else
! find the integral
  t1=c(np)/dble(np)
  t2=x-x0
  do i=np-1,1,-1
    t1=t1*t2+c(i)/dble(i)
  end do
  polynom=t1*t2
end if
return
end function
!EOC