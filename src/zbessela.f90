
! Copyright (C) 2008 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zbessela
! !INTERFACE:
subroutine zbessela(lmax,x,a)
! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum order of Bessel function (in,integer)
!   x    : real argument (in,real)
!   a    : array of returned values (out,real(0:lmax))
! !DESCRIPTION:
!   Computes variations of the spherical Bessel function, $a_l(x)=i^lj_l(ix)$,
!   for real argument $x$ and $l=0,1,\ldots,l_{\rm max}$. The recursion relation
!   $$ j_{l+1}(x)=\frac{2l+1}{x}j_l(x)+j_{l-1}(x) $$
!   is used upwards. For starting values there are
!   $$ a_0(x)=\frac{\sinh(x)}{x};\qquad a_1(x)=\frac{a_0(x)-\cosh(x)}{x} $$.
!   For $x\ll 1$ the asymtotic forms
!   $$ a_l(x)\approx\frac{(-x)^l}{(2l+1)!!} $$
!   are used.
!
! !REVISION HISTORY:
!   Created April 2008 from sbessel routine (Lars Nordstrom)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lmax
real(8), intent(in) :: x
real(8), intent(out) :: a(0:lmax)
! local variables
integer l
real(8) xi,a0,a1,at,t1,t2,xmin
if ((lmax.lt.0).or.(lmax.gt.50)) then
  write(*,*)
  write(*,'("Error(zbessela): lmax out of range : ",I8)') lmax
  write(*,*)
  stop
end if
if ((x.lt.0.d0).or.(x.gt.1.d8)) then
  write(*,*)
  write(*,'("Error(zbessela): x out of range : ",G18.10)') x
  write(*,*)
  stop
end if
xi=1.d0/x
xmin=1.d-7
select case (lmax)
 case(0)
   xmin=1.d-6
 case(1:2)
   xmin=1.d-4
 case(3:4)
   xmin=1.d-2
 case(5:)
   xmin=1.d0
end select
! treat x << 1
if (x.lt.xmin) then
  a(0)=1.d0
  t1=1.d0
  t2=1.d0
  do l=1,lmax
    t1=t1/dble(2*l+1)
    t2=-t2*x
    a(l)=t2*t1
  end do
  return
end if
! recurse up
a(0)=xi*(sinh(x))
if (lmax.eq.0) return
a(1)=xi*(a(0)-cosh(x))
if (lmax.eq.1) return
a0=a(0)
a1=a(1)
do l=2,lmax
  at=dble(2*l-1)*a1*xi+a0
  a0=a1
  a1=at
  a(l)=a1
end do
return
end subroutine
!EOC

