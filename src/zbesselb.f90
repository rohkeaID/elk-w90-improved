
! Copyright (C) 2008 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zbesselb
! !INTERFACE:
subroutine zbesselb(lmax,x,b)
! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum order of Bessel function (in,integer)
!   x    : real argument (in,real)
!   b    : array of returned values (out,real(0:lmax))
! !DESCRIPTION:
!   Computes variations of the spherical Bessel function
!   $b_l(x)=i^lh^{(1)}_l(ix)$, for real argument $x$ and
!   $l=0,1,\ldots,l_{\rm max}$. The recursion relation
!   $$ j_{l+1}(x)=\frac{2l+1}{x}j_l(x)+j_{l-1}(x) $$
!   is used upwards. For starting values there are
!   $$ b_0(x)=-\frac{e^{-x}}{x};\qquad b_1(x)=b_0(x)\left\{1+\frac{1}{x}
!    \right\}. $$
!   For $x\ll 1$ the asymtotic forms
!   $$ b_l(x)\approx\frac{-(2l-1)!!}{(-x)^{l+1}} $$
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
real(8), intent(out) :: b(0:lmax)
! local variables
integer l
real(8) xi,b0,b1,bt,t3,t4
if ((lmax.lt.0).or.(lmax.gt.50)) then
  write(*,*)
  write(*,'("Error(zbesselb): lmax out of range : ",I8)') lmax
  write(*,*)
  stop
end if
if ((x.lt.0.d0).or.(x.gt.1.d8)) then
  write(*,*)
  write(*,'("Error(zbesselb): x out of range : ",G18.10)') x
  write(*,*)
  stop
end if
xi=1.d0/x
! treat x << 1
if (x.lt.1.d-7) then
  b(0)=-xi
  t3=-1.d0
  t4=xi
  do l=1,lmax
    t3=t3*dble(2*l-1)
    t4=t4*xi
    b(l)=t4*t3
  end do
  return
end if
! recurse up
b(0)=-xi*exp(-x)
if (lmax.eq.0) return
b(1)=b(0)*(1.d0+xi)
if (lmax.eq.1) return
b0=b(0)
b1=b(1)
do l=2,lmax
  bt=dble(2*l-1)*b1*xi+b0
  b0=b1
  b1=bt
  b(l)=b1
end do
return
end subroutine
!EOC

