
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modrandom

! random number generator state
integer rndstate(0:5)
data rndstate / 1, 12, 234, 3456, 45678, 5678910 /

contains

!BOP
! !ROUTINE: randomu
! !INTERFACE:
real(8) function randomu()
! !DESCRIPTION:
!   Generates random numbers with a uniform distribution using the fifth-order
!   multiple recursive generator of P. L'Ecuyer, F. Blouin, and R. Coutre,
!   {\it ACM Trans. Modeling Comput. Simulation} {\bf 3}, 87 (1993). The
!   sequence of numbers $r_i$ is produced from
!   $$ x_i=(a_1 x_{i-1}+a_5 x_{i-5})\mod m $$
!   with $r_i=x_i/m$. The period is about $2^{155}$.
!
! !REVISION HISTORY:
!   Created January 2012 (JKD)
!EOP
!BOC
implicit none
! local variables
! parameters taken from the GNU Scientific Library (GSL)
integer(8), parameter :: a1=107374182, a5=104480, m=2147483647
integer i,i1,i5
data i / 0 /
i=modulo(i+1,6)
i1=modulo(i-1,6)
i5=modulo(i-5,6)
rndstate(i)=int(mod(a1*rndstate(i1)+a5*rndstate(i5),m))
randomu=dble(rndstate(i))/dble(m)
end function
!EOC

end module

