
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modrandom

! random number generator state
integer(8) rndstate(0:5)
data rndstate / 799047353, 1322018920, 1014372120, 1198189977, 832907020, &
                5678910 /

contains

!BOP
! !ROUTINE: randomu
! !INTERFACE:
real(8) function randomu()
! !DESCRIPTION:
!   Generates random numbers with a uniform distribution in the interval $[0,1]$
!   using the fifth-order multiple recursive generator of P. L'Ecuyer,
!   F. Blouin, and R. Coutre, {\it ACM Trans. Modeling Comput. Simulation}
!   {\bf 3}, 87 (1993). The sequence of numbers $r_i$ is produced from
!   $$ x_i=(a_1 x_{i-1}+a_5 x_{i-5})\mod m $$
!   with $r_i=x_i/m$. The period is about $2^{155}$.
!
! !REVISION HISTORY:
!   Created January 2012 (JKD)
!   Changed initial state, April 2017 (JKD)
!EOP
!BOC
implicit none
! local variables
! parameters taken from the GNU Scientific Library (GSL)
integer(8), parameter :: a1=107374182, a5=104480, m=2147483647
integer(8) i,i1,i5
data i / 0 /
i=modulo(i+1,6_8)
i1=modulo(i-1,6_8)
i5=modulo(i-5,6_8)
rndstate(i)=int(mod(a1*rndstate(i1)+a5*rndstate(i5),m))
randomu=dble(rndstate(i))/dble(m)

end function
!EOC

end module

