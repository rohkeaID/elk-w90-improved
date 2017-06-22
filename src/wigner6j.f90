
! Copyright (C) 2009 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: wigner6j
! !INTERFACE:
real(8) function wigner6j(j1,j2,j3,k1,k2,k3)
! !INPUT/OUTPUT PARAMETERS:
!   j1, j2, j3 : angular momentum quantum numbers (in,integer)
!   k1, k2, k3 : angular momentum quantum numbers (in,integer)
! !DESCRIPTION:
!   Returns the Wigner $6j$-symbol for integer arguments. This is computed using
!   the Racah formula:
!   \begin{align*}
!    &\left\{\begin{matrix} j_1 & j_2 & j_3 \\ k_1 & k_2 & k_3 \end{matrix}
!    \right\}=\sqrt{\Delta(j_1\,j_2\,j_3)\Delta(j_1\,k_2\,k_3)
!    \Delta(k_1\,j_2\,k_3)\Delta(k_1\,k_2\,j_3)}\,
!    \sum_n\frac{(-1)^n(n+1)!}{f(n)},
!   \end{align*}
!   where
!   \begin{align*}
!    f(n)=&(n-j_1-j_2-j_3)!\,(n-j_1-k_2-k_3)!\,(n-k_1-j_2-k_3)!\,
!    (n-k_1-k_2-j_3)! \\
!    &\times(j_1+j_2+k_1+k_2-n)!\,(j_2+j_3+k_2+k_3-n)!\,(j_1+j_3+k_1+k_3-n)!
!   \end{align*}
!   and
!   $$ \Delta(a,b,c)=\frac{(a+b-c)!\,(a-b+c)!\,(-a+b+c)!}{(a+b+c+1)!} $$
!   is the triangle coefficient, and where the sum is over all integers $n$ for
!   which the factorials in $f(n)$ have non-negative arguments. The Wigner-$6j$
!   function is zero unless each triad, $(j_1\,j_2\,j_3)$, $(j_1\,k_2\,k_3)$,
!   $(k_1\,j_2\,k_3)$ and $(k_1\,k_2\,j_3)$, satifies the triangle inequalities:
!   $$ |x-y|\le z \le x+y, $$
!   for triad $(x,y,z)$. See, for example, A. Messiah, {\it Quantum Mechanics},
!   Vol. 2., 1061-1066 (1962).
!
! !REVISION HISTORY:
!   Created August 2009 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: j1,j2,j3
integer, intent(in) :: k1,k2,k3
! local variables
integer n0,n1,n
real(8) sum,t1,t2
! external functions
real(8) triangle,factnm,factr
external triangle,factnm,factr
wigner6j=0.d0
if ((abs(j1-j2).gt.j3).or.((j1+j2).lt.j3)) return
if ((abs(j1-k2).gt.k3).or.((j1+k2).lt.k3)) return
if ((abs(k1-j2).gt.k3).or.((k1+j2).lt.k3)) return
if ((abs(k1-k2).gt.j3).or.((k1+k2).lt.j3)) return
if ((abs(j1).gt.50).or.(abs(j2).gt.50).or.(abs(j3).gt.50).or. &
    (abs(k1).gt.50).or.(abs(k2).gt.50).or.(abs(k3).gt.50)) then
  write(*,*)
  write(*,'("Error(wigner6j): arguments out of range :")')
  write(*,'(" j1, j2, j3 = ",3I8)') j1,j2,j3
  write(*,'(" k1, k2, k3 = ",3I8)') k1,k2,k3
  write(*,*)
  stop
end if
! range for summation
n0=max(j1+j2+j3,j1+k2+k3,k1+j2+k3,k1+k2+j3)
n1=min(j1+j2+k1+k2,j2+j3+k2+k3,j1+j3+k1+k3)
! Racah formula summation
sum=0.d0
do n=n0,n1
  t1=dble((-1)**n)*factr(n+1,n-(j1+j2+j3))
  t2=factnm(n-(j1+k2+k3),1)*factnm(n-(k1+j2+k3),1)*factnm(n-(k1+k2+j3),1)
  t2=t2*factnm(j1+j2+k1+k2-n,1)*factnm(j2+j3+k2+k3-n,1)*factnm(j1+j3+k1+k3-n,1)
  sum=sum+t1/t2
end do
! multiply by prefactor
wigner6j=sqrt(triangle(j1,j2,j3)*triangle(j1,k2,k3) &
             *triangle(k1,j2,k3)*triangle(k1,k2,j3))*sum
return
end function

real(8) function triangle(a,b,c)
implicit none
! arguments
integer, intent(in) :: a,b,c
! external functions
real(8) factnm,factr
external factnm,factr
triangle=factr(a+b-c,a+b+c+1)*factnm(a-b+c,1)*factnm(-a+b+c,1)
return
end function
!EOC

