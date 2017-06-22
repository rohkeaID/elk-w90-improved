
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: wigner3jf
! !INTERFACE:
real(8) function wigner3jf(j12,j22,j32,m12,m22,m32)
! !INPUT/OUTPUT PARAMETERS:
!   j12, j22, j32 : angular momentum quantum numbers times 2 (in,integer)
!   m12, m22, m32 : magnetic quantum numbers times 2 (in,integer)
! !DESCRIPTION:
!   Returns the Wigner $3j$-symbol for the case where the arguments may be
!   fractional, i.e. multiples of $\frac{1}{2}$. The input parameters to this
!   function are taken to be twice their actual values, which allows them to
!   remain integers. The formula used is identical to that in {\tt wigner3j}.
!
! !REVISION HISTORY:
!   Created January 2014 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: j12,j22,j32
integer, intent(in) :: m12,m22,m32
! local variables
integer jm1,jm2,jm3,n1,n2
integer l12,l22,l32,l42
integer k,k1,k2,l1,l2,l3
real(8) sgn,sum,t1
! external functions
real(8) factnm,factr
external factnm,factr
! check input variables
if ((j12.lt.0).or.(j22.lt.0).or.(j32.lt.0).or.(abs(m12).gt.j12).or. &
 (abs(m22).gt.j22).or.(abs(m32).gt.j32)) then
  write(*,*)
  write(*,'("Error(wigner3jf): invalid arguments :")')
  write(*,'("j12 = ",I8," j22 = ",I8," j32 = ",I8)') j12,j22,j32
  write(*,'("m12 = ",I8," m22 = ",I8," m32 = ",I8)') m12,m22,m32
  write(*,*)
  stop
end if
if ((j12.eq.0).and.(j22.eq.0).and.(j32.eq.0)) then
  wigner3jf=1.d0
  return
end if
if ((j12.gt.100).or.(j22.gt.100).or.(j32.gt.100)) then
  write(*,*)
  write(*,'("Error(wigner3jf): angular momenta out of range : ",3I8)') j12, &
   j22,j32
  write(*,*)
  stop
end if
jm1=j12+m12
jm2=j22+m22
jm3=j32+m32
if ((mod(jm1,2).ne.0).or.(mod(jm2,2).ne.0).or.(mod(jm3,2).ne.0)) then
  wigner3jf=0.d0
  return
end if
l12=j22-j12+j32
l22=j12-j22+j32
l32=j12+j22-j32
l42=j12+j22+j32
if ((mod(l12,2).ne.0).or.(mod(l22,2).ne.0).or.(mod(l32,2).ne.0).or. &
 (mod(l42,2).ne.0)) then
  wigner3jf=0.d0
  return
end if
l1=l12/2
l2=l22/2
l3=l32/2
if ((m12+m22+m32.ne.0).or.(l1.lt.0).or.(l2.lt.0).or.(l3.lt.0)) then
  wigner3jf=0.d0
  return
end if
n1=(j12-m12)/2
n2=(j22+m22)/2
k1=max(0,n1-l2,n2-l1)
k2=min(l3,n1,n2)
k=k1+(j22-j12+m32)/2
if (mod(k,2).ne.0) then
  sgn=-1.d0
else
  sgn=1.d0
end if
sum=0.d0
do k=k1,k2
  t1=sgn*factr(l1,l1-n2+k)*factr(l2,l2-n1+k)*factr(l3,l3-k)
  sum=sum+t1/(factnm(k,1)*factnm(n1-k,1)*factnm(n2-k,1))
  sgn=-sgn
end do
jm1=jm1/2
jm2=jm2/2
jm3=jm3/2
t1=factr(jm1,l1)*factr(jm2,l2)*factr(jm3,l3)
jm1=(j12-m12)/2
jm2=(j22-m22)/2
jm3=(j32-m32)/2
t1=t1*factr(jm3,1+l42/2)*factnm(jm1,1)*factnm(jm2,1)
wigner3jf=sum*sqrt(t1)
return
end function
!EOC

