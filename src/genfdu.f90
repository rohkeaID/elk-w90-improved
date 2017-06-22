
! Copyright (C) 2008 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genfdu
! !INTERFACE:
subroutine genfdu(i,u,j,f)
! !USES:
use moddftu
! !INPUT/OUTPUT PARAMETERS:
!   i : DFT+U entry (in,integer)
!   u : parameter U (inout,real)
!   j : parameter J (inout,real)
!   f : Slater parameters (inout,real)
! !DESCRIPTION:
!   Calculate the Slater parameters for DFT+$U$ calculation with different
!   approaches, see  {\it Phys. Rev. B} {\bf 80}, 035121 (2009). The relations
!   among Slater and Racah parameters are from E.U. Condon and G.H. Shortley,
!   {\it The Theory of Atomic Spectra},  The University Press, Cambridge (1935).
!
! !REVISION HISTORY:
!   Created July 2008 (Francesco Cricchio)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: i
real(8), intent(inout) :: u,j
real(8), intent(inout) :: f(0:2*lmaxdm)
! local variables
integer is,l,k,q
real(8) r1,r2
real(8) lambda,ufix
! automatic arrays
real(8) :: e(0:lmaxdm)
real(8) :: a(3,3),v1(3),v2(3)
! external functions
real(8) fyukawa,fyukawa0
external fyukawa,fyukawa0
is=idftu(1,i)
l=idftu(2,i)
! load input parameters to calculate Slater integrals
u=ujdu(1,i)
j=ujdu(2,i)
f(:)=fdu(:,i)
e(:)=edu(:,i)
lambda=lambdadu(i)
if (inpdftu.lt.4) then
! F(0) = U for any l-shell
  if (inpdftu.eq.1) f(0)=u
  select case(l)
  case(0)
! s electrons only f(0)=u
    if (inpdftu.eq.3) then
      f(0)=e(0)
      u=f(0)
    end if
  case(1)
! p electrons
    if (inpdftu.eq.1) then
! F(2) = 5.0 * J
      f(2)=5.d0*j
    else if (inpdftu.eq.3) then
! F(0) = E(0) + J= E(0) + 5/3 * E(1)
      f(0)=e(0)+(5.d0/3.d0)*e(1)
! F(2) = 5 * J = 25/3 * E1, Eq. 101
      f(2)=(25.d0/3.d0)*e(1)
    end if
  case(2)
! d electrons
    if (inpdftu.eq.1) then
! r1 = F(4)/F(2), see PRB 52, R5467 (1995)
      r1=0.625d0
      f(2)=(14.d0*j)/(1.d0+r1)
      f(4)=f(2)*r1
    else if (inpdftu.eq.3) then
! copy Racah parameters
      v1(1:3)=e(0:2)
! transformation matrix from Racah to Slater parameters
! obtained from inversion of Eq. 110-112, LN Notes 29-12-08
      a(1,1)=1.d0
      a(1,2)=1.4d0
      a(1,3)=0.d0
      a(2,1)=0.d0
      a(2,2)=0.1428571428571428d0
      a(2,3)=1.285714285714286d0
      a(3,1)=0.d0
      a(3,2)=2.8571428571428571d-2
      a(3,3)=-0.1428571428571428d0
! multiply transformation matrix by Racah parameters
      call r3mv(a,v1,v2)
! Slater parameters, Eq. 104-105, LN Notes 29-12-08
      f(0)=v2(1)
      f(2)=49.d0*v2(2)
      f(4)=441.d0*v2(3)
    end if
  case(3)
! f electrons
    if (inpdftu.eq.1) then
! r2 = F(6)/F(2), r1 = F(4)/F(2), see PRB 50, 16861 (1994)
      r1=451.d0/675.d0
      r2=1001.d0/2025.d0
      f(2)=6435.d0*j/(286.d0+195.d0*r1+250.d0*r2)
      f(4)=f(2)*r1
      f(6)=f(2)*r2
    else if (inpdftu.eq.3) then
! F(0) = E(0) + 9/7 * E(1) , Eq. 119, LN Notes 29-12-08
      f(0)=e(0)+(9.d0/7.d0)*e(1)
! copy Racah parameters
      v1(1:3)=e(1:3)
! transformation matrix from Racah to Slater parameters
! obtained from inversion of Eq. 120-122, LN Notes 29-12-08
      a(1,1)=2.3809523809523808d-2
      a(1,2)=3.404761904761904d0
      a(1,3)=0.2619047619047619d0
      a(2,1)=1.2987012987012984d-2
      a(2,2)=-1.688311688311688d0
      a(2,3)=5.1948051948051951d-2
      a(3,1)=2.1645021645021645d-3
      a(3,2)=7.5757575757575760d-2
      a(3,3)=-1.5151515151515152d-2
! multiply transformation matrix by Racah parameters
      call r3mv(a,v1,v2)
! Slater parameters, Eq. 115-117, LN Notes 29-12-08
      f(2)=225.d0*v2(1)
      f(4)=1089.d0*v2(2)
      f(6)=(184041.d0/25.d0)*v2(3)
    end if
  case default
    write(*,*)
    write(*,'("Error(genfdu): invalid l : ",I8)') l
    write(*,*)
    stop
  end select
else if (inpdftu.ge.4) then
! define energies for Slater parameters
  call energyfdu
! write energies for Slater parameters to a file
  call writeefdu
! calculate radial functions for Slater parameters
  call genfdufr
  if (inpdftu.eq.5) then
    ufix=udufix(i)
! calculate the lambda corresponding to udufix
! lambdadu0 is in/out and is initialized to 0 in readinput
    call findlambdadu(is,l,ufix,lambdadu0(i),lambda)
  end if
  do q=0,l
    k=2*q
    if (lambda.lt.1.d-2) then
! unscreened Slater parameters
      f(k)=fyukawa0(is,l,k)
    else
! screened Slater parameters
      f(k)=fyukawa(is,l,k,lambda)
    end if
  end do
end if
! calculate U and J from Slater integrals
if (inpdftu.ne.1) then
  u=f(0)
  select case(l)
  case(0)
    j=0.d0
  case(1)
! J = 1/5 * F(2)
    j=(1.d0/5.d0)*f(2)
  case(2)
! J = 1/14 * ( F(2) + F(4) ), Eq. 106, LN Notes 29-12-08
    j=(1.d0/14.d0)*(f(2)+f(4))
  case(3)
! J= 2/45 * F(2) + 1/33 * F(4) + 50/1287 * F(6), Eq. 118, LN Notes 29-12-08
    j=(2.d0/45.d0)*f(2)+(1.d0/33.d0)*f(4)+(50.d0/1287.d0)*f(6)
  case default
    write(*,*)
    write(*,'("Error(genfdu): invalid l : ",I8)') l
    write(*,*)
    stop
  end select
end if
! save calculated parameters
! (except Racah parameters that are provided only as input)
ujdu(1,i)=u
ujdu(2,i)=j
fdu(:,i)=f(:)
lambdadu(i)=lambda
return
end subroutine
!EOC

