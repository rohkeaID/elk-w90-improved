
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: radnucl
! !INTERFACE:
real(8) function radnucl(z)
! !INPUT/OUTPUT PARAMETERS:
!   z : atomic number (in,real)
! !DESCRIPTION:
!   Computes an approximate nuclear charge radius from the atomic number $Z$.
!   The nuclear mass number, $A$, is estimated using
!   $$ A=4.467\times 10^{-3}Z^2+2.163 Z-1.168, $$
!   [D. Andrae in {\it Relativistic Electronic Structure Theory - Fundamentals}
!   {\bf 11}, 203 (2002)], and the nuclear charge radius can be determined from
!   $$ r=\left(r_0+\frac{r_1}{A^{2/3}}+\frac{r_2}{A^{4/3}}\right)A^{1/3}, $$
!   where $r_0=0.9071$, $r_1=1.105$ and $r_2=-0.548$ [I. Angeli, {\it Atomic
!   Data and Nuclear Data Tables} {\bf 87}, 185 (2004)].
!
! !REVISION HISTORY:
!   Created October 2011 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: z
! local variables
! coefficients for computing mass number
real(8), parameter :: c2=4.467d-3, c1=2.163d0, c0=-1.168d0
! coefficients for computing charge radius (fm)
real(8), parameter :: r0=0.9071d0, r1=1.105d0, r2=-0.548d0
! Bohr radius in SI units (CODATA 2010)
real(8), parameter :: br_si=0.52917721092d-10
real(8) za,a,a13,a23,a43
za=abs(z)
! approximate nuclear mass number
if (za.le.1.d0) then
  a=1.d0
else
  a=c2*za**2+c1*za+c0
end if
! approximate nuclear charge radius
a13=a**(1.d0/3.d0)
a23=a13**2
a43=a13*a
radnucl=(r0+r1/a23+r2/a43)*a13
radnucl=radnucl*1.d-15/br_si
return
end function
!EOC

