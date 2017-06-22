
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: massnucl
! !INTERFACE:
real(8) function massnucl(z)
! !INPUT/OUTPUT PARAMETERS:
!   z : atomic number (in,real)
! !DESCRIPTION:
!   Computes an approximate nuclear mass from the atomic number $Z$. The nuclear
!   mass number, $A$, is first estimated using
!   $$ A=4.467\times 10^{-3}Z^2+2.163 Z-1.168, $$
!   [D. Andrae in {\it Relativistic Electronic Structure Theory - Fundamentals}
!   {\bf 11}, 203 (2002)]. Then the nuclear mass can be determined from:
!   $$ M=Z m_p+N m_n-\frac{B}{c^2}, $$
!   where $m_p$ is the proton mass, $m_n$ is the neutron mass and $B$ is the
!   nuclear binding energy. The binding energy is approximated by the
!   Weizs\"{a}cker formula:
!   $$ B=a_V A-a_S A^{2/3}-a_C Z^2 A^{-1/3}-a_{\rm sym}(Z-N)^2A^{-1}
!    +B_p+B_{\rm shell}. $$
!   See F. Yang and J. H. Hamilton in {\it Modern Atomic and Nuclear Physics},
!   Revised Edition 2010, for details on the quantities in this formula. In this
!   implementation, $B_p$ and $B_{\rm shell}$ are set to zero.
!
! !REVISION HISTORY:
!   Created February 2014 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: z
! local variables
! coefficients for computing mass number
real(8), parameter :: c2=4.467d-3, c1=2.163d0, c0=-1.168d0
! Weizsacker coefficients in MeV
real(8), parameter :: av=15.8d0, as=18.3d0, ac=0.72d0, asym=23.2d0
! proton and neutron masses in MeV/c^2 (CODATA 2010)
real(8), parameter :: mp=938.272046d0
real(8), parameter :: mn=939.565379d0
! atomic mass unit in MeV/c^2 (CODATA 2010)
real(8), parameter :: mu=931.494061d0
real(8) za,n,a,b
za=abs(z)
! approximate nuclear mass number
if (za.le.1.d0) then
  a=1.d0
else
  a=abs(c2*za**2+c1*za+c0)
end if
n=a-za
b=av*a-as*a**(2.d0/3.d0)-ac*(za**2)/a**(1.d0/3.d0)-asym*(za-n)**2/a
massnucl=(za*mp+n*mn-b)/mu
return
end function
!EOC

