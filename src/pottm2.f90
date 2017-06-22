
! Copyright (C) 2008 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: pottm2
! !INTERFACE:
subroutine pottm2(i,k1,p,vh,vx)
! !USES:
use moddftu
! !INPUT/OUTPUT PARAMETERS:
!   i  : DFT+U entry (in,integer)
!   k1 : k-index of tensor moment (in,integer)
!   p  : p-index of tensor moment  (in,integer)
!   vh : Hartree potential energy (out,real)
!   vx : exchange potential energy (out,real)
! !DESCRIPTION:
!   Calculates the DFT+$U$ Hartree and exchange potential energies for a 2-index
!   tensor moment component. See {\tt pottm3}.
!
! !REVISION HISTORY:
!   Created April 2008 (F. Cricchio and L. Nordstrom)
!   Modified, January 2014 (JKD)
!EOP
!BOC
implicit none
integer, intent(in) :: i
integer, intent(in) :: k1,p
real(8), intent(out) :: vh,vx
! local variables
integer l,k
real(8) nlk,t1,t2
! external functions
real(8) wigner3j,wigner6j,factnm
external wigner3j,wigner6j,factnm
l=idftu(2,i)
nlk=factnm(2*l,1)/sqrt(factnm(2*l-k1,1)*factnm(2*l+k1+1,1))
vh=0.d0
vx=0.d0
do k=0,2*l,2
  t1=0.5d0*(dble(2*l+1)*nlk*wigner3j(l,k,l,0,0,0))**2
  t2=0.5d0*dble((2*k1+1)*(-1)**k1)*wigner6j(l,l,k1,l,l,k)
  if (k.eq.k1) then
    if (p.eq.0) vh=t1*fdu(k1,i)
  end if
  vx=vx-t1*t2*fdu(k,i)
end do
return
end subroutine
!EOC

