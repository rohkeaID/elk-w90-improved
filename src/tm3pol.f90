
! Copyright (C) 2008 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tm3pol(l,k,p,r,w2,tm3p)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   l    : angular momentum (in,integer)
!   k    : k-index of tensor moment (in,integer)
!   p    : p-index of tensor moment (in,integer)
!   r    : r-index of tensor moment (in,integer)
!   w2   : modulus square of k-p-r tensor moment (in,real)
!   tm3p : polarisation (out,real)
! !DESCRIPTION:
!   Calculate the polarisation of each 3-index tensor component of the density
!   matrix, see {\it Phys. Rev. B} {\bf 80}, 035121 (2009).
!
! !REVISION HISTORY:
!   Created April 2008 (F. Cricchio and L. Nordstrom)
!EOP
!BOC
implicit none
! input variables
integer, intent(in) :: l,k,p,r
real(8), intent(in) :: w2
real(8), intent(out) :: tm3p
! local variables
integer g
real(8) nlk,t1
! external functions
real(8) wigner3j,factnm,factr
external wigner3j,factnm,factr
g=k+p+r
if (g.eq.0) then
  t1=sqrt(w2)
  tm3p=t1*(dble(2*(2*l+1))-t1)
else
  if (mod(g,2).eq.0) then
    t1=wigner3j(k,p,r,0,0,0)
  else
    t1=sqrt(factnm(g-2*p,1)*factnm(g-2*r,1)*factr(g-2*k,g+1))
    t1=t1*factnm(g,2)/(factnm(g-2*k,2)*factnm(g-2*p,2)*factnm(g-2*r,2))
  end if
  nlk=factnm(2*l,1)/sqrt(factnm(2*l-k,1)*factnm(2*l+k+1,1))
  tm3p=dble((2*k+1)*(2*r+1)*(2*l+1))*((t1*nlk)**2)*w2
end if
return
end subroutine
!EOC

