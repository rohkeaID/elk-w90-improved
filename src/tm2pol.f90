
! Copyright (C) 2008 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tm2pol(l,k,w2,tm2p)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   l    : angular momentum (in,integer)
!   k    : k-index of tensor moment (in,integer)
!   w2   : modulus square of k-p tensmom (in,real)
!   tm2p : polarisation (out,real)
! !DESCRIPTION:
!   Calculate the polarisation of each 2-index tensor component of the density
!   matrix, see {\it Phys. Rev. B} {\bf 80}, 035121 (2009).
!
! !REVISION HISTORY:
!   Created April 2008 (F. Cricchio and L. Nordstrom)
!EOP
!BOC
implicit none
! input variables
integer, intent(in) :: l,k
real(8),intent(in) :: w2
real(8),intent(out) :: tm2p
! local variables
real(8) nlk
! external functions
real(8) factnm
external factnm
nlk=factnm(2*l,1)/sqrt(factnm(2*l-k,1)*factnm(2*l+k+1,1))
tm2p=dble((2*k+1)*(2*l+1))*(nlk**2)*w2
return
end subroutine
!EOC

