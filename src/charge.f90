
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: charge
! !INTERFACE:
subroutine charge
! !USES:
use modmain
use modtest
! !DESCRIPTION:
!   Computes the muffin-tin, interstitial and total charges by integrating the
!   density.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
real(8) t1
! find the muffin-tin charges
call chargemt
! find the interstitial charge
t1=dot_product(rhoir(:),cfunir(:))
chgir=t1*omega/dble(ngtot)
! total calculated charge
chgcalc=chgmttot+chgir
! write total calculated charge to test file
call writetest(400,'calculated total charge',tol=1.d-6,rv=chgcalc)
return
end subroutine
!EOC

