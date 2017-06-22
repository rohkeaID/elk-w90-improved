
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rfinp
! !INTERFACE:
real(8) function rfinp(lrstp,rfmt1,rfir1,rfmt2,rfir2)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   lrstp : radial step length (in,integer)
!   rfmt1 : first function in real spherical harmonics for all muffin-tins
!           (in,real(lmmaxvr,nrmtmax,natmtot))
!   rfir1 : first real interstitial function in real-space (in,real(ngtot))
!   rfmt2 : second function in real spherical harmonics for all muffin-tins
!           (in,real(lmmaxvr,nrmtmax,natmtot))
!   rfir2 : second real interstitial function in real-space (in,real(ngtot))
! !DESCRIPTION:
!   Calculates the inner product of two real functions over the entire unit cell.
!   The input muffin-tin functions should have angular momentum cut-off
!   {\tt lmaxvr}. In the interstitial region, the integrand is multiplied with
!   the characteristic function, $\tilde{\Theta}({\bf r})$, to remove the
!   contribution from the muffin-tin. See routines {\tt rfmtinp} and
!   {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created July 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lrstp
real(8), intent(in) :: rfmt1(lmmaxvr,nrmtmax,natmtot),rfir1(ngtot)
real(8), intent(in) :: rfmt2(lmmaxvr,nrmtmax,natmtot),rfir2(ngtot)
! local variables
integer is,ias,ir
real(8) sum
! external functions
real(8) rfmtinp
external rfmtinp
sum=0.d0
! interstitial contribution
do ir=1,ngtot
  sum=sum+rfir1(ir)*rfir2(ir)*cfunir(ir)
end do
sum=sum*omega/dble(ngtot)
! muffin-tin contribution
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(is) REDUCTION(+:sum)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  sum=sum+rfmtinp(nrmt(is),nrmtinr(is),lrstp,rsp(:,is),r2sp(:,is), &
   rfmt1(:,:,ias),rfmt2(:,:,ias))
end do
!$OMP END DO
!$OMP END PARALLEL
rfinp=sum
return
end function
!EOC

