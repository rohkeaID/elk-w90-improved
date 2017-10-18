
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
integer is,ias
integer nr,nri,ir,i
real(8) t1
! automatic arrays
real(8) fr(nrmtmax)
! external functions
real(8) fintgt
external fintgt
! find the muffin-tin charges
chgmttot=0.d0
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
! extract the l=0 component from the muffin-tin density
  i=1
  do ir=1,nri
    fr(ir)=rhomt(i,ias)*r2sp(ir,is)
    i=i+lmmaxi
  end do
  do ir=nri+1,nr
    fr(ir)=rhomt(i,ias)*r2sp(ir,is)
    i=i+lmmaxo
  end do
  t1=fintgt(-1,nr,rsp(:,is),fr)
  chgmt(ias)=fourpi*y00*t1
  chgmttot=chgmttot+chgmt(ias)
end do
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

