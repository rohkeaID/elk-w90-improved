
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rhoinit
! !INTERFACE:
subroutine rhoinit
! !USES:
use modmain
use modomp
! !DESCRIPTION:
!   Initialises the crystal charge density. Inside the muffin-tins it is set to
!   the spherical atomic density. In the interstitial region it is taken to be
!   constant such that the total charge is correct. Requires that the atomic
!   densities have already been calculated.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ias,i
integer nr,nri,ir
real(8) t1
! set the muffin-tin density equal to the atomic plus excess density
t1=chgexs/omega
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  rhomt(:,ias)=0.d0
  i=1
  do ir=1,nri
    rhomt(i,ias)=(t1+rhosp(ir,is))/y00
    i=i+lmmaxi
  end do
  do ir=nri+1,nr
    rhomt(i,ias)=(t1+rhosp(ir,is))/y00
    i=i+lmmaxo
  end do
end do
! compute the muffin-tin charges
call chargemt
! add a constant to the density so the total charge is correct
t1=(chgtot-chgmttot)/omega
rhoir(:)=t1
t1=t1/y00
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  i=1
  do ir=1,nri
    rhomt(i,ias)=rhomt(i,ias)+t1
    i=i+lmmaxi
  end do
  do ir=nri+1,nr
    rhomt(i,ias)=rhomt(i,ias)+t1
    i=i+lmmaxo
  end do
end do
! compute the charges again
call charge
return
end subroutine
!EOC

