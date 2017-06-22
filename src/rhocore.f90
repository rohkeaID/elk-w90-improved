
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rhocore
! !INTERFACE:
subroutine rhocore
! !USES:
use modmain
! !DESCRIPTION:
!   Adds the core density and magnetisation to the muffin-tin functions. Also
!   computes the amount of leakage of core charge from the muffin-tin spheres
!   into the interstitial.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!   Fixed core moment direction, October 2012 (M. Meinert)
!EOP
!BOC
implicit none
! local variables
integer idm,ispn
integer is,ias,nr,ir
real(8) v(ndmag),sum,t1
! automatic arrays
real(8) fr(nrmtmax),gr(nrmtmax)
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  sum=0.d0
! loop over spin channels
  do ispn=1,nspncr
    do ir=1,nr
! add the core density to the muffin-tin density
      rhomt(1,ir,ias)=rhomt(1,ir,ias)+rhocr(ir,ias,ispn)/y00
      fr(ir)=rhocr(ir,ias,ispn)*r2sp(ir,is)
    end do
! compute the core charge inside the muffin-tins
    call fderiv(-1,nr,rsp(:,is),fr,gr)
    sum=sum+fourpi*gr(nr)
  end do
! core leakage charge
  chgcrlk(ias)=chgcr(is)-sum
! add to the magnetisation in the case of a spin-polarised core
  if (spincore) then
! compute the total moment in the muffin-tin
    do idm=1,ndmag
      do ir=1,nr
        fr(ir)=magmt(1,ir,ias,idm)*r2sp(ir,is)
      end do
      call fderiv(-1,nr,rsp(:,is),fr,gr)
      v(idm)=fourpi*y00*gr(nr)
    end do
! normalise
    if (ncmag) then
      t1=sqrt(v(1)**2+v(2)**2+v(3)**2)
    else
      t1=abs(v(1))
    end if
    if (t1.gt.1.d-10) v(:)=v(:)/t1
! add the core magnetisation to the total
    do ir=1,nr
      t1=abs((rhocr(ir,ias,1)-rhocr(ir,ias,2))/y00)
      magmt(1,ir,ias,:)=magmt(1,ir,ias,:)+t1*v(:)
    end do
  end if
end do
return
end subroutine
!EOC

