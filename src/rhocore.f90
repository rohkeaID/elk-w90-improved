
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
integer ispn,idm,is,ias
integer nr,nri,iro,ir,i
real(8) v(ndmag),sum,t1
! automatic arrays
real(8) fr(nrmtmax)
! external functions
real(8) fintgt
external fintgt
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  iro=nri+1
  sum=0.d0
! loop over spin channels
  do ispn=1,nspncr
! add the core density to the muffin-tin density
    i=1
    do ir=1,nri
      rhomt(i,ias)=rhomt(i,ias)+rhocr(ir,ias,ispn)/y00
      fr(ir)=rhocr(ir,ias,ispn)*r2sp(ir,is)
      i=i+lmmaxi
    end do
    do ir=iro,nr
      rhomt(i,ias)=rhomt(i,ias)+rhocr(ir,ias,ispn)/y00
      fr(ir)=rhocr(ir,ias,ispn)*r2sp(ir,is)
      i=i+lmmaxo
    end do
! compute the core charge inside the muffin-tins
    t1=fintgt(-1,nr,rsp(:,is),fr)
    sum=sum+fourpi*t1
  end do
! core leakage charge
  chgcrlk(ias)=chgcr(is)-sum
! add to the magnetisation in the case of a spin-polarised core
  if (spincore) then
! compute the total moment in the muffin-tin
    do idm=1,ndmag
      i=1
      do ir=1,nri
        fr(ir)=magmt(i,ias,idm)*r2sp(ir,is)
        i=i+lmmaxi
      end do
      do ir=iro,nr
        fr(ir)=magmt(i,ias,idm)*r2sp(ir,is)
        i=i+lmmaxo
      end do
      t1=fintgt(-1,nr,rsp(:,is),fr)
      v(idm)=fourpi*y00*t1
    end do
! normalise
    if (ncmag) then
      t1=sqrt(v(1)**2+v(2)**2+v(3)**2)
    else
      t1=abs(v(1))
    end if
    if (t1.gt.1.d-10) v(:)=v(:)/t1
! add the core magnetisation to the total
    i=1
    do ir=1,nri
      t1=abs((rhocr(ir,ias,1)-rhocr(ir,ias,2))/y00)
      magmt(i,ias,:)=magmt(i,ias,:)+t1*v(:)
      i=i+lmmaxi
    end do
    do ir=iro,nr
      t1=abs((rhocr(ir,ias,1)-rhocr(ir,ias,2))/y00)
      magmt(i,ias,:)=magmt(i,ias,:)+t1*v(:)
      i=i+lmmaxo
    end do
  end if
end do
return
end subroutine
!EOC

