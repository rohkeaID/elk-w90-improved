
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: moment
! !INTERFACE:
subroutine moment
! !USES:
use modmain
use modtest
! !DESCRIPTION:
!   Computes the muffin-tin, interstitial and total moments by integrating the
!   magnetisation.
!
! !REVISION HISTORY:
!   Created January 2005 (JKD)
!EOP
!BOC
implicit none
! local variables
integer idm,is,ias,ir
real(8) t1
! automatic arrays
real(8) fr(nrmtmax)
! external functions
real(8) fintgt,ddot
external fintgt,ddot
if (.not.spinpol) then
  mommt(:,:)=0.d0
  mommttot(:)=0.d0
  momir(:)=0.d0
  momtot(:)=0.d0
  return
end if
! find the muffin-tin moments
mommttot(:)=0.d0
do idm=1,ndmag
  do ias=1,natmtot
    is=idxis(ias)
    do ir=1,nrmt(is)
      fr(ir)=magmt(1,ir,ias,idm)*r2sp(ir,is)
    end do
    t1=fintgt(-1,nrmt(is),rsp(:,is),fr)
    mommt(idm,ias)=fourpi*y00*t1
    mommttot(idm)=mommttot(idm)+mommt(idm,ias)
  end do
end do
! find the interstitial moments
do idm=1,ndmag
  t1=ddot(ngtot,magir(:,idm),1,cfunir,1)
  momir(idm)=t1*omega/dble(ngtot)
end do
momtot(:)=mommttot(:)+momir(:)
! write total moment to test file
call writetest(450,'total moment',nv=ndmag,tol=2.d-2,rva=momtot)
! total moment magnitude
if (ncmag) then
  momtotm=sqrt(momtot(1)**2+momtot(2)**2+momtot(3)**2)
else
  momtotm=abs(momtot(1))
end if
return
end subroutine
!EOC

