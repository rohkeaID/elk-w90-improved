
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rfmtctof
! !INTERFACE:
subroutine rfmtctof(rfmt)
! !INPUT/OUTPUT PARAMETERS:
!   rfmt : real muffin-tin function (in,real(lmmaxvr,nrmtmax,natmtot))
! !DESCRIPTION:
!   Converts a real muffin-tin function from a coarse to a fine radial mesh by
!   using cubic spline interpolation. See {\tt rfinterp} and {\tt spline}.
!
! !REVISION HISTORY:
!   Created October 2003 (JKD)
!EOP
!BOC
use modmain
implicit none
! arguments
real(8), intent(inout) :: rfmt(lmmaxvr,nrmtmax,natmtot)
! local variables
integer ld,is,ias,ir,lm
ld=lmmaxvr*lradstp
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(is,ir,lm)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  do ir=1,nrmtinr(is),lradstp
    rfmt(lmmaxinr+1:lmmaxvr,ir,ias)=0.d0
  end do
! interpolate with a clamped spline
  do lm=1,lmmaxvr
    call rfinterp(nrcmt(is),rcmt(:,is),ld,rfmt(lm,1,ias),nrmt(is),rsp(:,is), &
     lmmaxvr,rfmt(lm,1,ias))
  end do
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine
!EOC

