
! Copyright (C) 2013 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine zfmtctof(zfmt)
use modmain
implicit none
! arguments
real(8), intent(inout) :: zfmt(2,lmmaxvr,nrmtmax,natmtot)
! local variables
integer ld1,ld2,is,ias,ir,lm
ld2=2*lmmaxvr
ld1=ld2*lradstp
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(is,ir,lm)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  do ir=1,nrmtinr(is),lradstp
    zfmt(:,lmmaxinr+1:lmmaxvr,ir,ias)=0.d0
  end do
! interpolate with a clamped spline
  do lm=1,lmmaxvr
! real part
    call rfinterp(nrcmt(is),rcmt(:,is),ld1,zfmt(1,lm,1,ias),nrmt(is), &
     rsp(:,is),ld2,zfmt(1,lm,1,ias))
! imaginary part
    call rfinterp(nrcmt(is),rcmt(:,is),ld1,zfmt(2,lm,1,ias),nrmt(is), &
     rsp(:,is),ld2,zfmt(2,lm,1,ias))
  end do
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine

