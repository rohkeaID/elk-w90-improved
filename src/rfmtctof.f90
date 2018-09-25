
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rfmtctof
! !INTERFACE:
subroutine rfmtctof(rfmt)
! !USES:
use modmain
use modomp
! !INPUT/OUTPUT PARAMETERS:
!   rfmt : real muffin-tin function (in,real(npmtmax,natmtot))
! !DESCRIPTION:
!   Converts a real muffin-tin function from a coarse to a fine radial mesh by
!   using cubic spline interpolation. See {\tt rfinterp} and {\tt spline}.
!
! !REVISION HISTORY:
!   Created October 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(inout) :: rfmt(npmtmax,natmtot)
! local variables
integer is,ias,lm,nthd
integer nr,nri,nro
integer iro,ir,npi,i
integer nrc,nrci,nrco
integer irco,irc,npci
! allocatable arrays
real(8), allocatable :: fi(:),fo(:),rfmt1(:)
if (lradstp.eq.1) return
call omp_hold(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(fi,fo,rfmt1,is,nr,nri,nro) &
!$OMP PRIVATE(iro,npi,nrc,nrci,nrco) &
!$OMP PRIVATE(irco,npci,lm,i,irc,ir) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ias=1,natmtot
  allocate(fi(nrcmtmax),fo(nrmtmax),rfmt1(npmtmax))
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  nro=nr-nri
  iro=nri+1
  npi=npmti(is)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  nrco=nrc-nrci
  irco=nrci+1
  npci=npcmti(is)
! interpolate up to lmaxi over entire muffin-tin
  do lm=1,lmmaxi
    i=lm
    do irc=1,nrci
      fi(irc)=rfmt(i,ias)
      i=i+lmmaxi
    end do
    do irc=irco,nrc
      fi(irc)=rfmt(i,ias)
      i=i+lmmaxo
    end do
    call rfinterp(nrc,rcmt(:,is),fi,nr,rsp(:,is),fo)
    i=lm
    do ir=1,nri
      rfmt1(i)=fo(ir)
      i=i+lmmaxi
    end do
    do ir=iro,nr
      rfmt1(i)=fo(ir)
      i=i+lmmaxo
    end do
  end do
! interpolate up to lmaxo on outer part of muffin-tin
  do lm=lmmaxi+1,lmmaxo
    i=npci+lm
    do irc=irco,nrc
      fi(irc)=rfmt(i,ias)
      i=i+lmmaxo
    end do
    call rfinterp(nrco,rcmt(irco,is),fi(irco),nro,rsp(iro,is),fo(iro))
    i=npi+lm
    do ir=iro,nr
      rfmt1(i)=fo(ir)
      i=i+lmmaxo
    end do
  end do
  call dcopy(npmt(is),rfmt1,1,rfmt(:,ias),1)
  deallocate(fi,fo,rfmt1)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
return
end subroutine
!EOC

