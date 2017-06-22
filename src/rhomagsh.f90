
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rhomagsh
! !INTERFACE:
subroutine rhomagsh
! !USES:
use modmain
! !DESCRIPTION:
!   Converts the muffin-tin density and magnetisation from spherical coordinates
!   to a spherical harmonic expansion. See {\tt rhomagk}.
!
! !REVISION HISTORY:
!   Created January 2009 (JKD)
!EOP
!BOC
implicit none
! local variables
integer idm,is,ias
integer nr,nri,nrc,nrci
! allocatable arrays
real(8), allocatable :: rfmt(:,:)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt,is,nr,nri) &
!$OMP PRIVATE(nrc,nrci,idm)
!$OMP DO
do ias=1,natmtot
  allocate(rfmt(lmmaxvr,nrcmtmax))
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmtinr(is)
  nrc=nrcmt(is)
  nrci=nrcmtinr(is)
! convert the density to spherical harmonics
  call rfcpy(nr,nri,rhomt(:,:,ias),rfmt)
  call rfsht(nrc,nrci,1,rfmt,lradstp,rhomt(:,:,ias))
! convert magnetisation to spherical harmonics
  if (spinpol) then
    do idm=1,ndmag
      call rfcpy(nr,nri,magmt(:,:,ias,idm),rfmt)
      call rfsht(nrc,nrci,1,rfmt,lradstp,magmt(:,:,ias,idm))
    end do
  end if
  deallocate(rfmt)
end do
!$OMP END DO
!$OMP END PARALLEL
return

contains

subroutine rfcpy(nr,nri,rfmt1,rfmt2)
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: rfmt1(lmmaxvr,nrmtmax)
real(8), intent(out) :: rfmt2(lmmaxvr,nrcmtmax)
! local variables
integer ir,irc
irc=0
do ir=1,nri,lradstp
  irc=irc+1
  call dcopy(lmmaxinr,rfmt1(:,ir),1,rfmt2(:,irc),1)
end do
do ir=nri+lradstp,nr,lradstp
  irc=irc+1
  call dcopy(lmmaxvr,rfmt1(:,ir),1,rfmt2(:,irc),1)
end do
return
end subroutine

end subroutine
!EOC

