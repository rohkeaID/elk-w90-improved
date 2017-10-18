
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
integer nrc,nrci,npc
! allocatable arrays
real(8), allocatable :: rfmt(:)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt,is,nrc,nrci,npc,idm)
!$OMP DO
do ias=1,natmtot
  allocate(rfmt(npcmtmax))
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
! convert the density to spherical harmonics
  rfmt(1:npc)=rhomt(1:npc,ias)
  call rfsht(nrc,nrci,rfmt,rhomt(:,ias))
! convert magnetisation to spherical harmonics
  do idm=1,ndmag
    rfmt(1:npc)=magmt(1:npc,ias,idm)
    call rfsht(nrc,nrci,rfmt,magmt(:,ias,idm))
  end do
  deallocate(rfmt)
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine
!EOC

