
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rhomagsh
! !INTERFACE:
subroutine rhomagsh
! !USES:
use modmain
use modomp
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
integer idm,is,ias,nthd
integer nrc,nrci,npc
! allocatable arrays
real(8), allocatable :: rfmt(:)
call omp_hold(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt,is,nrc,nrci,npc,idm) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ias=1,natmtot
  allocate(rfmt(npcmtmax))
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
! convert the density to spherical harmonics
  call dcopy(npc,rhomt(:,ias),1,rfmt,1)
  call rfsht(nrc,nrci,rfmt,rhomt(:,ias))
! convert magnetisation to spherical harmonics
  do idm=1,ndmag
    call dcopy(npc,magmt(:,ias,idm),1,rfmt,1)
    call rfsht(nrc,nrci,rfmt,magmt(:,ias,idm))
  end do
  deallocate(rfmt)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
return
end subroutine
!EOC

