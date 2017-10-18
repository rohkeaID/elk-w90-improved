
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine drhomagsh
use modmain
use modphonon
implicit none
! local variables
integer idm,is,ias
integer nrc,nrci,npc
! allocatable arrays
complex(8), allocatable :: zfmt(:)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zfmt,is,nrc,nrci) &
!$OMP PRIVATE(npc,idm)
!$OMP DO
do ias=1,natmtot
  allocate(zfmt(npcmtmax))
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
! convert the density derivative to spherical harmonics
  zfmt(1:npc)=drhomt(1:npc,ias)
  call zfsht(nrc,nrci,zfmt,drhomt(:,ias))
! convert the magnetisation derivative to spherical harmonics
  do idm=1,ndmag
    zfmt(1:npc)=dmagmt(1:npc,ias,idm)
    call zfsht(nrc,nrci,zfmt,dmagmt(:,ias,idm))
  end do
  deallocate(zfmt)
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine

