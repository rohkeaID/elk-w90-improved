
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genws
use modmain
use modomp
implicit none
! local variables
integer is,ias,nthd
integer nr,nri,nrc,nrci
! allocatable arrays
real(8), allocatable :: rfmt(:)
if (xcgrad.ne.4) return
! muffin-tin effective tau-DFT potential
call omp_hold(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt,is,nr,nri) &
!$OMP PRIVATE(nrc,nrci) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ias=1,natmtot
  allocate(rfmt(npcmtmax))
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
! convert to coarse radial mesh and spherical coordinates
  call rfmtftoc(nr,nri,wxcmt(:,ias),rfmt)
  call rbsht(nrc,nrci,rfmt,wsmt(:,ias))
  deallocate(rfmt)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
! interstitial tau-DFT potential
wsir(:)=wxcir(:)*cfunir(:)
return
end subroutine

