
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine genbs
use modmain
use modomp
implicit none
! local variables
integer idm,is,ia,ias,nthd
integer nr,nri,nrc,nrci,npc
real(8) cb,t1
! allocatable arrays
real(8), allocatable :: rfmt(:)
if (.not.spinpol) return
! coupling constant of the external field (g_e/4c)
cb=gfacte/(4.d0*solsc)
!------------------------------------!
!     muffin-tin Kohn-Sham field     !
!------------------------------------!
call omp_hold(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt,is,ia,nr,nri) &
!$OMP PRIVATE(nrc,nrci,npc,idm,t1) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ias=1,natmtot
  allocate(rfmt(npcmtmax))
  is=idxis(ias)
  ia=idxia(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
! exchange-correlation magnetic field in spherical coordinates
  do idm=1,ndmag
    call rfmtftoc(nr,nri,bxcmt(:,ias,idm),rfmt)
    call rbsht(nrc,nrci,rfmt,bsmt(:,ias,idm))
  end do
! add the external magnetic field
  t1=cb*(bfcmt(3,ia,is)+bfieldc(3))
  bsmt(1:npc,ias,ndmag)=bsmt(1:npc,ias,ndmag)+t1
  if (ncmag) then
    do idm=1,2
      t1=cb*(bfcmt(idm,ia,is)+bfieldc(idm))
      bsmt(1:npc,ias,idm)=bsmt(1:npc,ias,idm)+t1
    end do
  end if
  deallocate(rfmt)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
!-----------------------------------------------!
!     interstitial Kohn-Sham magnetic field     !
!-----------------------------------------------!
call omp_hold(ndmag,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(t1) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do idm=1,ndmag
  if (ncmag) then
    t1=cb*bfieldc(idm)
  else
    t1=cb*bfieldc(3)
  end if
! multiply by characteristic function
  bsir(:,idm)=(bxcir(:,idm)+t1)*cfunir(:)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
return
end subroutine

