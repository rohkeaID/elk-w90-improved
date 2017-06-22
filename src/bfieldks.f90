
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine bfieldks
use modmain
implicit none
! local variables
integer idm,is,ia,ias,nrc
real(8) cb,t1
if (.not.spinpol) return
! coupling constant of the external field (g_e/4c)
cb=gfacte/(4.d0*solsc)
!------------------------------------!
!     muffin-tin Kohn-Sham field     !
!------------------------------------!
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(is,ia,nrc,idm,t1)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  ia=idxia(ias)
  nrc=nrcmt(is)
! exchange-correlation magnetic field in spherical coordinates
  do idm=1,ndmag
    call rbsht(nrc,nrcmtinr(is),lradstp,bxcmt(:,:,ias,idm),1,bsmt(:,:,ias,idm))
  end do
! add the external magnetic field
  t1=cb*(bfcmt(3,ia,is)+bfieldc(3))
  bsmt(:,1:nrc,ias,ndmag)=bsmt(:,1:nrc,ias,ndmag)+t1
  if (ncmag) then
    do idm=1,2
      t1=cb*(bfcmt(idm,ia,is)+bfieldc(idm))
      bsmt(:,1:nrc,ias,idm)=bsmt(:,1:nrc,ias,idm)+t1
    end do
  end if
end do
!$OMP END DO
!$OMP END PARALLEL
!-----------------------------------------------!
!     interstitial Kohn-Sham magnetic field     !
!-----------------------------------------------!
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(t1)
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
return
end subroutine

