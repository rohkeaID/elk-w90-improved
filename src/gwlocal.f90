
! Copyright (C) 2017 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine gwlocal(vmt,vir,bmt,bir)
use modmain
implicit none
! arguments
real(8), intent(out) :: vmt(npcmtmax,natmtot),vir(ngtot)
real(8), intent(out) :: bmt(npcmtmax,natmtot,ndmag),bir(ngtot,ndmag)
! local variables
integer idm,is,ias
! automatic arrays
real(8) rfmt(npcmtmax)
do ias=1,natmtot
  is=idxis(ias)
! convert exchange-correlation potential to a coarse radial mesh
  call rfmtftoc(nrmt(is),nrmti(is),vxcmt(:,ias),rfmt)
! negate because V_xc should be removed from the self-energy
  rfmt(1:npcmt(is))=-rfmt(1:npcmt(is))
! convert to spherical coordinates
  call rbsht(nrcmt(is),nrcmti(is),rfmt,vmt(:,ias))
end do
! negate and multiply the interstitial V_xc by the characteristic function
vir(:)=-vxcir(:)*cfunir(:)
! do the same for B_xc in the spin-polarised case
if (spinpol) then
  do idm=1,ndmag
    do ias=1,natmtot
      is=idxis(ias)
      call rfmtftoc(nrmt(is),nrmti(is),bxcmt(:,ias,idm),rfmt)
      rfmt(1:npcmt(is))=-rfmt(1:npcmt(is))
      call rbsht(nrcmt(is),nrcmti(is),rfmt,bmt(:,ias,idm))
    end do
    bir(:,idm)=-bxcir(:,idm)*cfunir(:)
  end do
end if
return
end subroutine

