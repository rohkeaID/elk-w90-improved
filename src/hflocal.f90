
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hflocal(hyb,vmt,vir,bmt,bir)
use modmain
implicit none
! arguments
logical, intent(in) :: hyb
real(8), intent(out) :: vmt(npcmtmax,natmtot),vir(ngtot)
real(8), intent(out) :: bmt(npcmtmax,natmtot,ndmag),bir(ngtot,ndmag)
! local variables
integer idm,is,ias,np
! automatic arrays
real(8) rfmt1(npmtmax),rfmt2(npcmtmax)
! compute the Coulomb potential
call potcoul
! convert to spherical coordinates and store in output arrays
if (hyb) then
! hybrid functional case
  call potxc
  do ias=1,natmtot
    is=idxis(ias)
    np=npmt(is)
    rfmt1(1:np)=vclmt(1:np,ias)+vxcmt(1:np,ias)
    call rfmtftoc(nrmt(is),nrmti(is),rfmt1,rfmt2)
    call rbsht(nrcmt(is),nrcmti(is),rfmt2,vmt(:,ias))
  end do
  vir(:)=(vclir(:)+vxcir(:))*cfunir(:)
  if (spinpol) then
    do idm=1,ndmag
      do ias=1,natmtot
        is=idxis(ias)
        call rfmtftoc(nrmt(is),nrmti(is),bxcmt(:,ias,idm),rfmt1)
        call rbsht(nrcmt(is),nrcmti(is),rfmt1,bmt(:,ias,idm))
      end do
      bir(:,idm)=bxcir(:,idm)*cfunir(:)
    end do
  end if
else
! normal Hartree-Fock case
  do ias=1,natmtot
    is=idxis(ias)
    call rfmtftoc(nrmt(is),nrmti(is),vclmt(:,ias),rfmt1)
    call rbsht(nrcmt(is),nrcmti(is),rfmt1,vmt(:,ias))
  end do
  vir(:)=vclir(:)*cfunir(:)
end if
return
end subroutine

