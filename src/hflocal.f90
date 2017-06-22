
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hflocal(vmt,vir,bmt,bir)
use modmain
implicit none
! arguments
real(8), intent(out) :: vmt(lmmaxvr,nrcmtmax,natmtot),vir(ngtot)
real(8), intent(out) :: bmt(lmmaxvr,nrcmtmax,natmtot,ndmag),bir(ngtot,ndmag)
! local variables
integer idm,is,ias,ir,irc
! allocatable arrays
real(8), allocatable :: rfmt(:,:)
! compute the Coulomb potential
call potcoul
! generate the exchange-correlation potentials for hybrids
if (hybrid) call potxc
! convert to spherical coordinates and store in output arrays
if (hybrid) then
! hybrid functional case
  allocate(rfmt(lmmaxvr,nrcmtmax))
  do ias=1,natmtot
    is=idxis(ias)
    irc=0
    do ir=1,nrmt(is),lradstp
      irc=irc+1
      rfmt(:,irc)=vclmt(:,ir,ias)+vxcmt(:,ir,ias)
    end do
    call rbsht(nrcmt(is),nrcmtinr(is),1,rfmt,1,vmt(:,:,ias))
  end do
  deallocate(rfmt)
  vir(:)=(vclir(:)+vxcir(:))*cfunir(:)
  if (spinpol) then
    do idm=1,ndmag
      do ias=1,natmtot
        is=idxis(ias)
        call rbsht(nrcmt(is),nrcmtinr(is),lradstp,bxcmt(:,:,ias,idm),1, &
         bmt(:,:,ias,idm))
      end do
      bir(:,idm)=bxcir(:,idm)*cfunir(:)
    end do
  end if
else
! normal Hartree-Fock case
  do ias=1,natmtot
    is=idxis(ias)
    call rbsht(nrcmt(is),nrcmtinr(is),lradstp,vclmt(:,:,ias),1,vmt(:,:,ias))
  end do
  vir(:)=vclir(:)*cfunir(:)
end if
return
end subroutine

