
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rzfadd(za,zfmt,zfir,rfmt,rfir)
use modmain
implicit none
! arguments
complex(8), intent(in) :: za
complex(8), intent(in) :: zfmt(npcmtmax,natmtot),zfir(ngtot)
real(8), intent(inout) :: rfmt(npcmtmax,natmtot),rfir(ngtot)
! local variables
integer ias,is,npc
! add in muffin-tin region
do ias=1,natmtot
  is=idxis(ias)
  npc=npcmt(is)
  rfmt(1:npc,ias)=rfmt(1:npc,ias)+dble(za*zfmt(1:npc,ias))
end do
! add in interstitial region
rfir(:)=rfir(:)+dble(za*zfir(:))
return
end subroutine

