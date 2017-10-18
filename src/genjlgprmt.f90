
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genjlgprmt(lmax,ngp,gpc,ld,jlgprmt)
use modmain
implicit none
! arguments
integer, intent(in) :: lmax,ngp
real(8), intent(in) :: gpc(ngp)
integer, intent(in) :: ld
real(8), intent(out) :: jlgprmt(0:lmax,ld,nspecies)
! local variables
integer is,ig
real(8) t1
do is=1,nspecies
  do ig=1,ngp
    t1=gpc(ig)*rmt(is)
    call sbessel(lmax,t1,jlgprmt(:,ig,is))
  end do
end do
return
end subroutine

