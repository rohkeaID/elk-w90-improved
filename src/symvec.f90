
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine symvec(vc)
use modmain
implicit none
! arguments
real(8), intent(inout) :: vc(3)
! local variables
integer isym,ls
real(8) vs(3),v(3),t1
! make symmetric average
vs(:)=0.d0
do isym=1,nsymcrys
  ls=lsplsymc(isym)
  call r3mv(symlatc(:,:,ls),vc,v)
  vs(:)=vs(:)+v(:)
end do
! normalise
t1=1.d0/dble(nsymcrys)
vc(:)=t1*vs(:)
return
end subroutine

