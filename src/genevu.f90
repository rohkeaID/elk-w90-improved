
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genevu
use modmain
use modultra
implicit none
! local variables
integer ik0,ist
if (.not.ultracell) return

! loop over original k-points
do ik0=1,nkpt0
! loop over second-variational states
  do ist=1,nstsv
! solve the third-variational ultracell eigenvalue equation
    call eveqnu(ik0,ist)
  end do
end do


return
end subroutine

