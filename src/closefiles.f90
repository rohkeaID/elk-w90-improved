
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine closefiles
use modmain
implicit none
! local variables
logical opnd
integer iu
! loop over possible file units used by the code
do iu=20,1000
  inquire(unit=iu,opened=opnd)
  if (opnd) close(iu)
end do
return
end subroutine

