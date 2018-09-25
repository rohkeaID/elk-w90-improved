
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gentimes
use modmain
use modtddft
implicit none
! local variables
integer its
! number of time steps
ntimes=int(tstime/dtimes)+1
! generate the time-step array
if (allocated(times)) deallocate(times)
allocate(times(ntimes))
do its=1,ntimes
  times(its)=dble(its-1)*dtimes
end do
return
end subroutine

