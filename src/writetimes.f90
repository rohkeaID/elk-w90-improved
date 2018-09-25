
! Copyright (C) 2015 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writetimes
use modmain
use modtddft
implicit none
open(50,file='TIMESTEP.OUT',form='FORMATTED')
write(50,'(I8,G18.10)') itimes,times(itimes)
close(50)
return
end subroutine

