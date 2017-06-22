
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine potkst
use modmain
use modtddft
implicit none
! only adiabatic approximation currently available
call potks
return
end subroutine

